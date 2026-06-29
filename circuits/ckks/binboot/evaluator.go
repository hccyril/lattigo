// Package binboot 的评估器实现。
//
// 本文件实现论文 §3.1 Algorithm 2（BinBoot）和 §4.1 Algorithm 3（GateBoot）
// 的核心评估逻辑。两个算法均复用 Lattigo 标准自举流程的 ScaleDown、ModUp、
// CoeffsToSlots 和 SlotsToCoeffs 步骤，但在 EvalMod 阶段替换为论文提出的
// 三角函数 Chebyshev 逼近。
//
// 论文对应关系：
//   - Algorithm 2 (§3.1): BinBoot = StC ∘ Eval_fBinBoot ∘ ModRaise ∘ CtS
//   - Algorithm 3 (§4.1): GateBoot = StC ∘ Eval_fG ∘ ModRaise ∘ CtS(ct1+ct2)
//   - Theorem 1 (§3.2): BinBoot 二次噪声缩减
//   - Theorem 2 (§4.2): GateBoot 正确性
package binboot

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// ===================================================================
// 评估器结构体
// ===================================================================

// Evaluator 包装 Lattigo 标准 bootstrapping.Evaluator，在其基础上提供
// BinBoot（论文 Algorithm 2）和 GateBoot（论文 Algorithm 3）方法。
//
// 核心思路：复用标准自举的 ScaleDown/ModUp/CoeffsToSlots/SlotsToCoeffs，
// 仅替换 EvalMod 阶段的多项式逼近函数。
//
// 论文 §3.1: "The bootstrapping procedure is the same as the standard CKKS
// bootstrapping, except that the EvalMod step uses f_BinBoot instead of
// the sine/cosine approximation."
type Evaluator struct {
	// BTEvaluator 是标准 Lattigo 自举评估器
	BTEvaluator *bootstrapping.Evaluator

	// PolyEvaluator 用于 Chebyshev 多项式同态求值
	// 从 BTEvaluator.Mod1Evaluator.PolynomialEvaluator 获取
	PolyEvaluator *polynomial.Evaluator

	// BinBootPoly 是 f_BinBoot 的 Chebyshev 逼近多项式
	// 论文 §3.1: f_BinBoot(x) = (1 - cos(2πx)) / 2
	BinBootPoly bignum.Polynomial

	// GatePolys 预计算六个门函数的 Chebyshev 逼近
	// 论文 §4.1 Table 4
	GatePolys map[GateType]bignum.Polynomial

	// BootstrappingParameters 缓存自举参数
	BootstrappingParameters ckks.Parameters
}

// NewEvaluator 创建一个新的 BinBoot/GateBoot 评估器。
//
// 参数:
//   - btpEval: 已初始化的标准 bootstrapping.Evaluator
//   - degree: Chebyshev 多项式次数（论文实验中使用 30）
//   - K: 逼近区间 [-K, K]（论文中使用 4）
//
// 返回的评估器预计算了 BinBoot 和六个 GateBoot 的多项式逼近，
// 后续 Bootstrap/GateBoot 调用时无需重复计算。
//
// 论文 §3.1: "Eval_fBinBoot is the homomorphic evaluation of f_BinBoot
// via appropriate polynomial approximation."
// 论文 §4.1: "Step 4 consists in homomorphically evaluating a trigonometric
// function f_G that removes I and sends φ1+φ2 to G(φ1, φ2)."
func NewEvaluator(btpEval *bootstrapping.Evaluator, degree int, K float64) (*Evaluator, error) {

	// 从标准自举评估器获取多项式评估器
	// mod1.Evaluator 内部持有 polynomial.Evaluator
	polyEval := btpEval.Mod1Evaluator.PolynomialEvaluator

	// 【关键修复】获取 ScaleDown 的 MessageRatio（= Q0/Δ = 2^LogMessageRatio）
	// ScaleDown 会将消息归一化为 slot_values/MessageRatio，因此多项式需要
	// 逼近 f(MessageRatio * x) 以补偿归一化，确保 P(slot_values/MessageRatio + I) = f(slot_values)
	// 详见 NewBinBootPoly 和 NewGatePoly 的注释
	messageRatio := btpEval.Mod1Parameters.MessageRatio()

	// 预计算 BinBoot 多项式（包含 MessageRatio 补偿）
	binBootPoly := NewBinBootPoly(degree, K, messageRatio)

	// 预计算六个门函数多项式（包含 MessageRatio 补偿）
	gatePolys := make(map[GateType]bignum.Polynomial)
	for _, g := range []GateType{GateAND, GateOR, GateXOR, GateNAND, GateNOR, GateXNOR} {
		gatePolys[g] = NewGatePoly(g, degree, K, messageRatio)
	}

	return &Evaluator{
		BTEvaluator:             btpEval,
		PolyEvaluator:           polyEval,
		BinBootPoly:             binBootPoly,
		GatePolys:               gatePolys,
		BootstrappingParameters: btpEval.BootstrappingParameters,
	}, nil
}

// ===================================================================
// 论文 Algorithm 2: BinBoot
// ===================================================================

// Bootstrap 对单个密文执行 BinBoot 自举（论文 Algorithm 2）。
//
// 论文 §3.1 Algorithm 2:
//
//	Input: ct ∈ R_q^n (加密二进制消息 b/2 + ε)
//	1. CtS: CoeffsToSlots → 获取 slot 表示
//	2. ModRaise: 提升模数从 q 到 Q
//	3. StC: SlotsToCoeffs → 获取系数表示（此时消息为 (b+ε')/2 + I）
//	4. Eval_fBinBoot: 同态求值 f_BinBoot(x) = (1-cos(2πx))/2
//	Output: ct' ∈ R_Q^n (加密 b + ε''，其中 ε'' = O(ε²))
//
// 与标准 CKKS 自举的区别：
//   - EvalMod 阶段使用 f_BinBoot 替代 sine/cosine 逼近
//   - f_BinBoot 的导数在整数点处为零 → 二次噪声缩减（Theorem 1）
//   - 允许 q0 = 2Δ0（而非标准自举的 q0 >> Δ0）
//
// 【关键修复】ConjugateInvariant 环切换：
//
//	残差环为 ConjugateInvariant（N1=2^LogN），自举环为 Standard（N2=2^(LogN+1)）。
//	标准 Lattigo 自举的 ScaleDown/ModUp 操作在自举环（N2）上执行，
//	因此必须先通过 RealToComplexNew 将密文从残差环（N1）切换到自举环（N2），
//	自举完成后再通过 ComplexToRealNew 切换回残差环（N1）。
//	这与标准 EvaluateConjugateInvariant 的流程完全一致。
//
// 参数:
//   - ct: 待自举的密文（残差环 ConjugateInvariant），应加密二进制消息（值为 0 或 1/2）
//
// 返回: 自举后的密文（残差环），消息近似为 0 或 1（乘以 Δ0 后）
func (eval *Evaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	btpEval := eval.BTEvaluator

	// ====== 环切换：残差环(N1, ConjugateInvariant) → 自举环(N2, Standard) ======
	// 论文 Algorithm 2 的操作在自举环上执行，需要先将密文从 ConjugateInvariant
	// 环切换到 Standard 环。RealToComplexNew 通过 DomainSwitcher 完成此切换，
	// 将 N1 维度的实数密文扩展为 N2 维度的复数密文（虚部为零）。
	ctN2 := btpEval.RealToComplexNew(ct.CopyNew())

	// ====== 论文 Algorithm 2 Step 1-3: ScaleDown → ModUp → CoeffsToSlots ======
	// 这些步骤在自举环(N2)上执行，与标准自举流程完全一致

	// Step 1: ScaleDown — 将密文降至 level 0，scale 调整为 Q[0]/MessageRatio
	ctScaled, _, err := btpEval.ScaleDown(ctN2)
	if err != nil {
		return nil, fmt.Errorf("binboot: ScaleDown failed: %w", err)
	}

	// Step 2: ModUp — 将模数从 q 提升到 Q（含密钥切换）
	// 此步骤在自举环(N2)上执行，INTT 使用 N2 维度的 SubRing，不再有维度不匹配问题
	ctModUp, err := btpEval.ModUp(ctScaled)
	if err != nil {
		return nil, fmt.Errorf("binboot: ModUp failed: %w", err)
	}

	// Step 3: CoeffsToSlots — 同态编码，将系数表示转为 slot 表示
	ctReal, ctImag, err := btpEval.CoeffsToSlots(ctModUp)
	if err != nil {
		return nil, fmt.Errorf("binboot: CoeffsToSlots failed: %w", err)
	}

	// ====== 论文 Algorithm 2 Step 4: Eval_fBinBoot ======
	// 用论文三角函数 f_BinBoot(x)=(1-cos(2πx))/2 替代标准 sine/cosine 逼近
	// 论文 §3.1: "Eval_fBinBoot is the homomorphic evaluation of f_BinBoot"
	if ctReal, err = eval.evalBinBootPoly(ctReal); err != nil {
		return nil, fmt.Errorf("binboot: EvalBinBootPoly (real) failed: %w", err)
	}

	if ctImag != nil {
		if ctImag, err = eval.evalBinBootPoly(ctImag); err != nil {
			return nil, fmt.Errorf("binboot: EvalBinBootPoly (imag) failed: %w", err)
		}
	}

	// Step 5: SlotsToCoeffs — 同态解码，将 slot 表示转回系数表示
	ctOutN2, err := btpEval.SlotsToCoeffs(ctReal, ctImag)
	if err != nil {
		return nil, fmt.Errorf("binboot: SlotsToCoeffs failed: %w", err)
	}

	// ====== scale 补偿 ======
	// ComplexToRealNew 会将 scale 乘以 2，因此需要预先乘以 0.5 补偿
	// 这与标准 EvaluateConjugateInvariant 的处理一致（evaluator.go:494）
	ctOutN2.Scale = ctOutN2.Scale.Mul(rlwe.NewScale(1.0 / 2.0))

	// ====== 环切换：自举环(N2, Standard) → 残差环(N1, ConjugateInvariant) ======
	// 将自举后的密文从 Standard 环切换回 ConjugateInvariant 环
	ctOut := btpEval.ComplexToRealNew(ctOutN2)

	// 设置输出 scale 为残差参数的默认 scale
	ctOut.Scale = btpEval.ResidualParameters.DefaultScale()

	return ctOut, nil
}

// evalBinBootPoly 在密文上同态求值 f_BinBoot 的 Chebyshev 多项式。
//
// 这是论文 Algorithm 2 的 Step 4 核心实现。
// 与标准 mod1.EvaluateNew 的区别：
//   - 不使用 double-angle 公式（f_BinBoot 本身已包含 cos）
//   - 不使用 arcsine 逆函数
//   - 直接在 [-K, K] 上评估 Chebyshev 多项式
//
// 论文 §3.1: "Eval_fBinBoot is the homomorphic evaluation of f_BinBoot(x)=(1-cos(2x))/2
// via appropriate polynomial approximation."
//
// 论文 §3.2 Theorem 1: f_BinBoot 的导数在 x ∈ (1/2)Z 处为零，
// 因此 f_BinBoot(b/2 + I + ε/2) - b = O(ε²)，实现二次噪声缩减。
func (eval *Evaluator) evalBinBootPoly(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	poly := eval.BinBootPoly

	// 获取 Chebyshev 变量代换参数
	// Chebyshev 基要求将输入从 [A, B] 映射到 [-1, 1]
	// 变换: x' = scalar * x + constant
	//   scalar = 2 / (B - A)
	//   constant = (-A - B) / (B - A)
	// 对于对称区间 [-K, K]：scalar = 1/K, constant = 0
	scalar, constant := poly.ChangeOfBasis()

	ct := ctIn.CopyNew()

	// 【修复】正确的变量代换顺序：先乘 scalar，再 Rescale，最后加 constant
	// 之前代码先 Add 再 Mul，导致结果为 scalar*(ct+constant) 而非 scalar*ct+constant
	// 对于对称区间 constant=0 时两者等价，但为通用正确性修正执行顺序

	// Step 1: 乘以缩放因子 scalar — 将输入区间 [-K, K] 映射到 [-1, 1]
	// Mul by float64 scalar: scale 会乘以当前模数，需要后续 Rescale 恢复
	scalarFloat, _ := scalar.Float64()
	if err := eval.BTEvaluator.Mul(ct, scalarFloat, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev scaling: %w", err)
	}

	// Step 2: Rescale — 消耗一个层级，将 scale 恢复到合理范围
	if err := eval.BTEvaluator.Rescale(ct, ct); err != nil {
		return nil, fmt.Errorf("cannot rescale after Chebyshev scaling: %w", err)
	}

	// Step 3: 加上常数偏移 constant（对称区间时为 0，不影响结果）
	if err := eval.BTEvaluator.Add(ct, constant, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev offset: %w", err)
	}

	// 目标输出 scale: 自举参数的默认 scale
	targetScale := eval.BootstrappingParameters.DefaultScale()

	// 使用 polynomial.Evaluator.Evaluate 评估 Chebyshev 多项式
	// 论文 §3.1: 通过 Chebyshev 逼近同态求值 f_BinBoot
	ctOut, err := eval.PolyEvaluator.Evaluate(ct, poly, targetScale)
	if err != nil {
		return nil, fmt.Errorf("cannot evaluate BinBoot polynomial: %w", err)
	}

	return ctOut, nil
}

// ===================================================================
// 论文 Algorithm 3: GateBoot
// ===================================================================

// GateBoot 对两个密文执行 GateBoot 自举（论文 Algorithm 3），
// 在自举的同时完成二元门运算。
//
// 论文 §4.1 Algorithm 3:
//
//	Input: ct1, ct2 ∈ R_q^n (分别加密 b1/2 + ε1, b2/2 + ε2)
//	1. ct = ct1 + ct2 (加密 (b1+b2)/2 + ε1+ε2)
//	2. CtS ∘ ModRaise ∘ StC → 系数表示 (b1+b2+ε')/2 + I
//	   注意: GateBoot 中输入为 (b1+b2)/3 + I（因为 q0=3Δ0）
//	3. Eval_fG: 同态求值门函数 f_G(x)
//	Output: ct' ∈ R_Q^n (加密 G(b1, b2) + ε'')
//
// 论文 §4.1: "f_G((b1+b2)/3 + I) = G(b1, b2)"
// 三个关键点 x = 0, 1/3, 2/3 对应 (b1,b2) = (0,0), (1,0), (1,1)
//
// 论文 §4.2 Theorem 2: GateBoot 的正确性保证
//
// 【关键修复】ConjugateInvariant 环切换：
//
//	与 Bootstrap 相同，GateBoot 也需要在自举环(N2)上执行 ScaleDown/ModUp 等操作。
//	先在残差环(N1)上完成 ct1+ct2 加法，然后切换到自举环(N2)执行自举电路，
//	最后切换回残差环(N1)。
//
// 参数:
//   - ct1: 第一个输入密文（残差环 ConjugateInvariant，加密 b1/2 或 b1/3）
//   - ct2: 第二个输入密文（残差环 ConjugateInvariant，加密 b2/2 或 b2/3）
//   - gate: 门类型（GateNAND, GateAND 等）
//
// 返回: 自举后的密文（残差环），消息近似为 G(b1, b2)（0 或 1）
func (eval *Evaluator) GateBoot(ct1, ct2 *rlwe.Ciphertext, gate GateType) (*rlwe.Ciphertext, error) {

	// 检查门类型有效
	if gate < GateAND || gate > GateXNOR {
		return nil, fmt.Errorf("binboot: invalid gate type %d", gate)
	}

	btpEval := eval.BTEvaluator

	// ====== 论文 Algorithm 3 Step 0: ct = ct1 + ct2 ======
	// 论文 §4.1: "Input: ct1, ct2; Step 1: ct = ct1 + ct2"
	// 加法在残差环(N1, ConjugateInvariant)上执行，维度为 N1
	ctAdd := ct1.CopyNew()
	if err := eval.BTEvaluator.Add(ctAdd, ct2, ctAdd); err != nil {
		return nil, fmt.Errorf("binboot: cannot add ciphertexts: %w", err)
	}

	// ====== 环切换：残差环(N1, ConjugateInvariant) → 自举环(N2, Standard) ======
	// 将 ct1+ct2 的结果从 ConjugateInvariant 环切换到 Standard 环
	ctN2 := btpEval.RealToComplexNew(ctAdd)

	// ====== 论文 Algorithm 3 Step 1-3: ScaleDown → ModUp → CoeffsToSlots ======
	// 后续步骤在自举环(N2)上执行，与 BinBoot 相同

	// Step 1: ScaleDown
	ctScaled, _, err := btpEval.ScaleDown(ctN2)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot ScaleDown failed: %w", err)
	}

	// Step 2: ModUp — 在自举环(N2)上执行，维度匹配
	ctModUp, err := btpEval.ModUp(ctScaled)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot ModUp failed: %w", err)
	}

	// Step 3: CoeffsToSlots
	ctReal, ctImag, err := btpEval.CoeffsToSlots(ctModUp)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot CoeffsToSlots failed: %w", err)
	}

	// ====== 论文 Algorithm 3 Step 4: Eval_fG ======
	// 使用门函数三角多项式替换标准 EvalMod
	// 论文 §4.1: "Step 4 consists in homomorphically evaluating a trigonometric
	// function f_G that removes I and sends φ1+φ2 to G(φ1, φ2)."
	if ctReal, err = eval.evalGatePoly(ctReal, gate); err != nil {
		return nil, fmt.Errorf("binboot: GateBoot EvalGatePoly (real) failed: %w", err)
	}

	if ctImag != nil {
		if ctImag, err = eval.evalGatePoly(ctImag, gate); err != nil {
			return nil, fmt.Errorf("binboot: GateBoot EvalGatePoly (imag) failed: %w", err)
		}
	}

	// Step 5: SlotsToCoeffs
	ctOutN2, err := btpEval.SlotsToCoeffs(ctReal, ctImag)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot SlotsToCoeffs failed: %w", err)
	}

	// ====== scale 补偿 ======
	// ComplexToRealNew 会将 scale 乘以 2，预先乘以 0.5 补偿
	ctOutN2.Scale = ctOutN2.Scale.Mul(rlwe.NewScale(1.0 / 2.0))

	// ====== 环切换：自举环(N2, Standard) → 残差环(N1, ConjugateInvariant) ======
	ctOut := btpEval.ComplexToRealNew(ctOutN2)

	// 设置输出 scale 为残差参数的默认 scale
	ctOut.Scale = btpEval.ResidualParameters.DefaultScale()

	return ctOut, nil
}

// evalGatePoly 在密文上同态求值门函数 f_G 的 Chebyshev 多项式。
//
// 这是论文 Algorithm 3 的 Step 4 核心实现。
// 论文 §4.1 Table 4 定义了六个门的三角函数：
//   - f_G((b1+b2)/3 + I) = G(b1, b2)
//   - 三个关键点: x=0 → G(0,0), x=1/3 → G(1,0), x=2/3 → G(1,1)
//
// 与 evalBinBootPoly 的区别：
//   - 使用门函数多项式而非 BinBoot 多项式
//   - 门函数一般不是偶函数，不清理奇/偶系数
func (eval *Evaluator) evalGatePoly(ctIn *rlwe.Ciphertext, gate GateType) (*rlwe.Ciphertext, error) {

	poly, ok := eval.GatePolys[gate]
	if !ok {
		return nil, fmt.Errorf("binboot: gate polynomial for %s not found", gate.String())
	}

	// Chebyshev 变量代换（与 evalBinBootPoly 相同的正确顺序）
	scalar, constant := poly.ChangeOfBasis()

	ct := ctIn.CopyNew()

	// 【修复】正确的变量代换顺序：先乘 scalar，再 Rescale，最后加 constant
	// Step 1: 乘以缩放因子 scalar
	scalarFloat, _ := scalar.Float64()
	if err := eval.BTEvaluator.Mul(ct, scalarFloat, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev scaling: %w", err)
	}

	// Step 2: Rescale
	if err := eval.BTEvaluator.Rescale(ct, ct); err != nil {
		return nil, fmt.Errorf("cannot rescale after Chebyshev scaling: %w", err)
	}

	// Step 3: 加上常数偏移 constant
	if err := eval.BTEvaluator.Add(ct, constant, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev offset: %w", err)
	}

	// 评估门函数 Chebyshev 多项式
	targetScale := eval.BootstrappingParameters.DefaultScale()

	ctOut, err := eval.PolyEvaluator.Evaluate(ct, poly, targetScale)
	if err != nil {
		return nil, fmt.Errorf("cannot evaluate gate polynomial: %w", err)
	}

	return ctOut, nil
}

// ===================================================================
// 便捷方法
// ===================================================================

// GateBootNAND 是 GateBoot(ct1, ct2, GateNAND) 的便捷包装。
// 论文 §4.1: NAND 是最常用的门，因为它是功能完备的。
func (eval *Evaluator) GateBootNAND(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateNAND)
}

// GateBootAND 是 GateBoot(ct1, ct2, GateAND) 的便捷包装。
func (eval *Evaluator) GateBootAND(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateAND)
}

// GateBootOR 是 GateBoot(ct1, ct2, GateOR) 的便捷包装。
func (eval *Evaluator) GateBootOR(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateOR)
}

// GateBootXOR 是 GateBoot(ct1, ct2, GateXOR) 的便捷包装。
func (eval *Evaluator) GateBootXOR(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateXOR)
}

// GateBootNOR 是 GateBoot(ct1, ct2, GateNOR) 的便捷包装。
func (eval *Evaluator) GateBootNOR(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateNOR)
}

// GateBootXNOR 是 GateBoot(ct1, ct2, GateXNOR) 的便捷包装。
func (eval *Evaluator) GateBootXNOR(ct1, ct2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return eval.GateBoot(ct1, ct2, GateXNOR)
}

// ===================================================================
// 辅助：big.Float 转换
// ===================================================================

// bigFloatToFloat64 将 big.Float 转换为 float64，用于常数乘法。
func bigFloatToFloat64(x *big.Float) float64 {
	f, _ := x.Float64()
	return f
}
