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

	// 预计算 BinBoot 多项式
	binBootPoly := NewBinBootPoly(degree, K)

	// 预计算六个门函数多项式
	gatePolys := make(map[GateType]bignum.Polynomial)
	for _, g := range []GateType{GateAND, GateOR, GateXOR, GateNAND, GateNOR, GateXNOR} {
		gatePolys[g] = NewGatePoly(g, degree, K)
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
//   Input: ct ∈ R_q^n (加密二进制消息 b/2 + ε)
//   1. CtS: CoeffsToSlots → 获取 slot 表示
//   2. ModRaise: 提升模数从 q 到 Q
//   3. StC: SlotsToCoeffs → 获取系数表示（此时消息为 (b+ε')/2 + I）
//   4. Eval_fBinBoot: 同态求值 f_BinBoot(x) = (1-cos(2πx))/2
//   Output: ct' ∈ R_Q^n (加密 b + ε''，其中 ε'' = O(ε²))
//
// 与标准 CKKS 自举的区别：
//   - EvalMod 阶段使用 f_BinBoot 替代 sine/cosine 逼近
//   - f_BinBoot 的导数在整数点处为零 → 二次噪声缩减（Theorem 1）
//   - 允许 q0 = 2Δ0（而非标准自举的 q0 >> Δ0）
//
// 参数:
//   - ct: 待自举的密文，应加密二进制消息（值为 0 或 1/2）
//
// 返回: 自举后的密文，消息近似为 0 或 1（乘以 Δ0 后）
func (eval *Evaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	// 复用标准自举的 ScaleDown → ModUp → CoeffsToSlots
	// 这些步骤与论文 Algorithm 2 的 Step 1-3 对应
	btpEval := eval.BTEvaluator

	// Step 1: ScaleDown — 将密文降至 level 0，scale 调整为 Q[0]/MessageRatio
	ctScaled, _, err := btpEval.ScaleDown(ct.CopyNew())
	if err != nil {
		return nil, fmt.Errorf("binboot: ScaleDown failed: %w", err)
	}

	// Step 2: ModUp — 将模数从 q 提升到 Q（含密钥切换）
	ctModUp, err := btpEval.ModUp(ctScaled)
	if err != nil {
		return nil, fmt.Errorf("binboot: ModUp failed: %w", err)
	}

	// Step 3: CoeffsToSlots — 同态编码，将系数表示转为 slot 表示
	ctReal, ctImag, err := btpEval.CoeffsToSlots(ctModUp)
	if err != nil {
		return nil, fmt.Errorf("binboot: CoeffsToSlots failed: %w", err)
	}

	// Step 4: Eval_fBinBoot — 使用论文三角函数替换标准 EvalMod
	// 论文 §3.1: 用 f_BinBoot(x)=(1-cos(2πx))/2 替代 sine/cosine 逼近
	if ctReal, err = eval.evalBinBootPoly(ctReal); err != nil {
		return nil, fmt.Errorf("binboot: EvalBinBootPoly (real) failed: %w", err)
	}

	if ctImag != nil {
		if ctImag, err = eval.evalBinBootPoly(ctImag); err != nil {
			return nil, fmt.Errorf("binboot: EvalBinBootPoly (imag) failed: %w", err)
		}
	}

	// Step 5: SlotsToCoeffs — 同态解码，将 slot 表示转回系数表示
	ctOut, err := btpEval.SlotsToCoeffs(ctReal, ctImag)
	if err != nil {
		return nil, fmt.Errorf("binboot: SlotsToCoeffs failed: %w", err)
	}

	// 设置输出 scale 为默认 scale
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
	scalar, constant := poly.ChangeOfBasis()

	// 应用变量代换: ct' = scalar * ct + constant
	// 先做加法（常数项），再做乘法（缩放）
	ct := ctIn.CopyNew()

	// 加上常数偏移 constant
	if err := eval.BTEvaluator.Add(ct, constant, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev offset: %w", err)
	}

	// 乘以缩放因子 scalar
	// 注意: scalar 是 big.Float，需要转换为适当类型
	scalarFloat, _ := scalar.Float64()
	if err := eval.BTEvaluator.Mul(ct, scalarFloat, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev scaling: %w", err)
	}

	// Rescale 以保持 scale 在合理范围
	if err := eval.BTEvaluator.Rescale(ct, ct); err != nil {
		return nil, fmt.Errorf("cannot rescale after change of basis: %w", err)
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
//   Input: ct1, ct2 ∈ R_q^n (分别加密 b1/2 + ε1, b2/2 + ε2)
//   1. ct = ct1 + ct2 (加密 (b1+b2)/2 + ε1+ε2)
//   2. CtS ∘ ModRaise ∘ StC → 系数表示 (b1+b2+ε')/2 + I
//      注意: GateBoot 中输入为 (b1+b2)/3 + I（因为 q0=3Δ0）
//   3. Eval_fG: 同态求值门函数 f_G(x)
//   Output: ct' ∈ R_Q^n (加密 G(b1, b2) + ε'')
//
// 论文 §4.1: "f_G((b1+b2)/3 + I) = G(b1, b2)"
// 三个关键点 x = 0, 1/3, 2/3 对应 (b1,b2) = (0,0), (1,0), (1,1)
//
// 论文 §4.2 Theorem 2: GateBoot 的正确性保证
//
// 参数:
//   - ct1: 第一个输入密文（加密 b1/2 或 b1/3）
//   - ct2: 第二个输入密文（加密 b2/2 或 b2/3）
//   - gate: 门类型（GateNAND, GateAND 等）
//
// 返回: 自举后的密文，消息近似为 G(b1, b2)（0 或 1）
func (eval *Evaluator) GateBoot(ct1, ct2 *rlwe.Ciphertext, gate GateType) (*rlwe.Ciphertext, error) {

	// 检查门类型有效
	if gate < GateAND || gate > GateXNOR {
		return nil, fmt.Errorf("binboot: invalid gate type %d", gate)
	}

	// Step 0: ct = ct1 + ct2
	// 论文 §4.1 Algorithm 3: "Input: ct1, ct2; Step 1: ct = ct1 + ct2"
	ctAdd := ct1.CopyNew()
	if err := eval.BTEvaluator.Add(ctAdd, ct2, ctAdd); err != nil {
		return nil, fmt.Errorf("binboot: cannot add ciphertexts: %w", err)
	}

	// 后续步骤与 BinBoot 相同，但 EvalMod 使用门函数多项式
	btpEval := eval.BTEvaluator

	// Step 1: ScaleDown
	ctScaled, _, err := btpEval.ScaleDown(ctAdd)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot ScaleDown failed: %w", err)
	}

	// Step 2: ModUp
	ctModUp, err := btpEval.ModUp(ctScaled)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot ModUp failed: %w", err)
	}

	// Step 3: CoeffsToSlots
	ctReal, ctImag, err := btpEval.CoeffsToSlots(ctModUp)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot CoeffsToSlots failed: %w", err)
	}

	// Step 4: Eval_fG — 使用门函数三角多项式替换标准 EvalMod
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
	ctOut, err := btpEval.SlotsToCoeffs(ctReal, ctImag)
	if err != nil {
		return nil, fmt.Errorf("binboot: GateBoot SlotsToCoeffs failed: %w", err)
	}

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

	// Chebyshev 变量代换
	scalar, constant := poly.ChangeOfBasis()

	ct := ctIn.CopyNew()

	// 应用变量代换: ct' = scalar * ct + constant
	if err := eval.BTEvaluator.Add(ct, constant, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev offset: %w", err)
	}

	scalarFloat, _ := scalar.Float64()
	if err := eval.BTEvaluator.Mul(ct, scalarFloat, ct); err != nil {
		return nil, fmt.Errorf("cannot apply Chebyshev scaling: %w", err)
	}

	if err := eval.BTEvaluator.Rescale(ct, ct); err != nil {
		return nil, fmt.Errorf("cannot rescale after change of basis: %w", err)
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
