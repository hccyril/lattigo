// Package binboot 的多项式逼近生成器。
//
// 本文件实现论文 §3.1 和 §4.1 Table 4 中定义的三角函数的 Chebyshev 多项式逼近，
// 用于在 EvalMod 阶段替换标准 CKKS 自举的 sine/cosine 逼近。
//
// 论文对应关系：
//   - §3.1 公式: f_BinBoot(x) = (1 - cos(2πx)) / 2
//     性质: f(0+Z)=0, f(1/2+Z)=1，导数在 0、1/2 处为零 → 二次噪声缩减（Lemma 1）
//   - §4.1 Table 4: 六个对称二元门的三角函数 f_G(x)
//     性质: f_G((b1+b2)/3 + I) = G(b1, b2)，周期为 1
package binboot

import (
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// ===================================================================
// 论文 §3.1: BinBoot 三角函数
// f_BinBoot(x) = (1 - cos(2πx)) / 2
// ===================================================================

// fBinBoot 实现论文 §3.1 的核心三角函数 f_BinBoot(x) = (1-cos(2πx))/2。
//
// 论文原文（§3.1）:
//   "we choose the following: f_BinBoot(x) = (1/2)(1 - cos(2πx))"
//
// 性质（Lemma 1, §3.2）:
//   - f_BinBoot(b/2 + I + ε/2) - b = O(ε²)  （二次噪声缩减）
//   - f(0+Z) = 0, f(1/2+Z) = 1
//   - 导数在 x ∈ (1/2)Z 处为零，因此不放大噪声
func fBinBoot(x float64) float64 {
	return (1.0 - math.Cos(2.0*math.Pi*x)) / 2.0
}

// NewBinBootPoly 生成 f_BinBoot 的 Chebyshev 多项式逼近。
//
// 论文 §3.1: "Eval_fBinBoot is the homomorphic evaluation of f_BinBoot(x)=(1-cos(2x))/2
// via appropriate polynomial approximation."
//
// 参数:
//   - degree: Chebyshev 多项式次数（论文实验中使用 30）
//   - K: 逼近区间 [-K, K]，覆盖 ModRaise 引入的整数项 I 的范围
//
// 返回: bignum.Polynomial，可直接用于 circuits/ckks/polynomial.Evaluator.Evaluate
//
// 注意: 论文中 BinBoot 的输入经过 CtS 后为 (φ+ε+ε₂)/2 + I，
// 其中 I 是小整数，因此 K 取较小值（如 4）即可。
func NewBinBootPoly(degree int, K float64) bignum.Polynomial {
	// 使用 bignum.ChebyshevApproximation 在 [-K, K] 上逼近 f_BinBoot
	// 论文 §3.1 Figure 2 展示了该函数在 [-1, 1] 上的图像
	interval := bignum.Interval{
		A:     *bignum.NewFloat(-K, 128),
		B:     *bignum.NewFloat(K, 128),
		Nodes: degree,
	}

	poly := bignum.ChebyshevApproximation(func(x *bignum.Complex) *bignum.Complex {
		xf64, _ := x[0].Float64()
		return &bignum.Complex{
			newFloat(fBinBoot(xf64)),
			newFloat(0),
		}
	}, interval)

	// f_BinBoot 是偶函数（cos 是偶函数），标记 IsOdd=false 以优化评估
	poly.IsOdd = false
	poly.IsEven = true

	// 清零奇数次系数（偶函数的奇数次 Chebyshev 系数应为 0）
	for i := range poly.Coeffs {
		if i&1 == 1 {
			poly.Coeffs[i] = nil
		}
	}

	return poly
}

// ===================================================================
// 论文 §4.1 Table 4: GateBoot 六个门的三角函数
// ===================================================================

// fGateNAND 实现论文 Table 4 的 NAND 门函数:
//   f_NAND(x) = (2/3)(1 + sin(2πx + π/6))
//
// 验证（论文 §4.1）:
//   f_NAND(0) = (2/3)(1 + sin(π/6)) = (2/3)(1 + 0.5) = 1 = NAND(0,0)
//   f_NAND(1/3) = (2/3)(1 + sin(2π/3 + π/6)) = (2/3)(1 + sin(5π/6))
//               = (2/3)(1 + 0.5) = 1 = NAND(1,0)
//   f_NAND(2/3) = (2/3)(1 + sin(4π/3 + π/6)) = (2/3)(1 + sin(3π/2))
//               = (2/3)(1 - 1) = 0 = NAND(1,1)
func fGateNAND(x float64) float64 {
	return (2.0 / 3.0) * (1.0 + math.Sin(2.0*math.Pi*x+math.Pi/6.0))
}

// fGateAND 实现论文 Table 4 的 AND 门函数:
//   f_AND(x) = (1/3)(1 - 2*sin(2πx + π/6))
//
// 验证: AND(0,0)=1, AND(1,0)=0, AND(1,1)=1
//   f_AND(0) = (1/3)(1 - 2*0.5) = 0 ... 注意 AND = 1 - NAND
//   实际: f_AND(x) = 1 - f_NAND(x) = 1 - (2/3)(1+sin(...)) = (1/3) - (2/3)sin(...)
//   = (1/3)(1 - 2*sin(2πx + π/6))
func fGateAND(x float64) float64 {
	return (1.0 / 3.0) * (1.0 - 2.0*math.Sin(2.0*math.Pi*x+math.Pi/6.0))
}

// fGateOR 实现论文 Table 4 的 OR 门函数:
//   f_OR(x) = (1/3)(1 - cos(2πx)) + ...
//   论文 Table 4: f_OR(x) = (1/3)(1 - cos(2πx))  [需校验]
//
// 验证: OR(0,0)=0, OR(1,0)=1, OR(1,1)=1
//   f_OR(0) = (1/3)(1 - 1) = 0 ✓
//   f_OR(1/3) = (1/3)(1 - cos(2π/3)) = (1/3)(1 - (-0.5)) = 0.5 ...
//   需要参考论文原文：OR = (1/3)(1 - cos(2πx)) 不完全正确
//   实际论文 Table 4: f_OR(x) = (1/3)(1 - cos(2πx))，但需验证 2/3 处
//   f_OR(2/3) = (1/3)(1 - cos(4π/3)) = (1/3)(1 - (-0.5)) = 0.5
//   这不等于 1，因此论文 Table 4 的 OR 公式应为 (1/3)(2 - cos(2πx)) 或类似形式
//   根据论文 Figure 3 和 Table 4 的对称性，OR = 1 - NOR
func fGateOR(x float64) float64 {
	// OR = 1 - NOR, NOR = (1/3)(1 + 2cos(2πx))
	// OR = 1 - (1/3)(1 + 2cos(2πx)) = (2/3) - (2/3)cos(2πx) = (2/3)(1 - cos(2πx))
	return (2.0 / 3.0) * (1.0 - math.Cos(2.0*math.Pi*x))
}

// fGateXOR 实现论文 Table 4 的 XOR 门函数:
//   f_XOR(x) = (1/3)(1 + 2*sin(2πx - π/6))
//
// 验证: XOR(0,0)=0, XOR(1,0)=1, XOR(1,1)=0
//   f_XOR(0) = (1/3)(1 + 2*sin(-π/6)) = (1/3)(1 - 1) = 0 ✓
//   f_XOR(1/3) = (1/3)(1 + 2*sin(2π/3 - π/6)) = (1/3)(1 + 2*sin(π/2)) = (1/3)(3) = 1 ✓
//   f_XOR(2/3) = (1/3)(1 + 2*sin(4π/3 - π/6)) = (1/3)(1 + 2*sin(7π/6))
//             = (1/3)(1 + 2*(-0.5)) = 0 ✓
func fGateXOR(x float64) float64 {
	return (1.0 / 3.0) * (1.0 + 2.0*math.Sin(2.0*math.Pi*x-math.Pi/6.0))
}

// fGateNOR 实现论文 Table 4 的 NOR 门函数:
//   f_NOR(x) = (1/3)(1 + 2*cos(2πx))
//
// 验证: NOR(0,0)=1, NOR(1,0)=0, NOR(1,1)=0
//   f_NOR(0) = (1/3)(1 + 2) = 1 ✓
//   f_NOR(1/3) = (1/3)(1 + 2*cos(2π/3)) = (1/3)(1 - 1) = 0 ✓
//   f_NOR(2/3) = (1/3)(1 + 2*cos(4π/3)) = (1/3)(1 - 1) = 0 ✓
func fGateNOR(x float64) float64 {
	return (1.0 / 3.0) * (1.0 + 2.0*math.Cos(2.0*math.Pi*x))
}

// fGateXNOR 实现论文 Table 4 的 XNOR 门函数:
//   f_XNOR(x) = (2/3)(1 - sin(2πx - π/6))
//
// 验证: XNOR = 1 - XOR
//   f_XNOR(x) = 1 - (1/3)(1 + 2*sin(2πx - π/6))
//             = (2/3) - (2/3)sin(2πx - π/6) = (2/3)(1 - sin(2πx - π/6))
func fGateXNOR(x float64) float64 {
	return (2.0 / 3.0) * (1.0 - math.Sin(2.0*math.Pi*x-math.Pi/6.0))
}

// gateFunction 根据门类型返回对应的三角函数实现。
// 论文 §4.1 Table 4。
func gateFunction(g GateType) func(float64) float64 {
	switch g {
	case GateAND:
		return fGateAND
	case GateOR:
		return fGateOR
	case GateXOR:
		return fGateXOR
	case GateNAND:
		return fGateNAND
	case GateNOR:
		return fGateNOR
	case GateXNOR:
		return fGateXNOR
	default:
		panic("binboot: unknown gate type")
	}
}

// NewGatePoly 生成指定门的三角函数 f_G 的 Chebyshev 多项式逼近。
//
// 论文 §4.1: "Step 4 consists in homomorphically evaluating a trigonometric function f_G
// that removes I and sends φ1+φ2 to G(φ1, φ2)."
//
// 参数:
//   - gate: 门类型（GateNAND, GateAND 等）
//   - degree: Chebyshev 多项式次数
//   - K: 逼近区间 [-K, K]
//
// 返回: bignum.Polynomial，用于 EvalMod 阶段
//
// 注意: GateBoot 的输入为 (φ1+φ2+ε)/3 + I，因此逼近区间需覆盖 I 的范围。
// 三个关键点 x=0, 1/3, 2/3 对应 b1+b2 = 0, 1, 2。
func NewGatePoly(gate GateType, degree int, K float64) bignum.Polynomial {
	f := gateFunction(gate)

	interval := bignum.Interval{
		A:     *bignum.NewFloat(-K, 128),
		B:     *bignum.NewFloat(K, 128),
		Nodes: degree,
	}

	poly := bignum.ChebyshevApproximation(func(x *bignum.Complex) *bignum.Complex {
		xf64, _ := x[0].Float64()
		return &bignum.Complex{
			newFloat(f(xf64)),
			newFloat(0),
		}
	}, interval)

	return poly
}

// ===================================================================
// 辅助函数
// ===================================================================

// newFloat 创建一个指定值和精度的 big.Float（精度 128 位）。
func newFloat(v float64) *big.Float {
	return bignum.NewFloat(v, 128)
}

// ValidateGateFunction 在明文上验证门函数的正确性。
// 用于单元测试，确认 f_G(0), f_G(1/3), f_G(2/3) 与 G 的真值表一致。
//
// 论文 §4.1: "f_G((x1+x2)/3) = G(x1, x2) for all x1, x2 ∈ {0,1,2}"
func ValidateGateFunction(g GateType) bool {
	f := gateFunction(g)
	tol := 1e-10

	// 真值表: (b1, b2) -> b1+b2 -> expected G(b1,b2)
	cases := []struct {
		x, expected float64
	}{
		{0.0 / 3.0, gateTruthTable(g, 0, 0)},
		{1.0 / 3.0, gateTruthTable(g, 1, 0)},
		{2.0 / 3.0, gateTruthTable(g, 1, 1)},
	}

	for _, c := range cases {
		if math.Abs(f(c.x)-c.expected) > tol {
			return false
		}
	}
	return true
}

// gateTruthTable 返回门 G 的真值表值。
func gateTruthTable(g GateType, b1, b2 int) float64 {
	switch g {
	case GateAND:
		return float64(b1 * b2)
	case GateOR:
		return float64(b1 | b2)
	case GateXOR:
		return float64(b1 ^ b2)
	case GateNAND:
		return float64(1 - b1*b2)
	case GateNOR:
		return float64(1 - (b1 | b2))
	case GateXNOR:
		return float64(1 - (b1 ^ b2))
	default:
		return -1
	}
}

