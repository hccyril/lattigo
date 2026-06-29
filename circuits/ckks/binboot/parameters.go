// Package binboot 实现论文 "Bootstrapping Bits with CKKS" (BCKS24, 2024-767)
// 中的 BinBoot（Algorithm 2）与 GateBoot（Algorithm 3）两个核心算法。
//
// 本包在 Lattigo v6 现有 bootstrapping 框架基础上，替换 EvalMod 阶段的多项式逼近，
// 用论文提出的三角函数 f_BinBoot(x) = (1-cos(2πx))/2 和六个门的 f_G(x)（Table 4）
// 替代标准 CKKS 自举中的 sine/cosine 逼近，从而实现：
//   - BinBoot: q0 = 2Δ0，二次噪声缩减（Theorem 1）
//   - GateBoot: q0 = 3Δ0，自举同时完成二元门（Theorem 2）
package binboot

import (
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// ===================================================================
// 论文 §3.3 模数工程 & §5 Table 5: Param14 参数集
// ===================================================================
//
// 论文 Table 5 给出的 Param14 参数（低延迟）：
//   N = 2^14, (h, h~) = (256, 32), log2(QP) = 424, dnum = 13, depth = 2
//   模数链 (log2 q):
//     Base=32, StC=28, Mult=26×2, EvalMod=32×7, CtS=28×2
//   辅助模数 (log2 p): 32
//
// BinBoot 设置: q0 = 2*Δ0  （论文 §3.1: "Δ0 = q0/2"）
// GateBoot 设置: q0 = 3*Δ0 （论文 §4.1: "we may set q0 = 3*Δ0"）
//
// 由于 Lattigo 的 bootstrapping 框架要求底层 CKKS 参数与自举参数共享 Q 链，
// 这里我们构造与 Param14 模数链匹配的 ckks.ParametersLiteral 和
// bootstrapping.ParametersLiteral。

// Param14BinBootLiteral 对应论文 §5 Table 5 的 Param14，配置为 BinBoot 模式。
//
// 关键约束（论文 §3.1）:
//   - 基模 q0 = 2 * Δ0 （Base = 32, Δ0 = 31 bits，近似 2 倍关系）
//   - 环维度 N = 2^14（残差环 ConjugateInvariant）
//   - 自举环维度 N = 2^15（Standard，必须为 LogN+1，框架强制要求）
//   - 稀疏密钥 Hamming weight h~ = 32（EphemeralSecretWeight）
//   - 密集密钥 Hamming weight h = 256
//
// 【修复】LogQ 仅包含残差素数（Base + Mult 层级），自举电路素数（StC/EvalMod/CtS）
//
//	由 bootstrapping.NewParametersFromLiteral 根据因式分解深度自动追加。
//	之前误将电路素数写入 LogQ，导致模数链翻倍（LogQP≈910 而非≈460）。
//
// 【修复】Mod1Degree=32 提供 6 个 EvalMod 层级：
//   - Chebyshev 变量代换（Mul+Rescale）消耗 1 层
//   - 30 次多项式求值消耗 ceil(log2(31))=5 层
//   - 总计 6 层，与 Mod1Degree=32 的 depth=bits.Len64(32)=6 匹配
var Param14BinBootLiteral = struct {
	SchemeParams    ckks.ParametersLiteral
	BootstrapParams bootstrapping.ParametersLiteral
}{
	SchemeParams: ckks.ParametersLiteral{
		LogN: 14,
		// 残差模数链：仅包含 Base + Mult 层级
		// Base(32) + 2个Mult(45) = 3 个残差 Q 素数
		// 自举电路素数（StC/EvalMod/CtS）由框架自动追加
		LogQ: []int{
			32,     // Base: q0 ≈ 2*Δ0
			45, 45, // Mult: 非自举乘法层（供自举后运算使用）
		},
		LogP:            []int{32, 32, 32},
		Xs:              ring.Ternary{H: 256}, // h = 256（dense）
		LogDefaultScale: 31,                   // Δ0 ≈ q0/2
		RingType:        ring.ConjugateInvariant,
	},
	BootstrapParams: bootstrapping.ParametersLiteral{
		LogN:                  utils.Pointy(15), // 【修复】必须为残差 LogN+1（ConjugateInvariant 框架约束）
		EphemeralSecretWeight: utils.Pointy(32), // h~ = 32（sparse）
		LogMessageRatio:       utils.Pointy(1),  // MessageRatio=2，CosDiscrete 要求 dev>1
		Mod1Type:              0,                // CosDiscrete，将在 evaluator 中被自定义多项式覆盖
		Mod1Degree:            utils.Pointy(32), // 【修复】depth=6，为变量代换(1层)+多项式(5层)提供足够层级
		DoubleAngle:           utils.Pointy(0),  // BinBoot 不使用 double-angle
		K:                     utils.Pointy(4),
		EvalModLogScale:       utils.Pointy(32),
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{28}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{28}, {28}},
	},
}

// Param14GateBootLiteral 对应论文 §5 Table 5 的 Param14，配置为 GateBoot 模式。
//
// 关键约束（论文 §4.1）:
//   - 基模 q0 = 3 * Δ0 （Base=32, Δ0≈30 bits，近似 3 倍关系）
//   - 自举环维度 N = 2^15（Standard，必须为 LogN+1）
//   - 其余同 BinBoot
//
// 论文 §4.1: "we may set q0 = 3*Δ0, allowing for a small gain in overall modulus consumption"
// 三个等距点 0, 1/3, 2/3 对应 b1+b2 ∈ {0,1,2}。
//
// 【修复】同 Param14BinBootLiteral：
//   - LogQ 仅包含残差素数，电路素数由框架自动追加
//   - BootstrapParams.LogN = 15（LogN+1）
//   - Mod1Degree=32（提供 6 层 EvalMod）
var Param14GateBootLiteral = struct {
	SchemeParams    ckks.ParametersLiteral
	BootstrapParams bootstrapping.ParametersLiteral
}{
	SchemeParams: ckks.ParametersLiteral{
		LogN: 14,
		// 残差模数链：仅 Base + Mult 层级
		LogQ: []int{
			32,     // Base: q0 ≈ 3*Δ0
			45, 45, // Mult: 非自举乘法层
		},
		LogP:            []int{32, 32, 32},
		Xs:              ring.Ternary{H: 256},
		LogDefaultScale: 30, // Δ0 ≈ q0/3
		RingType:        ring.ConjugateInvariant,
	},
	BootstrapParams: bootstrapping.ParametersLiteral{
		LogN:                  utils.Pointy(15), // 【修复】必须为残差 LogN+1
		EphemeralSecretWeight: utils.Pointy(32),
		LogMessageRatio:       utils.Pointy(1), // MessageRatio=2，CosDiscrete 要求 dev>1
		Mod1Type:              0,
		Mod1Degree:            utils.Pointy(32), // 【修复】depth=6
		DoubleAngle:           utils.Pointy(0),
		K:                     utils.Pointy(4),
		EvalModLogScale:       utils.Pointy(32),
		SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{28}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{28}, {28}},
	},
}

// GateType 枚举论文 §4.1 Table 4 中的六个非平凡对称二元门。
type GateType int

const (
	GateAND  GateType = iota // f(x) = (1/3)(1 - 2*sin(2πx + π/6))
	GateOR                   // f(x) = (1/3)(1 - cos(2πx)) ... 见 Table 4
	GateXOR                  // f(x) = (1/3)(1 + 2*sin(2πx - π/6))
	GateNAND                 // f(x) = (2/3)(1 + sin(2πx + π/6))
	GateNOR                  // f(x) = (1/3)(1 + 2*cos(2πx))
	GateXNOR                 // f(x) = (2/3)(1 - sin(2πx - π/6))
)

// String 返回门类型的字符串表示，用于日志和测试输出。
func (g GateType) String() string {
	switch g {
	case GateAND:
		return "AND"
	case GateOR:
		return "OR"
	case GateXOR:
		return "XOR"
	case GateNAND:
		return "NAND"
	case GateNOR:
		return "NOR"
	case GateXNOR:
		return "XNOR"
	default:
		return "UNKNOWN"
	}
}
