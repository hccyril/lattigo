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
//   - 环维度 N = 2^14
//   - 稀疏密钥 Hamming weight h~ = 32（EphemeralSecretWeight）
//   - 密集密钥 Hamming weight h = 256
//   - dnum = 13（密钥切换 gadget rank）
//
// 注意：为适配 Lattigo 框架的素数生成约束，这里取 LogDefaultScale=31，
// Base prime=32，满足 q0 ≈ 2*Δ0 的比例关系（论文允许近似）。
var Param14BinBootLiteral = struct {
	SchemeParams    ckks.ParametersLiteral
	BootstrapParams bootstrapping.ParametersLiteral
}{
	SchemeParams: ckks.ParametersLiteral{
		LogN: 14,
		// 模数链：Base(32) + Mult(26×2) + StC(28) + CtS(28×2) + EvalMod(32×7)
		// 总计 13 个 Q 素数，对应 dnum=13
		LogQ: []int{
			32, // Base: q0 ≈ 2*Δ0
			26, 26, // Mult: 非自举乘法层
			28, // StC: Slots-to-Coeffs
			32, 32, 32, 32, 32, 32, 32, // EvalMod: 7 层
			28, 28, // CtS: Coeffs-to-Slots
		},
		LogP:            []int{32, 32, 32},
		Xs:              ring.Ternary{H: 256}, // h = 256（dense）
		LogDefaultScale: 31,                   // Δ0 ≈ q0/2
		RingType:        ring.ConjugateInvariant,
	},
	BootstrapParams: bootstrapping.ParametersLiteral{
		LogN:                        utils.Pointy(14),
		EphemeralSecretWeight:       utils.Pointy(32), // h~ = 32（sparse）
		LogMessageRatio:             utils.Pointy(1),  // q0/Δ0 = 2 → log2(2)=1
		Mod1Type:                    0,                // CosDiscrete，将在 evaluator 中被自定义多项式覆盖
		Mod1Degree:                  utils.Pointy(30),
		DoubleAngle:                 utils.Pointy(0), // BinBoot 不使用 double-angle
		K:                           utils.Pointy(4),
		EvalModLogScale:             utils.Pointy(32),
		SlotsToCoeffsFactorizationDepthAndLogScales:  [][]int{{28}},
		CoeffsToSlotsFactorizationDepthAndLogScales: [][]int{{28}, {28}},
	},
}

// Param14GateBootLiteral 对应论文 §5 Table 5 的 Param14，配置为 GateBoot 模式。
//
// 关键约束（论文 §4.1）:
//   - 基模 q0 = 3 * Δ0 （Base=32, Δ0≈30 bits，近似 3 倍关系）
//   - 其余同 BinBoot
//
// 论文 §4.1: "we may set q0 = 3*Δ0, allowing for a small gain in overall modulus consumption"
// 三个等距点 0, 1/3, 2/3 对应 b1+b2 ∈ {0,1,2}。
var Param14GateBootLiteral = struct {
	SchemeParams    ckks.ParametersLiteral
	BootstrapParams bootstrapping.ParametersLiteral
}{
	SchemeParams: ckks.ParametersLiteral{
		LogN: 14,
		LogQ: []int{
			32, // Base: q0 ≈ 3*Δ0
			26, 26, // Mult
			28, // StC
			32, 32, 32, 32, 32, 32, 32, // EvalMod
			28, 28, // CtS
		},
		LogP:            []int{32, 32, 32},
		Xs:              ring.Ternary{H: 256},
		LogDefaultScale: 30, // Δ0 ≈ q0/3
		RingType:        ring.ConjugateInvariant,
	},
	BootstrapParams: bootstrapping.ParametersLiteral{
		LogN:                        utils.Pointy(14),
		EphemeralSecretWeight:       utils.Pointy(32),
		LogMessageRatio:             utils.Pointy(1), // q0/Δ0 ≈ 3，取 log2 后近似
		Mod1Type:                    0,
		Mod1Degree:                  utils.Pointy(30),
		DoubleAngle:                 utils.Pointy(0),
		K:                           utils.Pointy(4),
		EvalModLogScale:             utils.Pointy(32),
		SlotsToCoeffsFactorizationDepthAndLogScales:  [][]int{{28}},
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
