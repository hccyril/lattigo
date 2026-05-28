// Package bootstrapping implements bootstrapping for fixed-point encrypted
// approximate homomorphic encryption over the complex/real numbers (CKKS scheme).
//
// ============================ 【中文文件说明】 ============================
// 本包实现 CKKS 同态加密方案的 Bootstrapping（自举）技术 [HEAAN2019]。
//
// 【什么是 Bootstrapping？】
// Bootstrapping 是一种"刷新"密文的技术。它将低 level 的密文转换为高 level 的密文，
// 同时重置噪声水平，使得可以继续执行同态运算。
//
// 【Bootstrapping 的 5 个步骤】（对应 [HEAAN2019] Algorithm 1）：
//   1. ScaleDown: 将密文缩放到 q/|m| 范围，使消息标准化
//   2. ModUp: 将模数从 q 扩展到 Q（整个模数链），使用密钥切换
//   3. CoeffsToSlots: 同态编码（DFT），系数表示 → 槽表示
//   4. EvalMod: 同态模约简（x mod 1），使用多项式近似（sin/cos 级数）
//   5. SlotsToCoeffs: 同态解码（IDFT），槽表示 → 系数表示
//
// 【数学原理】
//   解密: m' = c0 + c1·s = m·Δ + e（e 是小噪声）
//   Bootstrapping 核心：在密文上"同态地"执行解密+模约简操作
//   1. 同态解密: Enc(m·Δ + e) → CoeffsToSlots → 槽表示
//   2. 同态模约简: Enc(m·Δ + e) mod 1 ≈ Enc(m + e/Δ)，使用多项式近似
//   3. 同态编码: Enc(m + e/Δ) → SlotsToCoeffs → 系数表示
//   4. 结果: 新的密文 Enc(m·Δ' + e')，e' 是刷新后的噪声
// =====================================================================
package bootstrapping

import "github.com/tuneinsight/lattigo/v6/core/rlwe"

type Bootstrapper interface {

	// Bootstrap defines a method that takes a single Ciphertext as input and applies
	// an in place scheme-specific bootstrapping. The result is also returned.
	// An error should notably be returned if ct.Level() < Bootstrapper.MinimumInputLevel().
	Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

	// BootstrapMany defines a method that takes a slice of Ciphertexts as input and applies an
	// in place scheme-specific bootstrapping to each Ciphertext. The result is also returned.
	// An error should notably be returned if cts[i].Level() < Bootstrapper.MinimumInputLevel().
	BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error)

	// Depth is the number of levels consumed by the bootstrapping circuit.
	// This value is equivalent to params.MaxLevel() - OutputLevel().
	Depth() int

	// MinimumInputLevel defines the minimum level that the ciphertext
	// must be at when given to the bootstrapper.
	// For the centralized bootstrapping this value is usually zero.
	// For the collective bootstrapping it is given by the user-defined
	// security parameters
	MinimumInputLevel() int

	// OutputLevel defines the level that the ciphertext will be at
	// after the bootstrapping.
	// For the centralized bootstrapping this value is the maximum
	// level minus the depth of the bootstrapping circuit.
	// For the collective bootstrapping this value is usually the
	// maximum level.
	OutputLevel() int
}
