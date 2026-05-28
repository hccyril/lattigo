// Package ckks implements a RNS-accelerated version of the Homomorphic Encryption for Arithmetic for Approximate Numbers
// (HEAAN, a.k.a. CKKS) scheme. It provides approximate arithmetic over the complex numbers.package ckks
package ckks

// ============================ 【中文文件说明】 ============================
// 本包实现了 CKKS 同态加密方案（[CKKS2017] paper）。
// CKKS 是一种支持复数近似算术的全同态加密方案。
//
// 核心算法对应关系：
//   - Algorithm 1 (Encode/Decode): 编码/解码 → encoder.go
//   - Algorithm 2 (Encrypt): 加密 → 本文件的 NewEncryptor + rlwe.Encryptor
//   - Algorithm 3 (Decrypt): 解密 → 本文件的 NewDecryptor + rlwe.Decryptor
//   - Algorithm 4 (Multiply): 同态乘法 → evaluator.go
//
// 本文件提供便捷构造函数，封装底层 RLWE 层的对应函数。
// =========================================================================

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

// NewPlaintext allocates a new [rlwe.Plaintext].
//
// inputs:
//   - params: an [rlwe.ParameterProvider] interface
//   - level: the level of the plaintext
//
// output: a newly allocated [rlwe.Plaintext] at the specified level.
//
// Note: the user can update the field MetaData to set a specific scaling factor,
// plaintext dimensions (if applicable) or encoding domain, before encoding values
// on the created plaintext.
//
// 【中文说明】创建一个新的 CKKS 明文对象。
// 明文是编码后的消息在环 R_Q = Z_Q[X]/(X^N+1) 上的多项式表示。
// 【对应论文】[CKKS2017] — 明文是 Encode 算法的输出，表示为 m(X) ∈ R_Q
func NewPlaintext(params Parameters, level int) (pt *rlwe.Plaintext) {
	// 【步骤1】调用底层 RLWE 层的 NewPlaintext，分配一个空的明文多项式
	pt = rlwe.NewPlaintext(params, level)
	// 【步骤2】标记为 SIMD 批处理模式（使用槽编码）
	pt.IsBatched = true
	// 【步骤3】设置缩放因子 Scale = Δ = 2^LogDefaultScale，用于 CKKS 的定点数表示
	pt.Scale = params.DefaultScale()
	// 【步骤4】设置 SIMD 维度元数据，决定明文能编码多少数据
	pt.LogDimensions = params.LogMaxDimensions()
	return
}

// NewCiphertext allocates a new [rlwe.Ciphertext].
//
// inputs:
//   - params: an [rlwe.ParameterProvider] interface
//   - degree: the degree of the ciphertext
//   - level: the level of the Ciphertext
//
// output: a newly allocated [rlwe.Ciphertext] of the specified degree and level.
//
// 【中文说明】创建一个新的 CKKS 密文对象。
// 密文是 RLWE 密文，由 degree+1 个环多项式组成：
//   - degree=1: 标准密文 (c0, c1)
//   - degree=2: 乘法后密文 (c0, c1, c2)，需要重线性化
// 【对应论文】[CKKS2017] Algorithm 2 — 密文 ct = (c0, c1) ∈ R_Q²
func NewCiphertext(params Parameters, degree, level int) (ct *rlwe.Ciphertext) {
	// 【步骤1】调用底层 RLWE 层的 NewCiphertext，分配空的密文多项式数组
	ct = rlwe.NewCiphertext(params, degree, level)
	// 【步骤2】标记为 SIMD 批处理模式
	ct.IsBatched = true
	// 【步骤3】设置缩放因子 Scale = Δ = 2^LogDefaultScale
	ct.Scale = params.DefaultScale()
	// 【步骤4】设置 SIMD 维度元数据
	ct.LogDimensions = params.LogMaxDimensions()
	return
}

// NewEncryptor instantiates a new [rlwe.Encryptor].
//
// inputs:
//   - params: an [rlwe.ParameterProvider] interface
//   - key: *[rlwe.SecretKey] or *[rlwe.PublicKey]
//
// output: an [rlwe.Encryptor] instantiated with the provided key.
//
// 【中文说明】创建一个新的加密器，用于将明文加密为密文。
// key 可以是公钥（任何人都能加密）或私钥（仅密钥持有者能加密）。
// 【对应论文】[CKKS2017] Algorithm 2 (Encrypt)
// 加密公式: ct = Enc(pk, m) = (c0, c1)，其中 e0,e1 为高斯噪声多项式
func NewEncryptor(params Parameters, key rlwe.EncryptionKey) *rlwe.Encryptor {
	return rlwe.NewEncryptor(params, key)
}

// NewDecryptor instantiates a new [rlwe.Decryptor].
//
// inputs:
//   - params: an [rlwe.ParameterProvider] interface
//   - key: *[rlwe.SecretKey]
//
// output: an [rlwe.Decryptor] instantiated with the provided key.
//
// 【中文说明】创建一个新的解密器，用于将密文解密为明文。
// 解密必须使用私钥。
// 【对应论文】[CKKS2017] Algorithm 3 (Decrypt)
// 解密公式: m' = c0 + c1·s (mod Q)，其中 s 是私钥
func NewDecryptor(params Parameters, key *rlwe.SecretKey) *rlwe.Decryptor {
	return rlwe.NewDecryptor(params, key)
}

// NewKeyGenerator instantiates a new [rlwe.KeyGenerator].
//
// inputs:
//   - params: an [rlwe.ParameterProvider] interface
//
// output: an [rlwe.KeyGenerator].
//
// 【中文说明】创建一个新的密钥生成器。
// 可生成的密钥：私钥、公钥、重线性化密钥、Galois 密钥、评估密钥。
// 【对应论文】[CKKS2017] — 密钥生成基于 RLWE 问题的困难性
// 私钥 s 从低权重三元分布采样: s ∈ {-1, 0, 1}^N
func NewKeyGenerator(params Parameters) *rlwe.KeyGenerator {
	return rlwe.NewKeyGenerator(params)
}
