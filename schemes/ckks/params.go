package ckks

// ============================ 【中文文件说明】 ============================
// 本文件定义 CKKS 方案的参数结构体及其相关方法。
// 【对应论文】[CKKS2017] Section 5 (Parameter Selection)
//
// CKKS 参数决定了：
//   - 环的度 N = 2^LogN（影响安全性和槽数量）
//   - 模数链 Q = q₀ × q₁ × ... × qₗ（影响最大乘法深度）
//   - 缩放因子 Δ = 2^LogDefaultScale（影响数值精度）
//   - 精度模式 PREC64/PREC128（决定每次 Rescale 消耗的 level 数）
// =========================================================================

import (
	"encoding/json"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// PrecisionMode is a variable that defines how many primes (one
// per machine word) are required to store initial message values.
// This also sets how many primes are consumed per rescaling.
//
// There are currently two modes supported:
//   - PREC64 (one 64-bit word)
//   - PREC128 (two 64-bit words)
//
// PREC64 is the default mode and supports reference plaintext scaling
// factors of up to 2^{64}, while PREC128 scaling factors of up to 2^{128}.
//
// The PrecisionMode is chosen automatically based on the provided initial
// `LogDefaultScale` value provided by the user.
//
// 【中文说明】精度模式，定义了 CKKS 方案中明文缩放因子的精度位数。
// 该模式决定了存储初始消息值需要多少机器字，以及每次 Rescale 消耗多少个素数。
//   - PREC64: 64 位精度，每次 Rescale 消耗 1 个 level（素数）
//   - PREC128: 128 位精度，每次 Rescale 消耗 2 个 level（素数）
// 【对应论文】[CKKS2017] Section 3.3 — 缩放因子 Δ 控制数值精度
type PrecisionMode int

const (
	NTTFlag = true
	// PREC64: 64 位精度模式，默认模式，缩放因子最大支持 2^64
	PREC64  = PrecisionMode(0)
	// PREC128: 128 位精度模式，缩放因子最大支持 2^128，每次 Rescale 消耗 2 个 level
	PREC128 = PrecisionMode(1)
)

// ParametersLiteral is a literal representation of CKKS parameters.  It has public
// fields and is used to express unchecked user-defined parameters literally into
// Go programs. The [NewParametersFromLiteral] function is used to generate the actual
// checked parameters from the literal representation.
//
// Users must set the polynomial degree (in log_2, LogN) and the coefficient modulus, by either setting
// the Q and P fields to the desired moduli chain, or by setting the LogQ and LogP fields to
// the desired moduli sizes (in log_2). Users must also specify a default initial scale for the plaintexts.
//
// Optionally, users may specify the error variance (Sigma), the secrets' density (H), the ring
// type (RingType) and the number of slots (in log_2, LogSlots). If left unset, standard default values for
// these field are substituted at parameter creation (see NewParametersFromLiteral).
//
// 【中文说明】CKKS 参数的字面量表示，用户可以直接填写原始数值来定义参数。
// 填写后必须通过 NewParametersFromLiteral 函数进行验证，生成正式的 Parameters 对象。
//
// 必填字段:
//   - LogN: 对数多项式次数，N = 2^LogN，决定环的度和最大槽数量
//   - Q: 密文模数链（素数列表），每个素数增加一个乘法深度（或通过 LogQ 指定比特长度）
//   - P: 辅助模数链（用于密钥交换优化），不参与同态容量计算（或通过 LogP 指定比特长度）
//   - LogDefaultScale: 默认缩放因子的对数值，Δ = 2^LogDefaultScale
type ParametersLiteral struct {
	LogN            int
	LogNthRoot      int
	Q               []uint64
	P               []uint64
	LogQ            []int `json:",omitempty"`
	LogP            []int `json:",omitempty"`
	Xe              ring.DistributionParameters
	Xs              ring.DistributionParameters
	RingType        ring.Type
	LogDefaultScale int
	// 各字段中文说明见上方【中文说明】注释块
}

// GetRLWEParametersLiteral returns the [rlwe.ParametersLiteral] from the target [ckks.ParameterLiteral].
func (p ParametersLiteral) GetRLWEParametersLiteral() rlwe.ParametersLiteral {
	return rlwe.ParametersLiteral{
		LogN:         p.LogN,
		LogNthRoot:   p.LogNthRoot,
		Q:            p.Q,
		P:            p.P,
		LogQ:         p.LogQ,
		LogP:         p.LogP,
		Xe:           p.Xe,
		Xs:           p.Xs,
		RingType:     p.RingType,
		NTTFlag:      NTTFlag,
		DefaultScale: rlwe.NewScale(math.Exp2(float64(p.LogDefaultScale))),
	}
}

// Parameters represents a parameter set for the CKKS cryptosystem. Its fields are private and
// immutable. See [ParametersLiteral] for user-specified parameters.
//
// 【中文说明】经过验证的 CKKS 参数集。字段私有且不可变，确保参数一致性。
// 封装了底层 RLWE 参数，并添加了 CKKS 特有的缩放因子信息。
type Parameters struct {
	rlwe.Parameters
}

// NewParametersFromLiteral instantiate a set of CKKS parameters from a [ParametersLiteral] specification.
// It returns the empty parameters [Parameters]{} and a non-nil error if the specified parameters are invalid.
//
// If the LogSlots field is left unset, its value is set to LogN-1 for the Standard ring and to LogN for
// the conjugate-invariant ring.
//
// See [rlwe.NewParametersFromLiteral] for default values of the other optional fields.
//
// 【中文说明】从参数字面量创建经过验证的 CKKS 参数集。
// 该函数会验证所有参数的合法性，包括素数是否满足 NTT 条件、缩放因子是否有效等。
// 【对应论文】[CKKS2017] Section 5.1 — 参数选择原则
func NewParametersFromLiteral(pl ParametersLiteral) (Parameters, error) {
	rlweParams, err := rlwe.NewParametersFromLiteral(pl.GetRLWEParametersLiteral())
	if err != nil {
		return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: %w", err)
	}

	if pl.LogDefaultScale > 128 {
		return Parameters{}, fmt.Errorf("cannot NewParametersFromLiteral: LogDefaultScale=%d > 128 or < 0", pl.LogDefaultScale)
	}

	return Parameters{rlweParams}, nil
}

// StandardParameters returns the CKKS parameters corresponding to the receiver
// parameter set. If the receiver is already a standard parameter set
// (i.e., RingType==Standard), then the method returns the receiver.
func (p Parameters) StandardParameters() (pckks Parameters, err error) {
	if p.RingType() == ring.Standard {
		return p, nil
	}
	pckks = p
	pckks.Parameters, err = pckks.Parameters.StandardParameters()
	return
}

// ParametersLiteral returns the [ParametersLiteral] of the target [Parameters].
func (p Parameters) ParametersLiteral() (pLit ParametersLiteral) {
	return ParametersLiteral{
		LogN:            p.LogN(),
		LogNthRoot:      p.LogNthRoot(),
		Q:               p.Q(),
		P:               p.P(),
		Xe:              p.Xe(),
		Xs:              p.Xs(),
		RingType:        p.RingType(),
		LogDefaultScale: p.LogDefaultScale(),
	}
}

// GetRLWEParameters returns a pointer to the underlying RLWE parameters.
func (p Parameters) GetRLWEParameters() *rlwe.Parameters {
	return &p.Parameters
}

// MaxLevel returns the maximum ciphertext level
//
// 【中文说明】返回密文的最大层级。level = |Q| - 1，即模数链中素数的个数减 1。
// level 越高，可用的模数越多，能执行的乘法次数越多。每次 Rescale 消耗 level（1 或 2）。
func (p Parameters) MaxLevel() int {
	return p.QCount() - 1
}

// MaxDimensions returns the maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
//
// 【中文说明】返回单个明文多项式可以 SIMD 打包的最大矩阵维度。
// CKKS 使用 SIMD 技术将多个复数打包到一个明文多项式中：
//   - Standard 环: 最多 N/2 个槽（slot），即 1×(N/2) 的矩阵
//   - ConjugateInvariant 环: 最多 N 个槽，即 1×N 的矩阵
func (p Parameters) MaxDimensions() ring.Dimensions {
	switch p.RingType() {
	case ring.Standard:
		return ring.Dimensions{Rows: 1, Cols: p.N() >> 1}
	case ring.ConjugateInvariant:
		return ring.Dimensions{Rows: 1, Cols: p.N()}
	default:
		// Sanity check
		panic("cannot MaxDimensions: invalid ring type")
	}
}

// LogMaxDimensions returns the log2 of maximum dimension of the matrix that can be SIMD packed in a single plaintext polynomial.
func (p Parameters) LogMaxDimensions() ring.Dimensions {
	switch p.RingType() {
	case ring.Standard:
		return ring.Dimensions{Rows: 0, Cols: p.LogN() - 1}
	case ring.ConjugateInvariant:
		return ring.Dimensions{Rows: 0, Cols: p.LogN()}
	default:
		// Sanity check
		panic("cannot LogMaxDimensions: invalid ring type")
	}
}

// MaxSlots returns the total number of entries (slots) that a plaintext can store.
// This value is obtained by multiplying all dimensions from MaxDimensions.
func (p Parameters) MaxSlots() int {
	dims := p.MaxDimensions()
	return dims.Rows * dims.Cols
}

// LogMaxSlots returns the total number of entries (slots) that a plaintext can store.
// This value is obtained by summing all log dimensions from LogDimensions.
func (p Parameters) LogMaxSlots() int {
	dims := p.LogMaxDimensions()
	return dims.Rows + dims.Cols
}

// LogDefaultScale returns the log2 of the default plaintext
// scaling factor (rounded to the nearest integer).
//
// 【中文说明】返回默认缩放因子的对数值（以 2 为底，四舍五入到整数）。
// 缩放因子 Δ = 2^LogDefaultScale，用于 CKKS 的定点数表示。
// 编码时: m_scaled = round(m × Δ)，解密后: m = m_scaled / Δ
func (p Parameters) LogDefaultScale() int {
	return int(math.Round(math.Log2(p.DefaultScale().Float64())))
}

// EncodingPrecision returns the encoding precision in bits of the plaintext values which
// is max(53, log2(DefaultScale)).
func (p Parameters) EncodingPrecision() (prec uint) {
	if log2scale := math.Log2(p.DefaultScale().Float64()); log2scale <= 53 {
		prec = 53
	} else {
		prec = uint(log2scale)
	}

	return
}

// PrecisionMode returns the precision mode of the parameters.
// This value can be [ckks.PREC64] or [ckks.PREC128].
//
// 【中文说明】返回精度模式。根据 LogDefaultScale 自动选择：
//   - LogDefaultScale ≤ 64 → PREC64（每次 Rescale 消耗 1 个 level）
//   - LogDefaultScale > 64 → PREC128（每次 Rescale 消耗 2 个 level）
func (p Parameters) PrecisionMode() PrecisionMode {
	if p.LogDefaultScale() <= 64 {
		return PREC64
	}
	return PREC128
}

// LevelsConsumedPerRescaling returns the number of levels (i.e. primes)
// consumed per rescaling. This value is 1 if the precision mode is PREC64
// and is 2 if the precision mode is PREC128.
//
// 【中文说明】返回每次重缩放（Rescale）操作消耗的 level 数量。
//   - PREC64 模式: 消耗 1 个 level（1 个素数）
//   - PREC128 模式: 消耗 2 个 level（2 个素数）
// 这决定了模数链的有效乘法深度。
func (p Parameters) LevelsConsumedPerRescaling() int {
	switch p.PrecisionMode() {
	case PREC128:
		return 2
	default:
		return 1
	}
}

// GetOptimalScalingFactor returns a scaling factor b such that Rescale(a * b) = c
func (p Parameters) GetOptimalScalingFactor(a, c rlwe.Scale, level int) (b rlwe.Scale) {
	b = rlwe.NewScale(1)
	Q := p.Q()
	for i := 0; i < p.LevelsConsumedPerRescaling(); i++ {
		b = b.Mul(rlwe.NewScale(Q[level-i]))
	}
	return
}

// MaxDepth returns the maximum depth enabled by the parameters,
// which is obtained as p.MaxLevel() / p.LevelsConsumedPerRescaling().
//
// 【中文说明】返回最大乘法深度。
// 计算公式: MaxDepth = MaxLevel() / LevelsConsumedPerRescaling()
// 即模数链中可用于乘法操作的总次数（每次乘法后需要一次 Rescale）。
// 例如: 若有 2 个素数，PREC64 模式，则 MaxDepth = 1（只能做 1 次乘法+Rescale）
func (p Parameters) MaxDepth() int {
	return p.MaxLevel() / p.LevelsConsumedPerRescaling()
}

// LogQLvl returns the size of the modulus Q in bits at a specific level
func (p Parameters) LogQLvl(level int) int {
	tmp := p.QLvl(level)
	return tmp.BitLen()
}

// QLvl returns the product of the moduli at the given level as a [big.Int]
func (p Parameters) QLvl(level int) *big.Int {
	tmp := bignum.NewInt(1)
	for _, qi := range p.Q()[:level+1] {
		tmp.Mul(tmp, bignum.NewInt(qi))
	}
	return tmp
}

// GaloisElementForRotation returns the Galois element for generating the
// automorphism phi(k): X -> X^{5^k mod 2N} mod (X^{N} + 1), which acts as a
// cyclic rotation by k position to the left on batched plaintexts.
//
// Example:
// Recall that batched plaintexts are 2xN/2 matrices of the form [m, conjugate(m)]
// (the conjugate is implicitly ignored) thus given the following plaintext matrix:
//
// [a, b, c, d][conj(a), conj(b), conj(c), conj(d)]
//
// a rotation by k=3 will change the plaintext to:
//
// [d, a, b, c][conj(d), conj(a), conj(b), conj(c)]
//
// Providing a negative k will change direction of the cyclic rotation to the right.
//
// Note that when using the ConjugateInvariant variant of the scheme, the conjugate is
// dropped and the matrix becomes an 1xN matrix.
func (p Parameters) GaloisElementForRotation(k int) uint64 {
	return p.Parameters.GaloisElement(k)
}

// GaloisElementForComplexConjugation returns the Galois element for generating the
// automorphism X -> X^{-1 mod NthRoot} mod (X^{N} + 1). This automorphism
// acts as a swapping the rows of the plaintext algebra when the plaintext
// is batched.
//
// Example:
// Recall that batched plaintexts are 2xN/2 matrices of the form [m, conjugate(m)]
// (the conjugate is implicitly ignored) thus given the following plaintext matrix:
//
// [a, b, c, d][conj(a), conj(b), conj(c), conj(d)]
//
// the complex conjugation will return the following plaintext matrix:
//
// [conj(a), conj(b), conj(c), conj(d)][a, b, c, d]
//
// Note that when using the ConjugateInvariant variant of the scheme, the conjugate is
// dropped and this operation is not defined.
func (p Parameters) GaloisElementForComplexConjugation() uint64 {
	return p.Parameters.GaloisElementOrderTwoOrthogonalSubgroup()
}

// GaloisElementsForInnerSum returns the list of Galois elements necessary to apply the method
// `InnerSum` operation with parameters batch and n.
func (p Parameters) GaloisElementsForInnerSum(batch, n int) []uint64 {
	return rlwe.GaloisElementsForInnerSum(p, batch, n)
}

// GaloisElementsForReplicate returns the list of Galois elements necessary to perform the
// `Replicate` operation with parameters batch and n.
func (p Parameters) GaloisElementsForReplicate(batch, n int) []uint64 {
	return rlwe.GaloisElementsForReplicate(p, batch, n)
}

// GaloisElementsForTrace returns the list of Galois elements required for the for the Trace operation.
// Trace maps X -> sum((-1)^i * X^{i*n+1}) for 2^{LogN} <= i < N.
func (p Parameters) GaloisElementsForTrace(logN int) []uint64 {
	return rlwe.GaloisElementsForTrace(p, logN)
}

// Equal compares two sets of parameters for equality.
func (p Parameters) Equal(other *Parameters) bool {
	return p.Parameters.Equal(&other.Parameters)
}

// MarshalBinary returns a []byte representation of the parameter set.
// This representation corresponds to the one returned by MarshalJSON.
func (p Parameters) MarshalBinary() ([]byte, error) {
	return p.MarshalJSON()
}

// UnmarshalBinary decodes a []byte into a parameter set struct
func (p *Parameters) UnmarshalBinary(data []byte) (err error) {
	return p.UnmarshalJSON(data)
}

// MarshalJSON returns a JSON representation of this parameter set. See Marshal from the [encoding/json] package.
func (p Parameters) MarshalJSON() ([]byte, error) {
	return json.Marshal(p.ParametersLiteral())
}

// UnmarshalJSON reads a JSON representation of a parameter set into the receiver Parameter. See Unmarshal from the [encoding/json] package.
func (p *Parameters) UnmarshalJSON(data []byte) (err error) {
	var params ParametersLiteral
	if err = json.Unmarshal(data, &params); err != nil {
		return
	}
	*p, err = NewParametersFromLiteral(params)
	return
}

func (p *ParametersLiteral) UnmarshalJSON(b []byte) (err error) {
	var pl struct {
		LogN            int
		LogNthRoot      int
		Q               []uint64
		P               []uint64
		LogQ            []int
		LogP            []int
		Pow2Base        int
		Xe              map[string]interface{}
		Xs              map[string]interface{}
		RingType        ring.Type
		LogDefaultScale int
	}

	err = json.Unmarshal(b, &pl)
	if err != nil {
		return err
	}

	p.LogN = pl.LogN
	p.LogNthRoot = pl.LogNthRoot
	p.Q, p.P, p.LogQ, p.LogP = pl.Q, pl.P, pl.LogQ, pl.LogP
	if pl.Xs != nil {
		p.Xs, err = ring.ParametersFromMap(pl.Xs)
		if err != nil {
			return err
		}
	}
	if pl.Xe != nil {
		p.Xe, err = ring.ParametersFromMap(pl.Xe)
		if err != nil {
			return err
		}
	}
	p.RingType = pl.RingType
	p.LogDefaultScale = pl.LogDefaultScale
	return err
}
