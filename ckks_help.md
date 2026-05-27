# CKKS同态加密方案 - Lattigo库学习指南

## 目录

1. [概述](#1-概述)
2. [CKKS方案背景](#2-ckks方案背景)
3. [核心模块详解](#3-核心模块详解)
4. [编码(Encode/Decode)](#4-编码encodedecode)
5. [加密/解密(Encrypt/Decrypt)](#5-加密解密encryptdecrypt)
6. [同态运算(Evaluator)](#6-同态运算evaluator)
7. [重线性化(Relinearization)](#7-重线性化relinearization)
8. [重缩放(Rescaling)](#8-重缩放rescaling)
9. [Bootstrapping](#9-bootstrapping)
10. [参数设置](#10-参数设置)
11. [术语表](#11-术语表)

---

## 1. 概述

本指南帮助您学习Lattigo库中CKKS同态加密方案的核心模块。通过本指南，您将理解：
- CKKS方案的基本原理
- Lattigo库中的关键数据结构和方法
- 如何使用Lattigo进行同态计算
- Bootstrapping的工作机制

### 参考文献

本指南基于以下论文：

- **[CKKS2017]** Cheon, J. H., Kim, A., Kim, M., & Song, Y. (2017). *Homomorphic Encryption for Arithmetic of Approximate Numbers*. ASIACRYPT 2017.
- **[BV2014]** Brakerski, Z., & Vaikuntanathan, V. (2014). *Efficient Fully Homomorphic Encryption from (Standard) LWE*. FOCS 2014.
- **[BGV2014]** Brakerski, Z., Gentry, C., & Vaikuntanathan, V. (2014). *(Leveled) Fully Homomorphic Encryption without Bootstrapping*. ITCS 2014.
- **[HEAAN2019]** Cheon, J. H., et al. (2019). *Better Bootstrapping for Approximate Homomorphic Encryption*. SAC 2019.

---

## 2. CKKS方案背景

### 2.1 问题陈述

CKKS方案解决的问题是：在密文上直接进行复数的近似算术运算（包括加法、乘法）。

### 2.2 核心数学对象

**环 (Ring)**:
$$R = \mathbb{Z}_q[X]/(X^N + 1)$$

其中 $N = 2^n$ 是多项式次数，$q$ 是模数。

**消息编码**:
消息向量 $\mathbf{m} = (m_1, ..., m_{N/2}) \in \mathbb{C}^{N/2}$ 被编码为环 $R$ 上的多项式。

### 2.3 CKKS算法概述

```
┌─────────────────────────────────────────────────────────────┐
│                    CKKS 加密流程                            │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  消息 m ──┬──> Encode ──> 明文 pt ──> Encrypt ──> 密文 ct  │
│           │                          ↑                      │
│           │                          │                      │
│           └──> Scale(×Δ) ────────────┘                      │
│                                                             │
│  密文 ct ──┬──> Decrypt ──> 解密多项式 ──> Decode ──> 消息  │
│           │                                                    │
│           └──> Scale(÷Δ) ────────────┘                        │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

---

## 3. 核心模块详解

### 3.1 模块结构

Lattigo库中CKKS相关的核心文件位于以下目录：

```
schemes/ckks/           # CKKS方案核心实现
├── ckks.go            # 参数定义和常量
├── params.go          # 参数管理
├── encoder.go        # 编码/解码器
├── evaluator.go      # 同态运算评估器
├── bridge.go        # 域切换(标准/共轭不变)
├── linear_transformation.go  # 线性变换
├── precision.go     # 精度统计
└── scaling.go      # 缩放因子管理

circuits/ckks/bootstrapping/  # Bootstrapping实现
├── bootstrapping.go          # Bootstrapper接口
├── evaluator.go             # Bootstrapping评估器
├── parameters.go            # Bootstrap参数
├── parameters_literal.go    # Bootstrap参数字面量
├── keys.go                  # Bootstrap密钥
├── sk_bootstrapper.go       # 密钥切换Bootstrapper
└── default_parameters.go    # 默认参数集
```

### 3.2 参数模块 (params.go, ckks.go)

**文件路径**: `schemes/ckks/ckks.go`, `schemes/ckks/params.go`

**核心类型**:

```go
// ParametersLiteral - 用户定义的参数字面量
// 用于创建经过验证的参数集
type ParametersLiteral struct {
    LogN            int                    // 对数多项式次数 (N = 2^LogN)
    Q               []uint64               // 密文模数素数列表
    P               []uint64               // 辅助模数素数列表
    LogDefaultScale int                    // 默认缩放因子 (log2)
    RingType        ring.Type              // 环类型 (Standard/ConjugateInvariant)
    Xs              ring.DistributionParameters  // 密钥分布
    Xe              ring.DistributionParameters  // 错误分布
}

// Parameters - 经过验证的CKKS参数集
// 所有字段都是私有的和不可变的
type Parameters struct {
    rlwe.Parameters  // 嵌入的RLWE参数
}
```

**关键方法**:

| 方法名 | 功能描述 | 源文件 |
|--------|----------|--------|
| `NewParametersFromLiteral(pl)` | 从字面量创建参数集 | `ckks.go:90` |
| `MaxLevel()` | 返回最大密文级别 | `ckks.go:135` |
| `MaxDepth()` | 返回最大乘法深度 | `ckks.go:230` |
| `MaxSlots()` | 返回最大槽数量 | `ckks.go:167` |
| `LogDefaultScale()` | 返回对数缩放因子 | `ckks.go:181` |
| `PrecisionMode()` | 返回精度模式 (PREC64/PREC128) | `ckks.go:199` |

**使用示例**:

```go
// 创建参数字面量
paramsLit := ckks.ParametersLiteral{
    LogN: 14,              // N = 2^14 = 16384
    Q: []uint64{
        0x80000000080001,  // q0 ≈ 2^55
        0x2000000a0001,    // q1 ≈ 2^45
        0x2000000e0001,    // q2 ≈ 2^45
    },
    LogDefaultScale: 45,   // Δ = 2^45
    P: []uint64{
        0x80000000130001,  // p ≈ 2^55
        0x7fffffffe90001,   // p ≈ 2^55
    },
}

// 从字面量创建参数
params, err := ckks.NewParametersFromLiteral(paramsLit)
```

---

## 4. 编码(Encode/Decode)

### 4.1 功能描述

编码将复数向量转换为环 $R_q$ 上的多项式表示，解码执行反向操作。

### 4.2 算法原理

根据 **[CKKS2017]** Algorithm 1 (Encode):

```
输入: 消息向量 (m_1, ..., m_{N/2}) ∈ C^{N/2}
缩放因子 Δ

1. 特殊离散傅里叶变换 (DFT):
   对消息应用特殊的N/2点DFT，得到系数 (a_0, ..., a_{N/2-1})

2. 填充:
   将DFT结果放入多项式:
   a(X) = Σ_{i=0}^{N/2-1} a_i X^i + Σ_{i=1}^{N/2} conj(a_{N/2-i}) X^{N-i}

3. 缩放:
   返回 round(Δ · a(X)) mod q
```

### 4.3 源文件

**文件路径**: `schemes/ckks/encoder.go`

### 4.4 核心类型

```go
// Encoder - CKKS编码器
// 支持两种编码域:
// 1. Slots: 使用特殊DFT，支持SIMD并行和逐元素乘法
// 2. Coefficients: 直接编码，不支持SIMD
type Encoder struct {
    parameters Parameters  // 编码使用的参数
    prec uint              // 精度(位数)
    m int                  // 循环群阶
    rotGroup []int         // 旋转群元素
    roots interface{}      // DFT根预计算
}
```

### 4.5 关键方法

| 方法名 | 功能 | 对应论文算法 |
|--------|------|-------------|
| `NewEncoder(params)` | 创建编码器 | - |
| `Encode(values, pt)` | 将复数/浮点向量编码到明文 | Algorithm 1 |
| `Decode(pt, values)` | 从明文解码到复数/浮点向量 | Algorithm 1 |
| `Embed(values, meta, poly)` | 通用编码接口 | Algorithm 1 |
| `IFFT(values, logN)` | 特殊逆DFT | Section 3.2 |
| `FFT(values, logN)` | 特殊DFT | Section 3.2 |

**方法签名**:

```go
// Encode - 将FloatSlice编码到Plaintext
// values可以是:
//   - []complex128: 复数向量
//   - []float64: 实数向量
//   - []*big.Float: 高精度浮点向量
//   - []*bignum.Complex: 高精度复数向量
func (ecd Encoder) Encode(values interface{}, pt *rlwe.Plaintext) error

// Decode - 从Plaintext解码到FloatSlice
func (ecd Encoder) Decode(pt *rlwe.Plaintext, values interface{}) error
```

### 4.6 使用示例

```go
// 创建编码器
encoder := ckks.NewEncoder(params)

// 创建明文对象
plaintext := ckks.NewPlaintext(params, params.MaxLevel())

// 编码复数向量
values := []complex128{
    complex(1.0, 0.5),
    complex(2.0, 1.0),
    complex(3.0, 1.5),
    complex(4.0, 2.0),
}
if err := encoder.Encode(values, plaintext); err != nil {
    log.Fatal(err)
}

// 解码
decoded := make([]complex128, len(values))
if err := encoder.Decode(plaintext, decoded); err != nil {
    log.Fatal(err)
}
```

### 4.7 编码参数说明

| 参数 | 含义 | 影响 |
|------|------|------|
| `Level()` | 明文多项式的模数级别 | 影响精度和可用操作数 |
| `Scale` | 缩放因子 | 影响数值精度和噪声容忍度 |
| `LogDimensions` | 对数槽维度 | 控制编码的槽数量 |

---

## 5. 加密/解密(Encrypt/Decrypt)

### 5.1 功能描述

加密将明文转换为密文，解密执行反向操作。

### 5.2 算法原理

根据 **[CKKS2017]** Algorithm 2 (Encrypt) 和 Algorithm 3 (Decrypt):

**加密 (Encrypt)**:
```
输入: 公钥 pk = (b, a)，明文 m
输出: 密文 ct = (c_0, c_1)

1. 采样错误多项式 e_0, e_1 ~ χ_error (离散高斯)

2. 计算:
   c_0 = m + e_0 - a·s
   c_1 = e_1 - b·s
   (其中 s 是私钥)

返回 ct = (c_0, c_1)
```

**解密 (Decrypt)**:
```
输入: 私钥 s，密文 ct = (c_0, c_1)
输出: 消息 m'

m' = c_0 + c_1·s = m + e_0 + e_1·e  (包含噪声)
```

### 5.3 源文件

**文件路径**: `schemes/ckks/params.go`

### 5.4 核心函数

| 函数名 | 功能 | 源文件 |
|--------|------|--------|
| `NewEncryptor(params, key)` | 创建加密器 | `params.go:51` |
| `NewDecryptor(params, key)` | 创建解密器 | `params.go:62` |
| `NewPlaintext(params, level)` | 创建明文 | `params.go:20` |
| `NewCiphertext(params, degree, level)` | 创建密文 | `params.go:36` |

**类型签名**:

```go
// Plaintext - CKKS明文
// 基于rlwe.Plaintext, 添加CKKS特定元数据
type Plaintext = rlwe.Plaintext

// Ciphertext - CKKS密文
// 基于rlwe.Ciphertext
type Ciphertext = rlwe.Ciphertext

// Encryptor - 加密器
type Encryptor = rlwe.Encryptor

// Decryptor - 解密器
type Decryptor = rlwe.Decryptor
```

### 5.5 关键方法

加密器方法（在`rlwe.Encryptor`中）:

```go
// EncryptNew - 加密明文返回新密文
func (enc *Encryptor) EncryptNew(pt *Plaintext) (*Ciphertext, error)

// Encrypt - 加密明文到指定密文
func (enc *Encryptor) Encrypt(pt *Plaintext, ct *Ciphertext) error

// EncryptWithPrng - 使用指定PRNG加密
func (enc *Encryptor) EncryptWithPrng(pt *Plaintext, prng sampling.PRNG, ct *Ciphertext) error
```

解密器方法（在`rlwe.Decryptor`中）:

```go
// DecryptNew - 解密密文返回新明文
func (dec *Decryptor) DecryptNew(ct *Ciphertext) (*Plaintext, error)

// Decrypt - 解密密文到指定明文
func (dec *Decryptor) Decrypt(ct *Ciphertext, pt *Plaintext) error
```

### 5.6 使用示例

```go
// 生成密钥
keyGenerator := rlwe.NewKeyGenerator(params)
secretKey := keyGenerator.GenSecretKeyNew()        // 私钥
publicKey := keyGenerator.GenPublicKeyNew(secretKey) // 公钥

// 创建加密器/解密器
encryptor := ckks.NewEncryptor(params, publicKey)
decryptor := ckks.NewDecryptor(params, secretKey)

// 编码消息
encoder := ckks.NewEncoder(params)
plaintext := ckks.NewPlaintext(params, params.MaxLevel())
encoder.Encode(values, plaintext)

// 加密
ciphertext, _ := encryptor.EncryptNew(plaintext)

// 解密
decryptedPlain, _ := decryptor.DecryptNew(ciphertext)

// 解码
decodedValues := make([]complex128, len(values))
encoder.Decode(decryptedPlain, decodedValues)
```

---

## 6. 同态运算(Evaluator)

### 6.1 功能描述

评估器提供在密文上的同态运算操作，包括加法、减法、乘法、旋转等。

### 6.2 源文件

**文件路径**: `schemes/ckks/evaluator.go`

### 6.3 核心类型

```go
// Evaluator - CKKS同态运算评估器
// 组合了编码器和RLWE评估器
type Evaluator struct {
    *Encoder                // CKKS编码器
    *rlwe.Evaluator         // RLWE运算评估器
    pool *rlwe.BufferPool   // 内存缓冲池
}
```

### 6.4 关键方法

#### 6.4.1 加法和减法

```go
// Add - 同态加法
// ct_out = ct_a + ct_b 或 ct_a + scalar
func (eval Evaluator) Add(op0 *Ciphertext, op1 Operand, opOut *Ciphertext) error

// AddNew - 同态加法返回新密文
func (eval Evaluator) AddNew(op0 *Ciphertext, op1 Operand) (*Ciphertext, error)

// Sub - 同态减法
func (eval Evaluator) Sub(op0 *Ciphertext, op1 Operand, opOut *Ciphertext) error

// SubNew - 同态减法返回新密文
func (eval Evaluator) SubNew(op0 *Ciphertext, op1 Operand) (*Ciphertext, error)
```

#### 6.4.2 乘法

```go
// Mul - 同态乘法 (不含重线性化)
// ct_out = ct_a * ct_b (度变为2, 需要后续处理)
func (eval Evaluator) Mul(op0 *Ciphertext, op1 Operand, opOut *Ciphertext) error

// MulNew - 同态乘法返回新密文
func (eval Evaluator) MulNew(op0 *Ciphertext, op1 Operand) (*Ciphertext, error)

// MulRelin - 同态乘法 (含重线性化)
// ct_out = ct_a * ct_b (度保持为1)
func (eval Evaluator) MulRelin(op0 *Ciphertext, op1 Operand, opOut *Ciphertext) error

// MulRelinNew - 同态乘法(含重线性化)返回新密文
func (eval Evaluator) MulRelinNew(op0 *Ciphertext, op1 Operand) (*Ciphertext, error)
```

#### 6.4.3 重缩放

```go
// Rescale - 重缩放操作
// 将scale除以当前模数, 降低level
// 对应论文 Section 3.3
func (eval Evaluator) Rescale(op0, opOut *Ciphertext) error

// RescaleTo - 重缩放到指定scale
func (eval Evaluator) RescaleTo(op0 *Ciphertext, minScale Scale, opOut *Ciphertext) error

// ScaleUp - 向上缩放
func (eval Evaluator) ScaleUp(op0 *Ciphertext, scale Scale, opOut *Ciphertext) error
```

#### 6.4.4 旋转和共轭

```go
// Rotate - 循环旋转槽
// k > 0: 左旋k个位置
// k < 0: 右旋|k|个位置
func (eval Evaluator) Rotate(op0 *Ciphertext, k int, opOut *Ciphertext) error

// RotateNew - 旋转返回新密文
func (eval Evaluator) RotateNew(op0 *Ciphertext, k int) (*Ciphertext, error)

// Conjugate - 复数共轭
func (eval Evaluator) Conjugate(op0 *Ciphertext, opOut *Ciphertext) error

// ConjugateNew - 共轭返回新密文
func (eval Evaluator) ConjugateNew(op0 *Ciphertext) (*Ciphertext, error)
```

#### 6.4.5 其他操作

```go
// DropLevel - 降低level(不重缩放)
func (eval Evaluator) DropLevel(op0 *Ciphertext, levels int)

// InnerSum - 内和运算
// 对一批槽求和
func (eval Evaluator) InnerSum(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) error

// RotateAndAdd - 旋转并累加
func (eval Evaluator) RotateAndAdd(ctIn *Ciphertext, batchSize, n int, opOut *Ciphertext) error
```

### 6.5 乘法算法原理

根据 **[CKKS2017]** Algorithm 4 (Multiply):

```
输入: 密文 ct = (c_0, c_1), ct' = (c'_0, c'_1)
输出: 密文 ct'' = (d_0, d_1, d_2)

1. 乘法:
   d_0 = c_0 * c'_0
   d_1 = c_0 * c'_1 + c_1 * c'_0
   d_2 = c_1 * c'_1

2. 返回 ct'' = (d_0, d_1, d_2)
   注意: 度变为2, 需要重线性化或重缩放
```

### 6.6 使用示例

```go
// 创建评估器(需要评估密钥)
evaluator := ckks.NewEvaluator(params, evaluationKeys)

// 加法
ctAdd, _ := evaluator.AddNew(ctA, ctB)

// 常数乘法
ctScale, _ := evaluator.MulNew(ctA, complex(2.0, 0))

// 密文乘法
ctMul, _ := evaluator.MulNew(ctA, ctB)

// 重缩放(每次乘法后必须执行)
evaluator.Rescale(ctMul, ctMul)

// 旋转
ctRot, _ := evaluator.RotateNew(ctA, 1)  // 左旋1个位置
```

---

## 7. 重线性化(Relinearization)

### 7.1 功能描述

密文乘法后，多项式度从1变为2。重线性化将密文度降回1，控制密文大小。

### 7.2 算法原理

根据 **[BGV2014]** Section 3:

```
输入: 密文 ct = (c_0, c_1, c_2)，重线性化密钥 rlk
输出: 密文 ct' = (c'_0, c'_1)

1. 分解: 将 c_2 分解为基b下的系数
   c_2 = Σ_i c_2[i] * b^i

2. 密钥切换:
   使用rlk计算 c'_0 = c_0 + Σ_i c_2[i] * rlk[i].c_0
               c'_1 = c_1 + Σ_i c_2[i] * rlk[i].c_1

返回 ct' = (c'_0, c'_1)
```

### 7.3 源文件

**文件路径**: `core/rlwe/keygenerator.go`, `core/rlwe/keys.go`

### 7.4 密钥生成

```go
// GenRelinearizationKey - 生成重线性化密钥
func (kgen *KeyGenerator) GenRelinearizationKeyNew(sk *SecretKey) *RelinearizationKey

// GenEvaluationKeys - 生成完整评估密钥集
func (kgen *KeyGenerator) GenEvaluationKeysNew(params Parameters) *EvaluationKeys
```

### 7.5 使用示例

```go
// 生成评估密钥(包括重线性化密钥和Galois密钥)
keyGenerator := rlwe.NewKeyGenerator(params)
sk := keyGenerator.GenSecretKeyNew()
rlk := keyGenerator.GenRelinearizationKeyNew(sk)

// 创建带有重线性化密钥的评估器
evk := rlwe.NewMemEvaluationKeySet(rlk)
evaluator := ckks.NewEvaluator(params, evk)

// 使用MulRelin进行乘法(自动重线性化)
ctResult, _ := evaluator.MulRelinNew(ctA, ctB)
```

---

## 8. 重缩放(Rescaling)

### 8.1 功能描述

每次密文乘法后，scale会变为原来的平方倍。重缩放将scale除以当前模数，同时降低level，控制噪声增长。

### 8.2 算法原理

根据 **[CKKS2017]** Section 3.3:

```
输入: 密文 ct = (c_0, c_1)，当前level = i
缩放因子 Δ

输出: 密文 ct' = (c'_0, c'_1)

1. 缩放:
   对每个多项式系数除以 q_i (当前模数) 并取整
   c'_j = round(c_j / q_i)

2. 更新scale:
   Δ' = Δ / q_i

3. 降低level:
   返回 level = i - 1 的密文

注意: 重缩放会消耗一个level(模数)
```

### 8.3 缩放因子管理

```go
// Scale类型 - 表示缩放因子
type Scale struct {
    Value *big.Float  // 实际缩放值
}

// NewScale - 从数值创建缩放因子
func NewScale(v interface{}) Scale

// Mul - 缩放因子乘法
func (s Scale) Mul(t Scale) Scale

// Div - 缩放因子除法
func (s Scale) Div(t Scale) Scale

// Log2 - 返回log2(缩放因子)
func (s Scale) Log2() float64
```

### 8.4 精度模式

| 模式 | 描述 | 每个rescale消耗的level |
|------|------|------------------------|
| `PREC64` | 64位精度，默认 | 1 |
| `PREC128` | 128位精度 | 2 |

当`LogDefaultScale > 64`时，自动使用`PREC128`模式。

### 8.5 使用示例

```go
// 乘法后必须重缩放
ctMul, _ := evaluator.MulNew(ctA, ctB)
fmt.Printf("乘法后scale: 2^%.0f\n", ctMul.Scale.Log2())

evaluator.Rescale(ctMul, ctMul)
fmt.Printf("Rescale后scale: 2^%.0f, level: %d\n", ctMul.Scale.Log2(), ctMul.Level())

// 继续下一轮乘法
ctMul2, _ := evaluator.MulNew(ctMul, ctB)
evaluator.Rescale(ctMul2, ctMul2)
```

---

## 9. Bootstrapping

### 9.1 功能描述

Bootstrapping是一种重置同态计算能力的技术。通过在密文上执行解密操作，将level恢复到最大，同时清除大部分噪声。

### 9.2 为什么需要Bootstrapping

```
┌────────────────────────────────────────────────────────────┐
│                    密文容量消耗                              │
├────────────────────────────────────────────────────────────┤
│                                                             │
│  初始: [Level=Max, Scale=Δ, Noise=ε]                        │
│    │                                                         │
│    ├──> 加法 ──> [Level不变, Scale不变, Noise += ε']         │
│    │                                                         │
│    ├──> 乘法 ──> [Level-1, Scale/Δ, Noise *= Δ]              │
│    │                                                         │
│    └──> Rescale ──> [Level-1, Scale=q_i, Noise控制]          │
│                                                             │
│  问题: Level最终会耗尽, 此时无法再执行乘法                    │
│  解决: Bootstrapping                                          │
│                                                             │
└────────────────────────────────────────────────────────────┘
```

### 9.3 算法原理

根据 **[HEAAN2019]** 和 **[CKKS2017]**:

```
输入: 密文 ct 位于 level 0
输出: 密文 ct' 位于最大 level

Bootstrapping电路包含5个步骤:

┌──────────────────────────────────────────────────────────────┐
│ 1. ScaleDown                                                 │
│    - 将密文缩放到 q/|m| 的范围                                │
│    - 使消息标准化到固定区间                                   │
├──────────────────────────────────────────────────────────────┤
│ 2. ModUp                                                     │
│    - 将模数从 q 扩展到 Q (整个模数链)                          │
│    - 使用密钥切换技术实现                                     │
├──────────────────────────────────────────────────────────────┤
│ 3. CoeffsToSlots (同态编码)                                   │
│    - 使用特殊的同态DFT                                        │
│    - 将系数表示转换为槽表示                                    │
│    - 使得后续操作可以逐元素进行                                │
├──────────────────────────────────────────────────────────────┤
│ 4. EvalMod (同态模约简)                                       │
│    - 核心步骤: 计算 m mod 1                                    │
│    - 使用多项式近似 (Chebyshev/Cosine级数)                     │
│    - 配合Double-Angle公式减少多项式度                           │
├──────────────────────────────────────────────────────────────┤
│ 5. SlotsToCoeffs (同态解码)                                   │
│    - 使用特殊的同态IDFT                                        │
│    - 将槽表示转换回系数表示                                    │
└──────────────────────────────────────────────────────────────┘
```

### 9.4 源文件

**文件路径**: `circuits/ckks/bootstrapping/`

| 文件 | 功能 |
|------|------|
| `bootstrapping.go` | Bootstrapper主接口 |
| `evaluator.go` | Bootstrapping评估器实现 |
| `parameters.go` | Bootstrap参数管理 |
| `parameters_literal.go` | Bootstrap参数字面量 |
| `keys.go` | Bootstrap密钥生成 |
| `sk_bootstrapper.go` | SK切换Bootstrapper |

### 9.5 核心类型

```go
// Parameters - Bootstrapping参数
// 包含原始(残差)参数和Bootstrapping电路参数
type Parameters struct {
    ResidualParameters      ckks.Parameters  // 原始参数
    BootstrappingParameters ckks.Parameters  // Bootstrap电路参数
    SlotsToCoeffsParameters dft.MatrixLiteral  // 解码矩阵参数
    Mod1ParametersLiteral   mod1.ParametersLiteral  // 模约简参数
    CoeffsToSlotsParameters dft.MatrixLiteral  // 编码矩阵参数
    EphemeralSecretWeight   int  // 短暂密钥权重
}

// Evaluator - Bootstrapping评估器
type Evaluator struct {
    Parameters
    *ckks.Evaluator
    DFTEvaluator  *dft.Evaluator
    Mod1Evaluator *mod1.Evaluator
    *EvaluationKeys
    // ...
}
```

### 9.6 关键方法

```go
// NewEvaluator - 创建Bootstrapping评估器
func NewEvaluator(btpParams Parameters, evk *EvaluationKeys) (*Evaluator, error)

// Bootstrap - 对单个密文执行Bootstrapping
func (eval Evaluator) Bootstrap(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// BootstrapMany - 对多个密文执行Bootstrapping
func (eval Evaluator) BootstrapMany(cts []rlwe.Ciphertext) ([]rlwe.Ciphertext, error)

// Evaluate - 执行Bootstrapping电路
func (eval Evaluator) Evaluate(ctIn *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// ScaleDown - Bootstrapping第1步
// EvalMod - Bootstrapping第4步
// CoeffsToSlots - Bootstrapping第3步
// SlotsToCoeffs - Bootstrapping第5步
```

### 9.7 参数配置

```go
// ParametersLiteral - Bootstrap参数字面量
type ParametersLiteral struct {
    LogN  *int  // Bootstrap环度, 默认16
    
    LogSlots  *int  // 槽数量, 默认LogN-1
    
    // 编码/解码深度和缩放
    CoeffsToSlotsFactorizationDepthAndLogScales [][]int  // 默认: 4层, 每层2^56
    SlotsToCoeffsFactorizationDepthAndLogScales [][]int  // 默认: 3层, 每层2^39
    
    // 模约简参数
    EvalModLogScale *int  // EvalMod缩放, 默认60
    Mod1Type mod1.Type  // 近似类型: CosDiscrete/SinContinuous
    K *int  // 近似区间, 默认16
    Mod1Degree *int  // 多项式度, 默认30
    DoubleAngle *int  // Double-Angle迭代次数, 默认3
    
    // 安全参数
    EphemeralSecretWeight *int  // 短暂密钥权重, 默认32
    
    // 精度提升
    IterationsParameters *IterationsParameters  // 迭代参数
}

// IterationsParameters - 用于提升精度的迭代参数
type IterationsParameters struct {
    BootstrappingPrecision []float64  // 每轮精度, 如[]float64{16, 16}
    ReservedPrimeBitSize int  // 预留素数大小
}
```

### 9.8 使用示例

```go
// 1. 创建原始CKKS参数
schemeParams := ckks.ParametersLiteral{
    LogN: 14,
    LogQ: []int{50, 40, 40},
    LogDefaultScale: 40,
}
params, _ := ckks.NewParametersFromLiteral(schemeParams)

// 2. 创建Bootstrap参数
btpParamsLit := bootstrapping.ParametersLiteral{}
btpParams, _ := bootstrapping.NewParametersFromLiteral(params, btpParamsLit)

// 3. 生成Bootstrap密钥
sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()
btpKeys, _, _ := btpParams.GenEvaluationKeys(sk)

// 4. 创建Bootstrap评估器
btpEvaluator, _ := bootstrapping.NewEvaluator(btpParams, btpKeys)

// 5. 加密
encoder := ckks.NewEncoder(params)
encryptor := ckks.NewEncryptor(params, sk)
plaintext := ckks.NewPlaintext(params, 0)
encoder.Encode(values, plaintext)
ct, _ := encryptor.EncryptNew(plaintext)

// 6. 执行乘法耗尽level
ct1, _ := evaluator.MulNew(ct, ct)
evaluator.Rescale(ct1, ct1)
ct2, _ := evaluator.MulNew(ct1, ct1)
evaluator.Rescale(ct2, ct2)

// 7. Bootstrap恢复level
ctBootstrap, _ := btpEvaluator.Bootstrap(ct2)
fmt.Printf("Bootstrap后level: %d (恢复到%d)\n", 
    ctBootstrap.Level(), params.MaxLevel())

// 8. 继续计算
ct3, _ := evaluator.MulNew(ctBootstrap, ct)
```

---

## 10. 参数设置

### 10.1 参数选择原则

| 参数 | 建议值 | 说明 |
|------|--------|------|
| LogN | 13-16 | 安全边际随N增加而增加 |
| LogDefaultScale | 40-60 | 影响精度,应接近素数大小 |
| Q素数大小 | 45-60位 | 影响level消耗速度 |
| LogSlots | LogN-1 | 最大槽数量 |
| P素数大小 | 55-65位 | sqrt(#Qi)个素数 |

### 10.2 深度计算

```
最大乘法深度 = (|Q| - 1) / LevelsConsumedPerRescaling

其中:
- |Q|: Q中素数的个数
- LevelsConsumedPerRescaling: PREC64为1, PREC128为2
```

### 10.3 安全边际

Lattigo使用 **[USEC2021]** 的安全估计:

```go
// params.T() 返回估计的安全边际(比特数)
// 标准建议: >= 128比特
```

### 10.4 示例参数集

```go
// 小规模演示参数 (LogN=10)
smallParams := ckks.ParametersLiteral{
    LogN: 10,
    Q: []uint64{0x80000000080001, 0x2000000a0001},
    LogDefaultScale: 40,
    P: []uint64{0x80000000130001},
}

// 标准参数 (LogN=14)
standardParams := ckks.ParametersLiteral{
    LogN: 14,
    Q: []uint64{
        0x80000000080001,  // 55位
        0x2000000a0001,    // 45位
        0x2000000e0001,    // 45位
        0x2000001d0001,    // 45位
        0x1fffffcf0001,    // 45位
    },
    LogDefaultScale: 45,
    P: []uint64{
        0x80000000130001,
        0x7fffffffe90001,
    },
}

// 高精度参数 (LogN=15)
highPrecisionParams := ckks.ParametersLiteral{
    LogN: 15,
    LogQ: []int{60, 50, 50, 50, 50},  // 使用比特长度
    LogDefaultScale: 50,
}
```

---

## 11. 术语表

| 英文术语 | 中文术语 | 说明 |
|----------|----------|------|
| Ring (R) | 环 | $R = \mathbb{Z}_q[X]/(X^N+1)$ |
| Polynomial degree (N) | 多项式次数 | $N = 2^n$ |
| Ciphertext modulus (Q) | 密文模数 | 素数乘积 |
| Auxiliary modulus (P) | 辅助模数 | 用于密钥交换 |
| Scale (Δ) | 缩放因子 | 控制精度和噪声 |
| Slot | 槽 | SIMD并行单元 |
| Level | 层级 | 模数链中的位置 |
| Depth | 深度 | 可执行乘法次数 |
| Noise | 噪声 | 加密引入的误差 |
| Error | 错误 | 同noise |
| Relinearization | 重线性化 | 降低密文度 |
| Rescaling | 重缩放 | 消耗level控制噪声 |
| Bootstrapping | 自举 | 重置计算能力 |
| Gadget decomposition | 工具分解 | 密钥切换技术 |
| NTT | 快速数论变换 | 多项式乘法优化 |
| Galois key | Galois密钥 | 旋转操作密钥 |
| Evaluation key | 评估密钥 | 同公钥,用于同态运算 |
| Secret key (sk) | 私钥 | 解密用 |
| Public key (pk) | 公钥 | 加密用 |
| Plaintext (pt) | 明文 | 编码后的消息 |
| Ciphertext (ct) | 密文 | 加密后的密文 |
| Encoder | 编码器 | 消息→明文 |
| Decoder | 解码器 | 明文→消息 |
| Encryptor | 加密器 | 明文→密文 |
| Decryptor | 解密器 | 密文→明文 |
| Evaluator | 评估器 | 同态运算 |
| SIMD | 单指令多数据 | 并行计算技术 |
| DFT | 离散傅里叶变换 | 编码用 |
| RNS | 余数系统 | 模数表示 |
| Conjugate Invariant | 共轭不变 | 特殊环变体 |

---

## 附录A: Lattigo核心文件索引

### CKKS方案 (`schemes/ckks/`)

| 文件 | 关键类型/函数 | 行号 |
|------|--------------|------|
| `ckks.go` | `ParametersLiteral`, `Parameters`, `NewParametersFromLiteral` | 35-101 |
| `params.go` | `NewPlaintext`, `NewCiphertext`, `NewEncryptor`, `NewDecryptor` | 20-74 |
| `encoder.go` | `Encoder`, `Encode`, `Decode`, `IFFT`, `FFT` | 58-821 |
| `evaluator.go` | `Evaluator`, `Add`, `Sub`, `Mul`, `Rescale`, `Rotate` | 16-1329 |

### Bootstrapping (`circuits/ckks/bootstrapping/`)

| 文件 | 关键类型/函数 | 行号 |
|------|--------------|------|
| `parameters.go` | `Parameters`, `NewParametersFromLiteral` | 18-355 |
| `parameters_literal.go` | `ParametersLiteral`, `GetLogN`, `GetEvalMod1LogScale` | 125-458 |
| `evaluator.go` | `Evaluator`, `Bootstrap`, `Evaluate`, `ScaleDown`, `EvalMod` | 22-1016 |
| `keys.go` | `EvaluationKeys`, `GenEvaluationKeys` | 15-119 |

### RLWE核心 (`core/rlwe/`)

| 文件 | 关键类型/函数 | 说明 |
|------|--------------|------|
| `params.go` | `RLWEParameters`, 安全性估计 | RLWE参数管理 |
| `keygenerator.go` | `KeyGenerator`, `GenSecretKeyNew`, `GenRelinearizationKeyNew` | 密钥生成 |
| `encryptor.go` | `Encryptor`, `EncryptNew` | 加密 |
| `decryptor.go` | `Decryptor`, `DecryptNew` | 解密 |
| `evaluator.go` | `Evaluator`, 基础同态运算 | RLWE运算 |
| `ciphertext.go` | `Ciphertext`, 密文结构 | 密文表示 |

### Ring模块 (`ring/`)

| 文件 | 关键类型/函数 | 说明 |
|------|--------------|------|
| `ring.go` | `Ring`, 多项式环运算 | 环运算基础 |
| `ntt.go` | `NTT`, `INTT` | 数论变换 |

---

## 附录B: 关键公式索引

| 公式名称 | 数学表达式 | 位置 |
|----------|-----------|------|
| 环定义 | $R = \mathbb{Z}_q[X]/(X^N+1)$ | `ckks.go` |
| 缩放编码 | $m' = \Delta \cdot m$ | `encoder.go` |
| 加密 | $ct = (c_0, c_1) = (m+e_0-a\cdot s, e_1-b\cdot s)$ | `rlwe/encryptor.go` |
| 解密 | $m' = c_0 + c_1 \cdot s = m + error$ | `rlwe/decryptor.go` |
| 乘法 | $ct_{mul} = (d_0, d_1, d_2)$ | `ckks/evaluator.go` |
| Rescale | $ct' = (c_0/q_i, c_1/q_i), \Delta' = \Delta/q_i$ | `ckks/evaluator.go` |

---

*本指南由Claude编写, 基于Lattigo v6源代码和CKKS相关论文。*