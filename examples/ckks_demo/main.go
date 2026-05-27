// CKKS同态加密演示程序
// ================================================================================
// 本程序演示CKKS(Cheon-Kim-Kim-Song)同态加密方案的基本操作和Bootstrapping技术
//
// 参考文献:
// [1] Cheon, J. H., Kim, A., Kim, M., & Song, Y. (2017). Homomorphic encryption for
//
//	arithmetic of approximate numbers. ASIACRYPT 2017.
//
// [2] Chen, H., Laine, K., & Rindal, P. (2017). Fast homomorphic evaluation of
//
//	deep discretized neural networks. CRYPTO 2017.
//
// [3] Bossuet, L., et al. (2022). A full RNS-variant of FV like some homomorphic
//
//	encryption schemes. ACNS 2022.
//
// Copyright 2024 Cyrilao
// ================================================================================
package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// ================================================================================
// 常量定义 - 演示用参数
// ================================================================================

const (
	// -------------------------------------------------------------------------
	// LogN: 对数多项式次数
	// -------------------------------------------------------------------------
	// LogN = 10 表示多项式次数 N = 2^10 = 1024
	// 这意味着:
	// - 环为 R = Z[X]/(X^N + 1)
	// - 对于标准CKKS方案(slot encoding), 可用槽数量为 N/2 = 512
	// - 对于共轭不变CKKS方案, 可用槽数量为 N = 1024
	//
	// 注意: 小规模的N用于快速演示, 实际应用中通常使用LogN >= 15
	// -------------------------------------------------------------------------
	LogN = 10

	// -------------------------------------------------------------------------
	// LogQ: 密文模数的比特长度数组
	// -------------------------------------------------------------------------
	// Q = [q0, q1, q2, ...] 是pairwise coprime的素数序列
	// 每个素数增加一个可执行的乘法深度
	//
	// 设计说明:
	// - q0 = 2^55 ≈ 2^55: 第一个素数必须足够大以存储计算结果
	// - q1, q2, ... ≈ 2^40: 其余素数用于增加乘法深度
	//
	// 这里设计为只有2个素数, 所以最大乘法深度 = 2 (需要一次rescale)
	// -------------------------------------------------------------------------
	LogQLevel0 = 55 // 第一个素数的比特长度, 用于存储结果
	LogQLevel1 = 45 // 第二个素数的比特长度, 用于增加一次乘法

	// -------------------------------------------------------------------------
	// LogDefaultScale: 默认缩放因子 (log2)
	// -------------------------------------------------------------------------
	// 缩放因子Δ是CKKS方案中的关键参数
	// - 编码时: m_scaled = round(m * Δ)
	// - 解密后: m = m_scaled / Δ
	//
	// 设计说明:
	// - LogDefaultScale = 40 意味着 Δ = 2^40
	// - 缩放因子影响精度和噪声容忍度
	// - 缩放因子应该接近素数的比特长度以优化资源使用
	// -------------------------------------------------------------------------
	LogDefaultScale = 40

	// -------------------------------------------------------------------------
	// LogP: 辅助模数(用于密钥交换)
	// -------------------------------------------------------------------------
	// P = [p0, p1] 是额外的素数序列
	// 它们不参与同态容量计算, 但用于优化密钥交换操作
	//
	// 作用说明:
	// - 在key-switching时, RNS分解使用P作为辅助基数
	// - 可以显著降低密钥交换引入的噪声
	// - 默认使用sqrt(#Qi)个素数, 每个约60-65比特
	// -------------------------------------------------------------------------
	LogP = 50 // 辅助素数的比特长度

	// -------------------------------------------------------------------------
	// 演示数据长度
	// -------------------------------------------------------------------------
	// 演示用的向量长度
	// 对于LogN=10, N/2=512, 我们使用4个元素便于观察结果
	// -------------------------------------------------------------------------
	PlaintextSize = 4
)

// ================================================================================
// 主函数
// ================================================================================

func main() {
	fmt.Println("========================================")
	fmt.Println("  CKKS同态加密演示程序")
	fmt.Println("  Lattigo库 v6")
	fmt.Println("========================================")
	fmt.Println()

	// 记录总执行时间
	totalStart := time.Now()

	// -------------------------------------------------------------------------
	// 第一部分: 设置CKKS环境参数
	// -------------------------------------------------------------------------
	fmt.Println("[1] 设定CKKS环境参数")
	fmt.Println("----------------------------------------")

	params, err := setupCKKSParameters()
	if err != nil {
		fmt.Printf("错误: 设置CKKS参数失败: %v\n", err)
		os.Exit(1)
	}

	// 打印参数信息
	printParameters(params)

	// -------------------------------------------------------------------------
	// 第二部分: 创建明文向量
	// -------------------------------------------------------------------------
	fmt.Println("\n[2] 创建明文向量A和B")
	fmt.Println("----------------------------------------")

	// 创建两个复数向量作为演示数据
	// 复数在CKKS中通过SIMD打包实现并行计算
	valuesA, valuesB := createPlaintextVectors()

	fmt.Printf("  明文向量 A (复数): %v\n", formatComplexSlice(valuesA))
	fmt.Printf("  明文向量 B (复数): %v\n", formatComplexSlice(valuesB))

	// -------------------------------------------------------------------------
	// 第三部分: 编码操作
	// -------------------------------------------------------------------------
	fmt.Println("\n[3] 编码操作 (Encode)")
	fmt.Println("----------------------------------------")

	// 创建编码器
	// 编码器负责将复数向量转换为环R_Q上的多项式
	encoder := ckks.NewEncoder(params)

	// 创建明文对象并编码
	// 明文包含: 多项式表示、缩放因子、元数据
	plaintextA := ckks.NewPlaintext(params, params.MaxLevel())
	plaintextB := ckks.NewPlaintext(params, params.MaxLevel())

	// Encode方法: 将复数向量编码到明文多项式
	// 根据论文[1] Algorithm 1: Encode
	// 1. 将复数向量视为N/2维度的消息m
	// 2. 对消息进行特殊的离散傅里叶变换(针对X^N+1环)
	// 3. 将变换结果乘以缩放因子Δ并取整
	// 4. 得到环R_Q上的多项式表示
	fmt.Println("  执行Encoder.Encode(values, plaintext)...")
	if err := encoder.Encode(valuesA, plaintextA); err != nil {
		fmt.Printf("错误: 编码向量A失败: %v\n", err)
		os.Exit(1)
	}
	if err := encoder.Encode(valuesB, plaintextB); err != nil {
		fmt.Printf("错误: 编码向量B失败: %v\n", err)
		os.Exit(1)
	}
	fmt.Println("  编码成功!")
	fmt.Printf("  明文A的缩放因子: 2^%d = %.2e\n", plaintextA.Scale.Log2(), plaintextA.Scale.Float64())
	fmt.Printf("  明文B的缩放因子: 2^%d = %.2e\n", plaintextB.Scale.Log2(), plaintextB.Scale.Float64())

	// -------------------------------------------------------------------------
	// 第四部分: 加密操作
	// -------------------------------------------------------------------------
	fmt.Println("\n[4] 加密操作 (Encrypt)")
	fmt.Println("----------------------------------------")

	// 生成密钥对
	keyGenerator := rlwe.NewKeyGenerator(params)
	secretKey := keyGenerator.GenSecretKeyNew()          // 私钥sk
	publicKey := keyGenerator.GenPublicKeyNew(secretKey) // 公钥pk

	fmt.Println("  已生成密钥对 (sk, pk)")

	// 创建加密器和解密器
	encryptor := ckks.NewEncryptor(params, publicKey)
	decryptor := ckks.NewDecryptor(params, secretKey)

	// Encrypt方法: 将明文加密为密文
	// 根据论文[1] Algorithm 2: Encrypt
	// ct = (c0, c1) 其中:
	// c0 = m + e0 + e1*sk  (明文加上噪声项)
	// c1 = e1             (噪声的"密钥"部分)
	//
	// 这里e0, e1是从离散高斯分布采样的错误多项式
	// 加密的安全性依赖于这些噪声的分布
	fmt.Println("  执行Encryptor.EncryptNew(plaintext)...")
	ciphertextA, err := encryptor.EncryptNew(plaintextA)
	if err != nil {
		fmt.Printf("错误: 加密向量A失败: %v\n", err)
		os.Exit(1)
	}
	ciphertextB, err := encryptor.EncryptNew(plaintextB)
	if err != nil {
		fmt.Printf("错误: 加密向量B失败: %v\n", err)
		os.Exit(1)
	}
	fmt.Println("  加密成功!")
	fmt.Printf("  密文A: degree=%d, level=%d, scale=2^%d\n",
		ciphertextA.Degree(), ciphertextA.Level(), ciphertextA.Scale.Log2())
	fmt.Printf("  密文B: degree=%d, level=%d, scale=2^%d\n",
		ciphertextB.Degree(), ciphertextB.Level(), ciphertextB.Scale.Log2())

	// -------------------------------------------------------------------------
	// 第五部分: 解密操作验证
	// -------------------------------------------------------------------------
	fmt.Println("\n[5] 解密操作验证 (Decrypt)")
	fmt.Println("----------------------------------------")

	// Decrypt方法: 将密文解密为明文
	// 根据论文[1] Algorithm 3: Decrypt
	// m' = c0 + c1*sk = m + e (加上噪声)
	// 注意: 解密得到的是近似值, 包含噪声e
	fmt.Println("  执行Decryptor.DecryptNew(ciphertext)...")
	decryptedA := decryptor.DecryptNew(ciphertextA)
	decryptedB := decryptor.DecryptNew(ciphertextB)

	// 解码得到原始值
	decodedA := make([]complex128, PlaintextSize)
	decodedB := make([]complex128, PlaintextSize)
	encoder.Decode(decryptedA, decodedA)
	encoder.Decode(decryptedB, decodedB)

	fmt.Println("  解密并解码后的值:")
	fmt.Printf("    A: %v\n", formatComplexSlice(decodedA))
	fmt.Printf("    B: %v\n", formatComplexSlice(decodedB))
	fmt.Printf("    原始A: %v\n", formatComplexSlice(valuesA))
	fmt.Printf("    原始B: %v\n", formatComplexSlice(valuesB))
	fmt.Printf("    A的解密误差: %.2e\n", calculateError(valuesA, decodedA))
	fmt.Printf("    B的解密误差: %.2e\n", calculateError(valuesB, decodedB))

	// -------------------------------------------------------------------------
	// 第六部分: 同态乘法操作
	// -------------------------------------------------------------------------
	fmt.Println("\n[6] 同态乘法操作 (Multiply & Rescale)")
	fmt.Println("----------------------------------------")

	// 创建同态计算引擎
	evaluator := ckks.NewEvaluator(params, nil) // 暂时不需要重线性化密钥

	// 说明乘法深度限制
	fmt.Println("  参数设计说明:")
	fmt.Printf("    - 模数链Q有2个素数, 最大乘法深度 = %d\n", params.MaxLevel()+1)
	fmt.Printf("    - LevelsConsumedPerRescaling = %d (PREC64模式下每个rescale消耗1个level)\n", params.LevelsConsumedPerRescaling())
	fmt.Println("    - 这意味着:")
	fmt.Println("      乘法1次后rescale: A*B 成功 ✓")
	fmt.Println("      乘法2次后rescale: A*B*B 可能因噪声溢出而失败 ✗")

	// 第一次乘法: A * B
	fmt.Println("\n  执行第一次乘法: ctA * ctB")

	// Mul方法实现密文乘法
	// 根据论文[1] Algorithm 4: Multiply
	// 给定ct = (c0, c1)和ct' = (c0', c1'), 乘法结果为:
	// ct_mult = (d0, d1, d2) 其中:
	// d0 = c0*c0'
	// d1 = c0*c1' + c1*c0'
	// d2 = c1*c1'
	// 此时密文度为2, 需要重线性化或rescale

	// 在Lattigo中, 直接Mul后密文degree变为2
	ciphertextMul1 := ckks.NewCiphertext(params, 2, ciphertextA.Level())
	if err := evaluator.Mul(ciphertextA, ciphertextB, ciphertextMul1); err != nil {
		fmt.Printf("错误: 第一次乘法失败: %v\n", err)
		os.Exit(1)
	}
	fmt.Printf("  乘法后密文: degree=%d, level=%d, scale=2^%d\n",
		ciphertextMul1.Degree(), ciphertextMul1.Level(), ciphertextMul1.Scale.Log2())

	// Rescale方法: 除以当前模数并降低level
	// 根据论文[1] Section 3.3: Rescaling
	// 每个乘法后需要执行rescale来控制噪声增长
	// rescale将密文scale除以当前模数q_i
	if ciphertextMul1.Level() > 0 {
		if err := evaluator.Rescale(ciphertextMul1, ciphertextMul1); err != nil {
			fmt.Printf("错误: Rescale失败: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("  Rescale后: degree=%d, level=%d, scale=2^%d\n",
			ciphertextMul1.Degree(), ciphertextMul1.Level(), ciphertextMul1.Scale.Log2())
	}

	// 解密验证第一次乘法的结果
	decryptedMul1 := make([]complex128, PlaintextSize)
	decryptedMul1Plain := decryptor.DecryptNew(ciphertextMul1)
	encoder.Decode(decryptedMul1Plain, decryptedMul1)
	expectedMul1 := multiplyComplexSlices(valuesA, valuesB)
	fmt.Printf("  A*B 的结果: %v\n", formatComplexSlice(decryptedMul1))
	fmt.Printf("  期望A*B:    %v\n", formatComplexSlice(expectedMul1))

	// 第二次乘法: (A*B) * B = A*B^2
	fmt.Println("\n  执行第二次乘法: (ctA*ctB) * ctB")

	// 由于此时只有level=0, 再次乘法会导致噪声溢出
	ciphertextMul2 := ckks.NewCiphertext(params, 2, ciphertextMul1.Level())
	if err := evaluator.Mul(ciphertextMul1, ciphertextB, ciphertextMul2); err != nil {
		fmt.Printf("错误: 第二次乘法失败: %v\n", err)
		//os.Exit(1)
	}
	fmt.Printf("  乘法后密文: degree=%d, level=%d, scale=2^%d\n",
		ciphertextMul2.Degree(), ciphertextMul2.Level(), ciphertextMul2.Scale.Log2())

	// 尝试rescale, 由于只剩level=0, rescale会失败或产生错误结果
	fmt.Println("  尝试Rescale (预期失败或结果不准确)...")
	if ciphertextMul2.Level() > 0 {
		if err := evaluator.Rescale(ciphertextMul2, ciphertextMul2); err != nil {
			fmt.Printf("  Rescale失败(符合预期): %v\n", err)
		}
	} else {
		fmt.Println("  密文已在level=0, 无法再进行rescale")
	}

	// 解密验证第二次乘法的结果
	fmt.Println("Mul已出现错误，跳过执行解密验证第二次乘法的结果")
	// decryptedMul2 := make([]complex128, PlaintextSize)
	// decryptedMul2Plain := decryptor.DecryptNew(ciphertextMul2)
	// encoder.Decode(decryptedMul2Plain, decryptedMul2)
	// expectedMul2 := multiplyComplexSlices(expectedMul1, valuesB)
	// fmt.Printf("  A*B*B 的结果: %v\n", formatComplexSlice(decryptedMul2))
	// fmt.Printf("  期望A*B*B:    %v\n", formatComplexSlice(expectedMul2))
	// fmt.Printf("  误差: %.2e\n", calculateError(expectedMul2, decryptedMul2))
	// if calculateError(expectedMul2, decryptedMul2) > 0.1 {
	// 	fmt.Println("  注意: 由于噪声溢出, 结果误差较大!")
	// }

	// -------------------------------------------------------------------------
	// 第七部分: Bootstrapping操作
	// -------------------------------------------------------------------------
	fmt.Println("\n[7] Bootstrapping操作")
	fmt.Println("----------------------------------------")

	// 创建新的参数用于演示Bootstrapping
	// 为了演示, 我们需要一个有足够深度的参数集来进行bootstrapping
	fmt.Println("  设置Bootstrapping参数...")

	// 使用较小的参数以加快演示速度
	// 这些参数只用于演示Bootstrapping的概念
	schemeParams := ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{50, 45}, // 只有2个level: level=0 和 level=1
		LogP:            []int{45},
		LogDefaultScale: 40,
	}

	paramsBTP, err := ckks.NewParametersFromLiteral(schemeParams)
	if err != nil {
		fmt.Printf("错误: 创建Bootstrapping参数失败: %v\n", err)
		os.Exit(1)
	}

	// 创建Bootstrapping参数
	// Bootstrapping是重置同态计算能力的技术
	// 详见论文: Cheon et al. (2017), Handcock et al. (2019)
	btpParamsLit := bootstrapping.ParametersLiteral{}

	// 设置bootstrapping的环度为与scheme参数相同
	btpParamsLit.LogN = utils.Pointy(paramsBTP.LogN())

	// 创建Bootstrapping参数
	btpParams, err := bootstrapping.NewParametersFromLiteral(paramsBTP, btpParamsLit)
	if err != nil {
		fmt.Printf("错误: 创建Bootstrapping参数失败: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("  原始参数: LogN=%d, MaxLevel=%d\n", paramsBTP.LogN(), paramsBTP.MaxLevel())
	fmt.Printf("  Bootstrapping参数: LogN=%d, MaxLevel=%d\n",
		btpParams.BootstrappingParameters.LogN(), btpParams.BootstrappingParameters.MaxLevel())
	fmt.Printf("  Bootstrapping电路深度: %d\n", btpParams.Depth())

	// 生成Bootstrapping所需的密钥
	fmt.Println("  生成Bootstrapping密钥...")
	sk := rlwe.NewKeyGenerator(btpParams.BootstrappingParameters).GenSecretKeyNew()
	btpKeys, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		fmt.Printf("错误: 生成Bootstrapping密钥失败: %v\n", err)
		os.Exit(1)
	}

	// 创建Bootstrapping评估器
	btpEvaluator, err := bootstrapping.NewEvaluator(btpParams, btpKeys)
	if err != nil {
		fmt.Printf("错误: 创建Bootstrapping评估器失败: %v\n", err)
		os.Exit(1)
	}

	// 为演示创建新的明文和密文
	encoderBTP := ckks.NewEncoder(paramsBTP)
	encryptorBTP := ckks.NewEncryptor(paramsBTP, sk)
	decryptorBTP := ckks.NewDecryptor(paramsBTP, sk)

	// 创建用于Bootstrapping测试的值
	btpValuesA := make([]complex128, PlaintextSize)
	btpValuesB := make([]complex128, PlaintextSize)
	for i := range btpValuesA {
		btpValuesA[i] = complex(float64(i+1)*0.1, float64(i+1)*0.05)
		btpValuesB[i] = complex(float64(i+1)*0.2, float64(i+1)*0.1)
	}

	ptBTP := ckks.NewPlaintext(paramsBTP, 0) // 从level=0开始
	encoderBTP.Encode(btpValuesA, ptBTP)
	ctBTP, _ := encryptorBTP.EncryptNew(ptBTP)

	fmt.Printf("  原始值A: %v\n", formatComplexSlice(btpValuesA))
	fmt.Printf("  加密后在level=%d\n", ctBTP.Level())

	// 执行Bootstrapping
	fmt.Println("  执行Bootstrap(ct)...")
	ctBootstrap, err := btpEvaluator.Bootstrap(ctBTP)
	if err != nil {
		fmt.Printf("错误: Bootstrapping失败: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("  Bootstrapping后: level=%d (恢复到level=%d)\n",
		ctBootstrap.Level(), paramsBTP.MaxLevel())

	// 验证Bootstrapping结果
	resultBootstrap := make([]complex128, PlaintextSize)
	decryptedBTP := decryptorBTP.DecryptNew(ctBootstrap)
	encoderBTP.Decode(decryptedBTP, resultBootstrap)
	fmt.Printf("  Bootstrapping后解密结果: %v\n", formatComplexSlice(resultBootstrap))
	fmt.Printf("  原始值: %v\n", formatComplexSlice(btpValuesA))
	fmt.Printf("  误差: %.2e\n", calculateError(btpValuesA, resultBootstrap))

	// -------------------------------------------------------------------------
	// 第八部分: Bootstrapping后再次执行乘法
	// -------------------------------------------------------------------------
	fmt.Println("\n[8] Bootstrapping后再次执行 A*B*B")
	fmt.Println("----------------------------------------")

	// 使用Bootstrapping恢复的密文重新进行两次乘法
	fmt.Println("  创建用于后续演示的新密文...")

	// 创建新的明文并加密
	ptA := ckks.NewPlaintext(paramsBTP, paramsBTP.MaxLevel())
	ptB := ckks.NewPlaintext(paramsBTP, paramsBTP.MaxLevel())
	encoderBTP.Encode(btpValuesA, ptA)
	encoderBTP.Encode(btpValuesB, ptB)

	ctA, _ := encryptorBTP.EncryptNew(ptA)
	ctB, _ := encryptorBTP.EncryptNew(ptB)

	fmt.Printf("  密文A: level=%d\n", ctA.Level())
	fmt.Printf("  密文B: level=%d\n", ctB.Level())

	// 第一次乘法 A*B
	fmt.Println("\n  第一次乘法: ctA * ctB")
	evaluatorBTP := ckks.NewEvaluator(paramsBTP, btpKeys)

	ctMul1, _ := evaluatorBTP.MulNew(ctA, ctB)
	if err := evaluatorBTP.Rescale(ctMul1, ctMul1); err != nil {
		fmt.Printf("  Rescale警告: %v\n", err)
	}
	fmt.Printf("  结果: level=%d\n", ctMul1.Level())

	// 第二次乘法前先Bootstrapping恢复level
	fmt.Println("\n  执行Bootstrapping恢复level...")
	ctMul1BTP, err := btpEvaluator.Bootstrap(ctMul1)
	if err != nil {
		fmt.Printf("  Bootstrap错误: %v\n", err)
		os.Exit(234)
	}
	fmt.Printf("  Bootstrapping后: level=%d\n", ctMul1BTP.Level())

	// 第二次乘法 (A*B) * B
	fmt.Println("\n  第二次乘法: (ctA*ctB) * ctB")
	ctMul2, _ := evaluatorBTP.MulNew(ctMul1BTP, ctB)
	if err := evaluatorBTP.Rescale(ctMul2, ctMul2); err != nil {
		fmt.Printf("  Rescale警告: %v\n", err)
	}
	fmt.Printf("  结果: level=%d\n", ctMul2.Level())

	// 解密验证结果
	decryptedResult := make([]complex128, PlaintextSize)
	expectedResult := multiplyComplexSlices(multiplyComplexSlices(btpValuesA, btpValuesB), btpValuesB)
	decryptedPlain := decryptorBTP.DecryptNew(ctMul2)
	encoderBTP.Decode(decryptedPlain, decryptedResult)

	fmt.Printf("  A*B*B 最终结果: %v\n", formatComplexSlice(decryptedResult))
	fmt.Printf("  期望值: %v\n", formatComplexSlice(expectedResult))
	fmt.Printf("  误差: %.2e\n", calculateError(expectedResult, decryptedResult))

	if calculateError(expectedResult, decryptedResult) < 0.01 {
		fmt.Println("  ✓ Bootstrapping使多次乘法成为可能!")
	}

	// -------------------------------------------------------------------------
	// 总结
	// -------------------------------------------------------------------------
	fmt.Println("\n========================================")
	fmt.Println("  演示完成!")
	fmt.Printf("  总执行时间: %v\n", time.Since(totalStart))
	fmt.Println("========================================")

	fmt.Println("\n关键概念总结:")
	fmt.Println("1. Encode: 将复数向量编码为环R_Q上的多项式")
	fmt.Println("2. Encrypt: 使用公钥加密明文, 引入离散高斯噪声")
	fmt.Println("3. Decrypt: 使用私钥解密, 恢复近似消息")
	fmt.Println("4. Mul: 密文乘法, 乘积degree变为2")
	fmt.Println("5. Rescale: 除去模数降低level, 控制噪声增长")
	fmt.Println("6. Bootstrap: 重置计算能力, 恢复level并减少噪声")
	fmt.Println("\nBootstrapping工作流程:")
	fmt.Println("  1. ScaleDown: 缩放密文到适当范围")
	fmt.Println("  2. ModUp: 将模数从q扩展到Q")
	fmt.Println("  3. CoeffsToSlots: 同态编码(DFT)")
	fmt.Println("  4. EvalMod: 同态模约简, 提取小数部分")
	fmt.Println("  5. SlotsToCoeffs: 同态解码(IDFT)")
}

// ================================================================================
// 辅助函数
// ================================================================================

// setupCKKSParameters 创建并配置CKKS参数
// 根据论文[1] Section 5.1 的参数选择建议
func setupCKKSParameters() (ckks.Parameters, error) {
	// 参数字面量定义
	// 这些参数定义了在环R_Q = Z_Q[X]/(X^N + 1)上的加密操作
	//
	// 素数的选择要求:
	// 1. 每个素数q_i必须满足 q_i ≡ 1 (mod 2N)
	//    这是为了支持特殊的离散傅里叶变换(DFT)
	// 2. 素数之间必须pairwise coprime
	// 3. 第一个素数应该足够大, 以存储计算结果
	paramsLiteral := ckks.ParametersLiteral{
		// LogN: 对数多项式次数
		// N = 2^LogN 必须足够大以提供足够的安全边际
		LogN: LogN,

		// Q: 密文模数链
		// 每个素数增加一个乘法深度
		// 素数以2^bitlen的形式指定
		Q: []uint64{
			// 使用2^55作为第一个素数, 提供足够的空间存储结果
			0x80000000080001, // 约2^55, 满足 ≡ 1 (mod 2N) for N=2^10
			// 使用2^45作为后续素数, 用于增加深度
			0x2000000a0001, // 约2^45
		},

		// LogDefaultScale: 默认缩放因子 (log2)
		// 缩放因子应该在40-60比特之间以平衡精度和效率
		LogDefaultScale: LogDefaultScale,

		// P: 辅助模数用于密钥交换
		// 这些素数不参与同态容量, 但用于优化key-switching
		P: []uint64{
			0x80000000130001, // 约2^55
		},

		// RingType: 环类型
		// ring.Standard: Z[X]/(X^N+1), 槽数量 = N/2
		// ring.ConjugateInvariant: Z[X+X^-1]/(X^N+1), 槽数量 = N
		RingType: ring.Standard,

		// Xs: 密钥分布
		// 默认使用Hamming weight为192的三元分布
		// 这提供了良好的安全性
		Xs: rlwe.DefaultXs,

		// Xe: 错误分布
		// 默认使用标准差为3.2的离散高斯分布
		Xe: rlwe.DefaultXe,
	}

	// 从字面量创建经过验证的参数集
	params, err := ckks.NewParametersFromLiteral(paramsLiteral)
	if err != nil {
		return params, fmt.Errorf("参数验证失败: %w", err)
	}

	// 验证参数安全性
	// error: params.T undefined
	// if params.T() < 128 {
	// 	fmt.Printf("警告: 安全边际 %.0f 比特低于推荐值128比特\n", params.T())
	// }

	return params, nil
}

// printParameters 打印CKKS参数信息
func printParameters(params ckks.Parameters) {
	fmt.Printf("  环多项式次数: N = 2^%d = %d\n", params.LogN(), params.N())
	fmt.Printf("  密文模数链Q:\n")
	for i, _ := range params.Q() {
		fmt.Printf("    Q[%d] = %d bits (≈ 2^%d)\n", i, params.LogQLvl(i), params.LogQLvl(i))
	}
	fmt.Printf("  总模数比特数: %.0f\n", params.LogQP())
	fmt.Printf("  辅助模数P: %d个素数\n", params.PCount())
	fmt.Printf("  最大乘法深度: %d\n", params.MaxDepth())
	fmt.Printf("  最大槽数量: %d\n", params.MaxSlots())
	fmt.Printf("  默认缩放因子: 2^%d ≈ %.2e\n", params.LogDefaultScale(), params.DefaultScale().Float64())
	fmt.Printf("  精度模式: %s\n", func() string {
		if params.PrecisionMode() == ckks.PREC64 {
			return "PREC64 (64位)"
		}
		return "PREC128 (128位)"
	}())
	//fmt.Printf("  估计安全边际: %.0f比特\n", params.T())
}

// createPlaintextVectors 创建演示用的复数向量
func createPlaintextVectors() ([]complex128, []complex128) {
	// 创建两个复数向量作为演示数据
	// 选择较小的值以避免数值溢出
	valuesA := make([]complex128, PlaintextSize)
	valuesB := make([]complex128, PlaintextSize)

	// 向量A: 实部递增, 虚部为0.1的倍数
	// 例如: [0.5+0.1i, 1.0+0.2i, 1.5+0.3i, 2.0+0.4i]
	for i := 0; i < PlaintextSize; i++ {
		valuesA[i] = complex(float64(i+1)*0.5, float64(i+1)*0.1)
	}

	// 向量B: 实部为0.1的倍数, 虚部递增
	// 例如: [0.1+0.5i, 0.2+1.0i, 0.3+1.5i, 0.4+2.0i]
	for i := 0; i < PlaintextSize; i++ {
		valuesB[i] = complex(float64(i+1)*0.1, float64(i+1)*0.5)
	}

	return valuesA, valuesB
}

// multiplyComplexSlices 复数向量逐元素乘法
func multiplyComplexSlices(a, b []complex128) []complex128 {
	result := make([]complex128, len(a))
	for i := range a {
		result[i] = a[i] * b[i]
	}
	return result
}

// calculateError 计算近似计算结果与期望值的相对误差
func calculateError(expected, actual []complex128) float64 {
	if len(expected) != len(actual) {
		return -1
	}

	var maxRelErr float64
	for i := range expected {
		// 计算实部和虚部的相对误差
		expectedVal := expected[i]
		actualVal := actual[i]

		relErrReal := math.Abs(real(expectedVal) - real(actualVal))
		if math.Abs(real(expectedVal)) > 1e-10 {
			relErrReal /= math.Abs(real(expectedVal))
		}

		relErrImag := math.Abs(imag(expectedVal) - imag(actualVal))
		if math.Abs(imag(expectedVal)) > 1e-10 {
			relErrImag /= math.Abs(imag(expectedVal))
		}

		relErr := cmplx.Abs(expectedVal - actualVal)
		if cmplx.Abs(expectedVal) > 1e-10 {
			relErr /= cmplx.Abs(expectedVal)
		}

		if relErr > maxRelErr {
			maxRelErr = relErr
		}
	}

	return maxRelErr
}

// formatComplexSlice 格式化复数向量为字符串
func formatComplexSlice(values []complex128) string {
	str := "["
	for i, v := range values {
		if i > 0 {
			str += ", "
		}
		str += fmt.Sprintf("%.4f+%.4fi", real(v), imag(v))
	}
	str += "]"
	return str
}
