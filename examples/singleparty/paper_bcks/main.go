// Package main 实现论文 "Bootstrapping Bits with CKKS" (BCKS24) 的示例入口。
//
// 本程序演示 GateBoot 算法（论文 Algorithm 3）：
//   1. 初始化 Param14 参数集
//   2. 生成密钥和自举密钥
//   3. 加密两个随机二进制数组
//   4. 对每对密文执行 GateBootNAND（自举 + NAND 门）
//   5. 解密并验证结果与明文 NAND 一致
//
// 论文对应：§4.1 Algorithm 3、§4.2 Theorem 2、§5 Table 5 Param14
package main

import (
	"flag"
	"fmt"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/binboot"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

// flagShort 控制是否使用缩小参数以加快运行速度（不安全）。
var flagShort = flag.Bool("short", false, "run with smaller insecure parameters for speed")

func main() {
	flag.Parse()

	fmt.Println("=== BCKS24 GateBoot 示例 ===")
	fmt.Println("论文: Bootstrapping Bits with CKKS (2024-767)")
	fmt.Println()

	// ==============================
	// === 1) 参数初始化 ===
	// ==============================

	// 使用论文 §5 Table 5 的 Param14 参数集（GateBoot 模式）
	// 论文 §4.1: q0 = 3*Δ0
	lit := binboot.Param14GateBootLiteral

	// 如果 -short 标志，缩小 N 以加快速度
	if *flagShort {
		lit.SchemeParams.LogN = 12
		lit.BootstrapParams.LogN = utils.Pointy(13) // ConjugateInvariant 需要 LogN+1
	}

	fmt.Println("初始化残差参数 (ResidualParameters)...")

	// 从字面量创建残差 CKKS 参数
	params, err := ckks.NewParametersFromLiteral(lit.SchemeParams)
	if err != nil {
		panic(fmt.Errorf("cannot create residual parameters: %w", err))
	}

	fmt.Printf("  LogN=%d, LogSlots=%d, H=%d, LogQP=%f, levels=%d, LogScale=%d\n",
		params.LogN(),
		params.LogMaxSlots(),
		params.XsHammingWeight(),
		params.LogQP(),
		params.MaxLevel(),
		params.LogDefaultScale(),
	)

	// ==============================
	// === 2) 自举参数 ===
	// ==============================

	fmt.Println("初始化自举参数 (BootstrappingParameters)...")

	btpParams, err := bootstrapping.NewParametersFromLiteral(params, lit.BootstrapParams)
	if err != nil {
		panic(fmt.Errorf("cannot create bootstrapping parameters: %w", err))
	}

	fmt.Printf("  LogN=%d, LogQP=%f, levels=%d, EphemeralSecretWeight=%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.EphemeralSecretWeight,
	)

	// ==============================
	// === 3) 密钥生成 ===
	// ==============================

	fmt.Println("生成密钥...")

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println("生成自举评估密钥 (可能需要较长时间)...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(fmt.Errorf("cannot generate evaluation keys: %w", err))
	}
	fmt.Println("  完成")

	// ==============================
	// === 4) 创建评估器 ===
	// ==============================

	fmt.Println("创建标准自举评估器...")
	btpEval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		panic(fmt.Errorf("cannot create bootstrapping evaluator: %w", err))
	}

	fmt.Println("创建 BinBoot/GateBoot 评估器...")
	// 论文 §3.1: degree=30, K=4
	binbootEval, err := binboot.NewEvaluator(btpEval, 30, 4)
	if err != nil {
		panic(fmt.Errorf("cannot create binboot evaluator: %w", err))
	}

	// ==============================
	// === 5) 加密二进制数据 ===
	// ==============================

	slots := params.MaxSlots()
	if slots > 64 {
		slots = 64 // 限制测试规模
	}

	// 生成两个随机二进制数组
	b1 := GenBinaryVector(slots)
	b2 := GenBinaryVector(slots)

	// 计算明文 NAND 结果用于验证
	expected := make([]float64, slots)
	for i := 0; i < slots; i++ {
		expected[i] = float64(1 - b1[i]*b2[i]) // NAND = NOT(AND)
	}

	fmt.Printf("加密 %d 个二进制位...\n", slots)

	// 将二进制值编码为 CKKS 明文
	// GateBoot 模式下，消息编码为 b/3（因为 q0=3*Δ0）
	// 论文 §4.1: 输入为 (b1+b2)/3 + I
	values1 := make([]complex128, slots)
	values2 := make([]complex128, slots)
	for i := 0; i < slots; i++ {
		values1[i] = complex(float64(b1[i])/3.0, 0)
		values2[i] = complex(float64(b2[i])/3.0, 0)
	}

	pt1 := ckks.NewPlaintext(params, params.MaxLevel())
	if err := encoder.Encode(values1, pt1); err != nil {
		panic(err)
	}
	pt2 := ckks.NewPlaintext(params, params.MaxLevel())
	if err := encoder.Encode(values2, pt2); err != nil {
		panic(err)
	}

	ct1, err := encryptor.EncryptNew(pt1)
	if err != nil {
		panic(err)
	}
	ct2, err := encryptor.EncryptNew(pt2)
	if err != nil {
		panic(err)
	}

	// 将密文降至 level 0（自举要求）
	for ct1.Level() > 0 {
		btpEval.DropLevel(ct1, 1)
	}
	for ct2.Level() > 0 {
		btpEval.DropLevel(ct2, 1)
	}

	// ==============================
	// === 6) GateBoot NAND ===
	// ==============================

	fmt.Println("执行 GateBoot NAND（自举 + NAND 门）...")

	ctResult, err := binbootEval.GateBootNAND(ct1, ct2)
	if err != nil {
		panic(fmt.Errorf("GateBootNAND failed: %w", err))
	}

	fmt.Println("  完成")

	// ==============================
	// === 7) 解密验证 ===
	// ==============================

	fmt.Println("解密并验证结果...")

	ptResult := decryptor.DecryptNew(ctResult)
	resultValues := make([]complex128, slots)
	if err := encoder.Decode(ptResult, resultValues); err != nil {
		panic(err)
	}

	// 验证结果
	correct, total := 0, 0
	for i := 0; i < slots; i++ {
		val := real(resultValues[i])
		exp := expected[i]

		if ValidateBinaryTolerance(val, exp) {
			correct++
		}
		total++

		if i < 8 { // 打印前 8 个结果
			fmt.Printf("  [%d] b1=%d, b2=%d, NAND=%d, result=%.6f, expected=%.1f %s\n",
				i, b1[i], b2[i], int(exp), val, exp,
				statusStr(ValidateBinaryTolerance(val, exp)))
		}
	}

	fmt.Printf("\n结果: %d/%d 正确 (%.1f%%)\n", correct, total, float64(correct)/float64(total)*100)

	if correct == total {
		fmt.Println("✓ GateBoot NAND 验证通过！")
	} else {
		fmt.Println("✗ GateBoot NAND 验证失败（部分结果超出容差）")
	}
}

// statusStr 返回验证状态的字符串表示。
func statusStr(ok bool) string {
	if ok {
		return "✓"
	}
	return "✗"
}

// GenBinaryVector 生成一个长度为 n 的随机二进制向量（元素为 0 或 1）。
func GenBinaryVector(n int) []int {
	v := make([]int, n)
	for i := range v {
		v[i] = rand.Intn(2)
	}
	return v
}

// ValidateBinaryTolerance 验证解密值是否在二进制容差范围内。
//
// 论文 §4.1: GateBoot 输出应为 G(b1,b2) ∈ {0, 1}。
// 容差阈值 0.25: "1" > 0.75, "0" < 0.25。
func ValidateBinaryTolerance(actual, expected float64) bool {
	tol := 0.25
	return absFloat(actual-expected) < tol
}

func absFloat(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
