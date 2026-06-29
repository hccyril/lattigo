_按以下格式规则添加_
```
# 日期 标题 
## 问题
...
## 原因分析
...
## 修改过程
...
```

---

# 2026-06-28 BCKS24 论文核心实现完成

## 问题

实现论文 "Bootstrapping Bits with CKKS" (BCKS24, 2024-767) 中的 BinBoot（Algorithm 2）和 GateBoot（Algorithm 3）两个核心算法，基于 Lattigo v6 框架。

## 原因分析

论文的核心创新是用三角函数 `f(x)=(1-cos(2πx))/2` 替换标准 CKKS 自举中的 sine/cosine EvalMod 逼近，使 `q0 = 2Δ0`（BinBoot）或 `q0 = 3Δ0`（GateBoot），大幅节省模数预算。需要在 Lattigo v6 现有 bootstrapping 框架基础上，替换 EvalMod 阶段的多项式逼近。

## 修改过程

### 1. 创建 `circuits/ckks/binboot/parameters.go`
- 定义 `Param14BinBootLiteral` 和 `Param14GateBootLiteral` 参数集（论文 §5 Table 5）
- BinBoot: LogN=14, LogDefaultScale=31 (q0≈2Δ0), ConjugateInvariant
- GateBoot: LogN=14, LogDefaultScale=30 (q0≈3Δ0), ConjugateInvariant
- 定义 `GateType` 枚举（GateAND/OR/XOR/NAND/NOR/XNOR）及 String() 方法

### 2. 创建 `circuits/ckks/binboot/polynomials.go`
- 实现 `fBinBoot(x) = (1-cos(2πx))/2`（论文 §3.1）
- 实现六个门三角函数 `fGateNAND/AND/OR/XOR/NOR/XNOR`（论文 §4.1 Table 4）
- `NewBinBootPoly(degree, K)`: 生成 f_BinBoot 的 Chebyshev 逼近，标记偶函数并清零奇数次系数
- `NewGatePoly(gate, degree, K)`: 生成门函数的 Chebyshev 逼近
- `ValidateGateFunction(g)`: 明文验证门函数在 0, 1/3, 2/3 处与真值表一致
- 修复: 移除 `math/cmplx` 无用导入，修正 `bigFloat` 类型别名为 `*big.Float`

### 3. 创建 `circuits/ckks/binboot/evaluator.go`
- `Evaluator` 结构体包装 `bootstrapping.Evaluator`，预计算 BinBoot 和六个 GateBoot 多项式
- `NewEvaluator(btpEval, degree, K)`: 从标准自举评估器获取 polynomial.Evaluator，预计算所有多项式
- `Bootstrap(ct)`: 论文 Algorithm 2 — ScaleDown → ModUp → CoeffsToSlots → evalBinBootPoly → SlotsToCoeffs
- `GateBoot(ct1, ct2, gate)`: 论文 Algorithm 3 — Add → ScaleDown → ModUp → CoeffsToSlots → evalGatePoly → SlotsToCoeffs
- `evalBinBootPoly/evalGatePoly`: Chebyshev 变量代换 + polynomial.Evaluate
- 便捷方法: GateBootNAND/AND/OR/XOR/NOR/XNOR

### 4. 创建 `examples/singleparty/paper_bcks/main.go`
- 演示 GateBoot NAND 流程: 参数初始化 → 密钥生成 → 加密二进制数组 → GateBootNAND → 解密验证
- `GenBinaryVector(n)`: 随机二进制数组生成
- `ValidateBinaryTolerance(actual, expected)`: 容差 0.25 验证

### 5. 编译验证
- `go build ./circuits/ckks/binboot/` — 通过
- `go build ./examples/singleparty/paper_bcks/` — 通过
- `go vet ./circuits/ckks/binboot/ ./examples/singleparty/paper_bcks/` — 通过

### 修改文件清单
| 文件 | 操作 | 说明 |
|------|------|------|
| `circuits/ckks/binboot/parameters.go` | 新建 | Param14 参数集 + GateType 枚举 |
| `circuits/ckks/binboot/polynomials.go` | 新建 | 三角函数 + Chebyshev 逼近生成器 |
| `circuits/ckks/binboot/evaluator.go` | 新建 | BinBoot/GateBoot 评估器实现 |
| `examples/singleparty/paper_bcks/main.go` | 新建 | GateBoot NAND 示例入口 |
| `docs/todo.md` | 更新 | 标记所有阶段完成 |
| `docs/dev-logs.md` | 更新 | 本条开发日志 |
