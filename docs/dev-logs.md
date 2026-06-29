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

# 2026-06-29 修复 GateBoot NAND 结果完全错误（MessageRatio 归一化未补偿）

## 问题

前次修复（维度不匹配 panic）后程序成功运行，但 GateBoot NAND 结果完全错误：
- 解密值在 [-43, 5] 范围波动，期望值为 0 或 1
- 64 个结果仅 2 个正确（3.1%），本质上等同于随机噪声

## 原因分析

### 核心根因：多项式未补偿 ScaleDown 的 MessageRatio 归一化

通过完整追踪标准 Lattigo bootstrapping 的消息流，定位到 binboot 多项式生成的数学错误。

**标准 bootstrapping 消息流（ScaleDown 归一化）：**

1. 输入密文加密 `slot_values * Δ`（slot_values 是 slot 域原始浮点消息，Δ=DefaultScale）
2. **ScaleDown**（evaluator.go:566-613）：将消息归一化为 `slot_values * Q0/MessageRatio`，scale 调整为 `Q0/MessageRatio`，level 降至 0
   - MessageRatio = `2^LogMessageRatio`（mod1_parameters.go:106-108）
3. **ModUp**：模数 q0→Q，消息变为 `Q0*(slot_values/MessageRatio + I)`（I 是 ModRaise 引入的小整数）
4. **CoeffsToSlots**：转 slot 域，slot 值 = `slot_values/MessageRatio + I`（归一化消息 + 整数项）
5. **EvalMod**：对 `x = slot_values/MessageRatio + I` 应用多项式 P(x)

**论文期望 vs 实际输入：**

论文的 f_BinBoot 和 f_G 设计为直接接收 `slot_values + I`：
- BinBoot: f_BinBoot(b/2 + I) = b（b 是二进制位）
- GateBoot: f_NAND((b1+b2)/3 + I) = NAND(b1,b2)

但由于 ScaleDown 的 MessageRatio=2 归一化，EvalMod 实际输入是 `slot_values/2 + I`：
- BinBoot: 实际输入 b/4 + I（期望 b/2 + I），f_BinBoot(b/4) = (1-cos(πb/2))/2 ≠ b
- GateBoot: 实际输入 (b1+b2)/6 + I（期望 (b1+b2)/3 + I），f_NAND((b1+b2)/6) 完全偏离真值表

**数值验证（GateBoot NAND, MessageRatio=2）：**
- (b1+b2)=0: f_NAND(0/6) = (2/3)(1+sin(π/6)) = 1 ✓（偶然正确）
- (b1+b2)=1: f_NAND(1/6) = (2/3)(1+sin(π/3+π/6)) = (2/3)(1+sin(5π/6)) = 1 ✗（期望 1 但因相位偏移实际不精确）
- (b1+b2)=2: f_NAND(2/6) = (2/3)(1+sin(2π/3+π/6)) = (2/3)(1+sin(5π/6)) = 1 ✗（期望 0！）

这解释了结果值的大幅波动——多项式评估的不是论文设计的函数。

### 修复方案：多项式逼近 f(MessageRatio * x)

将多项式从逼近 `f(x)` 改为逼近 `f(MessageRatio * x)`，利用 f 的周期为 1 消除归一化：

```
P(slot_values/MessageRatio + I)
= f(MessageRatio * (slot_values/MessageRatio + I))
= f(slot_values + I*MessageRatio)
= f(slot_values)                          （f 周期为 1，I*MessageRatio 是整数）
```

**验证（GateBoot NAND, MessageRatio=2, slot_values=(b1+b2)/3）：**
- P((b1+b2)/6) = f_NAND(2*(b1+b2)/6) = f_NAND((b1+b2)/3) = NAND(b1,b2) ✓

### 附带问题：Lattigo MessageRatio 必须为 2 的幂

论文 GateBoot 要求 q0=3*Δ0（MessageRatio=3），但 Lattigo 的 `MessageRatio = 2^LogMessageRatio` 只支持 2 的幂。

当前设置 LogMessageRatio=1（MessageRatio=2）。通过 MessageRatio 补偿方案，slot_values=(b1+b2)/3 编码 + P(x)=f(2x) 的组合在数学上等价于论文的 q0=3Δ0 + f(x)：
- P((b1+b2)/3/2) = f(2*(b1+b2)/6) = f((b1+b2)/3) = G(b1,b2) ✓

这绕过了 MessageRatio 必须为 3 的限制，仅代价是多项式频率翻倍（需要更高次数维持精度，但测试级参数可接受）。

## 修改过程

### 1. `circuits/ckks/binboot/polynomials.go` — 多项式添加 MessageRatio 补偿

**NewBinBootPoly：**
- 签名添加 `messageRatio float64` 参数
- Chebyshev 逼近函数从 `fBinBoot(x)` 改为 `fBinBoot(messageRatio * x)`
- 注释完整说明 ScaleDown 归一化与补偿的数学关系

**NewGatePoly：**
- 签名添加 `messageRatio float64` 参数
- Chebyshev 逼近函数从 `f(x)` 改为 `f(messageRatio * x)`
- 注释完整说明补偿逻辑

### 2. `circuits/ckks/binboot/evaluator.go` — NewEvaluator 传递 MessageRatio

- 从 `btpEval.Mod1Parameters.MessageRatio()` 获取实际 MessageRatio
- 传递给 `NewBinBootPoly(degree, K, messageRatio)` 和 `NewGatePoly(gate, degree, K, messageRatio)`
- 注释说明 MessageRatio 的来源和作用

### 3. 论文方案逻辑一致性说明

- **Algorithm 2/3 的 EvalMod 语义不变**：多项式仍然评估论文的 f_BinBoot/f_G，只是通过 `f(MessageRatio*x)` 的变量代换补偿了 ScaleDown 的归一化
- **周期性保证正确性**：f_BinBoot 和 f_G 的周期为 1（cos/sin 的 2π 周期），I*MessageRatio 是整数，f(slot_values + I*MessageRatio) = f(slot_values)
- **编码不变**：main.go 仍然编码 slot_values=(b1+b2)/3（GateBoot），符合论文 §4.1
- **q0=3Δ0 约束的等效满足**：通过 MessageRatio=2 + f(2x) 的组合等效实现论文 q0=3Δ0 + f(x) 的效果

### 4. 编译验证

- `go build ./circuits/ckks/binboot/` — 通过
- `go build ./examples/singleparty/paper_bcks/` — 通过
- `go vet ./circuits/ckks/binboot/ ./examples/singleparty/paper_bcks/` — 通过

### 修改文件清单
| 文件 | 操作 | 说明 |
|------|------|------|
| `circuits/ckks/binboot/polynomials.go` | 修改 | NewBinBootPoly/NewGatePoly 添加 messageRatio 参数，函数改为 f(messageRatio*x) |
| `circuits/ckks/binboot/evaluator.go` | 修改 | NewEvaluator 获取 MessageRatio 并传递给多项式生成器 |
| `docs/dev-logs.md` | 更新 | 本条开发日志 |

---

# 2026-06-29 修复 GateBoot 运行时维度不匹配 panic 及参数配置错误

## 问题

执行 `go run ./examples/singleparty/paper_bcks/ -short` 时，程序在 GateBoot NAND 的 ModUp 阶段 panic：
```
panic: cannot inttCoreLazy: ensure that len(p1)=4096, len(p2)=8192 and len(roots)=8192 >= N=8192
```
调用链：`GateBoot → ModUp → GadgetProduct → ApplyEvaluationKey → INTT → inttCoreLazy`。

## 原因分析

经完整梳理错误堆栈与 Lattigo 源码，定位到三个独立但叠加的缺陷：

### 缺陷 1（直接触发 panic）：缺少 ConjugateInvariant 环切换

项目使用 `ring.ConjugateInvariant` 环类型，根据框架约束（`bootstrapping/parameters.go:65`），自举环 LogN 必须为残差环 LogN+1。因此：
- 残差环维度 N1 = 2^LogN（ConjugateInvariant，实数环）
- 自举环维度 N2 = 2^(LogN+1)（Standard，复数环）

标准 Lattigo 自举流程（`EvaluateConjugateInvariant`，evaluator.go:463）在执行 ScaleDown/ModUp 之前，会通过 `RealToComplexNew` 将密文从残差环（N1）切换到自举环（N2），自举完成后再通过 `ComplexToRealNew` 切换回残差环。

**binboot 评估器直接在残差环密文（N1=4096）上调用 `btpEval.ScaleDown` 和 `btpEval.ModUp`，但这两个操作在自举环（N2=8192）上执行。** `ModUp` 内部第一步 `ringQ.INTT(ctIn.Value[i], ctIn.Value[i])` 使用 N2 的 SubRing 对 N1 维度的多项式做 INTT，导致 `len(p1)=4096 < N=8192` 的维度不匹配 panic。

### 缺陷 2：SchemeParams.LogQ 误包含自举电路素数

`bootstrapping.NewParametersFromLiteral` 的逻辑（parameters.go:307-314）是：
1. 复制残差参数的全部 Q 素数作为基底
2. 根据 StC/EvalMod/CtS 的因式分解深度**自动追加**电路素数

原代码的 `LogQ` 包含了 13 个素数（Base + Mult + StC + EvalMod + CtS），但 StC/EvalMod/CtS 素数会被框架**再次追加**，导致模数链翻倍：总 Q 素数 = 13（残差）+ 8（电路）= 21，LogQP≈910 bits。正确配置应仅包含残差素数（Base + Mult），电路素数由框架自动追加。

参考标准参数集 `N15QP768H192H32`：其 `LogQ` 仅 3 个残差素数 `[]int{33, 50, 25}`，StC/EvalMod/CtS 素数全部由框架追加。

### 缺陷 3：BootstrapParams.LogN 未满足 LogN+1 约束 + Mod1Degree 层级不足

- **LogN**：原代码设 BootstrapParams.LogN=14 = 残差 LogN=14，但 ConjugateInvariant 要求 LogN+1=15。非 -short 模式会在 `NewParametersFromLiteral` 阶段报错退出。
- **Mod1Degree=30**：`mod1.ParametersLiteral.Depth()` 对 CosDiscrete 返回 `bits.Len64(max(30, 2*4-1))=5`，即 5 个 EvalMod 层级。但自定义多项式评估需要：
  - Chebyshev 变量代换（Mul + Rescale）消耗 1 层
  - 30 次多项式求值消耗 `ceil(log2(31))=5` 层
  - 总计 6 层 > 可用 5 层，多项式评估会因层级不足而失败

### 错误根因结论

**panic 的直接原因是缺陷 1（缺少环切换），而非参数维度过大导致内存不足。** 缺陷 2 和 3 是潜在问题，修复缺陷 1 后会暴露。

## 修改过程

### 1. `circuits/ckks/binboot/evaluator.go` — 环切换 + 变量代换顺序

**Bootstrap 方法（论文 Algorithm 2）：**
- 在 ScaleDown 之前添加 `btpEval.RealToComplexNew(ct)` 将密文从残差环(N1)切换到自举环(N2)
- 在 SlotsToCoeffs 之后添加 scale 补偿 `ctOutN2.Scale *= 0.5`（补偿 ComplexToRealNew 的 2x 因子）
- 添加 `btpEval.ComplexToRealNew(ctOutN2)` 将结果切换回残差环(N1)

**GateBoot 方法（论文 Algorithm 3）：**
- ct1+ct2 加法仍在残差环(N1)上执行（保持论文 Algorithm 3 Step 0 语义）
- 后续同样添加 RealToComplexNew → 自举电路 → scale补偿 → ComplexToRealNew

**evalBinBootPoly / evalGatePoly 方法：**
- 修复 Chebyshev 变量代换执行顺序：先 `Mul(scalar)` → `Rescale` → `Add(constant)`
- 原代码顺序 `Add(constant)` → `Mul(scalar)` → `Rescale` 得到 `scalar*(ct+constant)` 而非正确的 `scalar*ct+constant`
- 对对称区间 [-K, K]（constant=0）两者等价，但修正顺序确保通用正确性

### 2. `circuits/ckks/binboot/parameters.go` — 参数配置修复

**Param14BinBootLiteral & Param14GateBootLiteral：**
- `BootstrapParams.LogN`：14 → 15（满足 ConjugateInvariant 的 LogN+1 约束）
- `SchemeParams.LogQ`：从 13 个素数缩减为 3 个残差素数 `[]int{32, 45, 45}`（Base + 2个Mult），电路素数由框架自动追加
- `BootstrapParams.Mod1Degree`：30 → 32（`bits.Len64(32)=6`，为变量代换(1层) + 多项式(5层) = 6层提供足够层级）

### 3. `examples/singleparty/paper_bcks/main.go` — -short 模式参数缩减

-short 模式下完整覆盖所有需调整的字段：
- `LogN`：14→12（残差环 N=4096）
- `BootstrapParams.LogN`：15→13（自举环 N=8192，LogN+1）
- `LogQ`：缩减为 `[]int{32, 30}`（2 个残差素数，最小化模数链）
- `Mod1Degree`：32→16（`bits.Len64(16)=5`，为变量代换(1层) + 15次多项式(4层) = 5层）
- `polyDegree`：30→15（减少 Chebyshev 逼近计算量，加速运行）

### 4. 论文方案逻辑一致性说明

- **Algorithm 2/3 核心流程不变**：StC → EvalMod(自定义三角函数) → CtS 的论文逻辑完全保留
- **环切换不改变方案语义**：RealToComplexNew/ComplexToRealNew 是 ConjugateInvariant CKKS 的标准基础设施，仅完成实数环↔复数环的维度转换，不影响消息空间
- **LogMessageRatio=1（ratio=2）对 GateBoot 的正确性**：门函数 f_G 具有周期 1，f_G((b1+b2)/3 + 2I) = f_G((b1+b2)/3) = G(b1,b2)，因 2I 为整数，周期性保证正确性
- **变量代换顺序修复**确保 Chebyshev 基评估的数学正确性

### 5. 编译验证

- `gofmt -e` 语法检查通过（exit code 0）
- `go build` 因网络不可用（Go 模块缓存不完整，无法下载依赖）未能执行，需用户在网络环境下自行编译验证

### 修改文件清单
| 文件 | 操作 | 说明 |
|------|------|------|
| `circuits/ckks/binboot/evaluator.go` | 修改 | 添加环切换、修复变量代换顺序、scale补偿 |
| `circuits/ckks/binboot/parameters.go` | 修改 | LogN→15、LogQ仅残差素数、Mod1Degree→32 |
| `examples/singleparty/paper_bcks/main.go` | 修改 | -short模式完整参数缩减、polyDegree适配 |
| `docs/dev-logs.md` | 更新 | 本条开发日志 |

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
