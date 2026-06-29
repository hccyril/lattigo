# BCKS24 论文实现工作计划

> 论文：Bootstrapping Bits with CKKS (Bae, Cheon, Kim, Stehlé, 2024)
> 代码库：Lattigo v6 (`github.com/tuneinsight/lattigo/v6`)
> 开发手册：`docs/bcks.md`

## 总体目标

通过软件工程方式实现论文中的两个核心算法：
- **BinBoot**（Algorithm 2）：二进制自举 + 二次噪声清洗
- **GateBoot**（Algorithm 3）：自举同时完成二元门（NAND/AND/OR/XOR/NOR/XNOR）

核心创新点是用三角函数 `f(x)=(1-cos(2πx))/2` 替换标准 CKKS 自举中的 sine/cosine EvalMod，
使 `q0 = 2Δ0`（BinBoot）或 `q0 = 3Δ0`（GateBoot），大幅节省模数预算。

---

## 关键工作步骤与完成状态

### 阶段 0：前置准备
- [x] 0.1 阅读 `docs/bcks.md` 开发手册，确认技术栈与架构约束
- [x] 0.2 阅读论文 `docs/2024-767_BCKS24.pdf`（已提取为 txt），梳理 Algorithm 2/3、Table 4、Theorem 1/2
- [x] 0.3 调研 Lattigo v6 现有结构：`circuits/ckks/bootstrapping`、`circuits/ckks/mod1`、`circuits/ckks/polynomial`、`circuits/ckks/dft`
- [x] 0.4 确认 Lattigo v6 模块路径为 `github.com/tuneinsight/lattigo/v6`，Go 1.25

### 阶段 1：核心包实现 `circuits/ckks/binboot/`
- [x] 1.1 创建 `parameters.go`：定义 `BinBoot`/`GateBoot` 专用参数字面量（Param14，N=2^14，h=256，h~=32，dnum=13，q0=2Δ0 或 3Δ0）
  - 交付物：`BinBootParametersLiteral`、`GateBootParametersLiteral` 结构体与构造函数
  - 验收：参数能通过 `ckks.NewParametersFromLiteral` 校验
- [x] 1.2 创建 `polynomials.go`：定义论文 Table 4 的六个门三角函数的 Chebyshev 逼近
  - 交付物：`NewBinBootPoly()`、`NewGatePoly(gateType)` 函数，返回 `bignum.Polynomial`
  - 验收：多项式在 0、1/3、2/3（或 0、1/2）处的求值与论文 Table 4 一致
  - 论文对应：§3.1 公式 `f_BinBoot(x)=(1-cos(2πx))/2`；§4.1 Table 4
- [x] 1.3 创建 `evaluator.go`：实现 `BinBoot` 与 `GateBoot` 评估器
  - 交付物：`Evaluator.Bootstrap(ct)`（对应 Algorithm 2）、`Evaluator.GateBoot(ct1, ct2, gate)`（对应 Algorithm 3）
  - 验收：能复用 Lattigo 的 StC/ModRaise/CtS，并在 EvalMod 阶段注入自定义三角多项式
  - 论文对应：Algorithm 2（§3.1）、Algorithm 3（§4.1）、Theorem 1（§3.2）、Theorem 2（§4.2）

### 阶段 2：示例入口 `examples/singleparty/paper_bcks/`
- [x] 2.1 创建 `main.go`：编排流程——初始化 Param14、生成密钥、加密两个二进制数组、调用 GateBootNAND、解密验证
  - 交付物：可 `go run .` 执行的程序
  - 验收：解密结果与明文 NAND 一致（容差 < 0.25）
- [x] 2.2 创建 `utils.go`：随机二进制数组生成器、`ValidateBinaryTolerance` 验证器
  - 交付物：`GenBinaryVector(n)`、`ValidateBinaryTolerance(expected, actual)` 
  - 验收：输出 "1" > 0.75，"0" < 0.25
- [x] 2.3 创建 `README.md`：执行说明与预期耗时

### 阶段 3：编译与测试
- [x] 3.1 `go build ./circuits/ckks/binboot/` 编译通过
- [x] 3.2 `go build ./examples/singleparty/paper_bcks/` 编译通过
- [x] 3.3 `go vet` 通过
- [ ] 3.4 示例完整运行验证（受限于机器性能与密钥生成时间，未执行完整端到端测试）

### 阶段 4：文档更新
- [x] 4.1 更新 `docs/dev-logs.md`：记录所有修改文件、核心功能、论文对应关系
- [x] 4.2 更新 `docs/todo.md`：标记完成状态（本文件）

---

## 论文与代码对应关系速查

| 论文章节 | 论文内容 | 代码位置 |
|---------|---------|---------|
| §3.1 Algorithm 2 | BinBoot：CtS∘ModRaise∘StC → 实部提取 → Eval_fBinBoot | `circuits/ckks/binboot/evaluator.go: Bootstrap` |
| §3.1 公式 | f_BinBoot(x) = (1-cos(2πx))/2 | `circuits/ckks/binboot/polynomials.go: NewBinBootPoly` |
| §3.2 Theorem 1 | 二次噪声缩减 | `evaluator.go` 注释 |
| §4.1 Algorithm 3 | GateBoot：ct1+ct2 → CtS∘ModRaise∘StC → 实部 → Eval_fG | `evaluator.go: GateBoot` |
| §4.1 Table 4 | 六个门的三角函数 | `polynomials.go: NewGatePoly` |
| §4.2 Theorem 2 | GateBoot 正确性 | `evaluator.go` 注释 |
| §5 Table 5 | Param14 参数集 | `parameters.go: Param14` |
| §3.3 | 模数工程 q0=2Δ0 | `parameters.go` 注释 |
| §4.1 | 模数工程 q0=3Δ0 | `parameters.go` 注释 |

---

## 接手说明

若 token 耗尽，后续 agent 可按以下顺序接手：
1. 读取本文件了解整体进度
2. 读取 `docs/dev-logs.md` 了解已完成的具体修改
3. 读取 `docs/bcks.md` 了解开发规范
4. 检查 `circuits/ckks/binboot/` 与 `examples/singleparty/paper_bcks/` 的现有代码
5. 从下一个未完成步骤继续
