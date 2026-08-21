# 开发日志 / Dev Logs

本文件按时间倒序记录每次任务的执行过程与产出。最新一次在最上面。

---

# 2026-08-10 - 生成 CKKS Bootstrapping 论文学习笔记 (2020-1203.ipynb)

## 任务目标

为刚转密码学方向的博士生（计算机科学本科背景）生成一份中文学习笔记，配合论文和 Lattigo 源码，使其"只读这一份笔记就能全面理解论文和代码"。

- **论文**：`my_ckks_study/2020-1203.pdf` — Bossuat et al., *Efficient Bootstrapping for Approximate Homomorphic Encryption with Non-Sparse Keys*, Eurocrypt 2021（扩展版，46 页）。
- **代码**：本仓库 Lattigo v6，自举实现主要在 `circuits/ckks/bootstrapping/`。
- **输出**：`my_ckks_study/2020-1203.ipynb`（Jupyter notebook，中文 markdown + Python 代码）。

## 执行过程

### 1. 论文精读

- 用 `pypdf` 提取了 46 页论文全文到 `paper_text.txt`（中间产物，已清理），完整阅读。
- 论文核心：首个 **128-bit 安全**、不强制稀疏密钥、有精确失败概率评估的 CKKS 自举，三大贡献：
  - **贡献①（§3）** Algorithm 2 — 深度最优 + 无误差的同态多项式求值（消除 Full-RNS 尺度偏差）。
  - **贡献②（§4）** Algorithm 6 — 双重提升 BSGS 矩阵-向量乘（线性变换快 ~2×）。
  - **贡献③（§5-6）** 5 步自举电路的系统化参数化。

### 2. 源码映射（4 个并行探索 agent）

并行派出 4 个 Explore agent，把论文每个算法映射到具体的 Lattigo v6 代码位置：

| 探索目标 | 关键发现 | 代码位置 |
|---|---|---|
| Algorithm 2 (EvalRecurse) | v6 重命名为 **Paterson-Stockmeyer**，函数是 `recursePS`（不叫 EvalRecurse） | `circuits/common/polynomial/polynomial.go:109` |
| 尺度传播（核心创新） | 由 `simEvaluator` 接口实现，CKKS 实现含 baby/giant step | `circuits/ckks/polynomial/polynomial_evaluator_sim.go:52,67` |
| Algorithm 6 (双重提升) | `MultiplyByDiagMatrixBSGS`，两层 hoisting 都在 `lintrans_evaluator.go` | `circuits/common/lintrans/lintrans_evaluator.go:280` |
| EvalSine (Algorithm 7) | `EvaluateAndScaleNew`，含双角循环 + 可选 arcsin | `circuits/ckks/mod1/mod1_evaluator.go:31` |
| 5 步主入口 | `bootstrap` 方法，代码把 ModRaise 拆成 ScaleDown+ModUp | `circuits/ckks/bootstrapping/evaluator.go:552` |
| 参数表（Table 5） | 8 套默认参数，分 Sparse/Dense 两组 | `circuits/ckks/bootstrapping/default_parameters.go` |

注意发现：v6 默认开启了 2022/024 的稀疏密钥封装（`EphemeralSecretWeight=32`），这是相对原论文的增量优化。

### 3. 笔记生成

- 用 `nbformat` 写了一个 `build_notebook.py` 脚本，以可编辑形式组织 24 个 cell。
- 笔记分 **8 个部分**（序言 + 7 章正文），约 6.6 万字：
  - Part 0：序言 + 论文↔代码速查表
  - Part 1：数学/密码学前置（分圆环、RLWE、RNS、自举动机）
  - Part 2：Full-RNS CKKS 方案（论文 §2.1）
  - Part 3：自举 5 步总览（论文 §2.2, §5.1）
  - Part 4：贡献① Algorithm 2（论文 §3），逐行对照 `recursePS`
  - Part 5：贡献② 密钥切换 + 双重提升（论文 §4）
  - Part 6：5 步细节（论文 §5），含 EvalSine 完整代码解读
  - Part 7：参数选择与安全性（论文 §6），含 Equation 1 复现
  - Part 8：结果、总结、学习建议（论文 §7-8）

### 4. 代码 cell 调试（toy 化）

8 个 Python 代码 cell 用 **toy 小维度**演示数学内核，用户明确要求"接受更低安全等级和准确性，只要逻辑对应论文"。最终全部无错执行：

| cell | 演示内容 | 对应论文 |
|---|---|---|
| 1 | R_Q 上 negacyclic 多项式乘法 | §2.1 |
| 2 | Full-RNS rescale 尺度偏差量化 | §3 引子 |
| 3 | Chebyshev 多项式 + 乘积恒等式 | §3.1 |
| 4 | Algorithm 2 尺度传播消除加法误差 | §3.2, Table 1 |
| 5 | 双重提升 vs 朴素/单提升 BSGS 复杂度曲线 | §4.3, Table 2 |
| 6 | EvalSine 双角公式降低多项式度数 | §5.4 |
| 7 | Irwin-Hall 失败概率公式（Equation 1）toy 复现 | §6.2, Table 4 |
| 8 | 失败概率 & 安全性约束可视化 | §6.1-6.2 |

调试中遇到并解决的问题：
- 三引号转义：`r"""..."""` 内嵌 docstring 改用 `'''`。
- Algorithm 2 递归模拟的尺度闭合性：改为 baby-step 因子反推的直观演示。
- BSGS 成本比较的 `n1*n2=n` 守卫：改为 log-ratio 扫描。
- Irwin-Hall 公式大 H 溢出：按用户建议**toy 化**（小 h 直接精确算），验证公式正确性而非追求大参数精度。

### 5. 验证

- `jupyter nbconvert --execute` 全部 24 cell 无错通过（输出 484KB executed notebook，验证后清理）。
- 关键代码引用（5 个文件:行号）全部存在于笔记中。
- 反向提取 `build_notebook.py`，**round-trip 验证一致**：重跑脚本生成的 notebook 与原 notebook cell 内容完全相同。

## 最终产出

| 文件 | 说明 |
|---|---|
| `my_ckks_study/2020-1203.ipynb` | 学习笔记主文件（24 cell，约 6.6 万字） |
| `my_ckks_study/build_notebook.py` | 笔记生成脚本（可编辑、可重建，round-trip 一致） |
| `my_ckks_study/dev-logs.md` | 本开发日志 |

## 长效约定（本次确立）

1. 大内容任务拆成多章节。
2. 保留所有中间脚本（`build_notebook.py` 等）。
3. 每次任务完成后更新 `dev-logs.md`，节标题格式 `#日期-时间-任务名称`。
