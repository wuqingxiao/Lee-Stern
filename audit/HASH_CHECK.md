# HASH_CHECK：承诺哈希与盐的实现核对

日期：2026-11-01。

## 1. 核对对象

- 实现：`F:\codex-pj\ac3\output\lee_stern_revision_20261029\implementation\protocol.py`
- SHA-256：`90BCACBA74138DF183C5AB326745502A35D05F627C49343E8037104F1798F34C`
- 与 `F:\codex-pj\ac3\research\_journal_full_build\implementation\protocol.py` 完全相同。
- 冻结历史副本 `historical_evidence\frozen_implementation\protocol.py` 的哈希不同（313EDF...），但哈希与盐相关的代码行相同。

## 2. 已核实的实现细节

| 项 | 实现 | 位置 |
|---|---|---|
| 哈希函数 | `hashlib.sha256` | 第 63 行 |
| 上下文前缀 | `frame(b'COM-SHA256-SALT32-v1', ctx)`，ctx = frame('LeeStern-ExactC-v1', p, n, k, pack(C), pack(H), pack(s))，即完整公开实例，没有用摘要代替 | 第 63 行 |
| 单个承诺 | 复制前缀后 update(frame(bytes([tag]), *parts, salt)) | |
| 编码单射 | 每个字段带 8 字节长度前缀 | frame() |
| 域分离 | tag = 1/2/3 区分三类承诺 | |
| 盐 | `secrets.token_bytes(32)`，每轮 3 个新盐 | 第 89 行 |
| 长度检查 | 承诺与验证都要求盐长 = 32、摘要长 = 32 | |
| 比较 | `hmac.compare_digest` | 第 132 行 |
| b=2 检查 | counts(v) == C，即 v ∈ S_C | |
| 其他检查 | σ 是置换；向量分量在 [0, p) | |
| 状态 | 一次性：先置 used，再验证挑战；任何错误都中止会话 | 第 102–107 行 |
| 随机数 | 操作系统随机源，精确拒绝采样，无取模偏差 | |
| 测试 | `test_state_rules.py`：10 passed，8 subtests passed | |
| 轮数 | 基准测试用 219 轮；代码注释写明“219 是实验轮数，不是 128 位声明” | |

## 3. 结论

- 绑定：SHA-256 通用经典碰撞代价 2^128，满足 λ_hash = 128。任一类承诺出现碰撞都会破坏第 2 节的抽取（见 SOUNDNESS_P0.md），损失因子 ≤ 3。多目标不帮助找碰撞，因为生日界本身已是最优通用攻击。
- 隐藏：在随机预言模型下来自 256 位盐。暴力搜索盐或求原像代价 2^256。cm2/cm3 的消息是随机置换后的随机向量，熵很高。
- 通信：每轮承诺 3 × 32 = 96 字节；挑战 1 字节；响应中盐 2 × 32 = 64 字节；与哈希/盐相关的合计 161 字节/轮。

## 4. 限定与需改动处

- 这是研究实现，没有做常数时间处理，侧信道不在范围内。
- HVZK 模拟器用真实的哑承诺，没有编程随机预言。这对交互式 HVZK 足够，论文里要按此描述。
- 论文 §4.1（第 188 行附近）把固定 SHA-256 实例说成“用于性能测量”。P0 下它就是安全配置，这句措辞要改（只标记，未改论文）。