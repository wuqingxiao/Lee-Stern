# 最后一轮竞争攻击与最终参数判定 v7

**本轮已经运行完成，最终表已附。无需再次搜索。**

先读 FINAL_PARAMETER_DECISION.md；完整 C 与精确数值见 FINAL_PARAMETER_DECISION.json。三个固定候选最终均判“证据不足”，未获批为完整认证的经典 128 位参数。高昂攻击的上预算不能倒置为安全下界。

Python 3.10+，仅标准库。可独立检查保存结果：

```powershell
python verify_final.py
```

如需检查小规模执行及精确枚举：

```powershell
python test_round.py
```

六项测试包括独立域/表示枚举、8! 组成分割枚举、不对称条件矩、实际打包连接及 18,720 个满行秩矩阵/偏移组合。verify_final.py 复算 6 份最优证书、核查全部 35,802 行的最小值和覆盖状态；这是模型一致性检查，不是全局安全证明。

已有运行约 127.5 秒，输出 final_round/。若仅为独立复现原来的有限网格，可另选空目录：

```powershell
python run_round.py --out my_reproduction
```

时间是本机观察，其他设备可能不同；并集计数子进程有固定 4 秒限额，因此复现时可能增加超时缺口，不能把超时当作攻击排除。本轮现存结果没有超时，22 个子进程失败均为显式计数容量上限。

主要文件：

- FINAL_PARAMETER_DECISION.md/.csv/.json：最终三行判定及证据边界，JSON 包含完整 C。
- THEORY_AND_SCOPE.md：固定 E、完整并集、精确计数、条件概率链、费用和未覆盖项。
- fixed_domain.py：固定 E 域的精确 DP 和双射反排名。
- final_model.py：修复费用/布局的适配与不对称评估。
- execution.py、vendor/：冻结 v6.1 内核；来源绑定见 SOURCE_BINDINGS.json。
- final_round/points.csv、coverage.json、best.json、SUMMARY.json、VERIFICATION.json：全部有限点、覆盖缺口、6 份最优保存记录与复核结果。
- inputs/：同一三个完整 C、旧普通攻击记录、v6.1 记录和原文/条件认证账本。

离线包省略本地约百 MB 的重复中间 selected.json 和三个 *_fixed_selected.json。保留每一个有限点、每个域的状态和最优完整证书，足以复核交付表；中间文件可以通过同一有限网格重新生成。原始完整目录没有删除或覆盖。

内存上限是声明打包保留区 64 GiB；不是 Python RSS。目标规模未实际执行指数级攻击；schoolbook 位收费约定未作完整处理器编译认证。当前成果可以支持可复核的条件评估，不能直接支持具体 128 位安全声明。本轮到此结束，不自动扩大外部参数或继续实验。
