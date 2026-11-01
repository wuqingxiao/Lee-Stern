"""Produce a bounded, evidence-qualified verdict table from frozen results."""
from collections import Counter, defaultdict
import csv
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write(path,obj):Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')


def main():
    out=HERE/'final_round';cfgs=read(HERE/'inputs/candidates.json')
    summary=read(out/'SUMMARY.json');coverage=read(out/'coverage.json');best=read(out/'best.json')
    v61=read(HERE/'inputs/v61_selected.json')
    ordinary={x['candidate']:x['prior_price'] for x in read(HERE/'inputs/prior_ordinary_schedules.json')['selected']}
    baseline=read(HERE/'inputs/baseline_and_authentication.json')
    with (out/'points.csv').open(encoding='utf-8-sig',newline='') as f:points=list(csv.DictReader(f))
    rows=[];groups=defaultdict(list)
    for p in points:groups[(p['candidate'],p['family'],p['asymmetric'])].append(p)
    for cfg in cfgs:
        name=cfg['name'];prev=min((x for x in v61 if x['candidate']==name),key=lambda x:x['declared_bit_upper_log2'])
        new=next(x for x in best if x['candidate']==name and x['final_selection_objective']=='declared_bit_upper_log2')
        auth=baseline[name]['authentication'];core=auth['required_FCSD_assumption']
        # This checks an inherited conditional ledger, not the reduction theorem.
        noncore=sum((F(x['exact']) for x in auth['terms'].values()),F(0))
        assert noncore==F(auth['noncore_sum']['exact'])
        assert core['charge_cap']==str(2**166) and core['RO_query_cap']==str(2**166)
        assert F(core['maximum_success'])==F(1,4)
        assert F(auth['conditional_total_if_core_at_most_one_quarter']['exact'])==F(1,4)+noncore<F(1,2)
        assert auth['Q']+int(auth['wrapper_charges']['total'])<2**166
        # Frozen Q also includes simulator commitment queries; it is not just RB.
        assert auth['R']*2**128 <= auth['Q'] < auth['R']*(2**128+3)
        mode_min={}
        for x in points:
            if x['candidate']==name:
                key=x['family']+'/'+x['schedule'];value=float(x['declared_bit_upper_log2'])
                mode_min[key]=min(mode_min.get(key,float('inf')),value)
        observed=min(ordinary[name]['log2_charge'],prev['declared_bit_upper_log2'],new['declared_bit_upper_log2'])
        row=dict(candidate=name,recipe=cfg['primary_recipe'],p=cfg['p'],n=cfg['n'],k=cfg['k'],weight=cfg['weight'],
                 C_sha256=cfg['C_sha256'],C=cfg['C'],
                 original_shell_formula_reference_log2=baseline[name]['original_formula']['log2_reference_cost'],
                 original_shell_reference_is_actual_C_or_same_unit=False,
                 previous_ordinary_upper_log2=ordinary[name]['log2_charge'],
                 previous_ordinary_success_lower=ordinary[name]['success_lower'],
                 previous_ordinary_reserved_GiB=ordinary[name]['memory']/2**30,
                 previous_v61_RB_upper_log2=prev['declared_bit_upper_log2'],
                 v7_best_upper_log2=new['declared_bit_upper_log2'],
                 v7_best_heuristic_log2=new['heuristic_log2'],v7_best_family=new['family'],
                 v7_best_parameters=new['parameters'],v7_best_asymmetric=new['asymmetric'],
                 v7_best_success_lower=new['success_lower'],v7_declared_reserved_GiB=new['memory']['reserved_bytes']/2**30,
                 v7_minima_by_family_and_schedule=mode_min,
                 best_available_constructive_upper_log2=observed,
                 v7_changes_unified_minimum=observed<ordinary[name]['log2_charge'],
                 v7_reduces_previous_RB_upper_by_bits=max(0,prev['declared_bit_upper_log2']-new['declared_bit_upper_log2']),
                 coverage_statuses=dict(Counter(x['status'] for x in coverage if x['candidate']==name)),
                 arithmetic_ledger_under_inherited_premise_passes=True,
                 authentication_N=auth['N'],hash_bits=auth['h'],salt_bits=auth['rho'],learning_sessions_cap=auth['learning_sessions_cap'],
                 authentication_reduction_log2_charge=auth['log2_reduction_charge'],
                 required_FCSD_budget_log2=166,required_FCSD_success_upper='1/4',
                 general_FCSD_premise_established=False,physical_RSS_certified=False,whole_bit_compiler_certified=False,
                 cheap_attack_certificate_found_in_saved_rows=observed<=166,
                 verdict='证据不足',supports_complete_classical128=False,
                 reason='现有数值为指定攻击的费用上预算，不能推出一般FCSD时间—成功率上界；完整并集及其他竞争攻击仍有覆盖缺口。',
                 security_bits=None)
        assert not row['cheap_attack_certificate_found_in_saved_rows']
        rows.append(row)
    result=dict(round_status='FINISHED',date='2026-10-25',target='classical128 authentication work-factor target',
                scope='frozen three complete compositions; explicit finite schedules; 64 GiB declared packed reservations',
                points=summary['points'],statuses=summary['coverage_statuses'],
                outcome='三组均为证据不足；本轮结束，不自动继续扩大搜索。',
                risk_flag_is_not_security_proof=True,security_bits=None,candidates=rows)
    write(HERE/'FINAL_PARAMETER_DECISION.json',result)
    columns=['candidate','recipe','p','n','k','weight','C_sha256','previous_ordinary_upper_log2',
             'previous_v61_RB_upper_log2','v7_best_upper_log2','v7_best_heuristic_log2',
             'best_available_constructive_upper_log2','v7_declared_reserved_GiB','v7_best_success_lower',
             'required_FCSD_budget_log2','required_FCSD_success_upper','verdict','supports_complete_classical128','reason']
    with (HERE/'FINAL_PARAMETER_DECISION.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)
    text=['# 最终参数判定表：本轮已结束','','日期：2026-10-25。目标为经典认证工作量目标；固定现有三个完整 C，不再自动扩展外部参数。',
          '', '**最终判定：三组候选均为“证据不足”。没有一组在现有材料下批准为完整认证的经典 128 位参数；这也不是已经发现三组均可被攻破。**',
          '', '以下数字均为 log2(声明位收费上预算)，成功目标至少 1/2、声明打包保留区≤64 GiB。它们是可构造策略的上预算，不是安全位数，不是全体攻击的下界；不包含处理器/RSS认证。',
          '', '| 固定 C 的候选 | (p,n,k,w) | 已有普通 PGE/Wagner 上预算 | v6.1 RB 上预算 | 本轮最佳上预算 | 最终判定 |',
          '|---|---|---:|---:|---:|---|']
    for x in rows:
        text.append(f"| {x['recipe']} / `{x['candidate'].split('_')[-1]}` | ({x['p']},{x['n']},{x['k']},{x['weight']}) | {x['previous_ordinary_upper_log2']:.3f} | {x['previous_v61_RB_upper_log2']:.3f} | {x['v7_best_upper_log2']:.3f} | **证据不足** |")
    text+=['', '本轮没有降低上述已有普通 PGE/Wagner 策略给出的最低上预算；两组较重组成的 RB 上预算分别改善约 3.622 和 11.147 位。这只是当前策略/界更紧，不能转化为安全位数变化。',
           '', '## 本轮实际覆盖', '',
           '- 共 35,802 点：fixed_E 34,650 点、完整置换并集 full_union 1,152 点。',
           '- 300 个 fixed_E 域全部完成；36 个 full_union 域提案中 8 个完成，22 个达到计数上限，6 个超过当前缓存保留量。未完成者均列入覆盖缺口。',
           '- 半径 r∈{2,3,4,6,8}、自由位置 q∈{1,4,8,16,32}、ell∈{2,8,16,32}；full_union 使用另行限定的较小网格，详见 config.json。',
           '- L 与 u 分开选择，包含可行的 L=p^u 点、内存上限子表、2^24 子表及全域流式；驻留、部分分块、全域流式都收费。只搜索了三种有限组成移位，不是任意不对称组成的全局优化。',
           '- 有限不对称分割在本轮没有优于对称分割。fixed_E 与 full_union 分开计算，不能互相替代。',
           '', '## 最优新增记录', '',
           '| 候选 | 基础域 | ell,u,r,q,a | L | 根块 | 调度 | 声明保留 GiB | 成功下界 | 启发式 log2费用 |',
           '|---|---|---|---:|---:|---|---:|---|---:|']
    for x in rows:
        p=x['v7_best_parameters']
        text.append(f"| {x['recipe']} | {x['v7_best_family']} | {p['ell']},{p['u']},{p['radius']},{p['q']},{p['a']} | {p['L']} | {p['block_records']} | {p['schedule']} | {x['v7_declared_reserved_GiB']:.9f} | {x['v7_best_success_lower']} | {x['v7_best_heuristic_log2']:.3f} |")
    text+=['', '启发式列依赖 Poisson 存活、密钥均匀化和整周期代理首次命中，不能据此宣称攻击不可行。最优驻留记录几乎用满声明的 64 GiB；实际实现的对象/分配器开销可能越界，因此这里没有宣称可用 Python 在该物理内存内运行目标攻击。',
           '', '## 为什么没有判为通过', '',
           '1. **逻辑缺口**：昂贵攻击的费用上界不推出所有更便宜攻击不存在。本轮的新计数、概率下界和费用上预算仍属这类证据。',
           '2. **覆盖缺口**：full_union 的未计数域、更多组成/多层混合策略、公开代数和约束带来的强化，以及高级汉明/多重集/格恢复，没有被本轮全部排除。',
           '3. **协议接口缺口**：冻结的条件认证账本需要 FCSD 在预算和查询数均≤2^166 下成功率≤1/4。已有账本的算术在此假设下通过，但这一底层前提尚未由参数评估建立。',
           '', '256 轮、384 位哈希和 256 位盐的理论认证配置按最多 2^20 次学习会话核对；继承的归约费用约 2^164.331，非核心误差和约 2^-54.016，若核心成功率≤1/4，则总界<1/2。这是条件工作量目标，不是声称成功概率≤2^-128；该哈希配置也尚未在原冻结协议中实现。',
           '', '## 原文基线保留方式', '',
           '| 候选 | 原 Li–Wang 总 Lee 重量壳公式参考 log2费用 | 能否直接用于完整 C 的 128 位结论 |',
           '|---|---:|---|']
    for x in rows:text.append(f"| {x['recipe']} | {x['original_shell_formula_reference_log2']:.3f} | 不能；原公式参考不是本次固定 C 与统一位收费评估 |")
    text+=['', '## 可采用的论文结论与本轮终点', '',
           '可以写：在明确列出的固定组成、有限策略、声明收费和内存模型下，完成了上述竞争评估；本轮未构造出预算≤2^166且成功率≥1/2的保存策略。不能写：已排除预算≤2^166、成功率>1/4的所有攻击，或已获得具体经典 128 位安全参数。',
           '', '这张表是本轮最终交付，无需再跑一轮才能读到判定。三个 C 保留为条件安全与实现研究对象；当前不推荐其中任何一个作为已验证的 128 位参数。后续论文按“协议修复、条件安全证明与可复核评估”表述，具体 128 位声明暂不采用。',
           '', '完整 C、来源哈希和每点结果分别见 inputs/candidates.json、SOURCE_BINDINGS.json、final_round/points.csv；理论与实现边界见 THEORY_AND_SCOPE.md。']
    (HERE/'FINAL_PARAMETER_DECISION.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(json.dumps(dict(rows=len(rows),verdicts=[x['verdict'] for x in rows],points=len(points)),ensure_ascii=False))


if __name__=='__main__':main()
