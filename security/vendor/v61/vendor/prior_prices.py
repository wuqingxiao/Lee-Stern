"""Final, bounded assessment of TWO saved Lee--Stern candidates.

Combines evidence without turning named-resource upper budgets into hardness.
Adds actual-C restricted-set points and an explicit abstract bit-charge ledger.
Cannot manufacture a positive success bound or a lattice operation conversion.
"""
from fractions import Fraction as F
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent
MEMORY=64*2**30
NAMES=('wide700_n320','low370_n448')


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def lg(v):
    v=F(v)
    return None if not v else math.log2(v.numerator)-math.log2(v.denominator)
def module(name,path):
    path=Path(path);sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
def same(a,b):return all(a[k]==b[k] for k in ('p','n','k','C'))


def resource_prices(cfg,rank_bits,address_bits):
    """Explicit analytic charge model; not a certified CPU/bit compiler.

    Unknown instructions are rejected. Full vectors, not only field mul,
    are priced. Upper-count budgets remain constructive upper budgets.
    """
    b=(cfg['p']-1).bit_length();a=address_bits;B=max(rank_bits,a,b,cfg['n'].bit_length())
    loop=64*(a+B)
    base={
        'field_multiplications':8*b*b+4*a,
        'field_additions_subtractions':4*b+4*a,
        'field_inversions':16*b**3+6*a,
        'field_coordinate_comparisons':b+2*a,
        'histogram_counter_comparisons':cfg['n'].bit_length()+2*a,
        'histogram_counter_updates':2*cfg['n'].bit_length()+2*a,
        'Lee_integer_additions':2*cfg['n'].bit_length()+2*a,
        'vector_coordinate_writes':b+2*a,
        'preparation_coordinate_writes':b+2*a,
        'partition_index_writes':cfg['n'].bit_length()+2*a,
        'base_sort_record_moves':(cfg['n']+cfg['n']-cfg['k'])*b+2*a,
        'root_sort_record_moves':(cfg['n']+cfg['n']-cfg['k'])*b+2*a,
        'rank_bigint_multiplications':8*B**2+4*a,
        'rank_bigint_divisions':16*B**2+4*a,
        'rank_bigint_additions_subtractions':4*B+4*a,
        'rank_bigint_comparisons':2*B+4*a,
        'rank_tree_node_visits_upper':4*B+8*a,
        'expected_random_bits':1,
        'partition_uniform_integer_calls':8*cfg['n'].bit_length()+4*a,
        'sample_uniform_integer_calls':8*B+4*a,
        'offset_uniform_field_calls':8*b+4*a,
        'basis_uniform_field_calls':8*b+4*a,
        'PGE_pivot_tests':b+2*a,
        'sort_payload_bytes_read_upper':8+2*a,
        'sort_payload_bytes_written_upper':8+2*a,
        'public_integer_multiplications':8*B**2+4*a,
        'public_integer_additions_subtractions':4*B+4*a,
        'coefficient_table_cell_writes':B+2*a,
        'reservoir_uniform_integer_calls':8*B+4*a,
        'first_join_queries':B+4*a,
        'first_join_records':B+4*a,
        'root_matches':B+4*a,
        'streamed_leaf_vectors':cfg['n']*b+4*a,
        'generator_symbol_tests_upper':b+4*a,
    }
    return {k:v+(0 if k=='expected_random_bits' else loop) for k,v in base.items()}


def price_pge_schedule(cfg,row):
    rank=row['rank_operand_max_bits']
    caps={k:int(v) for k,v in row['total_caps'].items()}
    address=max(36,max(v.bit_length() for v in caps.values()))
    prices=resource_prices(cfg,rank,address)
    unknown=sorted(set(caps)-set(prices))
    if unknown:raise ValueError('unpriced resources: '+str(unknown))
    charge=sum(caps[k]*prices[k] for k in caps)
    return dict(log2_abstract_bit_charge_upper_model=lg(charge),exact_charge=str(charge),
        prices=prices,rank_operand_bits=rank,address_and_long_counter_bits=address,
        all_named_resource_classes=len(caps),unknown_resource_classes=unknown,
        success_lower=row['ensemble_success_lower'],payload=row['memory']['reserved_bytes'],
        first_hit_stops_no_later_than_full_scan=True,
        cost_scope='declared schoolbook bit-RAM CHARGE ESTIMATE on source counted schedule',
        instruction_coverage_and_schoolbook_bounds_formally_compiled=False,
        processor_time_or_whole_RSS_certified=False,
        security_lower_bound=False)


def prepare_restricted(cfg,comp):
    templates=[]
    for radius in (1,2,3):
        ells=(8,24);maxh=(cfg['k']+max(ells))//2
        states=comp.elite_symmetric(cfg['C'],cfg['p'],radius,maxh,96,keep=1)
        for ell in ells:
            h=(cfg['k']+ell)//2
            choices=[]
            for (length,weight),paths in states.items():
                if length==h and 0<weight<=96 and weight%2==0 and weight%cfg['p']:
                    score,D=paths[0]
                    if D[0]>=1:choices.append((score,weight,D))
            for score,weight,D in sorted(choices,key=lambda x:(-x[0],x[1]))[:2]:
                templates.append(dict(ell=ell,radius=radius,a=weight//2,q=1,left=D,right=D,
                    proposal_rule='top two split-probability symmetric DP types at each h/radius; not cost-optimized or exhaustive',score=score))
    return templates


def restricted(cfg):
    # Names are deliberately explicit: the two domains must not get mixed.
    folder=RESEARCH/'stern_restricted_set_scalable_v1'
    domain=module('domain',folder/'domain.py')
    module('moments',folder/'moments.py');module('list_costs',folder/'list_costs.py')
    distinct=module('final_distinct_analyzer',folder/'analyzer.py')
    tagged=module('final_fixed_mask_analyzer',RESEARCH/'stern_restricted_bjmm_v3/analyzer.py')
    comp=module('final_restricted_compositions',RESEARCH/'stern_restricted_bjmm_v3/compositions.py')
    templates=prepare_restricted(cfg,comp)
    rows=[];rejected={};domain_limits=[]
    for tid,t in enumerate(templates):
        h=(cfg['k']+t['ell'])//2
        try:
            B=domain.Domain(cfg['p'],h,t['a'],1,t['radius'],max_candidates=150000)
        except domain.CapacityLimit as exc:
            B=None;domain_limits.append(dict(template=tid,reason=str(exc)))
        for mode in ('full','amortized'):
            for u in range(1,t['ell']):
                if mode=='amortized' and 3*u<t['ell']:continue
                try:
                    row=tagged.assess(cfg,t['left'],t['right'],t['ell'],t['radius'],2,u,mode,
                        block_records=2**20,mask_policy='fixed',include_source_proxy=False)
                    row['variant']='ALGORITHM3_FIXED_MASK';row['template']=tid
                    row['complete_layout_64GiB_certified']=False
                    row['full_layout_missing_areas']=row['payload_excludes']
                    row['partial_payload_below_64GiB']=int(row['record_and_main_matrix_payload_bytes'])<=MEMORY
                    rows.append(row)
                except ValueError as exc:
                    label='fixed: '+str(exc);rejected[label]=rejected.get(label,0)+1
                if B is not None:
                    try:
                        row=distinct.assess(cfg,t['left'],t['right'],t['ell'],t['a'],1,t['radius'],u,mode,
                            block_records=2**20,domain=B)
                        row['variant']='DISTINCT_PERMUTATION_CLOSURE';row['template']=tid
                        row['main_packed_payload_below_64GiB']=int(row['listed_attack_packed_payload_bytes']['exact'])<=MEMORY
                        row['complete_layout_64GiB_certified']=False
                        rows.append(row)
                    except ValueError as exc:
                        label='distinct: '+str(exc);rejected[label]=rejected.get(label,0)+1
        print(json.dumps(dict(restricted_template=cfg['name'],done=tid+1,total=len(templates),rows=len(rows))),flush=True)
    half=[];distinct_fit=[]
    for r in rows:
        cert=r.get('theoretical_counted_budget') if r['variant']=='DISTINCT_PERMUTATION_CLOSURE' else r.get('constant_success_counted_resource_certificate')
        if cert is not None:half.append(r)
        if r['variant']=='DISTINCT_PERMUTATION_CLOSURE' and r['main_packed_payload_below_64GiB'] and cert is not None:distinct_fit.append(r)
    return dict(candidate=cfg,templates=templates,rows=rows,
        summary=dict(templates=len(templates),evaluated_points=len(rows),rejections=rejected,
            finite_counted_half_success_certificates=len(half),distinct_main_payload_and_half_success_points=len(distinct_fit),
            domain_capacity_limits=domain_limits,
            fixed_mask_is_not_tagged_all_mask_domain=True,
            duplicate_free_domain_not_replaced_by_tagged_cardinality=True,
            full_representation_survival_upper_bound_established=False,
            scopes='source fixed-mask variant and distinct-vector closure; fixed-C geodesic SUBFAMILY proves success lower only',
            complete_half_success_bit_cost=None,assessed_as_generic_128_safe=False))


def auth_curve(budget,gap):
    """Explicit time/success interface, with missing FCSD curve retained."""
    delta=F(gap['failure_probability_upper']['exact'])
    overhead=F(budget['noncore_sum']['exact'])
    return dict(adversary_charge_cap=str(2**128),RO_query_cap=str(2**128),
        learning_sessions=2**20,rounds=256,hash_bits=384,salt_bits=256,
        exact_reduction_charge=budget['reduction_charge_exact'],
        log2_reduction_charge=budget['log2_reduction_charge'],
        exact_reduction_queries=str(budget['Q']),
        retained_sufficient_resource_envelope=str(2**166),
        required_FCSD_success_upper='1/4',
        sufficient_impersonation_upper_IF_FCSD_assumption=str(F(1,4)+overhead),
        sufficient_impersonation_upper_below_half=(F(1,4)+overhead<F(1,2)),
        exact_noncore_overhead=str(overhead),log2_noncore_overhead=lg(overhead),
        time_success_curve_for_actual_C_established=False,
        target_half_success_cost_cannot_alone_bound_quarter_success_at_envelope=True,
        conditional_lattice_FCSD_lower_IF_SVP_success_3_over_4=str(max(F(0),F(3,4)-delta)),
        assumptions=['classical ROM and ideal independent randomness','declared resource ledger',
            'single user, strictly sequential active learning followed by one impersonation',
            'FCSD hardness assumption at the declared reduction resources'],
        authentication_128_certified=False,
        old_219_round_SHA256_performance_configuration_unchanged=True)


def report(records):
    flat=[];coverage=[]
    closed_file=HERE/'results/LMMT_CLOSED_SUMMARY.json'
    closed={r['name']:r for r in read(closed_file)} if closed_file.exists() else {}
    for r in records:
        c=r['candidate'];s=r['screen'];p=r['pge_wagner'];L=r['lattice']
        for method,status,why in (
            ('Li-Wang_LMMT','SPECIFIED_TILED_SCHEDULE_CLOSED_NOT_COMPETITIVE' if c['name'] in closed else 'NOT_CLOSED',
             'exact finite profile union + pointwise GL survival + counted cutoff closes a full-shell schedule; optimal/full sigma and competitive best cost unknown' if c['name'] in closed else 'joint success lower zero; per-proposal charges alone insufficient'),
            ('PGE_Wagner','FINITE_CONSTRUCTIVE_SCHEDULE_ONLY','full named vector priced in explicit model; not best-attack hardness lower bound'),
            ('Restricted_Balls','FINITE_SUBFAMILY_ASSESSMENT_ONLY','actual-C fixed-mask and distinct domains evaluated; subfamily survival is not whole-attack upper'),
            ('elementary_Hamming','BASELINE_ONLY','exact selection events; advanced Hamming ISD not assessed'),
            ('lattice_64GiB','HEURISTIC_RISK_SCREEN','GSA/projection/enumeration recovery and cycle/bit conversion uncalibrated'),
            ('lattice_unrestricted','PUBLISHED_MODEL_RISK','full-SVP planted separation proved; finite solver success and cost conversion remain assumptions'),
            ('multiset_PKP','NOT_COVERED','mutually-distinct-element source tables not applicable to this repeated C'),
            ('authentication','CONDITIONAL_INTERFACE_ONLY','noncore budget passes; FCSD time/success upper at reduction envelope unavailable')):
            coverage.append(dict(candidate=c['name'],method=method,status=status,reason=why))
        flat.append(dict(name=c['name'],p=c['p'],n=c['n'],k=c['k'],w=c['weight'],
            extra_C_solution_bound_log2=s['exact']['log2_extra_witness_probability_upper'],
            LMMT_optimistic_materialized_MODEL=s['joint_summary']['optimistic_model_min_log2'],
            LMMT_success_lower_positive_tiled_points=0,
            LMMT_exact_subevent_closed_schedule=c['name'] in closed,
            LMMT_closed_schedule_declared_bit_upper_log2=closed.get(c['name'],{}).get('best_declared_bit_charge_upper_log2'),
            PGE_Wagner_declared_bit_budget_log2=p['ordinary']['log2_abstract_bit_charge_upper_model'],
            public_sum_declared_bit_budget_log2=p['public_sum']['log2_abstract_bit_charge_upper_model'],
            restricted_evaluated_points=r['restricted_summary']['evaluated_points'],
            restricted_counted_half_success_points=r['restricted_summary']['finite_counted_half_success_certificates'],
            full_SVP_default_model_operations_log2=L['ordinary_full_SVP_costs'][0]['log2_model_operations'],
            full_SVP_payload_bytes_log2=L['ordinary_full_SVP_costs'][0]['payload']['log2_coordinate_payload_bytes_model'],
            log2_reduction_charge=r['authentication']['log2_reduction_charge'],
            unrestricted_classical128='DO_NOT_RECOMMEND_MODEL_RISK_AND_MISSING_COVERAGE',
            memory_limited_classical128='NOT_ESTABLISHED',
            complete_authentication128='NOT_ESTABLISHED',
            demonstrated_concrete_break=False))
    write(HERE/'results/decision_table.json',flat)
    for name,data in (('decision_table',flat),('coverage',coverage)):
        with (HERE/'results'/f'{name}.csv').open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    write(HERE/'results/coverage.json',coverage)
    lines=['# 两组现有候选的最终安全评估 v3','',
        '范围：只复核 wide700_n320 与 low370_n448；不新增C或扩大候选长度，不执行目标攻击或协议计时。', '',
        '**结论：两组均不批准为完整经典128位认证参数。** 无内存限制模型存在已发表格费用风险；64 GiB限定模型仍缺同一调度的完整成功/费用与攻击覆盖。此处不是实际破解证书。', '',
        '上一轮将这两组继续称为优先定案候选，遗漏了已保存的格风险评估；本轮纠正并重新绑定。', '',
        '## 统一判定表','',
        '|候选|p,n,k,w|PGE/Wagner声明bit账本上预算 log2|公开和全类型账本上预算 log2|全维SVP默认模型操作 log2|归约账本 log2|最终认证|',
        '|---|---|---:|---:|---:|---:|---|']
    for r in flat:
        lines.append(f"|{r['name']}|{r['p']},{r['n']},{r['k']},{r['w']}|{r['PGE_Wagner_declared_bit_budget_log2']:.3f}|{r['public_sum_declared_bit_budget_log2']:.3f}|{r['full_SVP_default_model_operations_log2']:.3f}|{r['log2_reduction_charge']:.6f}|不批准|")
    lines+=['','这些列不作跨单位取最小值。PGE/Wagner列是声明bit模型中的构造上预算，不是安全下界；SVP列没有假定模型操作等于位操作。','',
        '## 三项工作对应结果','',
        '1. **固定C完整攻击费用：指定调度闭合，竞争估计仍不足。** 实际C的PGE/Wagner构造调度与全部公开和类型费用已保留；全部命名资源定价，不只换算域乘法。该费用转换是明确声明的估计模型，不是对Python/CPU的形式化编译。原LMMT数值检查点仍没有正下界；本轮另用精确公共组成事件并集和逐矩阵GL存活下界闭合具体分块调度，未将原零下界改写为正，也不以这个慢构造的费用上预算证明安全。',
        '2. **统一攻击面表：已形成，包含未覆盖项。** 新计算同一C下的Restricted-Balls固定自由位置及去重向量域；保留其表示子族成功下界与完整攻击费用边界。额外重算64GiB格枚举风险及已发表高内存SVP模型；高级汉明及重复多重集PKP仍未得到足够竞争估计，不能写成全部攻击已覆盖。',
        '3. **协议预算：精确账本已代入，底层困难性前提未成立。** 256轮、SHA-384、256位盐、2^20学习会话、一次最终冒充；t_B的实际log2约164.330917。继续保留2^166充分接口。必须约束FCSD在该资源内的成功概率≤1/4，不能仅凭半成功费用超过某指数替代时间—成功率曲线。', '',
        '## 格风险的条件含义','',
        '普通公开格 Λ={(z,t):Hz=t*s mod p} 由H,s构造。既有精确球计数给出唯一最短目标的分离失败界；若精确SVP以概率q成功，则恢复合法C秘密概率至少q−δ，不假设两事件独立。',
        '默认全维SVP模型分别约2^110.132和2^147.508模型操作。BDGL16的渐近常数不构成目标维数的完整位费用定理；较慢分支是小维拟合的大维敏感性，不是可任意选择的安全预测。',
        '本轮保留0/8/16/24/32指数位的单位换算敏感性与q=3/4的明确条件场景。高内存路径超64GiB，不能用它宣称64GiB内已破解；也不能以本机内存不足排除一般经典攻击者。', '',
        '## 可用于论文的结论','',
        '可写：协议修复与条件三性质证明成立于列明假设；已有组成的有限攻击评估与同参数实现优化已形成；本次两组候选因格风险和未闭合的攻击覆盖，不作为最终经典128位认证参数。',
        '不可写：已证明攻击者找不到解、已覆盖所有可能攻击、两个候选都已被实际破解、固定C模型数字就是方案安全位数。', '',
        '## 交付与终点','',
        'results/decision_table.csv：逐候选判定；results/coverage.csv：逐方法覆盖；results/*：概率、费用、Restricted-Balls、格与预算证据；BIT_LEDGER.md：单位映射；SOURCE_BINDINGS.json：源码和输入绑定。',
        '本轮完成的是两组候选的最终评估和不批准判定。LMMT存活下界、完整统一位费用及全部适用攻击覆盖尚未全部补齐；这些缺口不因报告定稿而变成已完成。', '']
    if closed:
        lines+=['## 新增：精确组成子事件的可执行调度闭合','',
            '原样本仅用于选择一个公开、有限的组成事件集合；随后对该集合重新做精确整数计数。这里没有统计误差，也没有将抽样均值当成真概率。成本仍覆盖全部叶壳、全部连接和根，分析子事件不缩减实际攻击列表。', '',
            '|候选|公开组成事件数|可行布局|成功率下界|命名载荷GiB|声明bit账本上预算 log2|',
            '|---|---:|---:|---|---:|---:|']
        for name,r in closed.items():
            lines.append(f"|{name}|{r['public_profiles']}|{r['feasible_layouts']}|{r['success_lower']}|{r['payload_GiB']:.3f}|{r['best_declared_bit_charge_upper_log2']:.3f}|")
        lines+=['','该表是某个实际完整攻击调度的成功/命名费用证书；bit列仍是显式转换模型。它不证明一般攻击至少需要如此多工作，也不替代完整表示集合的sigma或竞争攻击优化。详见 LMMT_CLOSURE_PROOF.md。','']
    (HERE/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    return flat


def main():
    (HERE/'results').mkdir(parents=True,exist_ok=True)
    sourcepaths=[RESEARCH/'stern_classical128_assessment_v2/evaluate.py',
        RESEARCH/'stern_lee_parameter_decision_v1/assess.py',RESEARCH/'stern_classical128_assessment_v2/results.json',
        RESEARCH/'stern_concrete_budget_v1/budget.py',RESEARCH/'stern_lattice_64gib_v1/costs.py',
        RESEARCH/'stern_restricted_set_scalable_v1/analyzer.py',RESEARCH/'stern_restricted_set_scalable_v1/domain.py',
        RESEARCH/'stern_restricted_set_scalable_v1/moments.py',RESEARCH/'stern_restricted_set_scalable_v1/list_costs.py',
        RESEARCH/'stern_restricted_bjmm_v3/analyzer.py',RESEARCH/'stern_restricted_bjmm_v3/counting.py',
        RESEARCH/'stern_restricted_bjmm_v3/rank_model.py',RESEARCH/'stern_restricted_bjmm_v3/compositions.py',
        RESEARCH/'stern_classical128_assessment_v2/sources/reduction_pinned.py']
    bindings={str(p):sha(p) for p in sourcepaths}
    old=module('final_existing_assessor',RESEARCH/'stern_lee_parameter_decision_v1/assess.py')
    extended=module('final_existing_lattice',RESEARCH/'stern_classical128_assessment_v2/evaluate.py')
    oldrecords=read(RESEARCH/'stern_classical128_assessment_v2/results.json')
    records=[]
    for name in NAMES:
        started=time.perf_counter()
        detailpath=RESEARCH/'stern_formula_candidate_screen_v1/results/details'/f'{name}.json'
        detail=read(detailpath);cfg=detail['candidate'];bindings[str(detailpath)]=sha(detailpath)
        saved=next(r for r in oldrecords if r['candidate']['name']==name)
        if not same(cfg,saved['candidate']):raise ValueError('full C differs')
        write(HERE/'results'/f'{name}_C.json',cfg)
        # Actual-C probability and layout evidence is re-bound, never extrapolated.
        td=read(RESEARCH/'stern_formula_candidate_screen_v1/results'/f'{name}_tiled.json')
        if any(r['log2_single_ensemble_success_lower'] is not None for r in td['layouts']):
            raise ValueError('expected zero-lower gap no longer current; revise final decision')
        write(HERE/'results'/f'{name}_LMMT.json',dict(joint_evidence=detail['joint_evidence'],
            tiled_evidence=td['source'],same_C_verified=True,feasible_named_layouts=td['feasible_layouts'],
            model_summary=detail['joint_summary'],half_success_work=None,
            status='NO_POSITIVE_JOINT_SUCCESS_LOWER_ON_CHECKED_POINTS',
            cannot_invert_probability_upper_into_constructive_work=True))
        ds=read(detail['disjoint_evidence']['path'])
        ordinary=price_pge_schedule(cfg,ds['ordinary_winner'])
        parts=[price_pge_schedule(cfg,r) for r in ds['public_sum_components']]
        total=sum(int(r['exact_charge']) for r in parts)
        public=dict(log2_abstract_bit_charge_upper_model=lg(total),exact_charge=str(total),
            all_types=len(parts),components=parts,
            success_lower=ds['public_sum_aggregate']['ensemble_success_lower'],
            costs_summed_not_weighted=True,
            payload=ds['public_sum_aggregate']['memory_bytes'],
            security_lower_bound=False)
        pge=dict(ordinary=ordinary,public_sum=public,unit='declared bit-charge estimate, not certified real CPU bit count')
        write(HERE/'results'/f'{name}_PGE_Wagner.json',pge)
        rb=restricted(cfg);write(HERE/'results'/f'{name}_restricted.json',rb)
        # Recompute exact ordinary separation and revised authentication ledger.
        separation=extended.separation(cfg)
        if separation!=saved['exact_ordinary_separation']:raise ValueError('old separation changed')
        budget=extended.patched_budget(old,cfg)
        if budget['reduction_charge_exact']!=saved['authentication_budget']['reduction_charge_exact']:
            raise ValueError('authentication resource mismatch')
        lattice64=old.lattice_panel(cfg)
        write(HERE/'results'/f'{name}_lattice64.json',lattice64)
        auth=auth_curve(budget,separation);write(HERE/'results'/f'{name}_authentication.json',auth)
        record=dict(candidate=cfg,screen=detail,pge_wagner=pge,restricted_summary=rb['summary'],
            lattice=saved,lattice64=lattice64,authentication=auth,
            decision=dict(approved_classical128=False,actual_break_demonstrated=False,
                high_memory_model_risk=True,bounded_memory_final_status='NOT_ESTABLISHED',
                all_three_theoretical_gaps_closed=False),elapsed_seconds=time.perf_counter()-started)
        records.append(record);write(HERE/'results/records.json',records)
        print(json.dumps(dict(done=name,seconds=record['elapsed_seconds'],restricted=rb['summary'],
            ordinary_bit_charge=pge['ordinary']['log2_abstract_bit_charge_upper_model'],
            public_sum_bit_charge=pge['public_sum']['log2_abstract_bit_charge_upper_model'],
            final_approved=False)),flush=True)
    assert all(sha(p)==h for p,h in bindings.items())
    write(HERE/'SOURCE_BINDINGS.json',bindings)
    flat=report(records)
    ledger='''# 声明的完整命名资源bit账本

本轮逐项定价已有PGE/Wagner调度的全部命名资源向量，而非只给域乘法乘log²p。
域乘法、加减、逆、坐标比较、直方图、随机位、均匀采样、排序移动、读写字节、秩大整数运算、树访问、控制和初始化均保留。

具体价格由assess.py::resource_prices给出，每个结果保存实际价格、操作数位长及完整向量。对未知资源拒绝定价。
按学校式多精度模型：乘法O(B²)，除法O(B²)，加减比较O(B)；域元位长b、秩位长B与长循环计数/地址位长单独计算。
每类操作附显式控制费用64(a+B)，排序移动按原记录最大字段宽度收费，字节读写分别收费；一次性设置已经在源total_caps中。

这是明确的估计模型，不是Python/CPU指令覆盖或位复杂度编译证明。每类截止上预算收费的和仍是某构造调度的上预算，不能成为一般攻击安全下界。
公开和全类型成本相加，顺序复用内存取最大，只有成功概率按类型权重加权。

LMMT只有每提案期望费和零的数值成功下界，不能套此价格制造完整半成功费。
Restricted-Balls的来源包含组合数调用、DP/排名和不同集合；未经过同等覆盖证明的费用仍留在原命名单位，不假称已统一。
格节点周期/模型操作与位费之间无已校准等式，保留0/8/16/24/32指数位场景，不任意设为1。
'''
    (HERE/'BIT_LEDGER.md').write_text(ledger,encoding='utf-8')
    print(json.dumps(flat,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
