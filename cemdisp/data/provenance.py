"""多井输入数据来源口径登记与摘要构造。

本模块在原有“注入流体来源口径”基础上，扩展到单井几何、地面施工程序、
鞋口边界同步三个层面。这里仅记录来源与口径，不做数据校验，也不输出质量
代理值，避免把 provenance 与 validation / quality proxy 混用。
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule


@dataclass(frozen=True)
class FluidProvenance:
    """单种流体的来源口径说明。"""

    status: str
    note: str


@dataclass(frozen=True)
class SectionProvenance:
    """单个非流体输入段的来源口径说明。"""

    status: str
    note: str


@dataclass(frozen=True)
class WellProvenance:
    """单口井的聚合来源口径，覆盖流体、几何、程序和同步边界。"""

    well_name: str
    fluid: Mapping[str, FluidProvenance]
    geometry: SectionProvenance
    program: SectionProvenance
    sync: SectionProvenance

    @property
    def sections(self) -> Mapping[str, object]:
        """按固定段名返回聚合来源，供报告摘要统一遍历。"""

        return {
            "fluid": self.fluid,
            "geometry": self.geometry,
            "program": self.program,
            "sync": self.sync,
        }


_STATUS_LABELS = {
    "field": "现场符合",
    "partial": "部分符合",
    "proxy": "代理/暂定",
    "unknown": "未标定",
    "design": "设计值",
    "optimized_input": "优化输入",
    "model_assumption": "模型假设",
}


WELL_FLUID_PROVENANCE: dict[str, dict[str, FluidProvenance]] = {
    "呼101": {
        "平衡液": FluidProvenance("partial", "密度1.85现场抽取；流变为 0708 设计六速口径（2011111 §1.4.2，60℃ 六速 100/59/43/26/6/5 → PV41/YP9.2）；LEGACY(2026-08-29 前) 代理 0.030/3.0。"),
        "驱油隔离液": FluidProvenance("field", "密度与流变来自现场抽取资料。"),
        "领浆": FluidProvenance("field", "密度与幂律参数来自现场实验（主检 2011122.pdf W301-22094，89℃）；另一独立复检报告 2011121.doc（93℃）为双报告并存，建议论文引用前做双口径敏感性重算。"),
        "尾浆": FluidProvenance("field", "密度与幂律参数来自现场实验（主检 2011122.pdf W301-22094，89℃）；另一独立复检报告 2011121.doc（93℃）为双报告并存，建议论文引用前做双口径敏感性重算。"),
        "轻泥浆": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
        "中置液": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
        "井浆": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
    },
    "呼102": {
        "钻井液": FluidProvenance("field", "密度2.02g/cm³、PV=66mPa·s、YP=10Pa 现场实测（20211.doc s1.2，油基）。"),
        "替浆液": FluidProvenance("field", "密度2.02g/cm³、PV=66mPa·s、YP=10Pa 现场实测（替浆用油基钻井液，20211.doc）。"),
        "平衡液": FluidProvenance("field", "密度1.90g/cm³现场实测（20211.doc）；流变实测 20213.doc §1.3（1.90@60℃ 六速 102/60/37/19/3/2 → PV42/YP18）；LEGACY(2026-08-29 前) 邻井代理 0.025/1.5。"),
        "隔离液": FluidProvenance("field", "密度2.05g/cm³现场实测（20211.doc）；塑粘 PV30 实测（20216）；YP=8.0 保留但无实测（20216 未给），仅密度与 PV 实测。"),
        "领浆": FluidProvenance("field", "密度2.10g/cm³施工值（20211.doc），幂律 n=0.737/K=0.947 现场实测（20234.doc）。"),
        "尾管水泥浆": FluidProvenance("field", "密度2.10g/cm³施工值（20211.doc），幂律 n=0.737/K=0.947 现场实测（20234.doc）。"),
        "压塞液": FluidProvenance("partial", "密度2.10g/cm³现场实测（20211.doc）；流变无实测，复用钻井液（proxy）。"),
        "后置液": FluidProvenance("field", "密度1.90g/cm³现场实测（20211.doc）；流变实测 PV30/YP8（20216）。"),
        "冲洗液": FluidProvenance("proxy", "LEGACY(2026-08 前)邻井代理值（呼探1-002/呼103），保留用于敏感性 FLUSHER 序列验证，不作为呼102现场实录。"),
    },
    "呼103": {
        "钻井液": FluidProvenance("field", "密度1.98g/cm³、PV=54mPa·s、YP=10Pa 现场实测（20313.doc 1.3.1，65C 白油基）。"),
        "平衡液": FluidProvenance("field", "密度1.88g/cm³、PV=25mPa·s、YP=1.5Pa 现场实测（20314.doc 施工记录表）；角色=mud（轻泥浆）。"),
        "隔离液1": FluidProvenance("partial", "设计密度2.00g/cm³（20313.doc 7.1.2）；流变实测 POWER_LAW n=0.54/K=2.12（203111 表7，140↘93℃）；LEGACY(2026-08-29 前) Bingham 代理 PV35/YP8。"),
        "隔离液2": FluidProvenance("partial", "设计密度1.95g/cm³（20313.doc 7.1.2）；流变实测 POWER_LAW n=0.54/K=2.12（203111 表7，140↘93℃）；LEGACY(2026-08-29 前) Bingham 代理 PV35/YP8。"),
        "隔离液(实际)": FluidProvenance("field", "实际密度1.92g/cm³（20314.doc）；流变实测 POWER_LAW n=0.54/K=2.12（203111 表7）；LEGACY(2026-08-29 前) Bingham 代理 PV35/YP8。"),
        "领浆": FluidProvenance("field", "现场密度2.05g/cm³与幂律参数n=0.82,K=0.67，来自203111.docx实验报告（A级）。"),
        "中间浆": FluidProvenance("field", "密度2.05g/cm³、n=0.76、K=1.11来自实验报告（203111.docx）；体积按设计35m³（0708 多数口径 4 处：7.1.1 用量表总量 83=22+35+26、7.2 流程、实际 35；LEGACY(2026-08-29 前) 36 系浆柱表孤值）。"),
        "尾浆": FluidProvenance("field", "现场密度2.05g/cm³与幂律参数n=0.76,K=1.14，来自203111.docx实验报告（A级）。"),
        "压塞液": FluidProvenance("partial", "密度2.10g/cm³设计值（20313.doc 7.1.3）；流变复用替浆液。"),
        "中置液": FluidProvenance("partial", "密度1.95g/cm³设计值（20313.doc 7.1.3）；流变复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度1.98g/cm³按现场钻井液；流变复用固井时钻井液。"),
    },
    "呼1-003": {
        "钻井液": FluidProvenance("field", "密度1.95g/cm³现场实测；幂律 n=0.631/K=0.751 化验报告（129C↓93C，六速122/71/52/21/7/6）。"),
        "平衡液": FluidProvenance("partial", "密度1.75g/cm³现场实测；流变设计初稿 PV55/YP9.2（无实测，design_value）。"),
        "隔离液1": FluidProvenance("field", "密度2.05g/cm³现场化验；幂律 n=0.668/K=1.245 化验报告（与隔离液2同流变）。"),
        "隔离液2": FluidProvenance("field", "密度1.95g/cm³现场化验；幂律 n=0.668/K=1.245 化验报告。"),
        "领浆": FluidProvenance("field", "密度2.05g/cm³现场；幂律 n=0.597/K=1.622 化验报告。"),
        "尾浆": FluidProvenance("field", "密度1.95g/cm³现场化验实测；幂律 n=0.585/K=1.673 化验报告。"),
        "压塞液": FluidProvenance("partial", "密度1.95g/cm³设计值；流变设计 PV30/YP8（design_value）。"),
        "保护液": FluidProvenance("partial", "密度1.95g/cm³设计值；流变设计 PV30/YP8.2（design_value）。"),
        "基液": FluidProvenance("partial", "密度1.02g/cm³设计值；流变设计 PV30/YP8（design_value）。"),
        "替钻井液": FluidProvenance("field", "密度1.95g/cm³、塑粘62mPa·s 施工记录表替浆记录。"),
        "井浆": FluidProvenance("partial", "密度1.95g/cm³现场口径；流变沿用代理 PV30/YP9.3（无实测）。"),
    },
    "呼1-004": {
        "钻井液": FluidProvenance("field", "密度1.90g/cm³现场（固井前）；PV53/YP8.5 常温65C实测。"),
        "先导浆": FluidProvenance("field", "密度1.75g/cm³现场实测；流变 PV58/YP9.8 为优化输入（optimized_input，非实测）。"),
        "隔离液1": FluidProvenance("field", "密度1.95g/cm³现场实测；流变 PV58/YP9.8 为优化输入（化验幂律 n=0.611/K=0.736）。"),
        "隔离液2": FluidProvenance("field", "密度1.75g/cm³现场实测；流变 PV65/YP10 为优化输入（化验幂律同隔离液1）。"),
        "领浆": FluidProvenance("field", "密度1.93g/cm³现场实测；流变 PV170/YP13 为优化输入（化验幂律 n=0.853/K=0.746）。"),
        "尾浆": FluidProvenance("field", "密度1.90g/cm³现场实测；流变 PV180/YP14 为优化输入（化验幂律 n=0.869/K=0.669）。"),
        "压塞液": FluidProvenance("field", "密度1.70g/cm³实际后置液（施工记录表）。"),
        "保护液": FluidProvenance("field", "密度1.90g/cm³现场实测（设计/作业史）。"),
        "替钻井液": FluidProvenance("field", "密度1.90g/cm³、塑粘50mPa·s 施工记录表替浆记录。"),
        "基液": FluidProvenance("field", "密度1.02g/cm³设计值（与现场一致）。"),
        "井浆": FluidProvenance("field", "密度1.90g/cm³现场（替浆段井浆）。"),
        "冲洗液（FLUSHER）": FluidProvenance("model_assumption", "合成流体（ρ1880/5m³），现场无独立冲洗液；仅验证序列可表达，默认不注入。"),
    },
    "呼探1": {
        "钻井液": FluidProvenance("field", "密度2.08g/cm³施工实测（204121'施工前由 2.12 降至 2.08'，LEGACY 2.12 为下套管灌浆/设计口径）；PV=76mPa·s、YP=8Pa 实测（204121 洗井结束前，六速 152/86/62/36/6/5）；LEGACY(2026-08-29 前) 66/10 为 204111 设计口径。"),
        "平衡液(先导泥浆)": FluidProvenance("partial", "密度2.00g/cm³现场（204151.doc，旧名冲洗液/冲洗型隔离液；抗污染油基泥浆型）；流变复用施工泥浆实测 PV76/YP8（204121，2026-08-29 校准）。"),
        "驱油隔离液": FluidProvenance("field", "密度2.05g/cm³现场（204151.doc，旧名加重隔离液）；流变实测 POWER_LAW n=0.649/K=0.585（204131，2026-08-29 校准）。"),
        "领浆": FluidProvenance("field", "密度2.10g/cm³现场（204151.doc/204111.doc）；幂律 n=0.732/K=0.933 现场实测（204131 检测报告，93℃，2020-10-24 签发）；LEGACY(2026-08-29 前) proxy 呼102 0.737/0.947。"),
        "尾浆": FluidProvenance("field", "密度1.90g/cm³现场（204151.doc/204111.doc）；幂律 n=0.666/K=0.906 现场实测（204131，93℃，2020-10-24 签发）；LEGACY(2026-08-29 前) proxy 呼1-001 0.886/0.453。"),
        "后置液": FluidProvenance("partial", "密度2.05g/cm³现场（204151.doc 后置液2m³）；流变无实测，沿用 proxy（PV35/YP8 旧隔离液代理口径）。"),
        "替浆液1": FluidProvenance("partial", "密度1.95g/cm³现场（204151.doc 管串36m³）；流变无实测，复用钻井液口径（proxy；LEGACY 沿 204111 固井前 66/10 口径，204121 施工实测 76/8 未联动）。"),
        "替浆液2": FluidProvenance("partial", "密度2.08g/cm³现场（204151.doc 环空49m³）；流变无实测，复用钻井液口径（proxy；LEGACY 沿 204111 固井前 66/10 口径，204121 施工实测 76/8 未联动）。"),
    },
    "呼探1-002": {
        "平衡液": FluidProvenance("partial", "密度来自现场（设计 1.85；批混实测 1.87，数采平台 1.868–1.872），流变缺实测，使用代理值维持（2026-08-29 校准注记）。"),
        "隔离液": FluidProvenance("field", "密度来自现场；流变 POWER_LAW n=0.545/K=1.338 化验实测（2026-08-29 校准）；LEGACY(2026-08-29 前) Bingham 30/8（YP8 无 0708 直接出处）。"),
        "领浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "中间浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "尾浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "压塞液": FluidProvenance("partial", "密度1.90g/cm³已见现场作业史，流变仍复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变复用钻井液。"),
        "中置液": FluidProvenance("partial", "密度1.90g/cm³已见现场作业史，流变仍复用替浆液。"),
        "井浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
    },
    "呼探1-001": {
        "平衡液": FluidProvenance("partial", "密度来自现场，流变缺项按代理值补齐。"),
        "隔离液": FluidProvenance("field", "现场密度与独立流变参数，来自实验报告。"),
        "领浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "中间浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "尾浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "压塞液": FluidProvenance("partial", "流变复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变缺项按代理值补齐。"),
        "替钻井液": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "中置液": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "基液": FluidProvenance("partial", "密度独立，流变复用替浆液。"),
        "井浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
    },
}

_WELL_NAME_ALIASES = {
    "呼探1井": "呼探1",
    "呼探1-001井（HT1-001）": "呼探1-001",
    "呼探1-002井（HT1-002）": "呼探1-002",
    "呼1-003井（HT1-003）": "呼1-003",
    "呼1-004井（HT1-004）": "呼1-004",
}


_DEFAULT_GEOMETRY_PROVENANCE = SectionProvenance(
    "partial",
    "井深、鞋深和双径向字段已进入 WellSpec 口径；仍需逐井补齐原始资料追溯。",
)
_DEFAULT_PROGRAM_PROVENANCE = SectionProvenance(
    "partial",
    "施工步骤、体积、排量和阶段标签已进入 PumpingSchedule 口径；现场时序仍需逐井核对。",
)
_DEFAULT_SYNC_PROVENANCE = SectionProvenance(
    "partial",
    "鞋口延迟和边界同步字段已有登记口径；暂不在本任务接入 bridge/casing_flow。",
)


WELL_PROVENANCE: dict[str, WellProvenance] = {
    well_name: WellProvenance(
        well_name=well_name,
        fluid=fluid_provenance,
        geometry=_DEFAULT_GEOMETRY_PROVENANCE,
        program=_DEFAULT_PROGRAM_PROVENANCE,
        sync=_DEFAULT_SYNC_PROVENANCE,
    )
    for well_name, fluid_provenance in WELL_FLUID_PROVENANCE.items()
}


def _section_to_summary(section: SectionProvenance) -> dict[str, str]:
    """把非流体来源段转换为中文摘要字段。"""

    return {
        "来源口径": _STATUS_LABELS.get(section.status, section.status),
        "说明": section.note,
    }


def _fluid_registry_to_summary(registry: Mapping[str, FluidProvenance]) -> dict[str, object]:
    """把单井流体登记表转换为不依赖施工实例的摘要。"""

    details = [
        {
            "流体名称": fluid_name,
            "角色": "unknown",
            "来源口径": _STATUS_LABELS.get(provenance.status, provenance.status),
            "说明": provenance.note,
        }
        for fluid_name, provenance in registry.items()
    ]
    non_field = [item for item in details if item["来源口径"] != _STATUS_LABELS["field"]]
    return {
        "是否全部现场符合": len(non_field) == 0,
        "注入流体总数": len(details),
        "非现场完全符合流体数": len(non_field),
        "明细": details,
    }


def build_injected_fluid_provenance_summary(
    well_name: str,
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
) -> dict[str, object]:
    """返回施工注入流体的现场来源口径汇总。"""

    fluid_by_name = {fluid.name: fluid for fluid in fluids}
    canonical_well_name = _WELL_NAME_ALIASES.get(well_name, well_name)
    registered = WELL_FLUID_PROVENANCE.get(canonical_well_name, {})
    scheduled_names: list[str] = []
    for step in schedule.steps:
        if step.fluid_name not in scheduled_names:
            scheduled_names.append(step.fluid_name)

    details: list[dict[str, str]] = []
    for fluid_name in scheduled_names:
        fluid = fluid_by_name.get(fluid_name)
        provenance = registered.get(fluid_name, FluidProvenance("unknown", "未建立明确来源口径。"))
        details.append(
            {
                "流体名称": fluid_name,
                "角色": fluid.role.value if fluid is not None else "unknown",
                "来源口径": _STATUS_LABELS.get(provenance.status, provenance.status),
                "说明": provenance.note,
            }
        )

    all_field = all(item["来源口径"] == _STATUS_LABELS["field"] for item in details)
    non_field = [item for item in details if item["来源口径"] != _STATUS_LABELS["field"]]
    return {
        "是否全部现场符合": all_field,
        "注入流体总数": len(details),
        "非现场完全符合流体数": len(non_field),
        "明细": details,
    }


def format_injected_fluid_provenance_markdown(summary: Mapping[str, object]) -> list[str]:
    """把来源口径汇总渲染为 Markdown 摘要行。"""

    total_count = int(cast(int, summary["注入流体总数"]))
    non_field_count = int(cast(int, summary["非现场完全符合流体数"]))
    details = cast(list[dict[str, str]], summary["明细"])

    lines = [
        "## 注入流体现场符合性检查",
        "",
        f"- 是否全部现场符合：{'是' if bool(summary['是否全部现场符合']) else '否'}",
        f"- 注入流体总数：{total_count}",
        f"- 非现场完全符合流体数：{non_field_count}",
        "",
    ]
    for detail in details:
        lines.append(
            f"- {detail.get('流体名称', '未知流体')}：{detail.get('来源口径', '未标定')}；{detail.get('说明', '')}"
        )
    return lines


__all__ = [
    "FluidProvenance",
    "SectionProvenance",
    "WellProvenance",
    "WELL_FLUID_PROVENANCE",
    "WELL_PROVENANCE",
    "build_injected_fluid_provenance_summary",
    "format_injected_fluid_provenance_markdown",
]
