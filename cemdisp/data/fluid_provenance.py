"""注入流体现场来源口径登记与检查。

本模块把目前 six-well loaders 中已经核实过的流体来源口径显式登记出来，
用于在 runner 导出结果时同步给出“现场符合 / 部分符合 / 代理 / 未知”判断，
避免在流体数据尚不完整时对模型结果产生过度解释。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, cast

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule


@dataclass(frozen=True)
class FluidProvenance:
    """单种流体的来源口径说明。"""

    status: str
    note: str


_STATUS_LABELS = {
    "field": "现场符合",
    "partial": "部分符合",
    "proxy": "代理/暂定",
    "unknown": "未标定",
}


WELL_FLUID_PROVENANCE: dict[str, dict[str, FluidProvenance]] = {
    "呼101": {
        "平衡液": FluidProvenance("field", "密度与流变来自现场抽取资料。"),
        "驱油隔离液": FluidProvenance("field", "密度与流变来自现场抽取资料。"),
        "领浆": FluidProvenance("field", "密度与幂律参数来自现场实验/抽取资料。"),
        "尾浆": FluidProvenance("field", "密度与幂律参数来自现场实验/抽取资料。"),
        "轻泥浆": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
        "中置液": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
        "井浆": FluidProvenance("partial", "按现场阶段名称建模，但性质复用钻井液。"),
    },
    "呼102": {
        "冲洗液": FluidProvenance("proxy", "邻井代理值。"),
        "隔离液": FluidProvenance("partial", "密度 1.85g/cm³ 来自汇总表，B级；流变参数参考邻井。"),
        "尾管水泥浆": FluidProvenance("partial", "现场仅确认平均浆密度 2.10g/cm³；流变 n/K 仍沿用 legacy 占位值。"),
        "替浆液": FluidProvenance("field", "密度 2.02g/cm³ 来自 10042.xlsx，A级确认；流变参数为文献暂定。"),
    },
    "呼103": {
        "平衡液": FluidProvenance("partial", "密度来自现场，流变缺实测，按冲洗液体系代理。"),
        "驱油冲洗液": FluidProvenance("partial", "密度 1.02g/cm³ 来自现场，流变采用代理值。"),
        "隔离液": FluidProvenance("field", "密度 1.80g/cm³、流变 n=0.587/K=0.880 来自 203211.docx 实验报告，A级提取；冲洗效率 97.3%。"),
        "领浆": FluidProvenance("field", "密度 1.90g/cm³、流变 n=0.838/K=0.587 来自 203211.docx 实验报告，A级提取；稠化 493min，API失水 42mL，24h强度 11.5MPa。"),
        "中间浆": FluidProvenance("partial", "中间浆体积/流变含代理成分。"),
        "尾浆": FluidProvenance("field", "密度 1.90g/cm³、流变 n=0.766/K=1.093 来自 203211.docx 实验报告，A级提取；稠化 295min，API失水 35mL，24h强度 14.7MPa。"),
        "压塞液": FluidProvenance("partial", "压塞液流变复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变代理为固井时钻井液一致。"),
    },
    "呼探1": {
        "平衡液": FluidProvenance("partial", "密度来自现场，流变缺项按代理值补齐。"),
        "驱油隔离液": FluidProvenance("partial", "密度来自现场，流变缺项按代理值补齐。"),
        "领浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "中间浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "尾浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "压塞液": FluidProvenance("partial", "压塞液性质复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变参考类似井。"),
        "轻泥浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "中置液": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "井浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
    },
    "呼探1-002": {
        "平衡液": FluidProvenance("partial", "密度来自现场，流变缺实测，使用代理值。"),
        "隔离液": FluidProvenance("partial", "密度来自现场，流变缺实测，使用代理值。"),
        "领浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "中间浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "尾浆": FluidProvenance("field", "现场密度与幂律参数。"),
        "压塞液": FluidProvenance("partial", "密度 1.90g/cm³ 已见现场作业史，流变仍复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变复用钻井液。"),
        "中置液": FluidProvenance("partial", "密度 1.90g/cm³ 已见现场作业史，流变仍复用替浆液。"),
        "井浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
    },
    "呼探1-001": {
        "平衡液": FluidProvenance("field", "密度 1.75g/cm³ 来自技术总结，A级确认；流变参数来自化验报告。"),
        "隔离液": FluidProvenance("field", "现场密度与独立流变参数。"),
        "领浆": FluidProvenance("field", "密度 2.05g/cm³、流变 n=0.825/K=0.842 来自 HT1-001 化验报告，A级确认；稠化 432min，24h强度 4.5MPa。"),
        "中间浆": FluidProvenance("field", "密度 1.90g/cm³、流变 n=0.816/K=0.512 来自 HT1-001 化验报告，A级确认；稠化 302min，24h强度 15.2MPa。"),
        "尾浆": FluidProvenance("field", "密度 1.90g/cm³、流变 n=0.807/K=0.527 来自 HT1-001 化验报告，A级确认；稠化 223min，24h强度 17.7MPa。"),
        "压塞液": FluidProvenance("field", "密度 1.60g/cm³ 来自作业史，A级确认；体积 5m³。"),
        "替浆液": FluidProvenance("field", "密度 1.92g/cm³ 来自技术总结，A级确认；总量 88.8m³。"),
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
    "WELL_FLUID_PROVENANCE",
    "build_injected_fluid_provenance_summary",
    "format_injected_fluid_provenance_markdown",
]
