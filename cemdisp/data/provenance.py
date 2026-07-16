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
        "冲洗液": FluidProvenance("proxy", "邻井代理值（呼探1-002）。提取表显示平衡液密度1.88g/cm³，B级汇总表，待实测验证。"),
        "隔离液": FluidProvenance("partial", "原邻井代理值；提取表显示驱油隔离液密度1.85g/cm³、n=0.587、K=0.880，B级汇总表交叉校验，已升级。"),
        "尾管水泥浆": FluidProvenance("partial", "现场确认密度2.10g/cm³（主作业）；提取表显示尾浆密度1.90g/cm³、n=0.766、K=1.093（B级，呼103邻井类比）。当前loader沿用2.10g/cm³主作业值，待进一步核实。"),
        "替浆液": FluidProvenance("field", "密度2.02g/cm³来自10042.xlsx主作业记录（A级）；流变PV=80mPa·s、YP=15Pa为文献暂定值。"),
        "先导浆": FluidProvenance("partial", "提取表显示密度1.60g/cm³，体积52m³，B级汇总表。当前loader未单独定义，待回写。"),
        "领浆": FluidProvenance("partial", "提取表显示密度1.90g/cm³、n=0.838、K=0.587，B级（呼103邻井类比）。当前loader未单独定义，待回写。"),
        "尾浆": FluidProvenance("partial", "提取表显示密度1.90g/cm³、n=0.766、K=1.093，B级（呼103邻井类比）。当前loader未单独定义，待回写。"),
    },
    "呼103": {
        "平衡液": FluidProvenance("partial", "密度来自现场（203211.docx），流变缺实测，按冲洗液体系代理。"),
        "驱油冲洗液": FluidProvenance("partial", "密度来自现场，流变采用代理值。"),
        "隔离液": FluidProvenance("partial", "密度/流变并非全部现场实测。"),
        "领浆": FluidProvenance("field", "现场密度2.05g/cm³与幂律参数n=0.82,K=0.67，来自203211.docx实验报告（A级）。"),
        "中间浆": FluidProvenance("partial", "中间浆体积10m³为代理值；密度2.05g/cm³、n=0.76、K=1.11来自实验报告。"),
        "尾浆": FluidProvenance("field", "现场密度2.05g/cm³与幂律参数n=0.76,K=1.14，来自203211.docx实验报告（A级）。"),
        "压塞液": FluidProvenance("partial", "压塞液流变复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变代理为固井时钻井液一致。"),
    },
    "呼探1": {
        "平衡液": FluidProvenance("partial", "密度来自现场（设计1.50 vs 实际1.75，待统一），流变缺项按代理值补齐。"),
        "驱油隔离液": FluidProvenance("partial", "密度来自现场，流变缺项按代理值补齐。204131.doc含唯一明确化学配方。"),
        "领浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "中间浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "尾浆": FluidProvenance("field", "现场密度与幂律参数，来自实验报告。"),
        "压塞液": FluidProvenance("partial", "压塞液性质复用替浆液。"),
        "替浆液": FluidProvenance("partial", "密度按现场，流变参考类似井。"),
        "轻泥浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "中置液": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
        "井浆": FluidProvenance("partial", "现场阶段存在，但性质复用替浆液。"),
    },
    "呼探1-002": {
        "平衡液": FluidProvenance("partial", "密度来自现场，流变缺实测，使用代理值。"),
        "隔离液": FluidProvenance("partial", "密度来自现场，流变缺实测，使用代理值。"),
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
