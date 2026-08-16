# 控压固井顶替效率模型 - 项目 CLAUDE.md

> 补充全局 CLAUDE.md（`C:\Users\30525\.claude\CLAUDE.md`）之外的项目专属上下文。
> 代码规范、图表中文化、legacy 处理规则等见同目录 `AGENTS.md`，此处不重复。

## 项目定位

多井通用的固井顶替效率 1D-2D 耦合模型（D2DGA 框架），沉淀为可复用 `cemdisp` Python 包。
验证对象：呼图壁区块 8 口尾管井。

## 目录布局

### 代码（本仓库 `cement model/`，git 根）
- `cemdisp/` - 主包：`data/`(井/流体/泵注加载) + `transport1d/`(套管内1D) + `models2d/`(环空2D D2DGA) + `diagnostics/`(诊断/质量代理) + `validation/`(CBL对比/质量守恒) + `reporting/`(图表) + `runners/`(各井运行入口)
- `scripts/` - 分析/绘图/消融脚本
- `tests/` - 单元测试
- `docs/` - 与代码强耦合的技术文档（设计规格、调研，见下方索引）
- `hu101model/` `hu102model/` - legacy 原型，不随意改

### 论文（仓库外 `控压固井项目/论文构思及草稿/`）
- `石油学报版本/` - 当前主投版本（EI，~8000 字）
- `石油勘探与开发版本/` `石油钻采工艺版本/` - 另两版差异化论文
- `0708/` - 现场井资料原始数据；`数据提取/` - 数据提取中间产物
- `论文草稿.docx`、`《石油学报》论文格式模版.doc/txt`

## 文档存放约定（重要）

**新生成的 markdown 及独立文档一律存到 obsidian 仓库，统一用 obsidian 查看：**

```
D:\obsidian\obsidian-storage\1.科研\新疆油田控压固井项目\
├── 固井顶替效率改进\   模型解析、改进方向、因果分析、调研报告
├── 论文撰写\          论文框架、撰写指南
└── （根目录）          综合调研报告、目录索引
```

- `cement model/docs/` 仅放与代码强耦合的设计规格/调研（`superpowers/` 下 plans/specs/research）；
- 调研综述、模型解析、论文构思等**独立文档放 obsidian**；
- 两个源目录的现状索引见 obsidian `项目文档目录索引.md`。

## 环境与运行

- conda 环境：`shenjingwangluo`；可编辑安装：`python -m pip install -e .`
- 测试：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -v`
