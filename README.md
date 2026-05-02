# 固井顶替效率建模项目骨架

这个目录现在包含两类内容：

1. 既有的历史脚本/算例与资料；
2. 新建立的、可持续扩展的 Python 项目骨架。

新骨架的目标不是立即重写所有旧模型，而是先把“顶替效率模型综述 + 可安装代码结构 + 统一变量命名”稳定下来，后续再逐步接入经验模型和 2D 机理模型。

## 当前新增结构

```text
docs/
  displacement_efficiency_review.md
cemdisp/
  correlations/
  models2d/
  utils/
tests/
pyproject.toml
environment.yml
```

## 建议工作流

### 1. 激活环境

```powershell
conda activate shenjingwangluo
```

### 2. 本地可编辑安装

```powershell
python -m pip install -e .
```

### 3. 运行基础测试

```powershell
python -m unittest discover -s tests -v
```

## 代码设计原则

- `cemdisp.correlations`：用于放经验/拟合型顶替效率模型；
- `cemdisp.models2d`：用于放 Hele-Shaw / 2DGA 类机理模型接口；
- `cemdisp.geometry`：负责环空几何与偏心参数；
- `cemdisp.fluid`：负责流体流变数据结构；
- `cemdisp.utils.dimensionless`：放公共计算函数；
- `cemdisp.utils.validators`：放适用范围校验和警告。

## 与旧脚本的关系

当前目录中的 `环空六浆柱注入阶段模型(1)` 被视为 **legacy 原型资产**：

- 先保留，不直接改写；
- 后续逐步把可复用公式迁入新包；
- 新代码优先解决：统一单位、统一变量名、解耦数据文件与脚本。

## 下一步建议

1. 先把一到两个经验模型写实；
2. 再接入旧脚本中的偏心摩阻与环空几何公式；
3. 最后再扩展到 2DGA 数值求解。
