# 呼103 paper provenance

- Paper loader: `hu103` / `paper_loader_v1`
- Geometry source: field measured 139.7mm lower liner interval with explicit composite-liner context
- Pumping schedule source: confirmed design/actual schedule represented by legacy loader for paper v1
- Fluid source: field/model-ready fluid table represented by legacy loader
- Standoff source: design centralization proxy adjusted by clearance; model_assumption
- CBL source: interpreted CBL for 139.7mm interval 7338-7712m
- Include in CBL metrics: True

## Notes
- 复合尾管建模改进重点井；主CBL验证仅采用139.7mm段7338-7712m合格率0.1206。
- 168.3mm段CBL 0.0004和整段综合0.0605单独保留在说明中，不混入主验证窗。
- 当前求解器若按单外径运行，paper metadata必须保留等效/代理说明。
