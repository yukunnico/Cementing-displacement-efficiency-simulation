# 呼101 paper provenance

- Paper loader: `hu101` / `paper_loader_v1`
- Geometry source: field_measured + equivalent composite liner approximation
- Pumping schedule source: field_measured construction sequence with shoe-lag correction inherited from legacy loader
- Fluid source: field_measured densities and rheology from extracted field report
- Standoff source: model_assumption profile documented in loader notes
- CBL source: interpreted official CBL report 100312.PDF
- Include in CBL metrics: True

## Notes
- 主验证候选井；CBL合格率采用0.6277。
- 上部168.3mm+下部139.7mm复合尾管按面积守恒等效几何进入当前单外径求解器。
- 193.7mm回接段资料不进入本 paper loader 主验证口径。
