# 呼102 paper provenance

- Paper loader: `hu102` / `paper_loader_v1`
- Geometry source: field-rich source package; legacy reconstructed loader profile adopted for paper v1
- Pumping schedule source: legacy reconstructed schedule from field-rich records
- Fluid source: field/extracted densities and rheology as represented by legacy loader
- Standoff source: model_assumption because continuous field standoff is unavailable
- CBL source: interpreted official CBL report 100413.PDF
- Include in CBL metrics: True

## Notes
- 主验证候选井；CBL合格率采用0.6665。
- 设计文档指出旧loader需重构；paper v1明确标注当前仍采用legacy reconstructed口径。
- CBL评价窗与目标层段必须分别读取，不得混用。
