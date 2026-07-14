# 呼探1-001井（HT1-001） paper provenance

- Paper loader: `ht1_001` / `paper_loader_v1`
- Geometry source: legacy/proxy geometry pending extracted-data cleanup
- Pumping schedule source: legacy/proxy schedule pending cleaned schedule generation
- Fluid source: legacy/proxy fluid table pending source freeze
- Standoff source: model_assumption
- CBL source: pending frozen CBL evaluation window
- Include in CBL metrics: False

## Notes
- 设计文档要求先清洗schedule并冻结CBL窗口后再进入定量验证。
- paper v1允许输出应用性主结果，但不纳入CBL误差统计。
- 几何和施工口径均需在后续资料清洗后升级。
