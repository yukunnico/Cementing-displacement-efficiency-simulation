| well_id | well_name_cn | predicted_cbl_efficiency | measured_cbl_pass_rate | absolute_error | relative_error_percent | included_in_metrics | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hu101 | 呼101 | 0.650372 | 0.6277 | 0.022672 | 3.611986 | True | 主验证候选井；CBL合格率采用0.6277。; 上部168.3mm+下部139.7mm复合尾管按面积守恒等效几何进入当前单外径求解器。; 193.7mm回接段资料不进入本 paper loader 主验证口径。 |
| hu102 | 呼102 | 0.379941 | 0.6665 | 0.286559 | 42.994576 | True | 主验证候选井；CBL合格率采用0.6665。; 设计文档指出旧loader需重构；paper v1明确标注当前仍采用legacy reconstructed口径。; CBL评价窗与目标层段必须分别读取，不得混用。 |
| hu103 | 呼103 | 1.0 | 0.1206 | 0.8794 | 729.187396 | True | 复合尾管建模改进重点井；主CBL验证仅采用139.7mm段7338-7712m合格率0.1206。; 168.3mm段CBL 0.0004和整段综合0.0605单独保留在说明中，不混入主验证窗。; 当前求解器若按单外径运行，paper metadata必须保留等效/代理说明。 |
| ht1_001 | 呼探1-001井（HT1-001） | 0.892118 |  |  |  | False | 设计文档要求先清洗schedule并冻结CBL窗口后再进入定量验证。; paper v1允许输出应用性主结果，但不纳入CBL误差统计。; 几何和施工口径均需在后续资料清洗后升级。 |
