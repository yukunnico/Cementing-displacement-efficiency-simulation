# 论文数据代码已删除（c19f5a5，待论文代码重写）；其测试暂挂起收集，避免 ImportError 阻断全量 pytest。
# 重写论文 loader（cemdisp.data.loaders.paper / scripts.paper_data）后删除本文件即可恢复。
collect_ignore = ["test_cbl_metrics.py", "test_paper_loaders.py", "test_paper_result_schema.py"]
