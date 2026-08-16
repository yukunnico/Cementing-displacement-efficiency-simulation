# -*- coding: utf-8 -*-
import pandas as pd
p = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井四开油层尾管数据表-标准格式6.11最终.xlsx"
xl = pd.ExcelFile(p)
print("SHEETS:", xl.sheet_names)
for s in xl.sheet_names:
    df = xl.parse(s, header=None)
    print(f"\n===== SHEET: {s} shape={df.shape} =====")
    # print first 40 rows, all columns (truncated cells)
    with pd.option_context('display.max_columns', 20, 'display.width', 250, 'display.max_colwidth', 30):
        print(df.head(40).to_string())
