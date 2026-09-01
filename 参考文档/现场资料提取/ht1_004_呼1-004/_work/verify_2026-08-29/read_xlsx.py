import openpyxl, sys
sys.stdout.reconfigure(encoding='utf-8')
for f in [r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-003与HT1-004固井设计实际施工对比.xlsx",
          r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井四开油层尾管数据表-标准格式6.11最终.xlsx"]:
    print("="*100); print("FILE:", f.split("\\")[-1])
    wb = openpyxl.load_workbook(f, data_only=True)
    for ws in wb.worksheets:
        print("--- SHEET:", ws.title, " dims:", ws.dimensions)
        for row in ws.iter_rows():
            vals = [("" if c.value is None else str(c.value).replace("\n","⏎")) for c in row]
            if any(v.strip() for v in vals):
                print(" | ".join(vals).rstrip(" |"))
