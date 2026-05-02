import pandas as pd
import warnings
warnings.filterwarnings(" ignore\)
file_path = r\D:\\users\\desktop\\research\\控压固井项目\\0708\\导出\\sisyphus_output\\h103_spacer_search_results.xlsx\
xl = pd.ExcelFile(file_path)
for sheet in xl.sheet_names:
 print(\=\*60)
 print(\SHEET:\, sheet)
 print(\=\*60)
 df = pd.read_excel(file_path, sheet_name=sheet, header=None)
 for idx, row in df.iterrows():
 for col_idx, val in enumerate(row):
 if pd.notna(val):
 print(\ROW\, idx, \COL\, col_idx, \:\, val)
