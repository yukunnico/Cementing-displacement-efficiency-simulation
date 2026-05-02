import pandas as pd
df = pd.read_excel(r" D:\\users\\desktop\\research\\控压固井项目\\0708\\导出\\sisyphus_output\\h103_spacer_search_results.xlsx\, sheet_name=\Summary\, header=None)
search_terms = [" PV\, \YP\, \n\, \K\, \屈服值\, \塑性\, \密度\, \体积\, \m3\, \尾管\, \尾浆\]
for idx, row in df.iterrows():
    row_text = str(row.tolist())
    for term in search_terms:
        if term in row_text:
            print(" Row:\, idx, \Term:\, term, \Found in:\, row_text[:300])
