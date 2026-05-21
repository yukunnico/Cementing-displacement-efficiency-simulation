import json

with open('参考文档/呼1-003/提取数据/extracted_data.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

print('Searching for key tables...')
print('Total tables:', len(data['tables']))

# Search for tables containing specific keywords
found_tables = []
for idx, table in enumerate(data['tables']):
    table_text = ' '.join([' '.join(row) for row in table['data']])
    
    # Look for key data tables
    category = None
    if any(kw in table_text for kw in ['密度', '粘度', '塑性粘度', '动切力']):
        category = 'drilling_fluid'
    
    if any(kw in table_text for kw in ['稠化时间', 'API失水', '抗压强度', '水泥浆']):
        category = 'cement_slurry'
    
    if any(kw in table_text for kw in ['扶正器', '弹性扶正器', '居中度']):
        category = 'centralizer'
    
    if any(kw in table_text for kw in ['流变', 'n值', 'K值', '流变参数']):
        category = 'rheology'
    
    if category:
        found_tables.append((category, idx, table))

print('Found', len(found_tables), 'relevant tables')

# Save detailed results
with open('参考文档/呼1-003/提取数据/key_tables.json', 'w', encoding='utf-8') as f:
    output = []
    for category, idx, table in found_tables:
        output.append({
            'category': category,
            'table_index': idx,
            'rows': table['rows'],
            'cols': table['cols'],
            'data': table['data']
        })
    json.dump(output, f, ensure_ascii=False, indent=2)

print('Saved to key_tables.json')
