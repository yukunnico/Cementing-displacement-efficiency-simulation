#!/usr/bin/env python
"""Part 4: target_intervals, temperature_pressure, construction_events, image_index, source_trace"""

import csv, os

out_dir = r"C:\Users\katsura\Desktop\research\控压固井项目\现场资料提取\ht1_001_呼探1-001"
os.makedirs(out_dir, exist_ok=True)

def write_csv(filename, headers, rows):
    filepath = os.path.join(out_dir, filename)
    with open(filepath, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in rows:
            writer.writerow(row)
    print(f"  OK {filename}: {len(rows)} rows, {os.path.getsize(filepath)} bytes")

# ======= File 11: target_intervals.csv =======
target_rows = [
    ["ht1_001", "K1s_水层1", 6060.0, 6064.0, "水层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值2.76%"],
    ["ht1_001", "K1s_水层2", 6092.0, 6098.0, "水层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值1.57%"],
    ["ht1_001", "K1q_干层", 7355.0, 7409.0, "干层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值3.87%"],
    ["ht1_001", "K1q_差气层1", 7409.0, 7423.0, "差气层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值3.92%"],
    ["ht1_001", "K1q_差气层2", 7460.0, 7462.0, "差气层(录井解释)-高峰显示", "施工设计/section 2.3", "录井解释结论", "high", "全烃最大值44.07%高峰"],
    ["ht1_001", "K1q_油气同层", 7504.0, 7566.0, "油气同层(录井解释)", "施工设计/section 2.3", "录井解释结论", "high", "全烃最大值6.21%"],
    ["ht1_001", "J3k3_差气层3", 7624.0, 7632.0, "差气层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值2.99%"],
    ["ht1_001", "J3k3_气层", 7642.0, 7661.0, "气层(录井解释)", "施工设计/section 2.3", "录井解释结论", "high", "全烃最大值4.83%"],
    ["ht1_001", "J3k3_差气层4", 7661.0, 7685.0, "差气层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值3.45%"],
    ["ht1_001", "J3k2_干层", 7734.0, 7746.0, "干层(录井解释)", "施工设计/section 2.3", "录井解释结论", "medium", "全烃最大值1.89%"],
    ["ht1_001", "K1q_J3k_主要目的层段", 6828.0, 7746.0, "清水河组+喀拉扎组;主要目的层", "施工设计/section 1.1", "地质设计", "high", "目的层:白垩系清水河组(K1q)+侏罗系喀拉扎组(J3k2)"],
    ["ht1_001", "裸眼模拟段", 5670.0, 7746.0, "尾管固井模拟段裸眼部分", "施工设计/section 4.2", "井身结构定义", "high", "裸眼段5670-7746m"],
    ["ht1_001", "重叠段", 5460.159, 5670.0, "尾管与技术套管重叠段", "施工设计/section 4.2", "设计", "high", "重叠段长度209.841m"],
    ["ht1_001", "领浆封固段", 5185.7, 5900.0, "领浆封固段(含超返)", "施工设计/section 4.2/7.1.1", "设计", "high", "超返274.5m至5185.7m"],
    ["ht1_001", "中间浆封固段", 5900.0, 7000.0, "中间浆封固段", "施工设计/section 4.2", "设计", "high", "裸眼段占高1100m"],
    ["ht1_001", "尾浆封固段", 7000.0, 7746.0, "尾浆封固段", "施工设计/section 4.2", "设计", "high", "裸眼段占高746m"],
]
write_csv("target_intervals.csv",
    ["well_id","interval_name","top_md_m","bottom_md_m","purpose","source_file","source_location","data_type","confidence","notes"],
    target_rows)

# ======= File 12: temperature_pressure_profile.csv =======
tp_rows = [
    ["ht1_001", 7746.0, 150.0, None, "measured", "施工设计/section 2.5", "电测静止温度150C;梯度1.7428C/100m", "high", "井底静止温度(BHST)"],
    ["ht1_001", 7746.0, 127.5, None, "calculated", "施工设计/section 2.5", "循环温度=150*0.85=127.5C", "high", "井底循环温度(BHCT);温度系数0.85"],
    ["ht1_001", 5460.159, 110.0, None, "design", "施工设计/section 9", "领浆顶部温度110C", "high", "领浆设计顶部静止温度"],
    ["ht1_001", 5900.0, 118.0, None, "design", "施工设计/section 9", "中间浆顶部温度118C", "high", "中间浆设计顶部静止温度"],
    ["ht1_001", 7000.0, 137.0, None, "design", "施工设计/section 9", "尾浆顶部温度137C", "high", "尾浆设计顶部静止温度"],
    ["ht1_001", 7746.0, None, 149.5, "calculated", "施工设计/section 3.3", "地层压力系数1.936推算", "medium", "Pe=1.936*0.00981*7746=147.07+井口套压"],
    ["ht1_001", 7746.0, None, 147.0, "design", "化验报告/Tables 2-6", "水泥浆实验压力147MPa", "high", "实验压力条件"],
    ["ht1_001", 7746.0, None, None, "measured", "施工设计/section 3.3", "静态承压ECD=2.011g/cm3", "high", "承压测试ECD=2.011g/cm3;罐区累计液面下降0.9m3"],
    ["ht1_001", 7746.0, None, None, "measured", "施工设计/section 3.2", "最大ECD=1.981g/cm3(施工期间)", "high", "施工期间最大ECD;小于承压值2.011"],
]
write_csv("temperature_pressure_profile.csv",
    ["well_id","md_m","temperature_c","pressure_mpa","data_type","source_file","source_location","confidence","notes"],
    tp_rows)

# ======= File 13: construction_events.csv =======
events_rows = [
    ["ht1_001", None, 7531.0, "钻进井漏", "钻进至7531m渗漏;密度1.98;漏失当量2.038;降密度至1.96;桥堵12方(30%)", "钻进期间",
     "施工设计/section 1.6", "2025.7.27 13:50", "high", "第1次井漏;非固井期间"],
    ["ht1_001", None, 7585.0, "钻进井漏", "钻进至7585m渗漏;密度1.95;漏失当量2.0;泵堵漏浆7m3", "钻进期间",
     "施工设计/section 1.6", "2025.8.1 09:40", "high", "第2次井漏"],
    ["ht1_001", None, 7623.0, "钻进井漏", "钻进至7623m渗漏;密度1.95;漏失当量2.013;泵堵漏浆10m3", "钻进期间",
     "施工设计/section 1.6", "2025.8.3 09:50", "high", "第3次井漏"],
    ["ht1_001", None, 7660.4, "钻进井漏", "钻进至7660.4m井漏;密度1.95;漏失当量2.002;静堵1h;泵堵漏浆12m3;降密度至1.93", "钻进期间",
     "施工设计/section 1.6", "2025.8.4 21:50", "high", "第4次井漏;最严重"],
    ["ht1_001", 0, None, "施工开始", "固井施工开始(冲洗管线试压30MPa)", "固井施工",
     "施工小结/section 9", "2025.8.25 21:05", "high", ""],
    ["ht1_001", 30, None, "试压", "冲洗管线试压30MPa通过", "管线试压",
     "施工小结/section 9", "21:05-21:15", "high", "试压通过"],
    ["ht1_001", 65, None, "排量偏差", "先导浆实际排量1.4m3/min(设计1.0);施工小结记载实际排量偏高", "先导浆注入",
     "施工小结/section 9", "21:15", "medium", "设计排量1.0,实际1.4m3/min"],
    ["ht1_001", 170, None, "停泵控压", "停泵倒闸门/压塞控压6.4MPa", "倒闸门控压",
     "施工小结/section 9", "23:51-00:00", "high", "精细控压操作"],
    ["ht1_001", 274, None, "碰压成功", "碰压7.5->14.8MPa;放回水断流;到量碰压成功", "碰压",
     "施工小结/section 9", "约01:34", "high", "碰压成功"],
    ["ht1_001", 274, None, "循环排混浆", "拔出中心管以2.4m3/min排量循环洗井排混浆和前置液", "循环洗井",
     "施工小结/section 9", "施工结束后", "high", ""],
    ["ht1_001", None, None, "控压候凝", "施工结束后加回压7.98MPa;关井憋压9.19-12.77MPa候凝;至领浆:隔离液=8:2强度达3.5MPa后开井", "候凝",
     "施工设计/section 5.1", "施工后候凝期间", "high", "精细控压候凝方案"],
    ["ht1_001", None, None, "固井质量评价", "固井质量统计:合格", "质量评价",
     "施工小结/section 10", "测井后", "low", "仅有'合格'定性结论;无CBL量化数据"],
]
write_csv("construction_events.csv",
    ["well_id","event_time_min","event_md_m","event_type","description","related_stage","source_file","source_location","confidence","notes"],
    events_rows)

# ======= File 14: image_extraction_index.csv =======
img_rows = [
    ["ht1_001", "img_001", "20512/205121.jpg", "jpg", 1, "浮鞋出厂检测报告(克拉玛依金鑫科技)", False, False, False, "vision", "high", True, "设备出厂检验报告;非CBL/VDL"],
    ["ht1_001", "img_002", "20512/205122.jpg", "jpg", 1, "浮箍出厂检测报告(西部钻探工程院)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_003", "20512/205123.jpg", "jpg", 1, "浮箍正向压力试验报告(28MPa/300s)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_004", "20512/205124.jpg", "jpg", 1, "浮箍出厂检测+反向压力试验(35MPa)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_005", "20512/205125.jpg", "jpg", 1, "浮箍正向压力试验报告(28MPa)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_006", "20514/205141/2051411.jpg", "jpg", 1, "浮鞋+浮箍出厂检测报告(与20512重复)", False, False, False, "vision", "high", True, "20514/205141中JPG与20512完全相同"],
    ["ht1_001", "img_007", "20514/205141/2051412.jpg", "jpg", 1, "浮鞋/浮箍正向压力试验(与205122重复)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_008", "20514/205141/2051413.jpg", "jpg", 1, "浮箍正向压力试验(与205123重复)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_009", "20514/205141/2051414.jpg", "jpg", 1, "浮箍出厂检测+反向压力试验(与205124重复)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_010", "20514/205141/2051415.jpg", "jpg", 1, "浮箍正向压力试验(与205125重复)", False, False, False, "vision", "high", True, ""],
    ["ht1_001", "img_011", "施工设计.doc", "word_embedded", "section 1.2.3", "井身结构示意图", True, True, False, "not_extracted", "medium", True, "尾管固井相关;待单独提取"],
    ["ht1_001", "img_012", "施工设计.doc", "word_embedded", "section 1.4.3", "X/Y轴井径曲线(1m间隔)", True, True, False, "not_extracted", "medium", True, "含尾管段井径连续曲线"],
    ["ht1_001", "img_013", "施工设计.doc", "word_embedded", "section 1.4.4", "平均井径曲线(1m间隔)", True, True, False, "not_extracted", "medium", True, "含尾管段平均井径连续曲线"],
    ["ht1_001", "img_014", "施工设计.doc", "word_embedded", "section 1.4.5", "井身轨迹图", True, False, False, "not_extracted", "medium", True, ""],
    ["ht1_001", "img_015", "施工设计.doc", "word_embedded", "section 6.3", "居中度模拟图(裸眼段居中度80.4%)", True, True, False, "not_extracted", "medium", True, "重要;居中度模拟结果"],
    ["ht1_001", "img_016", "施工设计.doc", "word_embedded", "section 7.1.5", "浆柱结构图与施工动态图", True, True, False, "not_extracted", "medium", True, "可验证泵注程序"],
    ["ht1_001", "img_017", "施工设计.doc", "word_embedded", "section 8.2.3", "斯伦贝谢顶替效率模拟图(平均水泥充填率98%)", True, True, False, "not_extracted", "medium", True, "顶替效率模拟结果"],
    ["ht1_001", "img_018", "化验报告.docx", "word_embedded", "Fig 1-4", "领浆稠化曲线+温度高点+升降温停机实验", True, True, False, "not_extracted", "low", True, "水泥浆性能曲线"],
    ["ht1_001", "img_019", "化验报告.docx", "word_embedded", "Fig 5-6", "中间浆稠化曲线+温度高点实验", True, True, False, "not_extracted", "low", True, ""],
    ["ht1_001", "img_020", "化验报告.docx", "word_embedded", "Fig 7-8", "尾浆稠化曲线+温密双高点实验", True, True, False, "not_extracted", "low", True, ""],
    ["ht1_001", "img_021", "化验报告.docx", "word_embedded", "Fig 9-14", "相容性实验稠化曲线(7:3,7:2:1,7:1:2等)", True, True, False, "not_extracted", "low", True, "相容性实验"],
    ["ht1_001", "img_022", "化验报告.docx", "word_embedded", "Fig 16-17", "静胶凝实验曲线(中间浆118C;尾浆150C)", True, True, False, "not_extracted", "medium", True, "失重分析和候凝补压可用"],
    ["ht1_001", "img_023", "HT1-001井油层尾管报告.pdf", "pdf", None, "PDF化验报告(含稠化曲线/强度等)", True, True, False, "not_extracted", "medium", True, "PDF版化验报告;可能与docx互补"],
]
write_csv("image_extraction_index.csv",
    ["well_id","image_id","source_file","source_type","page_number_or_image_index","detected_content","is_liner_related","contains_useful_data","extracted_to_file","ocr_or_vision_method","confidence","needs_human_review","notes"],
    img_rows)

# ======= File 15: source_trace.csv =======
trace_rows = [
    ["ht1_001", "施工设计", "well_type", "直井", None, "field_measured", "施工设计", "section 1.1", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "liner_hanger_md", 5460.159, "m", "design", "施工设计", "section 4.2/6.1", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工小结", "liner_bottom_md", 7746.0, "m", "field_measured", "施工小结", "section 9", "word_text", "COM_text_extraction", "high", False, "实际下深"],
    ["ht1_001", "施工设计", "borehole_sizes", "241.3/215.9", "mm", "design", "施工设计", "section 1.2.1", "word_table", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "caliper_profile", "69 points @30m", "mm", "field_measured", "施工设计", "section 1.4.1", "word_table", "COM_text_extraction", "high", False, "30m间隔"],
    ["ht1_001", "施工设计", "inclination_profile", "69 points @30m", "deg", "field_measured", "施工设计", "section 1.4.1", "word_table", "COM_text_extraction", "high", False, "30m间隔"],
    ["ht1_001", "施工设计", "pumping_schedule_design", "详见表", None, "design", "施工设计", "section 7.2", "word_table", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工小结", "pumping_schedule_actual", "详见表", None, "field_measured", "施工小结", "section 9", "word_text", "COM_text_extraction", "high", False, "实际施工数据"],
    ["ht1_001", "施工设计", "casing_string_combo", "详见表", None, "design", "施工设计", "section 6.1", "word_table", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "套管数据表", "casing_detail_139.7", "详见表", None, "field_measured", "套管数据表xlsx", "Sheet 139.7", "excel_table", "openpyxl_read", "high", False, "三方核对版"],
    ["ht1_001", "套管数据表", "casing_detail_168.3", "详见表", None, "field_measured", "套管数据表xlsx", "Sheet 168.3", "excel_table", "openpyxl_read", "high", False, "三方核对版"],
    ["ht1_001", "施工设计", "centralizer_layout", "101只(每2根1只)", None, "design", "施工设计", "section 6.2", "word_table", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "mud_rheology_65C", "六速读数详表", None, "field_measured", "施工设计", "section 1.3.1", "word_table", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "化验报告", "fluid_rheology_128C", "六速读数详表", None, "field_measured", "化验报告docx", "Table 8", "word_table", "docx2txt", "high", False, ""],
    ["ht1_001", "化验报告", "cement_properties", "详见表4-6", None, "field_measured", "化验报告docx", "Tables 4-6", "word_table", "docx2txt", "high", False, ""],
    ["ht1_001", "化验报告", "spacer_properties", "详见表7", None, "field_measured", "化验报告docx", "Table 7", "word_table", "docx2txt", "high", False, ""],
    ["ht1_001", "工程量清单", "cement_additives", "详见表", None, "field_measured", "工程量清单doc", "全文", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "技术总结", "actual_construction", "详见表", None, "field_measured", "技术总结doc", "section 2/3", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工小结", "cbl_quality_result", "合格", None, "interpreted", "施工小结", "section 10", "word_text", "COM_text_extraction", "low", True, "仅有定性结论"],
    ["ht1_001", "施工设计", "target_intervals", "详见表", None, "interpreted", "施工设计", "section 2.3", "word_table", "COM_text_extraction", "medium", False, "录井解释"],
    ["ht1_001", "施工设计", "drilling_anomalies", "4次井漏", None, "field_measured", "施工设计", "section 1.6", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "safety_window", "1.936-2.011", "g/cm3", "field_measured", "施工设计", "section 3.3", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "formation_pressure_coeff", 1.936, "g/cm3", "field_measured", "施工设计", "section 3.3", "word_text", "COM_text_extraction", "high", False, ""],
    ["ht1_001", "施工设计", "bottomhole_static_temp", 150.0, "C", "field_measured", "施工设计", "section 2.5", "word_text", "COM_text_extraction", "high", False, ""],
]
write_csv("source_trace.csv",
    ["well_id","data_file","field_name","value","unit","data_type","source_file","source_location","source_type","extraction_method","confidence","needs_human_review","notes"],
    trace_rows)

print("Files 11-15 done.")