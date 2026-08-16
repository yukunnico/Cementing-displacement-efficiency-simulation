# -*- coding: utf-8 -*-
"""生成图片索引/来源追踪 CSV"""
import csv, os

OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003"
WELL = "ht1_003"
DSGN = "HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批).doc"
REC  = "HT1-003井168.3+139.7mm尾管固井施工记录表.doc"
HIST = "HT1-003井油层尾管固井作业史.doc"
LAB  = "化验报告/HT1-003 油层尾管 化验报告.docx"

def w(fname, header, rows):
    p = os.path.join(OUT, fname)
    with open(p, "w", encoding="utf-8-sig", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(header)
        wr.writerows(rows)
    print(f"wrote {fname}: {len(rows)} rows")

# ============ 14. image_extraction_index ============
w("image_extraction_index.csv",
  ["well_id","image_id","source_file","source_type","page_number_or_image_index","detected_content","is_liner_related","contains_useful_data","extracted_to_file","ocr_or_vision_method","confidence","needs_human_review","notes"],
  [
  [WELL,"cbl_15_3500","HT1-003_固井质量测井评价图_15-3500_20260622_完井.pdf","pdf_text_layer","1页","CBL评价图头：测量井段15-3500m；固井质量不合格；违反七条红线两条；水泥面5307.54m；人工井底7514.21m；套管程序；目的井段7442-7618m",1,1,"cbl_evaluation.csv","pymupdf_text","",0,"有文字层图头直接提取；曲线未定量"],
  [WELL,"cbl_3500_7595","HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf","pdf_text_layer","1页","CBL评价图头：测量井段3500-7595m；固井质量不合格；技术说明同前；深度刻度7595-5100m",1,1,"cbl_evaluation.csv","pymupdf_text","",0,"有文字层；曲线未定量"],
  [WELL,"cbl_ibc_10_7520","HT1-003井IBC固井质量与套损评价成果图_10-7520m.pdf","pdf_text_layer","137页","IBC固井质量与套损评价(斯伦贝谢FITS/USIT)曲线：GR/ECCE/RSAV/声幅AWMN-AWMX/变密度/USIT/套损；仅曲线标注文字无结论文字层",1,0,"","pymupdf_text","",1,"137页大图；文字层仅标注；曲线定量需人工"],
  [WELL,"plan_pdf","HT1-003井施工方案(1).pdf","scanned_pdf","12页","固井施工方案(扫描件)",1,0,"","tesseract_ocr_chi_sim","",1,"无文字层；低清扫描OCR质量差，未能可靠提取数值；施工流程与设计书一致未新增关键数据"],
  [WELL,"brief_pdf","HT1-003井尾管固井施工技术交底书.pdf","scanned_pdf","6页","固井施工技术交底(扫描件)：施工流程、控压当量2.00/2.02/2.04g/cm3、浆柱、异常处置",1,0,"","tesseract_ocr_chi_sim","",1,"OCR部分可读确认与设计一致；数值需人工复核"],
  [WELL,"third_party_pdf","化验报告/HT1-003井油层尾管水泥浆第三方检测报告.pdf","scanned_pdf","9页","水泥浆第三方检测报告(北京工程院复检)",1,0,"","tesseract_ocr_chi_sim","",1,"无文字层；低清扫描OCR严重失真数值不可靠；以化验报告.docx为准"],
  [WELL,"torque_pdf","80027队HTI-003井139.7和168.28套管扭矩数据.pdf","pdf_text_layer","61页","上扣扭矩控制汇总表：209根套管上扣记录(139.7mm BGT2约25000N.m等)，全部合格",1,1,"casing_liner_string.csv","pymupdf_text","",0,"209条上扣数据均已合格"],
  [WELL,"daq_ti","数据采集/80027 HT1-003 尾管替桨.xlsm","xlsm_curve","数据表2236行","实际替浆泵注时间序列：2026-05-27 16:56-20:07，总泵量91.01m3排量峰值1.41m3/min泵压峰值27MPa",1,1,"pumping_schedule.csv","openpyxl","",0,"核对施工记录表替浆91.9m3一致"],
  [WELL,"daq_peil","数据采集/80027HT1-003尾管固井配隔离液.xlsm","xlsm_curve","数据表1150行","实际配隔离液：16:07-17:44密度峰值2.30总量17.2m3；16:33密度2.127总量16m3(隔离液1)，17:23密度2.078总量8.4m3(隔离液2)",1,1,"fluid_properties.csv","openpyxl","",0,"实际隔离液密度2.05-2.13 / 1.95-2.08"],
  [WELL,"daq_ch","数据采集/80027 HT1-003 尾管 抽灰.xlsm","xlsm_curve","数据表2089行","实际抽灰(水泥混配)：密度-0.16-2.69(含清水)，流量峰值1.10m3/min；设备时钟2010年未校准",1,1,"fluid_properties.csv","openpyxl","",1,"时钟日期错误(2010-10-08)仅作水泥浆密度/排量核对"],
  [WELL,"lab_docx","化验报告/HT1-003 油层尾管 化验报告.docx","docx_tables","139行","油层尾管水泥浆实验报告：领/尾浆配方、稠化、失水、强度、流变六速、相容性、静胶凝",1,1,"fluid_properties.csv, rheometer_readings.csv","python-docx","",0,"本次流体性能主数据源"],
  ])

# ============ 15. source_trace ============
w("source_trace.csv",
  ["well_id","data_file","field_name","value","unit","data_type","source_file","source_location","source_type","extraction_method","confidence","needs_human_review","notes"],
  [
  [WELL,"","well_name_cn","呼1-003","","field_measured",DSGN,"封面","word_text","word_com","high",0,""],
  [WELL,"","alias","HT1-003","","field_measured",DSGN,"封面","word_text","word_com","high",0,""],
  [WELL,"","block","准噶尔盆地南缘冲断带霍玛吐背斜带呼西背斜","","field_measured",DSGN,"封面/前言","word_text","word_com","high",0,""],
  [WELL,"","field","呼西背斜呼探1井区","","field_measured","CBL评价图头 地区","图头","pdf_text","pymupdf","high",0,""],
  [WELL,"","well_type","直井","","field_measured",DSGN,"封面","word_text","word_com","high",0,""],
  [WELL,"","cementing_type","油层尾管精细控压固井(168.3+139.7mm复合尾管)","","field_measured",DSGN,"标题","word_text","word_com","high",0,""],
  [WELL,"","actual_td","7618","m","field_measured",DSGN,"前言","word_text","word_com","high",0,"完钻井深；设计井深7560m"],
  [WELL,"","hanger_md","5307.539","m","field_measured",DSGN,"前言","word_text","word_com","high",0,"悬挂器喇叭口位置"],
  [WELL,"","hanger_body_md","5315.439","m","field_measured",HIST,"八(2)工具数据","word_text","antiword","high",0,"尾管悬挂器下深"],
  [WELL,"","liner_shoe_md","7618","m","field_measured",DSGN,"六.6.1","word_text","word_com","high",0,"加长引鞋下深"],
  [WELL,"","openhole","5568-7618","m","field_measured",DSGN,"一.1.2.1","word_text","word_com","high",0,"273.1mm技套鞋5568m"],
  [WELL,"","borehole_nominal","241.3+215.9","mm","field_measured",DSGN,"一.1.2.1","word_text","word_com","high",0,"四开双径钻头"],
  [WELL,"","mud_density","1.95","g/cm3","field_measured",DSGN,"一.1.3.1","word_table","word_com","high",0,"油基钻井液"],
  [WELL,"","mud_pv_yp","51/10","mPa.s/Pa","field_measured",DSGN,"一.1.3.1","word_table","word_com","high",0,"塑粘51/动切10"],
  [WELL,"","mud_six_speed","122/71/52/21/7/6","deg","field_measured",DSGN,"一.1.3.1","word_table","word_com","high",0,"Φ600-Φ3"],
  [WELL,"","caliper_electrical","66点30m间隔","mm","field_measured",DSGN,"一.1.4.2","word_table","word_com","high",0,"电测井径5570-7618m；5568-7096平均247.3mm扩大率2.5%；7096-7618平均218.9mm扩大率1.4%"],
  [WELL,"","inclination_electrical","67点30m间隔","deg","field_measured",DSGN,"一.1.4.1","word_table","word_com","high",0,"井斜方位温度5560-7618m；最大井斜4.2度@7447m"],
  [WELL,"","centralizer_count","75+22=97","只","field_measured",DSGN,"六.6.2","word_table","word_com","high",0,"整体弹性扶正器"],
  [WELL,"","centralization_percent","83","%","interpreted",DSGN,"六.6.3","word_text","word_com","medium",0,"裸眼段居中度模拟"],
  [WELL,"","pumping_actual","91.9","m3","field_measured",REC,"替浆","word_table","antiword","high",0,"替浆量；DAQ总泵量91.01m3核对一致"],
  [WELL,"","bump_pressure","19-25","MPa","field_measured",REC,"碰压","word_table","antiword","high",0,""],
  [WELL,"","cement_volume","39+28=67","m3","field_measured",DSGN,"七.7.1.1","word_table","word_com","high",0,"领浆39+尾浆28；灰量157t(设计=实际)"],
  [WELL,"","lead_thickening","455","min/70Bc","field_measured",LAB,"表4 领浆配方与性能","docx_table","python-docx","high",0,""],
  [WELL,"","tail_thickening","256","min/70Bc","field_measured",LAB,"表5 尾浆配方与性能","docx_table","python-docx","high",0,""],
  [WELL,"","lead_strength_24h","3.9","MPa","field_measured",LAB,"表4","docx_table","python-docx","high",0,""],
  [WELL,"","tail_strength_24h","15.1","MPa","field_measured",LAB,"表5","docx_table","python-docx","high",0,""],
  [WELL,"","cbl_conclusion","不合格","","field_measured","CBL评价图头 固井质量","图头","pdf_text","pymupdf","high",0,"违反七条红线2条"],
  [WELL,"","cbl_pass_rate","缺失(无文字结论)","","pending_digitization","CBL评价图","","pdf_text","pymupdf","low",1,"需人工数字化曲线求合格率"],
  [WELL,"","third_party_report","低清扫描OCR不可靠","","missing","第三方检测报告.pdf","扫描件","scanned_pdf","tesseract_ocr","low",1,"以化验报告.docx为准"],
  ])
print("INDEX/TRACE DONE")
