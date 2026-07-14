#!/usr/bin/env python
"""Part 3: pumping_schedule, fluid_properties, rheometer_readings, cbl_evaluation"""

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

# ======= File 7: pumping_schedule.csv (actual data from 施工小结 section 9) =======
pump_rows = [
    ["ht1_001", 1, "管线试压", "shutdown", "冲洗管线试压", None, None, 0.0, 10.0, None, "pressure_test_30MPa",
     "measured", "施工小结/section 9", "2025.8.25/26 21:05-21:15;试压30MPa", "high", ""],
    ["ht1_001", 2, "先导浆", "wash", "先导浆(平衡液/冲洗液)", 40.0, 1.4, 10.0, 45.0, 1.75, "",
     "measured", "施工小结/section 9", "21:15-21:50;泵压17.9->20.4MPa;占高1171m;排量设计1.0实际1.4", "high", "设计排量1.0m3/min;实际1.4"],
    ["ht1_001", 3, "驱油隔离液", "spacer", "驱油隔离液(1.98g/cm3)", 20.0, 1.0, 50.0, 70.0, 1.98, "",
     "measured", "施工小结/section 9", "22:05-22:30;泵压14.9->15.2MPa;占高577.8m", "high", ""],
    ["ht1_001", 4, "领浆", "lead_cement", "领浆(2.05g/cm3)", 20.6, 1.0, 70.0, 93.0, 2.05, "",
     "measured", "施工小结/section 9", "22:30-22:53;泵压15.7->18MPa;封固5186-5900m", "high", ""],
    ["ht1_001", 5, "中间浆", "intermediate_cement", "中间浆(1.90g/cm3)", 28.7, 1.0, 93.0, 125.0, 1.90, "",
     "measured", "施工小结/section 9", "22:53-23:25;泵压13.8->15.6MPa;封固5900-7000m", "high", ""],
    ["ht1_001", 6, "尾浆", "tail_cement", "尾浆(1.90g/cm3)", 22.1, 1.0, 125.0, 156.0, 1.90, "",
     "measured", "施工小结/section 9", "23:25-23:51;泵压13-15.2MPa;封固7000-7746m", "high", "实际注22.1m3(设计21.4m3)"],
    ["ht1_001", 7, "倒闸门压塞", "shutdown", "停泵倒闸门/压塞控压", None, None, 156.0, 165.0, None, "pressure_control",
     "measured", "施工小结/section 9", "23:51-00:00;控压6.4MPa", "high", "倒闸门+控压等待"],
    ["ht1_001", 8, "压塞液", "displacement", "压塞液(后置液/1.98g/cm3)", 2.0, 0.6, 165.0, 168.0, 1.98, "plug_set",
     "measured", "施工小结/section 9", "00:00-00:03;泵压7MPa", "high", "压塞2m3"],
    ["ht1_001", 9, "替钻井液1", "displacement", "替钻井液(1.92g/cm3)", 25.0, 1.4, 168.0, 187.0, 1.92, "",
     "measured", "施工小结/section 9", "00:03-00:22;泵压21-19MPa", "high", "实际25m3;设计23m3(含12+12)"],
    ["ht1_001", 10, "中置液(保护液)", "displacement", "中置液(保护液/驱油隔离液)", 10.0, 1.0, 187.0, 198.0, 1.98, "",
     "measured", "施工小结/section 9", "00:22-00:33;泵压11.8-10.6MPa", "high", "中置液10m3"],
    ["ht1_001", 11, "基液", "displacement", "隔离液基液(1.02g/cm3)", 3.0, 1.0, 198.0, 200.0, 1.02, "",
     "measured", "施工小结/section 9", "00:33-00:35;泵压11MPa", "high", "基液3m3"],
    ["ht1_001", 12, "替钻井液2", "displacement", "替钻井液(1.92g/cm3)", 53.7, 0.9, 200.0, 259.0, 1.92, "bump_plug",
     "measured", "施工小结/section 9", "00:35-01:34;泵压11.5->20.3->7.5MPa;排量1.2-0.6;碰压7.5->14.8MPa放回水断流", "high", "碰压成功;总替量94.7m3(设计)"],
    ["ht1_001", 13, "循环洗井", "shutdown", "拔出中心管循环排混浆", None, 2.4, 259.0, None, None, "circulate_contaminated",
     "measured", "施工小结/section 9", "以2.4m3/min排量循环;排混浆和前置液", "high", "施工结束后循环洗井"],
]
write_csv("pumping_schedule.csv",
    ["well_id","step_index","stage_name","fluid_role","fluid_name","volume_m3","rate_m3_min","start_time_min","end_time_min","density_g_cm3","event_tag","data_type","source_file","source_location","confidence","notes"],
    pump_rows)

# ======= File 8: fluid_properties.csv =======
fluid_rows = [
    ["ht1_001", "钻井液(井浆@65C)", "mud", 1.92, "Power-Law", 51.0, 6.0, None, 0.18, 0.828, None, 65.0,
     "measured", "施工设计/section 1.3.1", "井口65C;六速:114/63/49/28/5/4", "high", "油基钻井液;油水比86:14;初切2.5Pa终切7.5Pa"],
    ["ht1_001", "钻井液(化验复测@128-93C)", "mud", 1.93, "Power-Law", None, None, None, 0.557, 0.689, None, 128.0,
     "measured", "化验报告/Table 8", "128->93C;六速:114/63/49/28/5/4", "high", "温度变化导致K值从0.18增至0.557;n从0.828降至0.689"],
    ["ht1_001", "先导浆(平衡液)", "wash", 1.75, None, None, None, None, None, None, None, None,
     "design", "施工设计/section 5.1", "设计要求:密度1.75;粘度<60s;PV<35;YP<6;含砂<0.3%", "medium", "未找到实测流变数据"],
    ["ht1_001", "驱油隔离液", "spacer", 1.98, "Power-Law", None, None, None, 1.338, 0.545, None, 93.0,
     "measured", "化验报告/Table 7,8", "128->93C;六速:114/80/66/43/10/6;冲洗效率97.5%", "high", "沉降稳定性1.97/1.98/1.99;配方:水+3.5%BCJ-300S+8%XZ-JS1L+40%XZ-YCXJ+1%XZ-XP1L+195%重晶石"],
    ["ht1_001", "领浆(2.05g/cm3)", "lead_cement", 2.05, "Power-Law", None, None, None, 0.876, 0.811, None, 128.0,
     "measured", "化验报告/Table 4,8", "128C;六速:>300/268/194/113/10/7;稠化427min;失水43ml", "high", "G+30%铁矿粉+35%石英砂+3%WG;24h强度3.7MPa;48h强度7.6MPa;7d强度30.8MPa"],
    ["ht1_001", "中间浆(1.90g/cm3)", "intermediate_cement", 1.90, "Power-Law", None, None, None, 0.504, 0.871, None, 128.0,
     "measured", "化验报告/Table 5,8", "128C;六速:>300/221/165/92/7/4;稠化314min;失水40ml", "high", "G+35%石英砂+4%WG+3%BCE-311S+4%蛭石;24h强度14.2MPa;7d强度31.3MPa"],
    ["ht1_001", "尾浆(1.90g/cm3)", "tail_cement", 1.90, "Power-Law", None, None, None, 0.453, 0.886, None, 128.0,
     "measured", "化验报告/Table 6,8", "128C;六速:>300/218/163/88/6/4;稠化208min;失水43ml", "high", "G+35%石英砂+4%WG+4%BCY-200S+3%BCE-311S+4%蛭石;24h强度15.1MPa;7d强度32.9MPa"],
    ["ht1_001", "压塞液/中置液", "displacement", 1.98, None, None, None, None, None, None, None, None,
     "design", "施工设计/section 7.1.3", "与驱油隔离液配方相同", "medium", "未单独实测"],
    ["ht1_001", "隔离液基液", "displacement", 1.02, None, None, None, None, None, None, None, None,
     "design", "施工设计/section 7.1.3", "驱油隔离液基液(水化)", "medium", ""],
]
write_csv("fluid_properties.csv",
    ["well_id","fluid_name","fluid_role","density_g_cm3","rheology_model","pv_mpa_s","yp_pa","yield_stress_pa","consistency_k_pa_sn","flow_index_n","viscosity_mpa_s","test_temperature_c","data_type","source_file","source_location","confidence","notes"],
    fluid_rows)

# ======= File 9: rheometer_readings.csv =======
rheo_rows = [
    ["ht1_001", "钻井液", "128->93", 114, 63, 49, 28, 5, 4, "化验报告/Table 8", "128->93C降温测试", "high", "油基钻井液;K=0.557 n=0.689"],
    ["ht1_001", "驱油隔离液", "128->93", 114, 80, 66, 43, 10, 6, "化验报告/Table 8", "128->93C降温测试", "high", "K=1.338 n=0.545"],
    ["ht1_001", "领浆", "128->93", 300, 268, 194, 113, 10, 7, "化验报告/Table 8", "128->93C降温测试", "high", "PHI600>300未完全读数;K=0.876 n=0.811"],
    ["ht1_001", "中间浆", "128->93", 300, 221, 165, 92, 7, 4, "化验报告/Table 8", "128->93C降温测试", "high", "PHI600>300未完全读数;K=0.504 n=0.871"],
    ["ht1_001", "尾浆", "128->93", 300, 218, 163, 88, 6, 4, "化验报告/Table 8", "128->93C降温测试", "high", "PHI600>300未完全读数;K=0.453 n=0.886"],
]
write_csv("rheometer_readings.csv",
    ["well_id","fluid_name","temperature_c","theta_600","theta_300","theta_200","theta_100","theta_6","theta_3","source_file","source_location","confidence","notes"],
    rheo_rows)

# ======= File 10: cbl_evaluation.csv =======
cbl_rows = [
    ["ht1_001", "呼探1-001", 5460.159, 7746.0, None, "合格", "施工小结记载'固井质量统计:合格'，无具体CBL振幅值或分段合格率;无CBL/VDL解释图或报告",
     False, "measured", "施工小结/section 10", "定性结论;非定量CBL数据", "low",
     "CBL资料严重缺失:施工小结仅记载'合格'二字;20512和20514中所有JPG文件均为大陆架浮箍浮鞋出厂检测报告;20511中为浅层测井数据(1103-3940m)非尾管段;无独立CBL/VDL报告或测井图"],
    ["ht1_001", "呼探1-001", None, None, None, None, "CBL资料缺失:需人工清查是否另有独立CBL测井报告或固井质量评价图;无定量CBL数据无法参与定量验证",
     False, "missing", "N/A", "N/A", "high",
     "重要缺失:建议咨询新疆油田勘探事业部或测井公司获取HT1-001尾管段CBL/VDL测井资料"],
]
write_csv("cbl_evaluation.csv",
    ["well_id","well_name_cn","cbl_top_md_m","cbl_bottom_md_m","cbl_pass_rate","cbl_quality_class","interpretation_summary","include_in_validation","data_type","source_file","source_location","confidence","notes"],
    cbl_rows)

print("Files 7-10 done.")