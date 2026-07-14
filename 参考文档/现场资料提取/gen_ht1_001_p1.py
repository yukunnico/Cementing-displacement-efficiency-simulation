#!/usr/bin/env python
"""Generate all 16 extraction CSV/MD files for HT1-001 well"""
import csv, os, sys

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

print("Generating HT1-001 extraction files...\n")

# ======= File 1: well_basic_info.csv =======
write_csv("well_basic_info.csv",
    ["well_id","well_name_cn","alias","block","field","well_type","cementing_type","data_source","confidence","notes"],
    [["ht1_001", "呼探1-001", "HT1-001, HT1-001井, 呼探1-001井",
      "准噶尔盆地南缘冲断带霍玛吐背斜带呼西背斜", "新疆油田（天山天然气项目部）",
      "直井", "油层尾管固井/精细控压尾管固井",
      "HT1-001井139.7+168.3mm尾管完井固井施工设计（审核）.doc",
      "high", "四开钻进;设计井深7620m,实际7746m;钻井液体系:油基;构造位置明确，资料齐全"]
    ])

# ======= File 2: well_geometry.csv =======
wgeom = [
    ["ht1_001", "liner_hanger_md", 5460.159, "m", None, None, "design", "施工设计", "section 4.2/6.1", "high", "悬挂器送入工具顶深;喇叭口位置5460.159m"],
    ["ht1_001", "liner_top_md", 5460.159, "m", None, None, "design", "施工设计", "section 6.1/4.2", "high", "尾管顶深=悬挂器位置"],
    ["ht1_001", "liner_bottom_md", 7746.0, "m", None, None, "measured", "施工小结", "section 9", "high", "实际下深7746m"],
    ["ht1_001", "liner_shoe_md", 7746.0, "m", None, None, "design", "施工设计", "section 6.1", "high", "加长引鞋底深=尾管底深"],
    ["ht1_001", "openhole_top_md", 5670.0, "m", None, None, "design", "施工设计", "section 1.2.1", "high", "裸眼段顶深=技术套管下深"],
    ["ht1_001", "openhole_bottom_md", 7746.0, "m", None, None, "measured", "施工小结", "section 1.2.1", "high", "裸眼段底深=完钻井深"],
    ["ht1_001", "simulation_top_md", 5460.159, "m", None, None, "design", "施工设计", "section 4.2", "high", "建议从悬挂器顶开始模拟"],
    ["ht1_001", "simulation_bottom_md", 7746.0, "m", None, None, "design", "施工设计", "section 4.2", "high", "建议到井底=尾管鞋"],
    ["ht1_001", "liner_outer_diameter_168", 168.3, "mm", 5460.159, 7174.941, "design", "施工设计", "section 6.1/6.4.1", "high", "上部复合尾管168.3mm段"],
    ["ht1_001", "liner_outer_diameter_139", 139.7, "mm", 7174.941, 7746.0, "design", "施工设计", "section 6.1/6.4.1", "high", "下部复合尾管139.7mm段"],
    ["ht1_001", "liner_inner_diameter_168", 138.76, "mm", None, None, "calculated", "施工设计", "section 6.4.1", "high", "168.3-2*14.27=139.76mm;壁厚14.27mm(D=168.3-2*14.27)"],
    ["ht1_001", "liner_inner_diameter_139", 107.94, "mm", None, None, "calculated", "施工设计", "section 6.4.1", "high", "139.7-2*15.88=107.94mm;壁厚15.88mm"],
    ["ht1_001", "borehole_nominal_diameter_upper", 241.3, "mm", 5670.0, 7441.0, "design", "施工设计", "section 1.2.1", "high", "上部裸眼钻头尺寸241.3mm"],
    ["ht1_001", "borehole_nominal_diameter_lower", 215.9, "mm", 7441.0, 7746.0, "design", "施工设计", "section 1.2.1", "high", "下部裸眼钻头尺寸215.9mm"],
    ["ht1_001", "drill_bit_size_upper", 241.3, "mm", None, None, "design", "施工设计", "section 1.2.1", "high", ""],
    ["ht1_001", "drill_bit_size_lower", 215.9, "mm", None, None, "design", "施工设计", "section 1.2.1", "high", ""],
    ["ht1_001", "float_collar_1_md", 7733.722, "m", None, None, "design", "施工设计", "section 6.1", "high", "浮箍1底深7733.722m"],
    ["ht1_001", "float_collar_2_md", 7710.781, "m", None, None, "design", "施工设计", "section 6.1", "high", "浮箍2底深7710.781m"],
    ["ht1_001", "ball_seat_md", 7642.938, "m", None, None, "design", "施工设计", "section 6.1", "high", "球座底深=阻位7642.938m"],
    ["ht1_001", "float_shoe_md", 7746.0, "m", None, None, "design", "施工设计", "section 6.1", "high", "加长引鞋底深"],
    ["ht1_001", "overlap_length", 209.841, "m", None, None, "design", "施工设计", "section 4.2", "high", "重叠段5460.159~5670m=209.841m"],
    ["ht1_001", "lower_plug_depth", 103.326, "m", None, None, "design", "施工设计", "section 4.2/6.1", "high", "下塞深度103.326m(浮鞋到球座)"],
    ["ht1_001", "average_caliper_upper", 247.83, "mm", 5670.0, 7441.0, "measured", "施工设计", "section 1.4.2", "high", "5670-7441m平均井径;扩大率2.71%"],
    ["ht1_001", "average_caliper_lower", 229.46, "mm", 7441.0, 7746.0, "measured", "施工设计", "section 1.4.2", "high", "7441-7746m平均井径;扩大率6.28%"],
    ["ht1_001", "max_inclination", 1.7665, "deg", 7665.0, None, "measured", "施工设计", "section 1.4.2", "high", "最大井斜1.7665度@7665m"],
    ["ht1_001", "casing_273_id", 245.42, "mm", 0.0, 5670.0, "calculated", "施工设计", "section 1.2.1", "medium", "273.1-2*13.84=245.42mm"],
    ["ht1_001", "casing_273_od", 273.1, "mm", 0.0, 5670.0, "design", "施工设计", "section 1.2.1", "high", "技术套管外径"],
    ["ht1_001", "liner_hanger_max_od", 235.0, "mm", None, None, "design", "施工设计", "section 5.1", "high", "悬挂器液缸外径235mm;本体外径187mm"],
    ["ht1_001", "casing_liner_separation_md", 7174.941, "m", None, None, "design", "施工设计", "section 6.1", "high", "168.3mm/139.7mm变扣位置"],
    ["ht1_001", "openhole_enlargement_upper", 2.71, "%", None, None, "measured", "施工设计", "section 1.4.2", "high", "5670-7441m扩大率"],
    ["ht1_001", "openhole_enlargement_lower", 6.28, "%", None, None, "measured", "施工设计", "section 1.4.2", "high", "7441-7746m扩大率"],
]
write_csv("well_geometry.csv",
    ["well_id","item","value","unit","top_md_m","bottom_md_m","data_type","source_file","source_location","confidence","notes"],
    wgeom)

print("Files 1-2 done.")