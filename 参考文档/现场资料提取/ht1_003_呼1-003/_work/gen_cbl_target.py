# -*- coding: utf-8 -*-
"""生成 CBL评价/目标层/温压/施工异常 CSV"""
import csv, os

OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003"
WELL = "ht1_003"
DSGN = "HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批).doc"
REC  = "HT1-003井168.3+139.7mm尾管固井施工记录表.doc"
SUM  = "HT1-003井168.3+139.7mm油层尾管固井总结.doc"
HIST = "HT1-003井油层尾管固井作业史.doc"
CBL  = "HT1-003_固井质量测井评价图_15-3500与_3500-7595_20260622_完井.pdf 图头"

def w(fname, header, rows):
    p = os.path.join(OUT, fname)
    with open(p, "w", encoding="utf-8-sig", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(header)
        wr.writerows(rows)
    print(f"wrote {fname}: {len(rows)} rows")

# ============ 10. cbl_evaluation ============
w("cbl_evaluation.csv",
  ["well_id","well_name_cn","cbl_top_md_m","cbl_bottom_md_m","cbl_pass_rate","cbl_quality_class","interpretation_summary","include_in_validation","data_type","source_description","confidence","notes"],
  [
  [WELL,"呼1-003",5307.54,7618,None,"不合格","固井质量测井评价图(中油测井新疆分公司，FITS装备，LEAD5.0软件)结论：固井质量不合格，违反七条红线规定两条：(1)上层套管鞋、尾管重合段及其以上25m环空范围内，固井水泥一界面胶结质量未达到连续胶结中等及以上；(2)油气水层段内固井水泥一、二界面无连续胶结中等及以上。",1,"field_measured","HT1-003_固井质量测井评价图_15-3500与_3500-7595_20260622_完井.pdf 图头","high","官方解释结论文字抄录。无数值合格率给出。测量井段15-7595m；电测井深7595m；15-5568m多层套管VDL第二界面不作评价"],
  [WELL,"呼1-003",15.0,5568.0,None,"不作评价(多层套管)","15-5568m为多层套管(回接193.68mm内技术套管)，VDL第二界面不作评价；一界面正常评价",1,"field_measured",CBL+" 技术说明","high","电测井深7595m"],
  [WELL,"呼1-003",5568.0,7618.0,None,"不合格(尾管段)","尾管段168.28/139.7mm 5316.036-7618m，水泥返至5307.540m；目的井段7442-7618m；固井质量不合格，油气水层段一、二界面无连续胶结中等及以上",1,"field_measured","HT1-003_固井质量测井评价图_3500-7595 图头","medium","整体结论文字抄录；未给出分井段合格率数值"],
  [WELL,"呼1-003",None,None,None,"待人工数字化","CBL定量合格率(分段胶结质量/合格率百分比)图上无文字结论，需人工数字化声幅/变密度曲线后方可定量；分段深度-质量评价留待人工复核",1,"pending_digitization","HT1-003_固井质量测井评价图_15-3500 / _3500-7595 / IBC评价成果图_10-7520m","low","本次环境无图像输入能力，未对曲线进行像素级定量；深度段范围按文件名(15-3500/3500-7595/10-7520m)填写"],
  ])

# ============ 11. target_intervals ============
w("target_intervals.csv",
  ["well_id","interval_name","top_md_m","bottom_md_m","purpose","source_description","data_type","confidence","notes"],
  [
  [WELL,"地质分层-连木沁组",None,5953,"地层划分",DSGN+" 二.2.1","field_measured","high","底深5953m厚661m"],
  [WELL,"地质分层-胜金口组",None,6076,"地层划分",DSGN+" 二.2.1","field_measured","high","底深6076m厚123m"],
  [WELL,"地质分层-呼图壁河组",None,6747,"地层划分",DSGN+" 二.2.1","field_measured","high","底深6747m厚671m"],
  [WELL,"地质分层-清水河组",None,7502,"地层划分",DSGN+" 二.2.1","field_measured","high","底深7502m厚755m"],
  [WELL,"地质分层-喀拉扎组三段",7586,None,"地层划分",DSGN+" 二.2.1","field_measured","high","顶7586m厚84m"],
  [WELL,"地质分层-喀拉扎组二段",7618,None,"地层划分",DSGN+" 二.2.1","field_measured","high","7618m未穿厚32m"],
  [WELL,"水层-K1q",6860,6880,"高压水层(压稳关键)",DSGN+" 二.2.3","field_measured","high","全烃1.4574%↑3.8333%；水层压力最高，<1.98g/cm3当量无法完全压稳"],
  [WELL,"干层-K1q",7110,7114,"干层",DSGN+" 二.2.3","field_measured","high","全烃0.7859%↑1.6560%"],
  [WELL,"气层-K1q",7418,7440,"气层(主力油层段)",DSGN+" 二.2.3","field_measured","high","全烃5.9631%↑9.3621%；主力油层段7418-7440m"],
  [WELL,"差气层-J3k3",7451,7463,"差气层",DSGN+" 二.2.3","field_measured","high","全烃4.2530%↑8.2618%"],
  [WELL,"差气层-J3k3",7474,7485,"差气层",DSGN+" 二.2.3","field_measured","high","全烃4.6407%↑6.3640%"],
  [WELL,"差气层-J3k3",7490,7498,"差气层",DSGN+" 二.2.3","field_measured","high","全烃4.5346%↑7.1647%"],
  [WELL,"气层-J3k3",7536,7554,"气层",DSGN+" 二.2.3","field_measured","high","全烃6.2362%↑12.5426%"],
  [WELL,"差气层-J3k3",7564,7576,"差气层",DSGN+" 二.2.3","field_measured","high","全烃3.428%↑20.1342%；施工记录表油气层底7576m"],
  [WELL,"主力油层段",6860,7522,"主力油层段",DSGN+" 二.2.3","interpreted","medium","6860-6880m(水层压力最高)、7418-7440m、7506-7522m"],
  [WELL,"目的井段",7442,7618,"CBL评价目的井段",CBL+" 目的井段","field_measured","high","CBL评价图目的井段7442-7618m；作业史油层顶界7418m/射孔底界7578m"],
  [WELL,"漏失层-喀拉扎组",7463,None,"钻进漏失层",DSGN+" 一.1.6","field_measured","high","7463.71m漏失漏速16.6m3/h漏失3.6m3"],
  ])

# ============ 12. temperature_pressure_profile ============
w("temperature_pressure_profile.csv",
  ["well_id","md_m","temperature_c","pressure_mpa","data_type","source_description","confidence","notes"],
  [
  [WELL,0,60,"","field_measured",DSGN+" 二.2.6.1","high","出口温度60C(井底循环洗井两周后)"],
  [WELL,5290,123,"","field_measured",DSGN+" 二.2.6.1","high","悬挂器温度123C(电测)；作业史悬挂器静止/循环123/96C"],
  [WELL,5305,152,"","field_measured",DSGN+" 九/作业史","high","实验温度5305-6500m静止温度152C循环129.2C；领浆温度系数0.85"],
  [WELL,6500,152,"","field_measured",DSGN+" 九","high","6500-7618m循环温度129.2C"],
  [WELL,7618,150,"","field_measured",HIST+" 八","high","作业史：井底静止/循环温度150/129C；设计电测152C"],
  [WELL,7463,None,None,"field_measured",DSGN+" 三.3.3","high","漏层7463m：压稳当量≥1.98施工上限当量≤2.042g/cm3"],
  [WELL,6880,None,None,"field_measured",DSGN+" 三.3.3","high","水层6880m：压稳当量≥1.98施工上限当量≤2.157g/cm3"],
  [WELL,7618,None,None,"field_measured",DSGN+" 三.3.3","high","井底7618m：压稳当量≥1.98施工上限当量≤2.04g/cm3；地层压力系数1.94-1.967；环空静压145.96MPa/管内145.72MPa(设计)；施工控制井底当量1.98-2.04g/cm3"],
  ])

# ============ 13. construction_events ============
w("construction_events.csv",
  ["well_id","event_time_min","event_md_m","event_type","description","related_stage","source_description","confidence","notes"],
  [
  [WELL,None,7463.71,"lost_circulation","钻进期间井漏：7463.71m喀拉扎组，排量27.5-28L/s泵压36-38MPa，密度2.00g/cm3，出口流量29.1%↓21.7%，漏失3.6m3漏速16.6m3/h；密度降至1.95并泵入堵漏剂解除","","field_measured",DSGN+" 一.1.6","high","钻进复杂"],
  [WELL,None,7100,"water_influx_circulation","下套管至7100m循环期间：顶通压力10MPa，顶通后瞬间降至4MPa，反推环空出水，液柱压力降低导致泵压下降","","field_measured",SUM+" (二)/作业史 六","high","水层出水，循环泵压远低于控压"],
  [WELL,None,200,"water_influx_running","下套管至200m顶通：出口顶替出6方低密度钻井液，环空液柱压力降低0.01g/cm3，水层压稳失衡出水；低密度钻井液排出后返浆正常","","field_measured",SUM+" (三)固井分析","medium","总结建议：顶通及灌浆均用重浆确保压稳"],
  [WELL,None,None,"pump_ball_seat","投球坐挂：2026-05-27 9:40投硫化铜球，9:46开泵送球(排量0.57m3/min泵压7MPa)，10:10球到位憋压12MPa悬重257→135t下放7.1m坐挂成功，18MPa憋通球座","","field_measured",SUM+"/作业史 八(4)","high","大陆架一体式平衡液缸尾管悬挂器，施工正常"],
  [WELL,None,None,"release","倒扣丢手：14:00停泵调整悬重160t，正转30圈基本无回转，上提180t继续上提0.7m悬重不变，下压30t固井","","field_measured",SUM+"/作业史 八(4)","high","丢手成功"],
  [WELL,None,None,"bump_pressure","碰压19-25MPa，放回水800L断流，放回水正常","","field_measured",REC+"/固井总结","high","正常碰压"],
  [WELL,None,None,"circulate_high_pressure","循环排混浆：排量2.4m3/min压力28.7MPa，控压5MPa，循环103方出口密度降低","","field_measured",SUM+"/作业史 七.6","high","碰压后循环排混浆"],
  [WELL,None,None,"shut_in_pressure","候凝：2026-05-28 4:00开始，环空井口加压12.9MPa加压54h，憋压候凝72h","","field_measured",REC+" 候凝","high","设计候凝至领浆:隔离液2=8:2混浆强度3.5MPa后开井"],
  [WELL,None,None,"mud_low_return","下套管至7100m时泵压异常低(3-4MPa控压6-7MPa)，判断水柱已形成；期间活动上提无法提开，出水循环至管内后继续下套管保证到底","","field_measured",SUM+" (三)固井分析","medium","水层出水导致环空液柱压力异常"],
  ])
print("CBL/TARGET/TP/EVENTS DONE")
