# -*- coding: utf-8 -*-
"""生成泵注程序/流体性能/六速读数 CSV"""
import csv, os

OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003"
WELL = "ht1_003"
DSGN = "HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批).doc"
REC  = "HT1-003井168.3+139.7mm尾管固井施工记录表.doc"
SUM  = "HT1-003井168.3+139.7mm油层尾管固井总结.doc"
HIST = "HT1-003井油层尾管固井作业史.doc"
LAB  = "化验报告/HT1-003 油层尾管 化验报告.docx"

def w(fname, header, rows):
    p = os.path.join(OUT, fname)
    with open(p, "w", encoding="utf-8-sig", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(header)
        wr.writerows(rows)
    print(f"wrote {fname}: {len(rows)} rows")

# ============ 7. pumping_schedule (实际+设计) ============
w("pumping_schedule.csv",
  ["well_id","step_index","stage_name","fluid_role","fluid_name","volume_m3","rate_m3_min","start_time_min","end_time_min","density_g_cm3","event_tag","data_type","source_description","confidence","notes"],
  [
  [WELL,0,"initial_condition","mud","油基钻井液","","","0","0",1.95,"","field_measured",DSGN+" 一.1.3.1","high","井内原钻井液密度1.95g/cm3，粘度62s，塑粘51mPa.s，动切10Pa"],
  [WELL,1,"flush_pipe","shutdown","冲洗管线试压","","","16.5","18.6","","pressure_test_30MPa","field_measured",SUM+"/作业史","high","2026-05-27 16:30-16:38 冲洗管线试压30MPa"],
  [WELL,2,"pump_lead_spacer_pre","wash","先导浆(低密前置液)",28,1.2,18.6,21.0,1.75,"","field_measured",SUM+"/作业史 施工流程","high","17:00注完，泵压14-17.8MPa；泥浆泵泵送"],
  [WELL,3,"pump_spacer1","spacer","低失水驱油隔离液1",16,1.0,21.0,21.25,2.05,"","field_measured",SUM+"/作业史","high","17:15注完，泵压19MPa；实际配隔离液xlsm密度~2.05-2.13总量16m3"],
  [WELL,4,"pump_spacer2","spacer","低失水驱油隔离液2",10,1.0,21.25,21.58,1.95,"","field_measured",SUM+"/作业史","high","17:35注完，泵压18-19MPa；实际配隔离液xlsm密度~2.08"],
  [WELL,5,"pump_lead","lead_cement","领浆",39,1.2,21.58,22.27,2.05,"","field_measured",SUM+"/作业史 七.4","high","18:16注完，泵压22MPa；灰量157t；稠化455min/70Bc"],
  [WELL,6,"pump_tail","tail_cement","尾浆",28,1.2,22.27,22.83,1.95,"","field_measured",SUM+"/作业史 七.4","high","18:50注完，泵压22MPa；稠化256min/70Bc"],
  [WELL,7,"displace","displacement","替浆(钻井液1.95+压塞液2+中置液14)",91.9,1.4,22.83,24.17,1.95,"","field_measured",SUM+"/作业史 七.5","high","18:50-20:10(约80min)排量1.4-0.8m3/min，泵压15-19-26-19MPa；实际替桨xlsm总泵量91.01m3排量峰值1.41m3/min；胶塞试通过钻具内容积61m3与灌水测试一致"],
  [WELL,8,"bump_plug","plug","碰压","","","24.17","","","bump_pressure_19-25MPa","field_measured",REC+"/固井总结","high","20:10碰压19-25MPa，放回水800L断流(正常)"],
  [WELL,9,"pull_center_tube","plug","拆水泥头拔中心杆","","","25.17","","","bump_pressure_14MPa","field_measured",SUM+"/作业史","high","21:10拆水泥头憋压14MPa拔出中心杆"],
  [WELL,10,"circulate_excess","displacement","大排量循环排混浆",200,2.4,25.67,27.17,1.95,"circulate_excess_mud","field_measured",SUM+"/作业史","high","21:40循环排混浆排量2.4m3/min压力28.7MPa控压5MPa循环103方出口密度降低"],
  [WELL,11,"shut_in","shutdown","候凝","","","","","","shut_in_72h","field_measured",REC+"/DSGN 7.2","high","2026-05-28 4:00开始候凝72h；环空井口加压12.9MPa加压54h憋压候凝；设计候凝至领浆:隔离液2=8:2强度3.5MPa"],
  [WELL,100,"design_pump_lead_spacer_pre","wash","先导浆(设计)",28,1.2,"","",1.75,"","design_value",DSGN+" 七.7.2","high","设计排量1.2m3/min泵压12.7-16.1MPa泥浆泵泵送"],
  [WELL,101,"design_pump_spacer1","spacer","低失水驱油隔离液1(设计)",16,1.2,"","",2.05,"","design_value",DSGN+" 七.7.2","high","设计排量1.2(实际1.0)"],
  [WELL,102,"design_pump_spacer2","spacer","低失水驱油隔离液2(设计)",10,1.2,"","",1.95,"","design_value",DSGN+" 七.7.2","high","设计排量1.2(实际1.0)"],
  [WELL,103,"design_pump_lead","lead_cement","领浆(设计)",39,1.2,"","",2.05,"","design_value",DSGN+" 七.7.2","high","设计排量1.2，40min，泵压17-15MPa"],
  [WELL,104,"design_pump_tail","tail_cement","尾浆(设计)",28,1.4,"","",1.95,"","design_value",DSGN+" 七.7.2","high","设计排量1.4-1.2，30min"],
  [WELL,105,"design_disp1","displacement","压塞液(设计)",2,1.0,"","",1.95,"","design_value",DSGN+" 七.7.2","high","排量1.0-1.4"],
  [WELL,106,"design_disp2","displacement","钻井液(设计)",25,1.6,"","",1.95,"","design_value",DSGN+" 七.7.2","high","固井计量罐"],
  [WELL,107,"design_disp3","displacement","保护液(设计)",14,1.4,"","",1.95,"","design_value",DSGN+" 七.7.2","high","批混车混配"],
  [WELL,108,"design_disp4","displacement","钻井液(设计)",8,1.2,"","",1.95,"","design_value",DSGN+" 七.7.2","high",""],
  [WELL,109,"design_disp5","displacement","钻井液(设计)",14,1.0,"","",1.95,"","design_value",DSGN+" 七.7.2","high",""],
  [WELL,110,"design_disp6","displacement","钻井液(设计)",16,0.8,"","",1.95,"","design_value",DSGN+" 七.7.2","high",""],
  [WELL,111,"design_disp7","displacement","钻井液(设计)",12.9,0.7,"","",1.95,"","design_value",DSGN+" 七.7.2","high","末段；理论顶替总量(含后置液)91.9m3"],
  [WELL,112,"design_circulate","displacement","循环排混浆(设计)",200,2.4,"","",1.95,"","design_value",DSGN+" 七.7.2","high","循环排混浆90min；碰压后加回压3.5-6.5MPa"],
  ])
print("  (注: 设计书14.4.2另有一版正注预案: 压塞液3m3密度2.20+基液1m3密度1.02，末段10.9m3；与7.2版略有差异，按7.2主施工流程记录)")

# ============ 8. fluid_properties ============
w("fluid_properties.csv",
  ["well_id","fluid_name","fluid_role","density_g_cm3","rheology_model","pv_mpa_s","yp_pa","yield_stress_pa","consistency_k_pa_sn","flow_index_n","viscosity_mpa_s","test_temperature_c","data_type","source_description","confidence","notes"],
  [
  [WELL,"油基钻井液(下套管前)","mud",1.95,"幂律+宾汉",51,10,"","",0.631,"",65,"field_measured",DSGN+" 一.1.3.1","high","粘度62s；HTHP失水1.2mL/1.5mm@150C；固含40%；氯离子22000mg/L；油水比85:15；破乳电压814；六速122/71/52/21/7/6；初终切4/7.5Pa"],
  [WELL,"油基钻井液(洗井结束前)","mud",1.95,"","",10,"","","","",62,"field_measured",REC+" 洗井","high","屈服值10Pa，密度1.95g/cm3，粘度62s"],
  [WELL,"油基钻井液(化验报告129C)","mud",1.95,"幂律","","","",0.751,0.631,"","field_measured",LAB+" 表7 流变性能","high","129C↓93C六速122/71/52/21/7/6，n=0.631,K=0.751"],
  [WELL,"先导浆","wash",1.75,"宾汉",55,9.2,"","","","","","design_value",DSGN+" 七/设计初稿","medium","低密先导浆28m3；要求粘度<55s塑粘≤30mPa.s"],
  [WELL,"低失水驱油隔离液1","spacer",2.05,"幂律","","","",1.245,0.668,"","field_measured",LAB+" 表7/表6","high","129C↓93C六速239/177/123/89/10/6；稳定性上2.05/中2.06/下2.06；冲洗效率97.9%；配方:水+3.5%悬浮剂+40%驱油冲洗剂+8%降失水剂+1%消泡剂+215%重晶石"],
  [WELL,"低失水驱油隔离液2","spacer",1.95,"幂律","","","",1.245,0.668,"","field_measured",LAB+" 表7/表6","high","与隔离液1同流变；稳定性上1.94/中1.95/下1.95；冲洗效率97.7%；配方+180%重晶石；实际配隔离液xlsm密度~2.08"],
  [WELL,"领浆","lead_cement",2.05,"幂律","","","",1.622,0.597,"","field_measured",LAB+" 表4/表7","high","稠化455min/70Bc(454@40Bc)；API失水37ml；游离液0；沉降0.01；24h强度3.9MPa/48h7.4/7d30.3；7d膨胀率0.06%；7d渗透率0.018e-3um2；顶部静胶凝34h10min；配方:G+35%铁矿粉+35%石英砂+1%BCJ-300S+3%BCE-311S+4%BCY-200S+(5%BCG-200L+4%BCD-210L+3.4%XZ-HN1G+0.3%XZ-XP1L)湿混+52%H2O"],
  [WELL,"尾浆","tail_cement",1.95,"幂律","","","",1.673,0.585,"","field_measured",LAB+" 表5/表7","high","稠化256min/70Bc(253@40Bc)；API失水40ml；游离液0；沉降0.01；24h强度15.1MPa/7d32.5；7d膨胀率0.08%；7d渗透率0.025e-3um2；顶部静胶凝394min底部227min；配方:G+20%铁矿粉+35%石英砂+1%BCJ-300S+3%BCE-311S+4%BCY-200S+(5%BCG-200L+3%BCD-210L+2.3%XZ-HN1G+0.3%XZ-XP1L)湿混+56%H2O"],
  [WELL,"压塞液","plug",1.95,"宾汉",30,8,"","","","","","design_value",DSGN+" 七","medium","驱油隔离液类，2m3"],
  [WELL,"保护液/中置液","displacement",1.95,"宾汉",30,8.2,"","","","","","design_value",DSGN+" 七","medium","14m3批混车混配"],
  [WELL,"基液","displacement",1.02,"宾汉",30,8,"","","","","","design_value",DSGN+" 七/14.4.2","medium","1m3密度1.02g/cm3"],
  [WELL,"替浆钻井液","displacement",1.95,"","",62,"","","","","","field_measured",REC+" 替浆","high","替浆介质钻井液密度1.95g/cm3塑粘62mPa.s(记录表)"],
  ])

# ============ 9. rheometer_readings ============
w("rheometer_readings.csv",
  ["well_id","fluid_name","temperature_c","theta_600","theta_300","theta_200","theta_100","theta_6","theta_3","source_description","confidence","notes"],
  [
  [WELL,"油基钻井液",65,122,71,52,21,7,6,DSGN+" 一.1.3.1","high","常温性能(取样2026.5.18)"],
  [WELL,"油基钻井液",30,162,98,67,42,11,9,DSGN+" 一.1.3.4 低温性能","high","05.16 AV81/PV68/YP13/初6终10"],
  [WELL,"油基钻井液",40,154,91,62,40,9,8,DSGN+" 一.1.3.4","high","AV77/PV65/YP12/初5终9"],
  [WELL,"油基钻井液",50,142,82,54,35,8,7,DSGN+" 一.1.3.4","high","AV71/PV62/YP9/初4终8"],
  [WELL,"油基钻井液(129C↓93C)",93,122,71,52,21,7,6,LAB+" 表7","high","n=0.631 K=0.751"],
  [WELL,"低失水驱油隔离液(129C↓93C)",93,239,177,123,89,10,6,LAB+" 表7","high","n=0.668 K=1.245"],
  [WELL,"领浆(129C↓93C)",93,207,121,93,60,15,9,LAB+" 表7","high","n=0.597 K=1.622"],
  [WELL,"尾浆(129C↓93C)",93,196,118,89,56,14,10,LAB+" 表7","high","n=0.585 K=1.673"],
  ])
print("PUMP/FLUID DONE")
