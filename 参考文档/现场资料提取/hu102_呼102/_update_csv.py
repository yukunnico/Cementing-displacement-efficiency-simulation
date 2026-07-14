import csv, os

out_dir = 'C:/Users/katsura/Desktop/research/控压固井项目/现场资料提取/hu102_呼102'
os.makedirs(out_dir, exist_ok=True)

# ===================================================================
# 1. well_basic_info.csv
# ===================================================================
with open(os.path.join(out_dir, 'well_basic_info.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','well_name_cn','alias','block','field','well_type','cementing_type','target_formation','drilling_rig','drilling_company','cementing_company','total_depth_m','cbl_qualified_rate_pct','cbl_evaluation','data_source','confidence','notes'])
    w.writerow([
        'hu102','呼102','呼102井','准噶尔盆地南缘呼图壁','呼图壁气田','定向井','尾管固井(139.7mm)',
        'J3k(侏罗系喀拉扎组)','90013','西部钻探','固井准东公司',
        '7735','66.65','不合格(低于70%)',
        '20211.doc(固井小结); 20212.doc(施工设计); 20214.doc(施工总结); 20216.doc(施工记录表); 100413.PDF(CBL评价图)',
        'high',
        '139.7mm油层尾管固井6821.895-7735m; CBL合格率66.65%, 不合格段7405-7480m及7502-7540m; 下部结构单流阀失效,施工方案临场调整; 施工日期2022-11-21'
    ])
print("1. well_basic_info.csv")

# ===================================================================
# 2. well_geometry.csv
# ===================================================================
with open(os.path.join(out_dir, 'well_geometry.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','item','value','unit','top_md_m','bottom_md_m','data_type','source_file','source_location','confidence','notes'])
    rows = [
        ['hu102','total_depth_md','7735','m',None,None,'measured','20211.doc','section 1.1','high','完钻井深'],
        ['hu102','liner_hanger_top_md','6821.895','m',None,None,'measured','20211.doc','pipe string','high','悬挂器坐挂顶深/送入工具'],
        ['hu102','liner_hanger_md','7665.516','m',None,None,'measured','20211.doc','section 1.1','high','悬挂器坐挂位置'],
        ['hu102','liner_hanger_od','149.2','mm',None,None,'measured','20211.doc','section 1.9','high','FHDS550型悬挂器'],
        ['hu102','liner_shoe_md','7735','m',None,None,'measured','20211.doc','section 1.1','high','尾管鞋深'],
        ['hu102','liner_top_md','6821.895','m',None,None,'measured','20211.doc','pipe string','high','送入工具顶部'],
        ['hu102','liner_bottom_md','7735','m',None,None,'measured','20211.doc','section 1.1','high','尾管底深'],
        ['hu102','openhole_top_md','6720','m',None,None,'measured','20211.doc','section 1.1','high','裸眼段顶深(设计)'],
        ['hu102','openhole_bottom_md','7735','m',None,None,'measured','20211.doc','section 1.1','high','裸眼段底深'],
        ['hu102','simulation_top_md','6820','m',None,None,'derived','20211.doc','derived','medium','建议模拟段顶深'],
        ['hu102','simulation_bottom_md','7735','m',None,None,'derived','20211.doc','derived','medium','建议模拟段底深'],
        ['hu102','liner_length_md','913.105','m',None,None,'measured','20211.doc','section 1.9','high','尾管总长'],
        ['hu102','liner_od','139.7','mm',None,None,'measured','20211.doc','section 1.1','high','尾管外径'],
        ['hu102','liner_id_heavy','108.1','mm',None,None,'measured','20211.doc','pipe specs','high','15.8mm壁厚段内径(底部50根)'],
        ['hu102','liner_id_light','111.16','mm',None,None,'measured','20211.doc','pipe specs','high','14.27mm壁厚段内径(上部36根)'],
        ['hu102','drill_bit_size','190.5','mm',None,None,'measured','20211.doc','section 1.1','high','钻头尺寸'],
        ['hu102','borehole_nominal_diameter','190.5','mm',None,None,'measured','20211.doc','section 1.8','high','名义井径'],
        ['hu102','float_collar_md','6826.74','m',None,None,'measured','20211.doc','construction description','high','实探喇叭口/浮箍位置'],
        ['hu102','float_collar_design_md','6823.095','m',None,None,'measured','20211.doc','section 1.9','high','设计喇叭口位置'],
        ['hu102','float_shoe_md','6823.095','m',None,None,'measured','20211.doc','construction summary','high','浮鞋深度'],
        ['hu102','upper_casing_od','219.1','mm',None,None,'measured','20211.doc','section 1.1','high','上层技术套管外径(悬挂器以上)'],
        ['hu102','cbl_measure_top_md','6840','m',None,None,'measured','100413.PDF','CBL header','high','CBL测量段顶深'],
        ['hu102','cbl_measure_bottom_md','7665','m',None,None,'measured','100413.PDF','CBL header','high','CBL测量段底深'],
        ['hu102','designed_bottom_plug','7585','m',None,None,'measured','20211.doc','section 6','high','设计下塞面'],
        ['hu102','actual_bottom_plug','7447','m',None,None,'measured','20211.doc','section 6','high','实探下塞面'],
        ['hu102','designed_plug_length','150','m',None,None,'measured','20211.doc','section 6','high','设计塞长'],
        ['hu102','max_inclination','11.49','deg',None,None,'measured','20211.doc','section 1.6','high','最大井斜@7490m'],
        ['hu102','bottomhole_static_temp','149','degC',None,None,'measured','20211.doc','section 1.3','high','井底静止温度(BHST)'],
        ['hu102','cement_test_temp','134','degC',None,None,'derived','20211.doc','section 1.3','high','水泥浆试验温度(0.9x温度系数)'],
        ['hu102','safe_window','0.066','g/cm3',None,None,'calculated','20211.doc','section 1.7','high','井下安全窗口'],
    ]
    for r in rows:
        w.writerow(r)
print("2. well_geometry.csv")

# ===================================================================
# 3. casing_liner_string.csv
# ===================================================================
with open(os.path.join(out_dir, 'casing_liner_string.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','component','top_md_m','bottom_md_m','od_mm','id_mm','grade','connection','wall_mm','quantity','cumulative_length_m','manufacturer','spacing_m','data_type','source_file','confidence','notes'])
    rows = [
        ['hu102','float_shoe','7735','7735','139.7','108.1','TP140V','TP-G2','15.8','1','0.961','大陆架',None,'measured','20211.doc','high','浮鞋'],
        ['hu102','liner_pipe','7734.039','7722.899','139.7','108.1','TP140V','TP-G2','15.8','1','11.14','天钢',None,'measured','20211.doc','high','套管第1根(底部)'],
        ['hu102','float_collar_1','7711.648','7711.313','139.7','108.1','TP140V','TP-G2','15.8','1','0.335','大陆架',None,'measured','20211.doc','high','浮箍#1'],
        ['hu102','liner_pipe','7711.313','7688.649','139.7','108.1','TP140V','TP-G2','15.8','1','22.664','天钢',None,'measured','20211.doc','high','第3根'],
        ['hu102','float_collar_2','7688.649','7688.314','139.7','108.1','TP140V','TP-G2','15.8','1','0.335','大陆架',None,'measured','20211.doc','high','浮箍#2'],
        ['hu102','liner_pipe','7688.314','7665.781','139.7','108.1','TP140V','TP-G2','15.8','2','22.533','天钢',None,'measured','20211.doc','high','第5-6根'],
        ['hu102','ball_seat','7665.781','7665.516','139.7','108.1','TP140V','TP-G2','15.8','1','0.265','大陆架',None,'measured','20211.doc','high','球座'],
        ['hu102','liner_pipe_heavy','7665.516','7575.695','139.7','108.1','TP140V','TP-G2','15.8','8','89.821','天钢',None,'measured','20211.doc','high','底部厚壁段(15.8mm), 50根总长563.205m'],
        ['hu102','marker_joint_1','7575.695','7574.365','139.7','108.1','TP140V','TP-G2','15.8','1','1.33','天钢',None,'measured','20211.doc','high','标节1'],
        ['hu102','liner_pipe_heavy2','7574.365','7382.662','139.7','108.1','TP140V','TP-G2','15.8','17','191.703','天钢',None,'measured','20211.doc','high','中部厚壁段'],
        ['hu102','marker_joint_2','7382.662','7381.333','139.7','108.1','TP140V','TP-G2','15.8','1','1.329','天钢',None,'measured','20211.doc','high','标节2'],
        ['hu102','liner_pipe_heavy3','7381.333','7099.652','139.7','108.1','TP140V','TP-G2','15.8','25','281.681','天钢',None,'measured','20211.doc','high','上部厚壁段'],
        ['hu102','liner_pipe_thin','7099.652','6826.015','139.7','111.16','TP140V','TP-G2','14.27','25','273.637','天钢',None,'measured','20211.doc','high','上部薄壁段(14.27mm), 36根总长273.637m'],
        ['hu102','liner_hanger','6826.015','6823.095','139.7',None,'TP140V','TP-G2',None,'1','2.92','大陆架',None,'measured','20211.doc','high','悬挂器(FHDS550)'],
        ['hu102','tieback_sleeve','6823.095','6821.895','139.7',None,'TP140V','TP-G2',None,'1','1.2','大陆架',None,'measured','20211.doc','high','回接筒'],
        ['hu102','sending_tool','6821.895','6820.475','127',None,None,'NC50',None,'1','1.42',None,None,'measured','20211.doc','high','送入工具'],
        ['hu102','cross_over','6820.475','6819.367',None,None,None,'DSTJ50x411',None,'1','1.108',None,None,'measured','20211.doc','high','变扣'],
        ['hu102','drill_pipe','6819.367','5101.247','114.3','97.18','S135','DSTJ50','8.56','180','1718.12',None,None,'measured','20211.doc','high','送入钻具'],
        ['hu102','cross_over_2','5101.247','5100.154',None,None,None,'DSTJ50xFHDS550',None,'1','0.93',None,None,'measured','20211.doc','high','变扣2'],
        ['hu102','drill_collar','5100.154','0','149.2',None,'V150','FHDS550','9.65','532','5101.251',None,None,'measured','20211.doc','high','钻铤(含钻余-1.097)'],
        # Centralizers
        ['hu102','centralizer','6820','7350',None,None,'弹性扶正器','整体式','190.5','23',None,None,'22','design','20211.doc','high','整体弹性扶正器,外径190.5mm,间距22m'],
        ['hu102','centralizer','7350','7560',None,None,'弹性扶正器','整体式','190.5','19',None,None,'11','design','20211.doc','high','整体弹性扶正器,间距19m'],
        ['hu102','centralizer','7560','7735',None,None,'弹性扶正器','整体式','190.5','8',None,None,'8','design','20211.doc','high','密集布置(8m间距),共55个扶正器'],
    ]
    for r in rows:
        w.writerow(r)
print("3. casing_liner_string.csv")

# ===================================================================
# 4. caliper_profile.csv (combined with inclination)
# ===================================================================
with open(os.path.join(out_dir, 'caliper_profile.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','md_m','caliper_mm','inclination_deg','azimuth_deg','data_type','source_file','confidence','notes'])
    caliper_data = [
        [7120,193.4,4.82,7.3],[7130,193.62,4.89,7.26],[7140,188.18,5.39,7.4],[7150,188.54,6.01,7.1],
        [7160,188.62,6.48,6.81],[7170,188.28,6.67,6.68],[7180,190.24,7.02,6.91],[7190,191.01,7.07,6.56],
        [7200,192.16,7.38,6.5],[7210,192.14,7.38,6.31],[7220,194.04,7.58,6.23],[7230,190.22,7.65,6.19],
        [7240,190,7.84,6.75],[7250,189.49,8.29,6.8],[7260,189.29,8.43,6.2],[7270,188.9,8.51,6.32],
        [7280,191.28,8.37,6.16],[7290,188.81,8.21,5.91],[7300,189.67,8.35,6.54],[7310,188.89,8.6,6.21],
        [7320,190.23,8.96,6.27],[7330,190.57,9.29,6.17],[7340,189.74,9.44,5.98],[7350,191.36,9.74,5.94],
        [7360,193.59,9.82,5.93],[7370,194.47,9.8,5.94],[7380,195.17,10.62,5.81],[7390,196.07,10.53,5.76],
        [7400,194.36,10.73,5.82],[7410,194.61,10.96,5.78],[7420,195.17,11.09,6.1],[7430,196.89,11.33,5.92],
        [7440,200.58,11.44,6.23],[7450,193.4,11.49,8.05],[7460,194.74,10.39,6.45],[7470,192.77,8.7,6.64],
        [7480,184.85,6.99,7.16],[7490,185.08,5.94,6.68],[7500,187.88,6.06,7.24],[7510,194.08,5.35,6.54],
        [7520,191.24,4.63,6.71],[7530,189.35,4.53,7.35],[7540,186.73,4.53,6.95],[7550,188.22,4.6,7.14],
        [7560,189.35,4.68,7.47],[7570,198.72,4.34,7.53],[7580,194.26,4.21,6.58],[7590,187.92,4.21,6.91],
        [7600,187.85,4.24,7.04],[7610,187.22,3.76,7.04],[7620,187.74,3.97,7.09],[7630,194.96,3.83,7.43],
        [7640,188.03,3.42,7.24],[7650,187.68,3.73,7.11],[7660,188.75,3.47,7],[7670,188.14,3.77,7.53],
        [7680,188.41,4.33,7.36],[7690,189.94,3.71,6.98],[7700,191.94,4.16,7.38],[7710,190.32,3.74,6.73],
        [7720,190.31,3.73,7.17],[7735,191.36,3.7,7.5],
    ]
    for d in caliper_data:
        w.writerow(['hu102', d[0], d[1], d[2], d[3], 'measured', '20211.doc(section 1.6)', 'high', '电测井径+井斜(10m间隔)'])
print("4. caliper_profile.csv")

# ===================================================================
# 5. inclination_profile.csv
# ===================================================================
with open(os.path.join(out_dir, 'inclination_profile.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','md_m','inclination_deg','azimuth_deg','data_type','source_file','confidence','notes'])
    for d in caliper_data:
        w.writerow(['hu102', d[0], d[2], d[3], 'measured', '20211.doc(section 1.6)', 'high', ''])
    w.writerow(['hu102', '6820', None, None, 'derived', '20211.doc', 'medium', '尾管悬挂器顶(无测斜)'])
    w.writerow(['hu102', '7490', '11.49', '8.05', 'measured', '20211.doc', 'high', '最大井斜点'])
    w.writerow(['hu102', '7735', '3.7', '7.5', 'measured', '20211.doc', 'high', '井底'])
print("5. inclination_profile.csv")

# ===================================================================
# 6. pumping_schedule.csv
# ===================================================================
with open(os.path.join(out_dir, 'pumping_schedule.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','stage','fluid_type','density_gcm3','volume_m3','pump_rate_m3_min','pressure_mpa','start_time','end_time','duration_min','cumulative_strokes','back_pressure_mpa','data_type','source_file','confidence','notes'])
    schedule = [
        ['hu102','1.前置液','平衡液','1.90','15','0.3-0.6',None,'12:00','12:43','43',None,'0.0','measured','20211.doc','high','原方案:平衡液'],
        ['hu102','异常','',None,None,None,None,'14:00',None,None,None,'4.2','event','20211.doc','high','管内外联通,下部结构失效;放回水不断流,单流阀失效'],
        ['hu102','2.隔离液(试)','隔离液','2.05','5','0.6',None,'14:12','14:17','5',None,None,'measured','20211.doc','high','泵压15.5MPa'],
        ['hu102','方案调整','',None,None,None,None,'17:35',None,None,None,None,'event','20211.doc','high','尾管内全部替入2.10g/cm3隔离液,施工结束立即关井,不留上塞'],
        ['hu102','3.隔离液','隔离液','2.05','15','0.6-0.8','13-15.7','17:56','17:57','1',None,'0.0','measured','20211.doc','high',''],
        ['hu102','4.领浆','领浆(水泥)','2.10','10','0.6-0.8','15.7','18:10','18:25','15',None,'0.0','measured','20211.doc','high','封固段6720-7300m'],
        ['hu102','5.尾浆','尾浆(水泥)','2.10','7','0.6-0.8','15.5','18:25','18:26','1',None,'0.0','measured','20211.doc','high','封固段7300-7735m'],
        ['hu102','6.下胶塞','',None,None,None,None,'18:26',None,None,None,'1.9','event','20211.doc','high','停泵下胶塞'],
        ['hu102','7.压塞液','压塞液','2.10','7','0.3-0.8','6-13','18:41','18:50','9',None,'0.0','measured','20211.doc','high',''],
        ['hu102','8.后置液','保护液','1.90','6','0.6-0.8','13','18:50','18:53','3',None,'0.0','measured','20211.doc','high',''],
        ['hu102','9.替浆','钻井液(油基)','2.02','72','0.9-0.46','17.2-9.9','18:57','20:25','88','3944_accum','0.0','measured','20211.doc','high','74m3替浆;排量逐步从0.9降至0.46; 19:30平衡液进环空; 19:44隔离液进环空; 20:02领浆进环空; 20:13胶塞重合; 20:18尾浆进环空'],
        ['hu102','10.关井','',None,None,None,None,'20:45',None,'20',None,'4.2','event','20211.doc','high','停泵关井套压4.2MPa'],
        ['hu102','11.拔中心管','',None,None,None,None,'21:00',None,None,None,'4.2','event','20211.doc','high','控压4.2MPa上提7m拔出中心管,上提悬重220t'],
        ['hu102','12.循环排混浆','钻井液','2.02','41','1.3','20','22日00:30',None,None,None,'0.0','measured','20211.doc','high','敞压循环排除混浆41方,见隔离液、水泥浆混浆'],
        ['hu102','13.控压循环','钻井液','2.02',None,'1.3','21.9','22日11:00',None,None,None,'2.2','measured','20211.doc','high','控压2.2MPa循环'],
        ['hu102','14.关井候凝','',None,None,None,None,'24日21:00',None,None,None,'5.8-3.7','event','20211.doc','high','加回压候凝,5.8降至3.7MPa'],
    ]
    for r in schedule:
        w.writerow(r)
print("6. pumping_schedule.csv")

# ===================================================================
# 7. fluid_properties.csv
# ===================================================================
with open(os.path.join(out_dir, 'fluid_properties.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','fluid_type','property','value','unit','test_temp_c','test_pressure_mpa','data_type','source_file','confidence','notes'])
    props = [
        # Drilling fluid (oil-based mud)
        ['hu102','钻井液(油基)','密度','2.02','g/cm3',None,None,'measured','20211.doc(s1.2)','high','柴油基泥浆'],
        ['hu102','钻井液(油基)','粘度','73','s',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','塑性粘度PV','66','mPa.s',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','屈服值YP','10','Pa',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','HTHP失水(150C)','2','mL','150',None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','HTHP泥饼(150C)','1.0','mm','150',None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','油水比','85:15','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','固相含量','39','v%',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','含砂','0.3','v%',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','氯离子','14000','mg/L',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','钙离子','8600','mg/L',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','初切','3.5','Pa',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','终切','7','Pa',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','破乳电压','650','V',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi600','142','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi300','82','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi200','58','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi100','34','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi6','5','',None,None,'measured','20211.doc(s1.2)','high',''],
        ['hu102','钻井液(油基)','phi3','4','',None,None,'measured','20211.doc(s1.2)','high',''],
        # Lead slurry (领浆) - lab test from 20234.doc
        ['hu102','领浆(水泥浆)','试验密度','2.20','g/cm3','126','150','measured','20234.doc','high','实验室配方试验'],
        ['hu102','领浆(水泥浆)','施工密度','2.10','g/cm3',None,None,'design','20211.doc','high','实际施工密度(封固6720-7300m)'],
        ['hu102','领浆(水泥浆)','配方','G+12%WG+75%磁铁矿粉+30%SiO2+3%SUP+8%ST900L+2%SXY-2+1.5%ST400S+3.8%ST300R+0.3%ST500L','',None,None,'measured','20234.doc','high','湿混'],
        ['hu102','领浆(水泥浆)','水灰比','52','%',None,None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','造浆率','59','%',None,None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','稠化时间','371','min','126','150','measured','20234.doc','high','初始20Bc-终止70Bc; 升温170min'],
        ['hu102','领浆(水泥浆)','phi300','163','','93',None,'measured','20234.doc','high','ZNN-D6B'],
        ['hu102','领浆(水泥浆)','phi200','128','','93',None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','phi100','106','','93',None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','phi6','10','','93',None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','phi3','6','','93',None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','流变模式','幂律模式','',None,None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','n','0.737','',None,None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','K','0.947','Pa.s^n',None,None,'measured','20234.doc','high',''],
        ['hu102','领浆(水泥浆)','API失水(30min)','42','mL','126','6.9','measured','20234.doc','high','CHANDLER4222'],
        ['hu102','领浆(水泥浆)','24h抗压强度','17.9','MPa','148','20.7','measured','20234.doc','high','4207D'],
        ['hu102','领浆(水泥浆)','游离液','0','%','93',None,'measured','20234.doc','high','120min静置'],
        ['hu102','领浆(水泥浆)','48h顶部强度','9.2','MPa','115','20.7','measured','20234.doc','high',''],
        # Tail slurry (尾浆)
        ['hu102','尾浆(水泥浆)','试验密度','2.20','g/cm3','126','150','measured','20234.doc','high','实验室配方试验'],
        ['hu102','尾浆(水泥浆)','施工密度','2.10','g/cm3',None,None,'design','20211.doc','high','实际施工密度(封固7300-7735m)'],
        ['hu102','尾浆(水泥浆)','配方','G+12%WG+75%磁铁矿粉+30%SiO2+3%SUP+8%ST900L+2%SXY-2+1.5%ST400S+2.4%ST300R+0.3%ST500L','',None,None,'measured','20234.doc','high','ST300R减至2.4%'],
        ['hu102','尾浆(水泥浆)','稠化时间','255','min','126','150','measured','20234.doc','high','初始19Bc-终止70Bc'],
        ['hu102','尾浆(水泥浆)','phi300','163','','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','phi200','128','','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','phi100','106','','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','phi6','10','','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','phi3','6','','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','n','0.737','',None,None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','K','0.947','Pa.s^n',None,None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','API失水(30min)','41','mL','126','6.9','measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','24h抗压强度','21.3','MPa','148','20.7','measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','游离液','0','%','93',None,'measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','温度高点稠化(132C)','218','min','132','150','measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','密度高点稠化(2.23)','227','min','126','150','measured','20234.doc','high',''],
        ['hu102','尾浆(水泥浆)','静胶凝发挥时间','260','min',None,None,'measured','20234.doc','high',''],
        # Contamination
        ['hu102','混浆(7:1:2)','稠化时间','>390','min',None,None,'measured','20234.doc','high','领浆:隔离液:先导浆=7:1:2, 未稠'],
        ['hu102','混浆(7:2:1)','稠化时间','>390','min',None,None,'measured','20234.doc','high','领浆:隔离液:先导浆=7:2:1, 未稠'],
        ['hu102','混浆(1:1:1:1)','稠化时间','>390','min',None,None,'measured','20234.doc','high','领浆:隔离液:先导浆:泥浆=1:1:1:1, 未稠'],
        # Spacer fluids (actual)
        ['hu102','前置液(平衡液)','密度','1.90','g/cm3',None,None,'measured','20211.doc','high',''],
        ['hu102','隔离液','密度','2.05','g/cm3',None,None,'measured','20211.doc','high',''],
        ['hu102','后置液(保护液)','密度','1.90','g/cm3',None,None,'measured','20211.doc','high',''],
        ['hu102','压塞液','密度','2.10','g/cm3',None,None,'measured','20211.doc','high',''],
    ]
    for r in props:
        w.writerow(r)
print("7. fluid_properties.csv")

# ===================================================================
# 8. cbl_evaluation.csv
# ===================================================================
with open(os.path.join(out_dir, 'cbl_evaluation.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','md_top_m','md_bottom_m','interval_length_m','bond_quality','cbl_amplitude_pct','bond_index','vd1_description','qualified','data_type','source_file','confidence','notes'])
    rows = [
        ['hu102','6840','7119.8','279.8','非评价段(双层套管)',None,None,'双层套管,二界面不做评价',None,'measured','100413.PDF','high','6840-7119.8m为双层套管段'],
        ['hu102','7119.8','7405','285.2','胶结良好-中等','<15',None,'CBL低值段,胶结合格','yes','measured','100413.PDF','medium','良好-中等胶结'],
        ['hu102','7405','7480','75','胶结差','高',None,'连续胶结差段,不合格','no','measured','100413.PDF','high','违反红线:油气水层段无连续中等以上胶结'],
        ['hu102','7480','7502','22','胶结中等','中等',None,'过渡段','yes','measured','100413.PDF','medium','中等胶结'],
        ['hu102','7502','7540','38','胶结差','高',None,'连续胶结差段(>25m),不合格','no','measured','100413.PDF','high','违反红线:连续>25m胶结差'],
        ['hu102','7540','7665','125','胶结良好-中等','<15',None,'CBL低值段,胶结合格','yes','measured','100413.PDF','medium','良好-中等胶结'],
        ['hu102','6840','7665','825','综合评价','66.65%合格率',None,'整体不合格(低于70%红线)','no','measured','100413.PDF','high','CBL合格率66.65%<70%,违反红线规定; 目的井段固井质量不合格'],
    ]
    for r in rows:
        w.writerow(r)
print("8. cbl_evaluation.csv")

# ===================================================================
# 9. temperature_pressure_profile.csv
# ===================================================================
with open(os.path.join(out_dir, 'temperature_pressure_profile.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','md_m','temperature_c','pressure_mpa','ecd_gcm3','data_type','source_file','confidence','notes'])
    rows = [
        ['hu102','7120','147.8',None,None,'measured','20211.doc(s1.3)','high','电测温度@7120m'],
        ['hu102','7735','149',None,None,'measured','20211.doc(s1.3)','high','井底静止温度(BHST)'],
        ['hu102',None,'134',None,None,'derived','20211.doc(s1.3)','high','水泥浆试验温度(0.9x温度系数)'],
        ['hu102','5300',None,None,'2.04','calculated','20211.doc(s1.7)','high','注入2.30g/cm3重浆30m3后井底ECD'],
        ['hu102','7422',None,None,'2.044','calculated','20211.doc(s1.7)','high','7422m ECD(目的层顶K1q)'],
        ['hu102','7735',None,'8.60',None,'measured','20211.doc(s1.7)','high','承压试验(钻井液2.02,承压8.60MPa,10min降至7.00MPa)'],
        ['hu102','7735',None,None,'2.11','calculated','20211.doc(s1.7)','high','承压当量密度(7.00MPa折算)'],
        ['hu102',None,None,None,'0.066','calculated','20211.doc(s1.7)','high','井下安全窗口(g/cm3)=2.11-2.044'],
        ['hu102','6820',None,None,'2.053-2.093-2.066','calculated','20215.xlsx(Sheet2)','high','ECD计算范围(固井段顶,排量1.3m3/min)'],
        ['hu102','7120',None,None,'2.055-2.094-2.068','calculated','20215.xlsx(Sheet2)','high','ECD计算范围'],
        ['hu102','7300',None,None,'2.057-2.094-2.069','calculated','20215.xlsx(Sheet2)','high','ECD计算范围'],
        ['hu102','7422',None,None,'2.058-2.094-2.069','calculated','20215.xlsx(Sheet2)','high','ECD计算范围(目的层)'],
        ['hu102','7526',None,None,'2.058-2.094-2.07','calculated','20215.xlsx(Sheet2)','high','ECD计算范围'],
        ['hu102','7735',None,None,'2.059-2.095-2.071','calculated','20215.xlsx(Sheet2)','high','ECD计算范围(井底)'],
        ['hu102','3000',None,None,'2.054','calculated','20211.doc(s9.2)','high','提钻前3000m注入2.20重浆70m3后井底当量'],
    ]
    for r in rows:
        w.writerow(r)
print("9. temperature_pressure_profile.csv")

# ===================================================================
# 10. construction_events.csv
# ===================================================================
with open(os.path.join(out_dir, 'construction_events.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','date','time','event_type','description','depth_m','flow_rate_m3_min','pressure_mpa','weight_t','duration','data_type','source_file','confidence','notes'])
    events = [
        ['hu102','2022-11-17','01:30-02:00','下尾管','下部结构入井',None,None,None,None,'30min','measured','20211.doc','high','下尾管开始'],
        ['hu102','2022-11-17','07:40','下尾管完成','尾管下完(913m)',None,None,None,None,None,'measured','20211.doc','high','尾管总长913m'],
        ['hu102','2022-11-17','08:10','循环','循环','913','1.2','4.9',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-17','09:15','接悬挂器','接悬挂器,60t,空悬30t',None,None,None,'60',None,'measured','20211.doc','high',''],
        ['hu102','2022-11-17','19:00','下钻','下钻至2982m','2982',None,None,None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-17','19:35','灌浆','灌浆,装胶芯,接顶驱','2982',None,None,None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-17','20:11','顶通','顶通','2982','0.04','3.4',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-18','01:30','循环重浆','循环排重浆','5405','0.65','5.3',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-18','11:31','下钻','下钻至5405m,控压1.5-1.6MPa','5405',None,'2.3(静止)',None,None,'measured','20211.doc','high','控压下钻'],
        ['hu102','2022-11-18','15:40','顶通循环','顶通循环,出口温度40C','5405','0.54','5.6',None,None,'measured','20211.doc','high','出口温度40C'],
        ['hu102','2022-11-19','00:00','下钻','下钻至7096m灌浆,接顶驱','7096',None,'2.3',None,None,'measured','20211.doc','high','控压2.3MPa'],
        ['hu102','2022-11-19','03:40','下钻到位','下钻到位,245t','7735',None,None,'245',None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','03:50','探底','多下放1.4m探到底,继续下放1m,下压10t','7735',None,None,None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','04:02','顶通','顶通','7735','0.06','4.37-4.23',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','08:47','循环','循环','7735','0.6','11.6',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','17:35','循环','循环','7735','1','16.74',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','18:40','循环','循环(最大排量)','7735','1.09','19-18.5',None,None,'measured','20211.doc','high','排量达到最大'],
        ['hu102','2022-11-19','18:50','停泵投球','停泵,控压2MPa,卸钻杆,投球',None,None,'2',None,None,'measured','20211.doc','high','控压2MPa'],
        ['hu102','2022-11-19','20:32','送球到位','开泵送球到位','7735','0.18','4.3',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','20:35','憋压','憋压至11.3MPa,关井未圈压',None,None,'11.3',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','20:40','稳压','稳压11.3MPa,5min',None,None,'11.3',None,'5min','measured','20211.doc','high','稳压成功'],
        ['hu102','2022-11-19','20:43','坐挂','悬重由244t下放至238t后悬重不变',None,None,None,'244-238',None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','20:44','憋压','憋压至12MPa',None,None,'12',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','21:25','坐挂成功','由245t下放至207t(下压13t)回缩距3.04m',None,None,None,'245-207',None,'measured','20211.doc','high','坐挂成功:下压13t,回缩距3.04m'],
        ['hu102','2022-11-19','21:30','憋通球座','憋通球座','7735',None,'17.7',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-19','22:00','循环','循环','7735','1.1','17.5',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-20','00:00','循环','循环','7735','1.1','19.5',None,None,'measured','20211.doc','high','排量1.1,压力升至19.5MPa'],
        ['hu102','2022-11-20','09:30','循环','循环(低排量)','7735','0.55(30冲)','8.6-9.9',None,None,'measured','20211.doc','high','降低排量循环'],
        ['hu102','2022-11-20','12:07','计量演练','循环+罐面计量演练','7735','0.9(50冲)','16.1',None,None,'measured','20211.doc','high','正计量20方,反16方,差值4方:液面平稳不漏'],
        ['hu102','2022-11-20','15:41','循环压力上涨','循环压力上涨','7735','0.55(30冲)','12.1-10.18',None,None,'measured','20211.doc','high','出现明显压力上涨(憋堵征兆)'],
        ['hu102','2022-11-20','18:15','计量演练','循环+罐面计量演练','7735','1(55冲)','19.37-18.3',None,None,'measured','20211.doc','high','正计量42.5方,反36.5方,差值6方'],
        ['hu102','2022-11-20','18:45','停泵观察','停泵观察出口断流,回留5方','7735',None,None,None,'30min','measured','20211.doc','high',''],
        ['hu102','2022-11-20','18:50','倒扣丢手','悬重210t处正转40r,丢手成功','6826.33',None,None,'210-223',None,'measured','20211.doc','high','实探喇叭口6826.33m,下压14t;自由行程1.2m'],
        ['hu102','2022-11-20','19:24','发现渗漏','循环时发现渗漏(0.5/10min)','7735','1(55冲)','18.5',None,None,'measured','20211.doc','high','发现渗漏0.5方/10min'],
        ['hu102','2022-11-21','03:30','循环','循环(稳定)','7735','0.82(45冲)','16.4-16.17',None,None,'measured','20211.doc','high','稳定循环'],
        ['hu102','2022-11-21','11:07','循环','施工前最后一次循环','7735','0.9(50冲)','17.18',None,None,'measured','20211.doc','high','施工前循环'],
        ['hu102','2022-11-21','12:00','注前置液','注平衡液15m3',None,'0.3-0.6',None,None,'43min','measured','20211.doc','high','密度1.90g/cm3'],
        ['hu102','2022-11-21','14:00','异常-工具失效','管内外联通,下部结构失效',None,None,'4.2(套压)',None,None,'measured','20211.doc','high','泵入平衡液后放回水不断流,下部结构单流阀失效'],
        ['hu102','2022-11-21','17:35','方案调整','开会修改方案:尾管内全部替入2.10隔离液,不留上塞',None,None,None,None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','17:56','注隔离液','注隔离液15m3','7735','0.6-0.8','13-15.7',None,None,'measured','20211.doc','high','密度2.05g/cm3'],
        ['hu102','2022-11-21','18:10','注领浆','注领浆10m3','7735','0.6-0.8','15.7',None,'15min','measured','20211.doc','high','密度2.10g/cm3'],
        ['hu102','2022-11-21','18:25','注尾浆','注尾浆7m3','7735','0.6-0.8','15.5',None,None,'measured','20211.doc','high','密度2.10g/cm3'],
        ['hu102','2022-11-21','18:26','下胶塞','停泵下胶塞,控压1.9MPa',None,None,'1.9',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','18:41','注压塞液','注压塞液7m3','7735','0.3-0.8','6-10-13',None,'9min','measured','20211.doc','high','密度2.10g/cm3'],
        ['hu102','2022-11-21','18:50','注后置液','注保护液6m3','7735','0.6-0.8','13',None,'3min','measured','20211.doc','high','密度1.90g/cm3'],
        ['hu102','2022-11-21','18:57','替浆开始','替浆(钻井液)',None,'0.9','18.2',None,None,'measured','20211.doc','high','密度2.02g/cm3'],
        ['hu102','2022-11-21','19:30','平衡液进环空','替浆30m3,平衡液进环空',None,'0.9','17.2-17',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','19:44','隔离液进环空','替浆44m3,隔离液进环空',None,'0.9','15.7',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','20:02','领浆进环空','替浆59m3,领浆进环空',None,'0.73','13.5',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','20:13','胶塞重合','胶塞重合,剩余6方停泵',None,'0.64','12.96',None,None,'measured','20211.doc','high','泵冲累计3610,胶塞重合'],
        ['hu102','2022-11-21','20:18','尾浆进环空','替浆69m3,尾浆进环空',None,'0.64','13',None,None,'measured','20211.doc','high',''],
        ['hu102','2022-11-21','20:25','到量停泵','替浆72m3到量停泵',None,'0.46','9.9',None,None,'measured','20211.doc','high','泵冲累计3944,到量'],
        ['hu102','2022-11-21','20:45','关井','停泵关井套压4.2MPa,倒反循环管线',None,None,'4.2(套压)',None,'15min','measured','20211.doc','high',''],
        ['hu102','2022-11-21','21:00','拔中心管','控压4.2MPa上提7m拔出中心管',None,None,None,'220(悬重)',None,'measured','20211.doc','high',''],
        ['hu102','2022-11-22','00:30','循环排混浆','敞压循环排除混浆41方',None,'1.3','20',None,None,'measured','20211.doc','high','见隔离液、水泥浆混浆'],
        ['hu102','2022-11-22','11:00','控压循环','控压循环',None,'1.3','21.9',None,None,'measured','20211.doc','high','控压2.2MPa'],
        ['hu102','2022-11-24','21:00','关井候凝','加回压候凝5.8降至3.7MPa',None,None,None,None,None,'measured','20211.doc','high','候凝完成'],
        # Anomalies summary
        ['hu102','2022-11-17~20',None,'异常-生产组织','车辆核酸问题导致运灰/液车辆延后出车,大雪导致车辆误车',None,None,None,None,None,'measured','20211.doc','high','多辆车无连续3天核酸;5辆车信息有误;大雪地面泥泞'],
        ['hu102','2022-11-20',None,'异常-憋堵漏失','循环压力缓慢上涨,长时间循环导致憋堵加剧,发生漏失',None,None,None,None,None,'measured','20211.doc','high','产生吃入/回吐并存的呼吸效应'],
        ['hu102','2022-11-21','14:00','异常-工具失效','下部结构单流阀失效,放回水不断流',None,None,None,None,None,'measured','20211.doc','high','导致方案临场调整:不留上塞,带水泥头起钻循环排混浆'],
    ]
    for r in events:
        w.writerow(r)
print("10. construction_events.csv")

# ===================================================================
# 11. image_extraction_index.csv
# ===================================================================
with open(os.path.join(out_dir, 'image_extraction_index.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','file_path','file_name','file_type','category','description','extraction_status','notes'])
    images = [
        ['hu102','2/202/2022/20221/202211.jpg','202211.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段1','待OCR/视觉分析','深度段待确认'],
        ['hu102','2/202/2022/20221/202212.jpg','202212.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段2','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20221/202213.jpg','202213.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段3','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20221/202214.jpg','202214.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段4','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20221/202215.jpg','202215.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段5','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20221/202216.jpg','202216.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图-段6','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202221.jpg','202221.jpg','JPG(888x1857)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202222.jpg','202222.jpg','JPG(1920x910)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202223.jpg','202223.jpg','JPG(1520x720)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202224.jpg','202224.jpg','JPG(1920x864)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202225.jpg','202225.jpg','JPG(1280x1706)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202226.jpg','202226.jpg','JPG(1706x1280)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202227.jpg','202227.jpg','JPG(1706x1280)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/202228.jpg','202228.jpg','JPG(1706x1280)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','2/202/2022/20222/2022291.jpg','2022291.jpg','JPG(888x1857)','CBL/VDL测井图','CBL/VDL测井图','待OCR/视觉分析',''],
        ['hu102','1/1004/10041/100413.PDF','100413.PDF','PDF(1页)','固井质量评价图','呼102井CBL评价图(6840-7665m)','部分提取(文字层)','含完整井段CBL评价结论:合格率66.65%'],
        ['hu102','1/1004/10041/100422.PDF','100422.PDF','PDF(1页)','RCD固井质量成果图','RCD固井质量成果图-页1','待OCR(扫描图像)',''],
        ['hu102','1/1004/10041/100423.PDF','100423.PDF','PDF(1页)','RCD固井质量成果图','RCD固井质量成果图-页2','待OCR(扫描图像)',''],
        ['hu102','1/1004/10041/100424.PDF','100424.PDF','PDF(1页)','RCD固井质量成果图','RCD固井质量成果图-页3','待OCR(扫描图像)',''],
        ['hu102','1/1004/10041/100411.PDF','100411.PDF','PDF(1页,23MB)','固井质量评价图','CBL评价图(60-5449.6m)','部分提取(文字层)','193.7mm段,非尾管段'],
        ['hu102','1/1004/10041/100412.PDF','100412.PDF','PDF(1页)','固井质量评价图','CBL评价图(5449.6-6248m)','部分提取(文字层)','过渡段'],
        ['hu102','1/1004/10041/100419.PDF','100419.PDF','PDF(1页)','固井质量评价图','CBL评价图(5461-7085m)','部分提取(文字层)','顶部段'],
        ['hu102','2/202/2024/20242.jpg','20242.jpg','JPG','施工监测图','施工过程监测图','待视觉分析',''],
        ['hu102','2/202/2024/20243.jpg','20243.jpg','JPG','施工监测图','施工过程监测图','待视觉分析',''],
        ['hu102','2/202/2024/20244.jpg','20244.jpg','JPG','施工监测图','施工过程监测图','待视觉分析',''],
        ['hu102','2/202/2024/20245.png','20245.png','PNG','施工监测图','施工过程监测图','待视觉分析',''],
        ['hu102','2/202/2024/20246.png','20246.png','PNG','施工监测图','施工过程监测图','待视觉分析',''],
    ]
    for r in images:
        w.writerow(r)
print("11. image_extraction_index.csv")

# ===================================================================
# 12. source_trace.csv
# ===================================================================
with open(os.path.join(out_dir, 'source_trace.csv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['well_id','source_path','file_name','file_type','relevant_section','extraction_method','text_quality','cbl_relevant','notes'])
    sources = [
        ['hu102','2/202/2021/20211.doc','20211.doc','.doc(Word97-2003)','139.7mm尾管固井小结','Word COM Content.Text','good(中文字符正确)','yes','核心文档:井身结构/钻井液/井径+井斜/泵注/施工事件/扶正器/钻具'],
        ['hu102','2/202/2021/20211.wps','20211.wps','.wps','同20211.doc','未提取','未提取','yes','WPS格式,与.doc内容相同'],
        ['hu102','2/202/2021/20212.doc','20212.doc','.doc','139.7mm施工设计','Word COM Content.Text','good','yes','管柱/泵注程序/流体/地层/86根管柱明细'],
        ['hu102','2/202/2021/20214.doc','20214.doc','.doc','139.7mm施工总结','Word COM Content.Text','good','yes','施工方案/事件时间线/固井质量/探塞'],
        ['hu102','2/202/2021/20215.xlsx','20215.xlsx','.xlsx','套管单根数据(86根)','openpyxl','部分(中文sheet名乱码)','yes','86根套管单根数据明细; ECD计算; 管柱组合'],
        ['hu102','2/202/2021/20216.doc','20216.doc','.doc','施工记录表','Word COM Content.Text','good','yes','流体性能/施工参数记录'],
        ['hu102','2/202/2022/20221/*.jpg','202211-202216.jpg','.jpg(6张)','CBL/VDL测井图','待视觉分析(API不可用)','不可直接提取','yes','CBL测井图像(6张),1280x1706'],
        ['hu102','2/202/2022/20222/*.jpg','202221-2022291.jpg','.jpg(8张+1doc)','CBL/VDL测井图','待视觉分析','不可直接提取','yes','CBL测井图像(8张),多尺寸'],
        ['hu102','2/202/2022/20223.doc','20223.doc','.doc','193.7mm回接技术要求','antiword','good','no','非139.7mm尾管:193.7mm回接固井'],
        ['hu102','2/202/2022/20224.doc','20224.doc','.doc','193.7mm回接记录表','antiword','good','no','非139.7mm尾管'],
        ['hu102','2/202/2022/20225.pdf','20225.pdf','.pdf','193.7mm水泥浆试验','pdfplumber','不可直接提取','no','非139.7mm尾管'],
        ['hu102','2/202/2022/20229.xlsx','20229.xlsx','.xlsx','多管柱数据','openpyxl','部分','partial','含168.3mm+193.7mm+139.7mm,仅取139.7mm部分'],
        ['hu102','2/202/2023/20234.doc','20234.doc','.doc','水泥浆试验报告(领浆+尾浆)','Word COM Content.Text','good','yes','领浆/尾浆配方/稠化/流变/强度/失水/混浆污染试验'],
        ['hu102','2/202/2024/20241.xls','20241.xls','.xls','工具/监测数据','xlrd','部分','partial','193.7mm相关工具数据居多,139.7mm部分有限'],
        ['hu102','2/202/2024/20242-20246.*','.jpg/.png','施工监测图(5张)','待视觉分析','不可直接提取','yes','施工过程压力/排量监测图'],
        ['hu102','2/202/2025/20251-202593.*','各种格式','365.1mm技术套管','未提取(非尾管)','N/A','no','非139.7mm尾管:365.1mm技套资料'],
        ['hu102','2/202/2026/20261-202697.*','各种格式','473.1mm导管/表层','未提取(非尾管)','N/A','no','非139.7mm尾管:导管/表层资料'],
        ['hu102','1/1004/10041/100413.PDF','100413.PDF','.pdf(1页)','CBL评价图(6840-7665m)','pdfplumber字符提取','partial(文字层含关键结论)','yes','CBL合格率66.65%,不合格段7405-7480m+7502-7540m'],
        ['hu102','1/1004/10041/100411.PDF','100411.PDF','.pdf(1页)','CBL评价图(60-5449.6m)','pdfplumber','partial','no','193.7mm段CBL,非尾管段'],
        ['hu102','1/1004/10041/100412.PDF','100412.PDF','.pdf(1页)','CBL评价图(5449.6-6248m)','pdfplumber','partial','no','过渡段CBL,非尾管段'],
        ['hu102','1/1004/10041/100419.PDF','100419.PDF','.pdf(1页)','CBL评价图(5461-7085m)','pdfplumber','partial','partial','含上部尾管段CBL'],
        ['hu102','1/1004/10041/100422-100424.PDF','PDF(3页)','RCD固井质量成果图','pdfplumber(扫描,无文字)','不可直接提取','yes','扫描图像,无文字层,需OCR'],
    ]
    for r in sources:
        w.writerow(r)
print("12. source_trace.csv")

print("\nAll 12 CSV files written successfully!")
