import sys, csv
sys.stdout.reconfigure(encoding='utf-8')
base = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004"
rows = list(csv.DictReader(open(base+r"\caliper_profile.csv", encoding="utf-8-sig")))
print("caliper rows:", len(rows))
segs=[]
for i,r in enumerate(rows):
    top=float(r["md_m"]); bot=float(rows[i+1]["md_m"]) if i+1<len(rows) else 7660.0
    segs.append((top,bot,float(r["caliper_mm"]),r["notes"][:40]))
for a,b in [(5578,5630),(5630,5660),(5690,5720),(5990,6020),(6560,6600),(6980,7010),(7340,7378.05),(7378.05,7380),(7380,7400),(7460,7490),(7490,7520),(7520,7550),(7580,7610),(7610,7660)]:
    ss=[s for s in segs if s[0]>=a-0.01 and s[1]<=b+0.01]
    tot_mm=sum((s[1]-s[0])*s[2] for s in ss); tot_m=sum(s[1]-s[0] for s in ss)
    print(f"  {a}-{b}: {len(ss)} seg avg={tot_mm/tot_m if tot_m else float('nan'):.2f}  vals={sorted(set(s[2] for s in ss))}")
print("  305 rows:", [(s[0],s[1],s[2]) for s in segs if s[2]>=299])
print("  5578 vicinity:", [(s[0],s[1],s[2]) for s in segs if 5550<=s[0]<=5640])
print("  last 3:", segs[-3:])
print("  count>=300:", len([s for s in segs if s[2]>=299]))
# avg 5578-7521 and 7521-7660 per CSV
u=[s for s in segs if s[0]>=5578-0.01 and s[1]<=7521+0.01]
t1=sum((s[1]-s[0])*s[2] for s in u)/sum(s[1]-s[0] for s in u)
l=[s for s in segs if s[0]>=7521-0.01 and s[1]<=7660+0.01]
t2=sum((s[1]-s[0])*s[2] for s in l)/sum(s[1]-s[0] for s in l)
print(f"  CSV avg 5578-7521 = {t1:.2f} (n={len(u)}); 7521-7660 = {t2:.2f} (n={len(l)})")
print("  5578-5582 from table avg:", sum((min(s1,5582)-s0)*c for s0,s1,c in segs if s0<5582 and s1>5578)/sum(min(s1,5582)-max(s0,5578) for s0,s1,c in segs if s0<5582 and s1>5578))
inc = list(csv.DictReader(open(base+r"\inclination_profile.csv", encoding="utf-8-sig")))
print("incl rows:", len(inc))
pts=[(float(r["md_m"]), float(r["inclination_deg"])) for r in inc]
for md in (0,5600,5900,6410,7100,7340,7460,7520,7660):
    print(f"  md={md}:", [p for p in pts if abs(p[0]-md)<0.01])
print("  max:", max(pts,key=lambda p:p[1]))
