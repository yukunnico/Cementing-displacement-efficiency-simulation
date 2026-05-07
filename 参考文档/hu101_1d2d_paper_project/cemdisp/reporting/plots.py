
from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

def plot_all(result, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    m=result.metrics
    plt.figure(figsize=(8,5))
    plt.plot(m['time_min'],m['effective_efficiency'],label='全井段')
    plt.plot(m['time_min'],m['cbl_eval_interval_efficiency'],label='CBL评价段')
    plt.plot(m['time_min'],m['target_interval_efficiency'],label='目的层段')
    plt.xlabel('时间 / min'); plt.ylabel('顶替效率'); plt.title('呼101 1D-2D耦合顶替效率时程')
    plt.grid(alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(output_dir/'顶替效率时程.png',dpi=220); plt.close()

    p=result.depth_profiles.sort_values('井深_m')
    plt.figure(figsize=(7,8))
    plt.plot(p['周向平均水泥体积分数'],p['井深_m'],label='周向平均')
    plt.plot(p['宽边水泥体积分数'],p['井深_m'],label='宽边')
    plt.plot(p['窄边水泥体积分数'],p['井深_m'],label='窄边')
    plt.gca().invert_yaxis(); plt.xlabel('水泥体积分数'); plt.ylabel('井深 / m')
    plt.title('最终水泥体积分数深度剖面'); plt.grid(alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(output_dir/'最终水泥体积分数深度剖面.png',dpi=220); plt.close()

    plt.figure(figsize=(7,8))
    c=result.lead_field+result.tail_field
    extent=[result.geom['phi'][0],result.geom['phi'][-1],result.geom['md'][0],result.geom['md'][-1]]
    plt.imshow(c.T,aspect='auto',origin='lower',extent=extent)
    plt.colorbar(label='水泥体积分数')
    plt.xlabel('归一化方位：0=宽边，1=窄边'); plt.ylabel('井深 / m')
    plt.title('最终环空二维水泥浓度场')
    plt.tight_layout(); plt.savefig(output_dir/'最终环空二维水泥浓度场.png',dpi=220); plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(m['time_min'],m['front_wide_md_m'],label='宽边前沿')
    plt.plot(m['time_min'],m['front_narrow_md_m'],label='窄边前沿')
    plt.gca().invert_yaxis(); plt.xlabel('时间 / min'); plt.ylabel('前沿井深 / m')
    plt.title('宽边与窄边水泥前沿推进'); plt.grid(alpha=.25); plt.legend(); plt.tight_layout(); plt.savefig(output_dir/'宽窄边前沿推进.png',dpi=220); plt.close()

    seg=result.segment_efficiency
    plt.figure(figsize=(9,5))
    plt.bar(seg['井段_m'],seg['平均水泥体积分数'])
    plt.xticks(rotation=45,ha='right'); plt.ylabel('平均水泥体积分数'); plt.title('分段最终顶替效率')
    plt.tight_layout(); plt.savefig(output_dir/'分段最终顶替效率.png',dpi=220); plt.close()
