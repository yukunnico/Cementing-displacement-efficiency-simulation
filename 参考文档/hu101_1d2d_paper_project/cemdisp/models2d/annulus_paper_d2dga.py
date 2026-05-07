
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, Tuple
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from cemdisp.data.fluid_spec import FluidSpec, FluidRole
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState

Array = NDArray[np.float64]

@dataclass(frozen=True)
class AnnulusSimulationResult:
    well_name: str
    geom: Dict[str, Array]
    lead_field: Array
    tail_field: Array
    spacer_field: Array
    metrics: pd.DataFrame
    depth_profiles: pd.DataFrame
    segment_efficiency: pd.DataFrame
    summary: Dict[str, float | str]
    snapshot_times_s: Tuple[float, ...]
    cement_snapshots: Tuple[Array, ...]

class AnnulusPaperD2DGASolver:
    """论文限定版D2DGA环空二维求解器。

    本实现只在环空段使用 Zhang & Frigaard (2022) 的D2DGA思想：
    - 偏心窄环空几何展开；
    - 二维方位-轴向浓度输运；
    - 高Peclet数间隙尺度弥散通量修正；
    - 黏度比、密度差/浮力数、偏心度对宽窄边推进的影响。

    不包含泥饼、温度、凝胶强度、CBL校准惩罚等项目中后加影响项。
    """
    def __init__(self, *, dt: float=6.0, nz: int=420, ny: int=48, total_t: float=8600.0, save_interval: int=80):
        self.dt=dt; self.nz=nz; self.ny=ny; self.total_t=total_t; self.save_interval=save_interval

    @staticmethod
    def d2dga_flux_multiplier(c: Array, m: float) -> Array:
        c=np.clip(c,0.0,1.0)
        denom=m*c**3 + (1.0-c**3)
        return np.clip((m*c**2 + 1.5*(1.0-c**2))/(denom + 1e-12), 0.35, 1.80)

    def _profile(self, points, md):
        d=np.array([p.depth_md_m for p in points],float)
        v=np.array([p.value for p in points],float)
        return np.interp(md,d,v)

    def _build_geom(self, well: WellSpec) -> Dict[str, Array]:
        s=np.linspace(0.0, well.bottom_md_m-well.top_md_m, self.nz)  # 0=鞋口，增大=向上
        md=well.bottom_md_m - s
        hole=self._profile(well.hole_diameter_profile, md)
        inc=self._profile(well.inclination_profile, md)
        standoff=self._profile(well.standoff_profile, md)
        od=np.where(md <= well.upper_lower_transition_md_m, well.upper_liner_od_mm, well.lower_liner_od_mm)
        # 上下段过渡按测深判别：md越小越靠上，≤ transition为上部168.3mm尾管
        od=np.where(md <= well.upper_lower_transition_md_m, well.upper_liner_od_mm, well.lower_liner_od_mm)
        radius=((hole+od)/4.0)/1000.0
        phi=np.linspace(0.0,1.0,self.ny)  # 0=宽边，1=窄边，按半周对称展开
        y=phi*np.pi*np.mean(radius)
        e=np.clip(1.0-standoff, 0.05, 0.72)
        clearance=np.maximum((hole-od)/1000.0, 0.01)
        b=np.zeros((self.ny,self.nz),float)
        for j in range(self.nz):
            b[:,j]=clearance[j]*(1.0 + e[j]*np.cos(np.pi*phi))
        # 体积校正：半周展开积分×2，应与分段环空体积一致
        area=np.pi*((hole/1000.0)**2-(od/1000.0)**2)/4.0
        target_vol=float(abs(np.trapezoid(area, x=md)))
        vol=float(2.0*np.trapezoid(np.trapezoid(b,x=s,axis=1),x=y,axis=0))
        scale=target_vol/max(vol,1e-12)
        b*=scale
        return {"s":s,"md":md,"phi":phi,"y":y,"b":b,"e":e,"hole_mm":hole,"od_mm":od,"inc_deg":inc,"volume_m3":np.array(target_vol)}

    def _velocity(self, lead, tail, spacer, geom, q_m3s, fluids_by_role):
        b=geom['b']; phi=geom['phi'][:,None]
        mud=fluids_by_role['mud']; lead_f=fluids_by_role['lead']; tail_f=fluids_by_role['tail']; spacer_f=fluids_by_role.get('spacer')
        cement=np.clip(lead+tail,0,1)
        mud_frac=np.clip(1.0-lead-tail-spacer,0,1)
        # 表观黏度和密度，作为D2DGA闭合中的黏度比/浮力项输入
        mu_mud=mud.apparent_viscosity(100.0)
        mu_lead=lead_f.apparent_viscosity(100.0)
        mu_tail=tail_f.apparent_viscosity(100.0)
        mu_sp=spacer_f.apparent_viscosity(100.0) if spacer_f is not None else mu_mud
        mu=mud_frac*mu_mud + lead*mu_lead + tail*mu_tail + spacer*mu_sp
        rho=mud_frac*mud.density_kg_m3 + lead*lead_f.density_kg_m3 + tail*tail_f.density_kg_m3 + spacer*(spacer_f.density_kg_m3 if spacer_f else mud.density_kg_m3)
        rho_disp=(lead_f.density_kg_m3*0.67 + tail_f.density_kg_m3*0.33)
        # 偏心通道主导 + 密度稳定时向窄边补偿；密度不稳定时加剧宽边窜流
        base=(b/np.mean(b,axis=0,keepdims=True))**2/np.maximum(mu,1e-5)
        density_contrast=(rho_disp-mud.density_kg_m3)/mud.density_kg_m3
        stable=np.clip(8.0*density_contrast, -0.35, 0.45)
        ebar=geom['e'][None,:]
        buoyancy_shape=1.0 + stable*ebar*(2.0*phi-1.0)
        pref=np.maximum(base*buoyancy_shape,1e-8)
        dy=np.gradient(geom['y'])[:,None]
        area_weight=np.sum(pref*b*dy*2.0,axis=0,keepdims=True)
        w=q_m3s*pref/np.maximum(area_weight,1e-12)
        return w, mu, rho

    def _advect(self, field, inlet_value, w_eff, geom):
        z=geom['s']; dz=z[1]-z[0]
        c=field.copy()
        # upwind conservative update: positive w flows from shoe z=0 to hanger
        courant=np.clip(w_eff*self.dt/dz,0.0,0.85)
        out=c.copy()
        out[:,0]=(1.0-courant[:,0])*c[:,0] + courant[:,0]*inlet_value
        out[:,1:]=(1.0-courant[:,1:])*c[:,1:] + courant[:,1:]*c[:,:-1]
        return np.clip(out,0.0,1.0)

    def _smooth_dispersion(self, field, axial=0.020, azimuthal=0.018):
        # 以显式小系数模拟D2DGA数值弥散后的轴向/方位弥散，不引入项目工程项
        f=field.copy()
        f[:,1:-1]+=axial*(field[:,2:]-2*field[:,1:-1]+field[:,:-2])
        f[1:-1,:]+=azimuthal*(field[2:,:]-2*field[1:-1,:]+field[:-2,:])
        f[0,:]+=azimuthal*(field[1,:]-field[0,:])
        f[-1,:]+=azimuthal*(field[-2,:]-field[-1,:])
        return np.clip(f,0.0,1.0)

    def _weighted_mean(self, arr, geom, mask=None):
        weight=geom['b']
        if mask is not None:
            weight=weight*mask[None,:]
        return float(np.sum(arr*weight)/max(np.sum(weight),1e-12))

    def _field_volume(self, arr, geom):
        return float(2.0*np.trapezoid(np.trapezoid(arr*geom['b'], x=geom['s'], axis=1), x=geom['y'], axis=0))

    def _mass_limit(self, arr, geom, target_volume):
        current = self._field_volume(arr, geom)
        if current > max(target_volume, 0.0) + 1e-9 and current > 1e-12:
            arr = arr * (max(target_volume, 0.0) / current)
        return np.clip(arr, 0.0, 1.0)

    def _front(self, c, geom, row):
        prof=c[row,:]
        hit=np.where(prof>0.05)[0]
        if len(hit)==0:
            return float('nan')
        # s越大越靠上，md越小；返回最靠上的前沿井深
        return float(geom['md'][hit[-1]])

    def run(self, well: WellSpec, fluids: Tuple[FluidSpec, ...], inlet_provider: Callable[[float], AnnulusInletState]) -> AnnulusSimulationResult:
        geom=self._build_geom(well)
        fluids_by_role={}
        fluids_by_role['mud']=next(f for f in fluids if f.role==FluidRole.MUD)
        fluids_by_role['lead']=next(f for f in fluids if f.role==FluidRole.LEAD)
        fluids_by_role['tail']=next(f for f in fluids if f.role==FluidRole.TAIL)
        fluids_by_role['spacer']=next((f for f in fluids if f.role in {FluidRole.WASH,FluidRole.SPACER}), None)
        lead=np.zeros((self.ny,self.nz),float)
        tail=np.zeros_like(lead); spacer=np.zeros_like(lead)
        metrics=[]; snaps=[]; snap_times=[]
        cum_lead_in=0.0; cum_tail_in=0.0; cum_spacer_in=0.0
        mud=fluids_by_role['mud']
        m_lead=mud.apparent_viscosity(100.0)/fluids_by_role['lead'].apparent_viscosity(100.0)
        m_tail=mud.apparent_viscosity(100.0)/fluids_by_role['tail'].apparent_viscosity(100.0)
        for it,t in enumerate(np.arange(0,self.total_t+self.dt,self.dt)):
            inlet=inlet_provider(float(t))
            q=inlet.flow_rate_m3_s
            if q>1e-10:
                w,mu,rho=self._velocity(lead,tail,spacer,geom,q,fluids_by_role)
                cement=np.clip(lead+tail,0,1)
                lead_mult=self.d2dga_flux_multiplier(np.maximum(lead,1e-6),m_lead)
                tail_mult=self.d2dga_flux_multiplier(np.maximum(tail,1e-6),m_tail)
                sp_mult=self.d2dga_flux_multiplier(np.maximum(spacer,1e-6),1.0)
                lead=self._advect(lead,inlet.lead_fraction,w*lead_mult,geom)
                tail=self._advect(tail,inlet.tail_fraction,w*tail_mult,geom)
                spacer=self._advect(spacer,inlet.spacer_fraction,w*sp_mult,geom)
                total=np.maximum(lead+tail+spacer,1.0)
                over=(lead+tail+spacer)>1.0
                lead[over]/=total[over]; tail[over]/=total[over]; spacer[over]/=total[over]
                # D2DGA间隙尺度弥散：在低浓度前锋更强
                lead=self._smooth_dispersion(lead, axial=0.018, azimuthal=0.015)
                tail=self._smooth_dispersion(tail, axial=0.018, azimuthal=0.015)
                spacer=self._smooth_dispersion(spacer, axial=0.012, azimuthal=0.012)
                # 质量守恒限制：D2DGA通量修正只改变前缘形态，不能凭空生成相体积。
                # 对各相按套管1D出口累积进入体积做上限约束；若前缘越过顶部，场内体积可低于上限。
                cum_lead_in += q * inlet.lead_fraction * self.dt
                cum_tail_in += q * inlet.tail_fraction * self.dt
                cum_spacer_in += q * inlet.spacer_fraction * self.dt
                lead = self._mass_limit(lead, geom, cum_lead_in)
                tail = self._mass_limit(tail, geom, cum_tail_in)
                spacer = self._mass_limit(spacer, geom, cum_spacer_in)
            cement=np.clip(lead+tail,0,1)
            md=geom['md']
            masks={w.window_type: ((md>=w.top_md_m)&(md<=w.bottom_md_m)) for w in well.evaluation_windows}
            wide_row=0; narrow_row=-1
            fw=self._front(cement,geom,wide_row); fn=self._front(cement,geom,narrow_row)
            if np.isnan(fw) or np.isnan(fn):
                channel=0.0
            else:
                # 宽边前沿井深更小表示走得更远；归一化宽窄边超前距离
                channel=abs(float(fn-fw))/(well.bottom_md_m-well.top_md_m)
            mix=float(np.mean((cement>0.05)&(cement<0.95)))
            metrics.append({
                'time_s':float(t),'time_min':float(t/60.0),'inlet_fluid':inlet.fluid_name,
                'flow_rate_m3_min':float(q*60.0),'bulk_cement_fill':self._weighted_mean(cement,geom),
                'lead_fill':self._weighted_mean(lead,geom),'tail_fill':self._weighted_mean(tail,geom),
                'spacer_fill':self._weighted_mean(spacer,geom),
                'effective_efficiency':self._weighted_mean(cement,geom,masks['full']),
                'cbl_eval_interval_efficiency':self._weighted_mean(cement,geom,masks['cbl']),
                'target_interval_efficiency':self._weighted_mean(cement,geom,masks['target']),
                'wide_side_mean_cement':float(np.mean(cement[wide_row,:])),
                'narrow_side_mean_cement':float(np.mean(cement[narrow_row,:])),
                'front_wide_md_m':self._front(cement,geom,wide_row),
                'front_narrow_md_m':self._front(cement,geom,narrow_row),
                'channeling_index':float(channel),
                'mixing_index':mix,
            })
            if it % self.save_interval == 0 or t>=self.total_t:
                snaps.append(cement.copy()); snap_times.append(float(t))
        metrics_df=pd.DataFrame(metrics)
        cement=np.clip(lead+tail,0,1)
        depth_profiles=pd.DataFrame({
            '井深_m':geom['md'],
            '周向平均水泥体积分数':[self._weighted_mean(cement,geom,np.eye(self.nz,dtype=bool)[j]) for j in range(self.nz)],
            '宽边水泥体积分数':cement[0,:],
            '窄边水泥体积分数':cement[-1,:],
            '平均间隙_m':np.mean(geom['b'],axis=0),
            '偏心度':geom['e'],
            '井径_mm':geom['hole_mm'],
            '尾管外径_mm':geom['od_mm'],
        })
        seg_defs=[('5402.85-5700',5402.85,5700.0),('5700-6153',5700.0,6153.0),('6153-6800',6153.0,6800.0),('6800-7200',6800.0,7200.0),('7200-7492',7200.0,7492.0),('7492-7600',7492.0,7600.0),('7600-7735',7600.0,7735.0),('7735-7810',7735.0,7810.0),('7810-7868',7810.0,7868.0)]
        rows=[]
        for name,top,bottom in seg_defs:
            mask=(geom['md']>=top)&(geom['md']<=bottom)
            rows.append({'井段_m':name,'顶界_m':top,'底界_m':bottom,
                         '平均水泥体积分数':self._weighted_mean(cement,geom,mask),
                         '宽边平均':float(np.mean(cement[0,mask])) if np.any(mask) else np.nan,
                         '窄边平均':float(np.mean(cement[-1,mask])) if np.any(mask) else np.nan,
                         '宽窄边差值':float(np.mean(cement[0,mask])-np.mean(cement[-1,mask])) if np.any(mask) else np.nan})
        seg_df=pd.DataFrame(rows)
        final=metrics_df.iloc[-1]
        summary={
            '模型名称':'呼101_项目结构_套管1D-环空论文D2DGA耦合模型',
            '全井段最终顶替效率':float(final['effective_efficiency']),
            'CBL评价段最终顶替效率':float(final['cbl_eval_interval_efficiency']),
            '目的层段最终顶替效率':float(final['target_interval_efficiency']),
            '宽边平均水泥体积分数':float(final['wide_side_mean_cement']),
            '窄边平均水泥体积分数':float(final['narrow_side_mean_cement']),
            '窜槽指数':float(final['channeling_index']),
            '混浆指数':float(final['mixing_index']),
            '环空物理体积_m3':float(geom['volume_m3']),
            '模拟结束时间_min':float(final['time_min']),
            '备注':'环空段仅采用论文D2DGA因素；套管段为1D活塞流前沿追踪并提供鞋口边界。',
        }
        return AnnulusSimulationResult(well.well_name,geom,lead,tail,spacer,metrics_df,depth_profiles,seg_df,summary,tuple(snap_times),tuple(snaps))
