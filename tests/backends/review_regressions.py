#!/usr/bin/env python3
"""On-demand propagation checks for production-review source/material fixes.

Run CPU fixtures alone with --cpu, or compare both backends with --cpu and --gpu.
No tests are registered with pytest/CTest. Existing suites remain unchanged.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

from compare_waveforms import run, compare, GRID as COARSE_GRID, RECEIVERS

# Avoid the growing amplitudes observed in the coarser trial fixtures.
GRID = COARSE_GRID.replace('h=200', 'h=100')


def sac(path, values, dt):
    floats=np.full(70,-12345,dtype='<f4'); ints=np.full(40,-12345,dtype='<i4')
    floats[0]=dt;floats[5]=0;floats[6]=(len(values)-1)*dt
    ints[:7]=[2026,275,1,2,3,234,6];ints[9]=len(values);ints[15]=1;ints[35]=1
    path.write_bytes(floats.tobytes()+ints.tobytes()+b' '*192+np.asarray(values,dtype='<f4').tobytes())


def cases(root):
    common='time t=0.8\nsupergrid gp=8\nfileio path=output\ndeveloper cfl=0.8\n'
    block='block vp=6000 vs=3464 rho=2700\n'
    position='source x=1637 y=1643 z=1000 '
    histories=[]
    for component in range(6):
        phase=np.linspace(0,np.pi,128)
        histories.append((component+1)*np.sin(phase)**2+.1*component*np.sin(2*phase))
    dt=.0025
    (root/'force.dat').write_text(f'0 {dt} 128\n'+'\n'.join(map(str,histories[0]))+'\n')
    for component,values in zip(('x','y','z'),histories):sac(root/f'forces.{component}',values*1e10,dt)
    for component,values in zip(('xx','xy','xz','yy','yz','zz'),histories):sac(root/f'moments.{component}',values*1e14,dt)
    inputs={}
    for filtered in (False,True):
        filter_line='prefilter type=lowpass order=4 passes=2 fc2=12\n' if filtered else ''
        for name,source in (
            ('discrete-force',f'fx=1e10 fy=2e10 fz=-5e9 dfile={root}/force.dat'),
            ('sac-forces',f'sacbasedisp={root}/forces'),
            ('sac-moments',f'sacbase={root}/moments')):
            inputs[f'{name}-filter{int(filtered)}']=GRID+common+filter_line+block+position+source+'\n'+RECEIVERS
    # Both sides of the Cartesian interface, including fractional locations.
    refinement=GRID.replace('h=200','h=100')+common.replace('gp=4','gp=8')
    refinement+='refinement zmax=2000\ndeveloper ctol=1e-10 cmaxit=400 crelax=0.2\n'+block
    coarse_h=float(re.search(r'h=([0-9.]+)',GRID)[1])
    fine_h=coarse_h/2
    interface=2000.
    fine_nz=round(interface/fine_h)+1
    depths=[(kc-1+.46)*fine_h for kc in (fine_nz-3,fine_nz-2,fine_nz-1,fine_nz-6)]
    depths += [interface+(kc-1+.47)*coarse_h for kc in (1,2,3,6)]
    depths += [interface+epsilon for epsilon in (-1e-7,0.,1e-7)]
    for depth in depths:
        for kind,source in (('force','fx=1e10 fy=2e10 fz=-5e9'),('moment','mxy=1e15 mxz=2e15 myy=5e14')):
            name=f'interface-{kind}-z{depth:.17g}'
            if name in inputs:
                raise ValueError(f'Duplicate interface source depth: {name}')
            inputs[name]=refinement+f'source x=1637 y=1643 z={depth} '+source+' type=Gaussian freq=20 t0=0.05\n'+RECEIVERS
    inputs['low-vp-vs']=GRID+common+'block vp=3540 vs=3000 rho=2700\n'+position+'mxy=1e15 type=Gaussian freq=20 t0=0.05\n'+RECEIVERS
    # Isotropic stiffness represented through the anisotropic operator.
    rho=2700; mu=rho*3464**2;lam=rho*6000**2-2*mu
    diagonal={'c11','c22','c33'};shear={'c44','c55','c66'};offdiag={'c12','c13','c23'}
    coefficients=[f'c{i}{j}={lam+2*mu if f"c{i}{j}" in diagonal else mu if f"c{i}{j}" in shear else lam if f"c{i}{j}" in offdiag else 0}'
                  for i in range(1,7) for j in range(i,7)]
    inputs['anisotropic-isotropic-limit']=GRID+common+'anisotropy\nablock rho=2700 '+' '.join(coefficients)+'\n'+position+'mxy=1e15 type=Gaussian freq=20 t0=0.05\n'+RECEIVERS
    inputs['isotropic-limit-reference']=GRID+common+block+position+'mxy=1e15 type=Gaussian freq=20 t0=0.05\n'+RECEIVERS
    inputs['event-default-paths']=GRID+common+'event name=first\n'+block+position+'mxy=1e15 type=Gaussian freq=20 t0=0.05 event=first\n'+RECEIVERS
    topo='topography input=gaussian zmax=2000 order=4 gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500\n'
    for width in range(1,6):
        for attenuation in (False,True):
            material='block vp=6000 vs=3464 rho=2700 vpgrad=0.1 vsgrad=0.05 rhograd=0.02'
            material+=' qp=200 qs=100\n' if attenuation else '\n'
            inputs[f'curvi-extrapolate{width}-att{int(attenuation)}']=GRID.replace('h=200','h=100').rstrip()+f' extrapolate={width}\n'+common.replace('gp=4','gp=8')
            inputs[f'curvi-extrapolate{width}-att{int(attenuation)}']+=topo+'refinement zmax=800\n'+('attenuation nmech=1\n' if attenuation else '')+material+position+'mxy=1e15 type=Gaussian freq=20 t0=0.05\n'+RECEIVERS
    return inputs



def interface_selection(name, log):
    depth=float(name.split('-z')[1])
    grids={int(g):(float(h),int(nz)) for g,h,nz in
           re.findall(r'^\s*(\d+)\s+(\S+)\s+\d+\s+\d+\s+(\d+)\s+\d+\s+Cartesian\s*$',log.read_text(),re.MULTILINE)}
    if set(grids)!={0,1}:raise ValueError('Unexpected actual refinement grids')
    fine_h,fine_nz=grids[1];coarse_h,coarse_nz=grids[0]
    interface=(fine_nz-1)*fine_h
    if interface!=2000 or coarse_h!=2*fine_h:raise ValueError('Actual grids differ from source sweep')
    g=0 if depth>=interface else 1
    h,nz=grids[g];zmin=interface if g==0 else 0.
    kc=max(1,min(int(np.floor((depth-zmin)/h+1)),nz-1))
    return {'grid':g,'kc':kc,'Nz':nz,'h':h,'depth':depth}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu',type=Path,required=True)
    parser.add_argument('--gpu',type=Path)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--tasks',type=int,default=2)
    parser.add_argument('--case',action='append',help='Run only matching case prefixes')
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    report={}
    coverage={}
    for name,text in cases(root).items():
        if args.case and not any(name.startswith(prefix) for prefix in args.case):continue
        # Interface receivers observe meaningful arrivals on both grids.
        coordinates=[('station',(1800,1800,1600)),('edge',(1800,1800,2800))]
        native=run(args.cpu.resolve(strict=True),root/name/'cpu',text,args.tasks,False,coordinates)
        if name=='event-default-paths':native=root/name/'cpu/receivers.h5'
        if name.startswith('interface-'):
            coverage[name]=interface_selection(name,root/name/'cpu/run.log')
        if args.gpu:
            gpu=run(args.gpu.resolve(strict=True),root/name/'gpu',text,args.tasks,True,coordinates)
            if name=='event-default-paths':gpu=root/name/'gpu/receivers.h5'
            if name.startswith('interface-') and interface_selection(name,root/name/'gpu/run.log')!=coverage[name]:
                raise ValueError('Backend source grid/stencil selections differ')
            report[name]=compare(native,gpu,2e-5,1e-12,arrival_after=.15)
        else:
            report[name]=compare(native,native,2e-5,1e-12,arrival_after=.15)
        # These fixed sources should produce metre-scale displacements.
        # A generous ceiling prevents matching unstable traces from passing.
        if max(item['reference_peak'] for item in report[name].values()) > 1e3:
            raise ValueError(f'Unphysical displacement in fixed regression fixture: {name}')
        (root/'interface-coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
        (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PASS:',name,flush=True)
    # Check the independent isotropic limit for each backend that was executed.
    if {'isotropic-limit-reference','anisotropic-isotropic-limit'} <= set(report):
        for backend in ('cpu','gpu') if args.gpu else ('cpu',):
            compare(root/'isotropic-limit-reference'/backend/'output/receivers.h5',
                    root/'anisotropic-isotropic-limit'/backend/'output/receivers.h5',2e-5,1e-12,arrival_after=.15)
            print('PASS: anisotropic isotropic limit on',backend,flush=True)


if __name__=='__main__':main()
