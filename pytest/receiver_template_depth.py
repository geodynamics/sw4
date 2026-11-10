#!/usr/bin/env python3
"""Check USEZVALUE against inline absolute/relative receivers over topography.

Run inside an authorized allocation. The optional unsupported-baseline control
records the old Gaussian-topography failure or missing absolute-depth behavior.
A flat-terrain inline control must still run successfully.
"""
import argparse
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np


INPUT = """grid h=200 x=4000 y=4000 z=3000 extrapolate=1
time steps=6
fileio path=.
supergrid gp=4
topography input=gaussian zmax=2000 order=4 gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500
block vp=6000 vs=3464 rho=2700
source x=1600 y=1600 z=200 mxy=1e18 t0=0.01 freq=16.6667 type=Gaussian
rec x=1800 y=1800 z=300 file=s_absolute usgsformat=0 sacformat=0 hdf5format=1 hdf5file=inline.h5 writeEvery=2 downSample=2 nsew=0
rec x=1800 y=1800 depth=300 file=s_relative usgsformat=0 sacformat=0 hdf5format=1 hdf5file=inline.h5 writeEvery=2 downSample=2 nsew=0
rechdf5 infile={template} outfile=receiver.h5 writeEvery=2 downSample=2
"""


def run(args, binary, root, template, input_text=INPUT, expect_topography_failure=False):
    root.mkdir()
    (root/'case.in').write_text(input_text.format(template=template))
    command=[args.launcher,'--overlap','--exact','--nodes=1','--ntasks=2',
             '--ntasks-per-node=2','--cpus-per-task=4']
    if args.job_id:
        command.append('--jobid='+args.job_id)
    if args.gpu:
        command+=['--gpus-per-task=1','--gpu-bind=closest']
    command += [str(binary),str(root/'case.in')]
    env=dict(os.environ,OMP_NUM_THREADS='1')
    result=subprocess.run(command,cwd=root,env=env,text=True,stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,timeout=120)
    (root/'sw4.log').write_text(result.stdout)
    if expect_topography_failure and result.returncode != 0:
        assert 'GaussianHill::inverse__grid_mapping: poor convergence' in result.stdout, result.stdout[-8000:]
        assert '-1.00000000e+38' in result.stdout and 'Time step ' not in result.stdout, result.stdout[-8000:]
        return None
    assert result.returncode==0 and 'HDF5-DIAG' not in result.stdout,result.stdout[-8000:]
    return root


def station_equal(left,right):
    assert set(left)==set(right)
    for field in left:
        a,b=left[field][...],right[field][...]
        if np.issubdtype(a.dtype,np.number):
            np.testing.assert_allclose(a,b,rtol=2e-12,atol=1e-15,equal_nan=True,
                                       err_msg=left.name+'/'+field)
        else:
            np.testing.assert_array_equal(a,b)


def check(root, unsupported=False):
    with h5py.File(root/'inline.h5','r') as inline,h5py.File(root/'receiver.h5','r') as template:
        station_equal(inline['s_relative'],template['s_relative'])
        assert any(not np.allclose(inline['s_absolute'][component][...],
                                   inline['s_relative'][component][...],
                                   rtol=2e-12,atol=1e-15) for component in 'XYZ'), \
            'Depth fixture does not distinguish absolute and relative station placement'
        if unsupported:
            differences=[]
            for component in 'XYZ':
                a=inline['s_absolute'][component][...]
                b=template['s_absolute'][component][...]
                differences.append(not np.allclose(a,b,rtol=2e-12,atol=1e-15))
            assert any(differences),'Baseline did not reproduce the missing absolute-depth behavior'
        else:
            station_equal(inline['s_absolute'],template['s_absolute'])
            peak=max(float(np.max(np.abs(template[name][component][...])))
                     for name in ('s_absolute','s_relative') for component in 'XYZ')
            assert peak>0


def elevated_absolute(args, binary):
    """Both representations must retain a receiver below terrain at negative z."""
    template = args.work_dir / 'elevated-stations.h5'
    with h5py.File(template, 'w') as stream:
        for name, z, use_z in [('s_absolute', -50, 1), ('s_relative', 50, 0)]:
            group = stream.create_group(name)
            group['STX,STY,STZ'] = np.asarray([2000, 2000, z], dtype=np.float64)
            group['ISNSEW'] = np.asarray([0], dtype=np.int32)
            group['USEZVALUE'] = np.asarray([use_z], dtype=np.int32)
    deck = INPUT.replace('x=1800 y=1800 z=300', 'x=2000 y=2000 z=-50')
    deck = deck.replace('x=1800 y=1800 depth=300', 'x=2000 y=2000 depth=50')
    root = run(args, binary, args.work_dir/'elevated', template, input_text=deck)
    with h5py.File(root/'inline.h5', 'r') as inline, h5py.File(root/'receiver.h5', 'r') as receivers:
        for name in ('s_absolute', 's_relative'):
            assert name in inline and name in receivers, 'Valid elevated receiver omitted: ' + name
            station_equal(inline[name], receivers[name])
        for component in 'XYZ':
            absolute = receivers['s_absolute'][component][...]
            relative = receivers['s_relative'][component][...]
            assert np.all(np.isfinite(absolute))
            np.testing.assert_allclose(absolute, relative, rtol=2e-12, atol=1e-15)
        assert any(np.max(np.abs(receivers['s_absolute'][component][...])) > 0 for component in 'XYZ')
    print('PASS: absolute z=-50 and relative depth=50 retain identical nonzero motion below terrain z=-100')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sw4',required=True,type=Path)
    parser.add_argument('--unsupported-baseline',type=Path)
    parser.add_argument('--work-dir',required=True,type=Path)
    parser.add_argument('--job-id')
    parser.add_argument('--gpu',action='store_true')
    parser.add_argument('--launcher',default='srun')
    args=parser.parse_args()
    args.work_dir=args.work_dir.resolve()
    args.work_dir.mkdir(parents=True,exist_ok=False)
    template=args.work_dir/'stations.h5'
    with h5py.File(template,'w') as stream:
        for name,use_z in [('s_absolute',1),('s_relative',0)]:
            group=stream.create_group(name)
            group['STX,STY,STZ']=np.asarray([1800,1800,300],dtype=np.float64)
            group['ISNSEW']=np.asarray([0],dtype=np.int32)
            group['USEZVALUE']=np.asarray([use_z],dtype=np.int32)
    if args.unsupported_baseline:
        baseline=args.unsupported_baseline.resolve(strict=True)
        inline_only="\n".join(line for line in INPUT.splitlines()
                               if 'file=s_relative ' not in line and not line.startswith(('rechdf5 ','topography ')))+"\n"
        control=run(args,baseline,args.work_dir/'baseline-flat-control',template,input_text=inline_only)
        with h5py.File(control/'inline.h5','r') as stream:
            assert any(np.max(np.abs(stream['s_absolute'][component][...]))>0 for component in 'XYZ')
        failed=run(args,baseline,args.work_dir/'baseline',template,expect_topography_failure=True)
        if failed is not None:
            check(failed,unsupported=True)
        print('PASS: baseline reproduces unsupported depth behavior; flat absolute inline control works')
    check(run(args,args.sw4.resolve(strict=True),args.work_dir/'candidate',template))
    print('PASS: template absolute/relative depths match inline receivers over Gaussian topography')
    elevated_absolute(args, args.sw4.resolve(strict=True))


if __name__=='__main__':
    main()
