#!/usr/bin/env python3
"""On-demand SAC/rechdf5/USGS waveform and metadata equivalence checks.

Run inside an allocation. No tests are added to pytest or CTest. SAC and
receiver HDF5 store float32 samples; USGS text retains the solver precision.
Downsampled HDF5 is compared to the corresponding full-rate SAC/text samples.
"""
import argparse
from datetime import datetime, timedelta
import json
import os
from pathlib import Path
import re
import subprocess

import h5py
import numpy as np
from validation import check_solver_log


def read_sac(path):
    data = path.read_bytes()
    for endian in ('<', '>'):
        floats = np.frombuffer(data, dtype=endian+'f4', count=70)
        ints = np.frombuffer(data, dtype=endian+'i4', count=40, offset=280)
        if ints[6] == 6 and 0 < ints[9] < 10000000 and len(data) == 632 + 4*ints[9]:
            component = data[440+20*8:440+21*8].decode('ascii').strip()
            values = np.frombuffer(data, dtype=endian+'f4', count=ints[9], offset=632).astype(float)
            utc = datetime(int(ints[0]), 1, 1) + timedelta(days=int(ints[1])-1,
                    hours=int(ints[2]), minutes=int(ints[3]), seconds=int(ints[4]),
                    milliseconds=int(ints[5]))
            return component, values, floats, utc
    raise ValueError(f'Invalid SAC binary header/payload: {path}')


def check(directory, mode, orientation, downsample):
    columns = {
        'displacement': ('EW', 'NS', 'UP') if orientation else ('X', 'Y', 'Z'),
        'velocity': ('Vew', 'Vns', 'Vup') if orientation else ('Vx', 'Vy', 'Vz'),
        'div': ('Div',), 'curl': ('Curlx', 'Curly', 'Curlz'),
        'strains': ('Uxx', 'Uyy', 'Uzz', 'Uxy', 'Uxz', 'Uyz'),
        'displacementgradient': ('DUXDX','DUXDY','DUXDZ','DUYDX','DUYDY','DUYDZ','DUZDX','DUZDY','DUZDZ'),
    }[mode]
    output = directory/'output'
    text = np.loadtxt(output/'ascii.txt')
    if text.ndim != 2 or text.shape[1] != 1+len(columns) or not np.all(np.isfinite(text)):
        raise ValueError('Invalid USGS text shape/data')
    sac = {}
    for path in output.glob('ascii.*'):
        if path.suffix != '.txt':
            component, values, header, utc = read_sac(path)
            sac[component] = values, header, utc
    errors = {}
    with h5py.File(output/'receivers.h5') as stream:
        group = stream['station']
        dt = float(np.asarray(stream['DELTA']).reshape(-1)[0])
        datetime_value = stream.attrs['DATETIME']
        if mode in ('displacement', 'velocity'):
            unit = stream.attrs['UNIT']
            unit = unit.decode() if isinstance(unit, bytes) else str(unit)
            expected_unit = 'm/s' if mode == 'velocity' else 'm'
            if unit != expected_unit or f'({expected_unit})' not in (output/'ascii.txt').read_text():
                raise ValueError(f'Output units differ: {unit}, expected {expected_unit}')
        hdf_utc = datetime.fromisoformat(datetime_value.decode() if isinstance(datetime_value, bytes) else str(datetime_value))
        date = re.search(r'^# Date: UTC\s+(\S+)', (output/'ascii.txt').read_text(), re.MULTILINE)
        text_utc = datetime.strptime(date[1], '%m/%d/%Y:%H:%M:%S.%f')
        if text_utc != hdf_utc:
            raise ValueError(f'USGS/HDF5 UTC differs: {text_utc}, {hdf_utc}')
        npts = int(np.asarray(group['NPTS']).reshape(-1)[0])
        np.testing.assert_allclose(np.asarray(stream['ORIGINTIME']).reshape(-1)[0], .05,
                                   atol=1e-8, rtol=2e-7)
        expected_indices = np.arange(0, text.shape[0], downsample)
        if npts != len(expected_indices):
            raise ValueError(f'HDF5 NPTS does not describe decimated SAC/text: {npts}')
        for index, component in enumerate(columns):
            values, header, sac_utc = sac[component]
            hdf = np.asarray(group[component][:npts], dtype=float)
            if len(values) != len(text) or not np.all(np.isfinite(values)) or not np.all(np.isfinite(hdf)):
                raise ValueError('SAC/text lengths or finite values differ')
            peak = float(np.max(np.abs(text[:,index+1])))
            floor = 1e-12 + 2e-6*peak
            delta_text = float(np.max(np.abs(values-text[:,index+1])))
            delta_hdf = float(np.max(np.abs(hdf-values[expected_indices])))
            if max(delta_text, delta_hdf) > floor:
                raise ValueError(f'{mode}/{component} differs: text={delta_text}, HDF5={delta_hdf}, allowed={floor}')
            # SAC DELTA/B and HDF5 DELTA have float32 precision; USGS time is %e.
            relative_time = float(header[5]) + np.arange(len(values))*float(header[0])
            expected_time = (sac_utc-hdf_utc).total_seconds() + relative_time
            np.testing.assert_allclose(text[:,0], expected_time, rtol=8e-7, atol=1e-8)
            np.testing.assert_allclose(dt, float(header[0])*downsample, rtol=2e-7, atol=1e-10)
            np.testing.assert_allclose(header[6], relative_time[-1], rtol=2e-7, atol=1e-8)
            np.testing.assert_allclose((sac_utc-hdf_utc).total_seconds() + float(header[7]), .05,
                                       rtol=2e-7, atol=1e-8)
            if abs((sac_utc-hdf_utc).total_seconds()) >= .001:
                raise ValueError(f'SAC/HDF5 UTC differs: {sac_utc}, {hdf_utc}')
            for label, value in [('CMPAZ',header[57]),('CMPINC',header[58])]:
                np.testing.assert_allclose(np.asarray(group[component+label]).reshape(-1)[0], value, atol=1e-5, rtol=0)
            errors[component] = {'text_sac_difference': delta_text, 'hdf_sac_difference':delta_hdf, 'reference_peak':peak}
        if max(item['reference_peak'] for item in errors.values()) <= 1e-10:
            raise ValueError('No meaningful receiver signal')
        if int(np.asarray(group['ISNSEW']).reshape(-1)[0]) != orientation:
            raise ValueError('Receiver orientation metadata differs')
    return errors


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sw4',type=Path,required=True)
    p.add_argument('--backend',choices=('OPENMP','CUDA','HIP'),default='OPENMP')
    p.add_argument('--work-dir',type=Path,required=True)
    p.add_argument('--tasks',type=int,default=2)
    p.add_argument('--reader',type=Path,help='Optional native sw4_receiver_check executable')
    args=p.parse_args()
    exe=args.sw4.resolve(strict=True); root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    report={}
    cases=[(mode,orientation,downsample,topo) for mode in ('displacement','velocity')
           for orientation in (0,1) for downsample in (1,3) for topo in (False,True)]
    cases += [(mode,0,1,False) for mode in ('div','curl','strains','displacementgradient')]
    cases = [(*case, '234567') for case in cases]
    cases += [('displacement',0,1,False,fraction) for fraction in ('000123','999999')]
    for mode,orientation,downsample,topo,fraction in cases:
        name=f'{mode}-nsew{orientation}-ds{downsample}-topo{int(topo)}-utc{fraction}'
        directory=root/name;directory.mkdir(exist_ok=True)
        with h5py.File(directory/'stations.h5','w') as stream:
            group=stream.create_group('station');group['STX,STY,STZ']=[1800.,1800.,200.]
            group['ISNSEW']=np.array([orientation],dtype=np.int32);group['USEZVALUE']=np.array([1],dtype=np.int32)
        text='grid h=200 x=4000 y=4000 z=4000 lat=37 lon=-122 az=27\n'
        text+=f'time t=0.5 utcstart=10/02/2026:01:02:03.{fraction}\nsupergrid gp=4\nfileio path=output\ndeveloper cfl=0.8\n'
        if topo:
            text+='topography input=gaussian zmax=2000 order=4 gaussianAmp=100 gaussianXc=2000 gaussianYc=2000 gaussianLx=1500 gaussianLy=1500\n'
        text+='block vp=6000 vs=3464 rho=2700\nsource x=1600 y=1600 z=600 mxy=1e15 t0=0.05 freq=20 type=Gaussian\n'
        text+=f'rec x=1800 y=1800 z=200 file=ascii sta=station nsew={orientation} variables={mode} sacformat=1 usgsformat=1 writeEvery=9\n'
        text+=f'rechdf5 infile=stations.h5 outfile=receivers.h5 variables={mode} downsample={downsample} writeEvery=9\n'
        (directory/'case.in').write_text(text)
        command=['srun','--exclusive','--exact','-N','1','-n',str(args.tasks),'-c','4']
        command += ['--gres=none'] if args.backend=='OPENMP' else ['--gpus-per-task=1','--gpu-bind=single:1']
        env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
        env.update(OMP_NUM_THREADS='1',MPICH_GPU_SUPPORT_ENABLED='0')
        with (directory/'run.log').open('w') as log:
            result=subprocess.run(command+[str(exe),'case.in'],cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
        if result.returncode:raise RuntimeError(f'Solver failed: {directory}/run.log')
        check_solver_log(directory/'run.log')
        report[name]=check(directory,mode,orientation,downsample)
        if args.reader and mode=='displacement' and orientation==0 and downsample==1:
            # The native reader executable runs on CPUs even for GPU-written files.
            reader_command=['srun','--exclusive','--exact','--gres=none','-N','1','-n',str(args.tasks),'-c','4',
                            str(args.reader.resolve(strict=True)),str(directory/'case.in'),str(directory/'output')]
            with (directory/'reader.log').open('w') as log:
                reader=subprocess.run(reader_command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
            check_solver_log(directory/'reader.log')
            if reader.returncode or 'PASS:' not in (directory/'reader.log').read_text():
                raise RuntimeError(f'Reader round trip failed: {directory}/reader.log')
        (root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PASS:',name,flush=True)


if __name__=='__main__':main()
