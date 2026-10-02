#!/usr/bin/env python3
"""Explicit real-reader checks for malformed histories, endian and instrument basis."""
import argparse
from datetime import datetime, timedelta
import os
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np
from receiver_formats import read_sac


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reader',type=Path,required=True)
    parser.add_argument('--fixture',type=Path,required=True,help='Flat displacement/XYZ format case directory')
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args()
    fixture=args.fixture.resolve(strict=True);root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=True)
    env={k:v for k,v in os.environ.items() if not k.startswith('SLURM_') or k=='SLURM_JOB_ID'}
    env.update(OMP_NUM_THREADS='1')
    cases=('sac-missing','sac-truncated','sac-nan','sac-wrong-quantity','sac-different-dt',
           'hdf-missing','hdf-wrong-unit','hdf-nan','hdf-count-shape','hdf-short-component',
           'hdf-long-string','hdf-wrong-basis','hdf-zero-dt','hdf-wrong-downsample',
           'usgs-incomplete-row','usgs-wrong-quantity','sac-big-endian','instrument13','instrument45',
           'legacy-grid','ignore-utc')
    for name in cases:
        work=root/name;work.mkdir(exist_ok=True)
        output=work/'output';shutil.copytree(fixture/'output',output,dirs_exist_ok=True)
        (work/'case.in').write_text((fixture/'case.in').read_text())
        x=output/'ascii.x'
        if name.startswith('sac-'):
            data=bytearray(x.read_bytes())
            if name=='sac-missing':x.unlink()
            elif name=='sac-truncated':x.write_bytes(data[:634])
            elif name=='sac-nan':data[632:636]=np.float32(np.nan).tobytes();x.write_bytes(data)
            elif name=='sac-wrong-quantity':data[440+20*8:440+21*8]=b'Vx      ';x.write_bytes(data)
            elif name=='sac-different-dt':data[:4]=np.float32(.02).tobytes();x.write_bytes(data)
            elif name=='sac-big-endian':
                for path in output.glob('ascii.[xyz]'):
                    data=path.read_bytes()
                    path.write_bytes(np.frombuffer(data[:440],dtype='<u4').byteswap().tobytes()+data[440:632]+
                        np.frombuffer(data[632:],dtype='<u4').byteswap().tobytes())
        elif name.startswith('hdf-'):
            with h5py.File(output/'receivers.h5','r+') as f:
                g=f['station']
                if name=='hdf-missing':del g['Z']
                elif name=='hdf-wrong-unit':f.attrs['UNIT']='m/s'
                elif name=='hdf-nan':g['X'][2]=np.nan
                elif name=='hdf-count-shape':del g['NPTS'];g['NPTS']=[49,49]
                elif name=='hdf-short-component':
                    data=g['Z'][:int(g['NPTS'][0])-1];del g['Z'];g['Z']=data
                elif name=='hdf-long-string':f.attrs['DATETIME']='x'*256
                elif name=='hdf-wrong-basis':g['ISNSEW'][...]=1
                elif name=='hdf-zero-dt':f['DELTA'][...]=0
                elif name=='hdf-wrong-downsample':f['DOWNSAMPLE'][...]=3
        elif name.startswith('usgs-'):
            path=output/'ascii.txt';text=path.read_text()
            if name=='usgs-incomplete-row':
                lines=text.splitlines();lines[-1]=' '.join(lines[-1].split()[:-1]);text='\n'.join(lines)+'\n'
            else:text=text.replace('X displacement (m)','X velocity (m/s)')
            path.write_text(text)
        elif name=='legacy-grid':
            for path in output.glob('ascii.[xyz]'):
                data=bytearray(path.read_bytes())
                data[440+17*8:440+18*8]=b' '*8
                path.write_bytes(data)
        elif name=='ignore-utc':
            # Change absolute references by a day, leaving relative sample
            # times unchanged. Remove SAC's sub-millisecond reference residue
            # from B/E/O so its relative time starts at zero like text/HDF5.
            for path in output.glob('ascii.[xyz]'):
                data=bytearray(path.read_bytes())
                ints=np.frombuffer(data,dtype='<i4',count=40,offset=280)
                ints[1]+=1
                floats=np.frombuffer(data,dtype='<f4',count=70)
                floats[[5,6,7]]-=.000567
                path.write_bytes(data)
            with h5py.File(output/'receivers.h5','r+') as f:
                value=f.attrs['DATETIME']
                if isinstance(value,bytes): value=value.decode()
                f.attrs['DATETIME']=(datetime.fromisoformat(value)+timedelta(days=1)).isoformat()
            path=output/'ascii.txt'
            path.write_text(path.read_text().replace('10/02/2026:','10/03/2026:'))
        else:
            angle=float(name.removeprefix('instrument'))
            traces=[read_sac(output/f'ascii.{c}') for c in 'xyz']
            header=traces[0][2];xx,yy,zz=[item[1] for item in traces]
            north=header[40]*xx+header[41]*yy;east=header[42]*xx+header[43]*yy
            a=np.deg2rad(angle)
            signals=(np.cos(a)*north+np.sin(a)*east,-np.sin(a)*north+np.cos(a)*east,-zz)
            for channel,(signal,az,inc) in zip('xyz',zip(signals,(angle,angle+90,0),(90,90,0))):
                path=output/f'ascii.{channel}';data=bytearray(path.read_bytes())
                data[440+17*8:440+18*8]=b' '*8
                data[57*4:59*4]=np.asarray([az,inc],dtype='<f4').tobytes()
                data[632:]=np.asarray(signal,dtype='<f4').tobytes();path.write_bytes(data)
        command=['srun','--exclusive','--exact','--gres=none','-N','1','-n','2','-c','1',
                 str(args.reader.resolve(strict=True)),str(work/'case.in'),str(output),'displacement','0']
        if name in ('legacy-grid','ignore-utc'):
            command+=['1','grid' if name=='legacy-grid' else 'ignore-utc']
        with (work/'reader.log').open('w') as log:
            result=subprocess.run(command,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=120)
        log=(work/'reader.log').read_text()
        positive=name in ('sac-big-endian','instrument13','instrument45','legacy-grid','ignore-utc')
        if positive:
            valid=result.returncode==0 and 'PASS:' in log
        else:
            valid=result.returncode!=0 and any(message in log for message in
                ('SAC receiver read failed','Receiver HDF5 read failed','Fatal input error:'))
            valid=valid and 'Segmentation fault' not in log
        if not valid:raise RuntimeError(f'Unexpected reader result: {work}/reader.log\n{log[-3000:]}')
        print('PASS:',name,'accepted' if positive else 'rejected',flush=True)


if __name__=='__main__':main()
