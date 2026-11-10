#!/usr/bin/env python3
"""Use real MPI/HDF5 to check independent receiver metadata groups.

Four ranks form two simulation groups. The actual reader processes different
8- and 2,053-station files; receiver construction is a stub with a real
within-group collective. This checks distribution/isolation, not solver physics.
"""
import argparse
import os
from pathlib import Path
import subprocess

import h5py
import numpy as np

from receiver_template_metadata import STUBS, extract


MAIN=r'''
int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  try {
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    if(world_size!=4 || argc!=3) throw runtime_error("requires four MPI ranks and two templates");
    int group=world_rank/2;
    MPI_Comm_split(MPI_COMM_WORLD,group,world_rank,&station_group);
    MPI_Comm_rank(station_group,&rank_id);
    EW ew; ew.topo=false;
    vector<vector<TimeSeries *> > stations(1);
    readStationHDF5(&ew,argv[group+1],"output.h5",137,3,
                   TimeSeries::Velocity,0,&stations,100,100,false,
                   false,false,0,0,true,true,true,0,false,1);
    int count=group==0 ? 8 : 2053;
    if((int)stations[0].size()!=count) throw runtime_error("incorrect group receiver count");
    for(int k=0;k<count;++k) {
      ostringstream name;
      name << (group==0 ? 'a' : 'b') << setfill('0') << setw(6) << k;
      TimeSeries *s=stations[0][k];
      if(s->name!=name.str() || s->x!=k%99 || s->y!=group+1 || s->z!=0)
        throw runtime_error("receiver metadata crossed simulation groups");
    }
    int total_opens,total_visits;
    MPI_Allreduce(&opens,&total_opens,1,MPI_INT,MPI_SUM,station_group);
    MPI_Allreduce(&visits,&total_visits,1,MPI_INT,MPI_SUM,station_group);
    if(total_opens!=1 || total_visits!=1) throw runtime_error("template not read exactly once per group");
    for(TimeSeries *s:stations[0]) if(s->ts0==s) delete s->handle;
    for(TimeSeries *s:stations[0]) delete s;
    MPI_Comm_free(&station_group);
    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank==0) cout << "PASS: two independent MPI groups retain their own receiver metadata; one read per group\n";
  } catch(const exception &error) {
    cerr << error.what() << endl;
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  MPI_Finalize();
}
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hdf5-prefix',type=Path,required=True)
    parser.add_argument('--mpi-prefix',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--job-id')
    parser.add_argument('--cxx',help='MPI C++ wrapper; default: MPI prefix/bin/mpicxx')
    args=parser.parse_args()
    root=args.work_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    stubs=STUBS.replace('static const MPI_Comm station_group = (MPI_Comm)17;',
                        'static MPI_Comm station_group;')
    for macro in ['#define MPI_Bcast test_bcast','#define MPI_Abort test_abort',
                  '#define MPI_Wtime test_time']:
        stubs=stubs.replace(macro,'')
    old='mode(m), every(e), down(d), event(ev), topo(t), xyz(component) {}'
    new='''mode(m), every(e), down(d), event(ev), topo(t), xyz(component) {
      int one=1,total;
      MPI_Allreduce(&one,&total,1,MPI_INT,MPI_SUM,station_group);
      if(total!=2) throw runtime_error("receiver construction used the wrong group");
    }'''
    assert old in stubs
    stubs=stubs.replace(old,new)
    unit=root/'groups.cpp'
    unit.write_text(stubs+extract((args.source_root/'src/readhdf5.C').read_text())+MAIN)
    program=root/'groups'
    command=[args.cxx or str(args.mpi_prefix/'bin/mpicxx'),'-std=c++11','-DUSE_HDF5',
             '-I'+str(args.hdf5_prefix/'include'),'-I'+str(args.mpi_prefix/'include'),
             str(unit),'-L'+str(args.hdf5_prefix/'lib'),
             '-Wl,-rpath,'+str(args.hdf5_prefix/'lib'),'-lhdf5',
             '-Wl,-rpath,'+str(args.mpi_prefix/'lib'),'-o',str(program)]
    subprocess.run(command,check=True)
    files=[]
    for group,count in enumerate((8,2053)):
        filename=root/('group-'+str(group)+'.h5');files.append(filename)
        with h5py.File(filename,'w') as stream:
            for k in range(count):
                station=stream.create_group(('a' if group==0 else 'b')+'{:06d}'.format(k))
                station['STX,STY,STZ']=np.asarray([k%99,group+1,0],dtype=np.float64)
    launch=['srun','--overlap','--exact','--nodes=1','--ntasks=4',
            '--ntasks-per-node=4','--cpus-per-task=2']
    if args.job_id: launch.append('--jobid='+args.job_id)
    launch+=[str(program),*[str(path) for path in files]]
    result=subprocess.run(launch,cwd=root,env=dict(os.environ,OMP_NUM_THREADS='1'),
                          text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                          timeout=120)
    (root/'mpi.log').write_text(result.stdout)
    assert result.returncode==0,result.stdout[-8000:]
    print(result.stdout)


if __name__=='__main__': main()
