#!/usr/bin/env python3
"""Link a test-only main against an existing backend's actual SW4 objects.

Use the same compiler environment as the production build. Generated main,
commands and native results belong in a new task-owned work directory. Execute
only in an authorized allocation; this program does not submit a job.
"""
from pathlib import Path
import argparse,json,os,re,shlex,subprocess,sys
def execute(command,cwd):
 print(shlex.join([str(x) for x in command]),flush=True)
 subprocess.run(command,cwd=cwd,check=True)
def remap(tokens,mapping):return [mapping.get(x,x) for x in tokens]
def build(source, dest, build_dir=None, make_optdir=None, revision=None, python=sys.executable):
 source=source.resolve();dest=dest.resolve();dest.mkdir(parents=True,exist_ok=False)
 cuda = bool(make_optdir) or 'CUDA_DEFINES' in (build_dir/'CMakeFiles/sw4.dir/flags.make').read_text()
 e={'backend':'CUDA' if cuda else 'OPENMP'}
 main=dest/('native_main.cu' if e['backend']=='CUDA' else 'native_main.C');obj=dest/'native_main.o';dobj=dest/'native_link.o';binary=dest/'objectives-native'
 execute([python,str(source/'pytest/receiver_objective_numerics.py'),'--source-root',str(source),'--emit-native-main',str(main)],source)
 if make_optdir:
  result=subprocess.run(['make','-n','-W','src/main.C','sw4','optdir='+make_optdir,'GIT_VERSION='+revision],cwd=source,text=True,capture_output=True,check=True)
  (dest/'make-commands.txt').write_text(result.stdout);(dest/'make-query-stderr.log').write_text(result.stderr)
  lines=result.stdout.splitlines();compileline=next(x for x in lines if ' -c ../src/main.C' in x)
  linkline=next(x for x in lines if ' -dlink ' in x)
  finalline=next(x for x in lines if ' -o sw4 main.o file_link.o ' in x)
  cwd=source/make_optdir
  command=shlex.split(compileline.split(';',1)[1]);command=remap(command,{'../src/main.C':str(main)})+['-o',str(obj)]
  execute(command,cwd)
  command=remap(shlex.split(linkline.split(';',1)[1]),{'main.o':str(obj),'file_link.o':str(dobj)})
  execute(command,cwd)
  command=remap(shlex.split(finalline.split(';',1)[1]),{'main.o':str(obj),'file_link.o':str(dobj),'sw4':str(binary)})
  execute(command,cwd)
 else:
  cwd=build_dir.resolve();target=cwd/'CMakeFiles/sw4.dir'
  flags={}
  for line in (target/'flags.make').read_text().splitlines():
   if ' = ' in line and not line.startswith('#'):
    key,value=line.split(' = ',1);flags[key]=value
  text=(target/'build.make').read_text().splitlines()
  oldobj='CMakeFiles/sw4.dir/'+('gpu_sources/main.cu.o' if e['backend']=='CUDA' else 'src/main.C.o')
  line=next(x.strip() for x in text if x.startswith('\t') and ' -c ' in x and ' -o '+oldobj in x)
  for key,value in flags.items():line=line.replace('$('+key+')',value)
  assert '$(' not in line,line
  command=shlex.split(line)
  cindex=command.index('-c');command[cindex+1]=str(main)
  command=remap(command,{oldobj:str(obj),oldobj+'.d':str(dest/'native_main.d')})
  execute(command,cwd)
  mapping={oldobj:str(obj),'CMakeFiles/sw4.dir/cmake_device_link.o':str(dobj),'bin/sw4':str(binary)}
  if e['backend']=='CUDA':
   response=dest/'deviceObjects.rsp';response.write_text((target/'deviceObjects1.rsp').read_text().replace(oldobj,str(obj)))
   command=remap(shlex.split((target/'dlink.txt').read_text()),dict(mapping,**{'CMakeFiles/sw4.dir/deviceObjects1.rsp':str(response)}))
   execute(command,cwd)
  execute(remap(shlex.split((target/'link.txt').read_text()),mapping),cwd)
 return binary
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
 p.add_argument('--work-dir',type=Path,required=True)
 p.add_argument('--build-dir',type=Path,help='Existing CMake build with production SW4 objects')
 p.add_argument('--make-optdir',help='Existing RAJA make output directory name')
 p.add_argument('--revision',help='Version metadata for an archived RAJA source tree')
 a=p.parse_args()
 assert bool(a.build_dir)!=bool(a.make_optdir),'select exactly one existing build route'
 revision=a.revision
 if a.make_optdir and not revision:revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source_root,text=True).strip()
 print(build(a.source_root,a.work_dir,a.build_dir,a.make_optdir,revision))

if __name__=='__main__':main()
