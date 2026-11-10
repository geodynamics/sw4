#!/usr/bin/env python3
"""Check complete production receiver objectives against numerical/analytic oracles.

Only construction/storage and EW diagnostic access are replaced. The full selected
misfit, shift correlation, maximum search, adjoint-source and interpolant methods
are compiled unchanged. This checks receiver data operators, not the PDE adjoint.
Use --emit-native-main to exercise real TimeSeries objects with native backend
objects and allocators; run that generated executable in an authorized allocation.
"""
from pathlib import Path
import re,subprocess,json,argparse,hashlib,math

def end(text,start):
 a=text.index('{',start);level=0
 for i in range(a,len(text)):
  if text[i]=='{':level+=1
  elif text[i]=='}':
   level-=1
   if not level:return i+1
 raise ValueError('unbalanced method')
def selected(root,gpu,cxx="c++"):
 command=[cxx,'-E','-P','-x','c++']+(['-DSW4_USE_RAJA'] if gpu else [])+['-']
 src=(root/'src/TimeSeries.C').read_text();src=re.sub(r'^\s*#include[^\n]*','',src,flags=re.M)
 src=subprocess.check_output(command,input=src,text=True)
 names=['misfit','shiftfunc','compute_maxshift','misfit2','getwgh','getwgh5','set_window','window_weight']
 methods=[]
 for name in names:
  for m in re.finditer(r'(?:float_sw4|void) TimeSeries::'+name+r'\s*\([^;]*?\)\s*(?:const\s*)?\{',src):methods.append(src[m.start():end(src,m.start())])
 if sum('void TimeSeries::set_window' in m for m in methods)==1:
  h=(root/'src/TimeSeries.h').read_text();h=re.sub(r'^\s*#include[^\n]*','',h,flags=re.M);h=subprocess.check_output(command,input=h,text=True)
  for m in re.finditer(r'void set_window\s*\([^;]*?\)\s*\{',h):methods.append(h[m.start():end(h,m.start())].replace('void set_window','void TimeSeries::set_window',1))
 for m in re.finditer(r'static float_sw4 receiver_window_weight\([^;]*?\)\s*\{',src):methods.insert(0,src[m.start():end(src,m.start())])
 return '\n'.join(methods)
STUB=r'''
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
using namespace std;
using float_sw4=double;
struct EW {void get_utc(int*){};};
class TimeSeries {
public:
 bool m_myPoint=true,m_usgsFormat=true,m_hdf5Format=false,m_use_win=false,m_use_x=true,m_use_y=true,m_use_z=true;
 double m_t0=0,m_shift=0,m_dt=.01,m_winL=-1e38,m_winR=1e38,m_winL2=-1e38,m_winR2=-1e38;
 int mLastTimeStep=0,mAllocatedSize=0,m_nComp=3,m_misfit_scaling=1,m_event=0,m_global_event=0,m_utc[7]{};
 string m_staName="probe";EW dummy;EW* m_ew=&dummy;
 vector<vector<double>> data;vector<vector<float>> fdata;double* mRecordedSol[3];float* mRecordedFloats[3];
 TimeSeries(int n,double start,double dt,bool as_float=false):m_usgsFormat(!as_float),m_t0(start),m_dt(dt){allocateRecordingArrays(n-1,start,dt);}
 void allocateRecordingArrays(int last,double start,double dt){mLastTimeStep=last;mAllocatedSize=last+1;m_shift=start-m_t0;m_dt=dt;data.assign(3,vector<double>(last+1,0));fdata.assign(3,vector<float>(last+1,0));for(int q=0;q<3;++q){mRecordedSol[q]=data[q].data();mRecordedFloats[q]=fdata[q].data();}}
 double** getRecordingArray(){return mRecordedSol;}float** getRecordingArrayFloats(){return mRecordedFloats;}
 void set_window(double,double);void set_window(double,double,double,double);double window_weight(double)const;
 double misfit(TimeSeries&,TimeSeries*,double&,double&,double&);double misfit2(TimeSeries&,TimeSeries*);
 double compute_maxshift(TimeSeries&);
 void shiftfunc(TimeSeries&,double,double&,double&,double&,double** =nullptr);
 void getwgh(double,double[6],double[6],double[6]);void getwgh5(double,double[6],double[6],double[6]);
 void put(int q,int i,double value){data[q][i]=value;fdata[q][i]=value;}
 double get(int q,int i)const{return m_usgsFormat?data[q][i]:fdata[q][i];}
};
'''
MAIN = r'''
double wave(double t,int q){double x=(t-3.1)/.45;return (1+.2*q)*exp(-x*x)*(cos((4+q)*x)+.13*sin(2*x));}
void windows(TimeSeries&u,int w){if(w==1)u.set_window(2.2,3.9);if(w==2)u.set_window(2.1,2.9,3.2,4.1);if(w==3)u.set_window(2.0,3.5,3.0,4.2);}
int main(){cout<<setprecision(17);
 for(int w=0;w<4;++w)for(int scaling=0;scaling<2;++scaling)for(int floats=0;floats<2;++floats)for(int mask=0;mask<2;++mask){
  TimeSeries u(301,1.5,.01,floats),o(801,0,.01,floats),d(301,1.5,.01,false);u.m_misfit_scaling=scaling;
  windows(u,w); u.m_use_x = mask==0;
  for(int q=0;q<3;++q){for(int i=0;i<=o.mLastTimeStep;++i)o.put(q,i,wave(i*.01,q));for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,1.13*wave(1.5+i*.01-.137,q)+.01*sin(.1*i+q));}
  double g,h,h1;double value=u.misfit(o,&d,g,h,h1);
  double e=1e-5,a,b,c,ag,ah,ah1,bg,bh,bh1;
  o.m_shift=e;a=u.misfit(o,nullptr,ag,ah,ah1);o.m_shift=-e;b=u.misfit(o,nullptr,bg,bh,bh1);o.m_shift=0;
  double fdg=(a-b)/(2*e),fdh=(ag-bg)/(2*e);
  double adj=0;for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)adj+=d.data[q][i]*sin(.13*i+.7*q);
  double pe=floats?1e-3:1e-5;
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,u.get(q,i)+pe*sin(.13*i+.7*q));
  a=u.misfit(o,nullptr,ag,ah,ah1);
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,u.get(q,i)-2*pe*sin(.13*i+.7*q));
  b=u.misfit(o,nullptr,bg,bh,bh1);
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,u.get(q,i)+pe*sin(.13*i+.7*q));
  double adjfd=(a-b)/(2*pe);
  double f,df,ddf,fp,dfp,ddfp,fm,dfm,ddfm;
  u.shiftfunc(o,.137,f,df,ddf);u.shiftfunc(o,.137+e,fp,dfp,ddfp);u.shiftfunc(o,.137-e,fm,dfm,ddfm);
  double tt=u.misfit2(o,&d),ttadj=0;
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)ttadj+=d.data[q][i]*sin(.13*i+.7*q);
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,u.get(q,i)+pe*sin(.13*i+.7*q));
  a=u.misfit2(o,nullptr);
  for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,u.get(q,i)-2*pe*sin(.13*i+.7*q));
  b=u.misfit2(o,nullptr);
  cout<<"CASE "<<w<<" "<<scaling<<" "<<floats<<" "<<mask<<" "<<value<<" "<<g<<" "<<fdg<<" "<<h<<" "<<fdh<<" "<<adj<<" "<<adjfd<<" "<<f<<" "<<df<<" "<<(fp-fm)/(2*e)<<" "<<ddf<<" "<<(dfp-dfm)/(2*e)<<" "<<tt<<" "<<ttadj<<" "<<(a-b)/(2*pe)<<"\n";
 }

 // Independent closed-form quotient-rule oracle: linear data interpolate exactly.
 {TimeSeries u(41,1,.01),o(501,0,.01),d(41,1,.01);u.set_window(-100,100);u.m_use_y=u.m_use_z=false;
 double A=0,A1=0,A2=0,B=0,B1=0,B2=0;
 for(int i=0;i<=o.mLastTimeStep;++i)o.put(0,i,2+.7*i*.01);
 for(int i=0;i<=u.mLastTimeStep;++i){double t=1+i*.01,obs=2+.7*t,synt=1.2+.2*t,r=obs-synt;u.put(0,i,synt);A+=.5*r*r;A1-=r*.7;A2+=.7*.7;B+=obs*obs;B1-=2*obs*.7;B2+=2*.7*.7;}
 double g,h,h1,J=u.misfit(o,&d,g,h,h1),jg=A1/B-A*B1/(B*B),jh=A2/B-2*A1*B1/(B*B)-A*B2/(B*B)+2*A*B1*B1/(B*B*B);
 if(abs(J-A/B)>1e-12 || abs(g-jg)>1e-9 || abs(h-jh)>1e-8 || abs(h1-A2/B)>1e-8) throw runtime_error("linear normalized objective derivatives disagree with closed form");
 u.m_use_x=false;u.m_use_y=true;
 for(int q=0;q<3;++q){for(int i=0;i<=o.mLastTimeStep;++i)o.put(q,i,q==1?1:0);for(int i=0;i<=u.mLastTimeStep;++i)u.put(q,i,q==1?1:0);}
 double f,df,ddf;u.shiftfunc(o,.123,f,df,ddf);if(abs(f-41)>1e-8)throw runtime_error("excluding X also excluded valid Y correlation");
 cout<<"PASS analytic quotient-rule and component-selection controls\n";}

 // Independent quartic interpolation oracle at the last six-/five-point transition.
 {TimeSeries u(10,1.015,.1),o(21,0,.1),d(10,1.015,.1);u.set_window(-100,100);u.m_misfit_scaling=0;u.m_use_y=u.m_use_z=false;
 auto polynomial=[](double t){return 1+t+t*t+.1*t*t*t+.03*t*t*t*t;};
 auto derivative=[](double t){return 1+2*t+.3*t*t+.12*t*t*t;};
 for(int i=0;i<=o.mLastTimeStep;++i)o.put(0,i,polynomial(i*.1));
 double expected_curvature=0;
 for(int i=0;i<=u.mLastTimeStep;++i){double t=1.015+i*.1;u.put(0,i,polynomial(t));expected_curvature+=derivative(t)*derivative(t);}
 double g,h,h1,J=u.misfit(o,&d,g,h,h1);
 if(J>1e-20 || abs(g)>1e-9 || abs(h-expected_curvature)>1e-8*expected_curvature)throw runtime_error("endpoint quartic interpolant differs from exact polynomial");
 cout<<"PASS analytic endpoint interpolation and derivatives\n";}
 // Exact-size endpoint buffers expose both float and double stencil overreads.
 for(int floats=0;floats<2;++floats){TimeSeries u(301,5.505,.01,floats),o(801,0,.01,floats),d(301,5.505,.01);double g,h,h1,f,df,ddf;u.set_window(5,9);u.misfit(o,&d,g,h,h1);u.shiftfunc(o,.003,f,df,ddf);double tt=u.misfit2(o,&d);for(int q=0;q<3;++q)for(int i=0;i<=u.mLastTimeStep;++i)if(!isfinite(d.data[q][i]))throw runtime_error("zero-curvature travel-time source is nonfinite");if(!isfinite(tt))throw runtime_error("zero-data travel-time score is nonfinite");}
 cout<<"PASS endpoint storage and zero-curvature controls\n";
}
'''

def run(root, gpu, work, compiler, asan):
    work.mkdir(parents=True, exist_ok=False)
    methods=selected(root,gpu,compiler)
    unit=work/'objectives.cpp';unit.write_text(STUB+methods+MAIN)
    command=[compiler,'-std=c++11','-O1',str(unit),'-o',str(work/'objectives')]
    if asan:command+=['-fsanitize=address','-fno-omit-frame-pointer']
    subprocess.run(command,check=True)
    result=subprocess.run([str(work/'objectives')],capture_output=True,text=True)
    (work/'stdout.log').write_text(result.stdout);(work/'stderr.log').write_text(result.stderr)
    assert result.returncode==0,(result.returncode,result.stdout[-4000:],result.stderr[-4000:])
    rows=[[float(x) for x in line[5:].split()] for line in result.stdout.splitlines() if line.startswith('CASE ')]
    assert len(rows)==32
    comparisons=[('L2 shift derivative',5,6),('L2 shift curvature',7,8),('L2 data adjoint',9,10),('correlation derivative',12,13),('correlation curvature',14,15),('travel-time data adjoint',17,18)]
    for row in rows:
        assert all(math.isfinite(x) for x in row)
        for name,i,j in comparisons:
            tol=5e-3 if row[2] and name.endswith('adjoint') else 2e-5
            assert abs(row[i]-row[j])<=tol*max(.001,abs(row[i]),abs(row[j])),(name,row[:4],row[i],row[j])
    (work/'results.json').write_text(json.dumps(dict(source=str(root),gpu_path=gpu,source_sha256=hashlib.sha256((root/'src/TimeSeries.C').read_bytes()).hexdigest(),method_sha256=hashlib.sha256(methods.encode()).hexdigest(),rows=rows,asan=asan),indent=2)+'\n')
    print('PASS complete production objectives, derivatives and data adjoints:',root,'GPU path' if gpu else 'native path')
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--work-dir',type=Path)
    p.add_argument('--cpu-reference',type=Path)
    p.add_argument('--cxx',default='c++')
    p.add_argument('--raja',action='store_true',help='Select unified GPU implementation')
    p.add_argument('--asan',action='store_true')
    p.add_argument('--emit-native-main',type=Path,help='Generate a native test main for linking with this backend')
    args=p.parse_args();root=args.source_root.resolve()
    if args.emit_native_main:
        source=(root/'src/main.C').read_text()
        header=(Path(__file__).resolve().parent/'receiver_objective_native.h').resolve()
        source=source.replace('using namespace std;','using namespace std;\n#include '+json.dumps(str(header)),1)
        source,n=re.subn(r'simulation\.solve\s*\([^;]+;', 'receiver_objective_checks(simulation);',source,flags=re.S)
        assert n in (1,2),n
        args.emit_native_main.write_text(source)
        return
    assert args.work_dir,'--work-dir is required for numerical checks'
    rows=run(root,args.raja,args.work_dir/'candidate',args.cxx,args.asan)
    if args.cpu_reference:
        reference=run(args.cpu_reference.resolve(),False,args.work_dir/'cpu-reference',args.cxx,args.asan)
        for a,b in zip(rows,reference):
            assert a[:4]==b[:4]
            for i in [4,5,7,9,11,12,14,16,17]:
                assert abs(a[i]-b[i])<=2e-5*max(.001,abs(a[i]),abs(b[i])),('CPU parity',a[:4],i,a[i],b[i])
        print('PASS CPU reference parity')

if __name__=='__main__':main()
