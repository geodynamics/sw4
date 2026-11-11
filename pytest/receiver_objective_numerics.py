#!/usr/bin/env python3
"""Check receiver windows, objectives and data adjoints.

Run windows and numerics locally against extracted production methods.
Use native-build with an existing HDF5-enabled SW4 build and the same
compiler environment to link a driver against its production objects.
Run the resulting objectives-native binary in an authorized allocation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shlex
import subprocess
import tempfile


def end(text, start):
    opening = text.index("{", start)
    level = 0
    for index in range(opening, len(text)):
        if text[index] == "{":
            level += 1
        elif text[index] == "}":
            level -= 1
            if level == 0:
                return index + 1
    raise ValueError("unterminated native function")


def preprocess(path, compiler, gpu):
    command = [compiler, "-E", "-P", "-x", "c++"]
    if gpu:
        command.append("-DSW4_USE_RAJA")
    source = re.sub(r"^\s*#include[^\n]*", "", path.read_text(), flags=re.M)
    return subprocess.check_output(command + ["-"], input=source, text=True)


def selected(root,gpu,cxx="c++"):
 src=preprocess(root / "src/TimeSeries.C", cxx, gpu)
 names=['misfit','shiftfunc','compute_maxshift','misfit2','getwgh','getwgh5','set_window','window_weight']
 methods=[]
 for name in names:
  for m in re.finditer(r'(?:float_sw4|void) TimeSeries::'+name+r'\s*\([^;]*?\)\s*(?:const\s*)?\{',src):methods.append(src[m.start():end(src,m.start())])
 if sum('void TimeSeries::set_window' in m for m in methods)==1:
  h=preprocess(root / "src/TimeSeries.h", cxx, gpu)
  for m in re.finditer(r'void set_window\s*\([^;]*?\)\s*\{',h):methods.append(h[m.start():end(h,m.start())].replace('void set_window','void TimeSeries::set_window',1))
 for m in re.finditer(r'static float_sw4 receiver_window_weight\([^;]*?\)\s*\{',src):methods.insert(0,src[m.start():end(src,m.start())])
 return '\n'.join(methods)
STUB=r'''
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
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
 double misfit(TimeSeries&,TimeSeries*,double&,double&,double&);double misfit2(TimeSeries&,TimeSeries*,bool* = nullptr);
 double compute_maxshift(TimeSeries&,bool&);
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
  bool valid=false;
  double tt=u.misfit2(o,&d,&valid),ttadj=0;
  if(!valid) throw runtime_error("isolated matrix peak was excluded");
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
 // A linear observation has nonzero slope and zero or near-zero curvature.
 for(double dt : {.125,.1}) {
  double start=dt==.125 ? 3 : 3.017;
  TimeSeries u(25,start,dt),o(97,0,dt),d(25,start,dt);
  u.m_use_y = u.m_use_z = false;
  u.put(0,0,1);
  for(int i=0;i<=o.mLastTimeStep;++i)o.put(0,i,3*(i-40));
  double f,df,ddf;
  u.shiftfunc(o,2,f,df,ddf);
  if(abs(df-3/dt)>1e-9 || abs(ddf)>1e-9)
    throw runtime_error("invalid linear correlation fixture");
  bool resolved=true;
  for(int q=0;q<3;++q) for(int i=0;i<=u.mLastTimeStep;++i) d.put(q,i,1);
  double score=u.misfit2(o,&d,&resolved);
  if(resolved || score!=0)
    throw runtime_error("unresolved window contributed to the objective");
  for(int q=0;q<3;++q)
    for(int i=0;i<=u.mLastTimeStep;++i)
      if(d.data[q][i]!=0) throw runtime_error("unresolved window left adjoint forcing");
 }
 cout<<"PASS exclusion of zero and near-zero curvature windows\n";

 // A coarse neighbor can lie beyond a Ricker sidelobe trough even though
 // Newton refinement from the best sample resolves the translated peak.
 for(double width : {.05,.1,.2}) for(double delay : {.03,.137,.73,1.3}) {
  if(width==.05 && delay==.73) continue; // Not a valid prior CPU peak.
  TimeSeries u(801,0,.01),o(801,0,.01),d(801,0,.01);
  u.m_use_y=u.m_use_z=false;
  for(int i=0;i<=u.mLastTimeStep;++i) {
   double a=(i*.01-4-delay)/width,b=(i*.01-4)/width;
   u.put(0,i,(1-2*a*a)*exp(-a*a));
   o.put(0,i,(1-2*b*b)*exp(-b*b));
  }
  bool valid=false;
  double score=u.misfit2(o,&d,&valid);
  if(!valid || abs(sqrt(2*score)-delay)>2e-5)
   throw runtime_error("resolved Ricker peak differs from CPU translation control");
 }
 cout<<"PASS resolved Ricker peaks across frequency and coarse-lag phase\n";

 // Analytic weak peak: its curvature scales with amplitude, its lag does not.
 for(int exponent : {-10,-50}) {
  TimeSeries u(41,4,.125),o(97,0,.125),d(41,4,.125);
  u.m_use_y = u.m_use_z = false;
  u.put(0,0,1);
  for(int i=0;i<=o.mLastTimeStep;++i) {
   double x=i-(4+.1875)/.125;
   o.put(0,i,ldexp(1.,exponent)*(64-x*x));
  }
  bool valid=false;
  double score=u.misfit2(o,&d,&valid);
  if(!valid || abs(score-.5*.1875*.1875)>1e-12 || abs(d.data[0][1]+.1875*.125)>1e-12)
    throw runtime_error("resolvable weak peak disagrees with analytic objective/adjoint");
 }
 {TimeSeries u(41,4,.125),o(97,0,.125),d(41,4,.125);
  u.m_use_y=u.m_use_z=false;
  u.put(0,0,1);
  for(int i=0;i<=o.mLastTimeStep;++i) {
    const double periodic[4]={1,1,-1,-1};
    o.put(0,i,periodic[i%4]);
  }
  bool valid=true;
  double score=u.misfit2(o,&d,&valid);
  if(valid || score!=0) throw runtime_error("ambiguous periodic peaks were accepted");
 }
 cout<<"PASS weak isolated peak and competing peak controls\n";

}
'''

def check_numerics(root, gpu, work, compiler, asan):
    work.mkdir(parents=True, exist_ok=False)
    methods=selected(root,gpu,compiler)
    unit=work/'objectives.cpp';unit.write_text(STUB+methods+MAIN)
    command=[compiler,'-std=c++11','-O1',str(unit),'-o',str(work/'objectives')]
    if asan:command+=['-fsanitize=address,float-cast-overflow',
                      '-fno-sanitize-recover=float-cast-overflow','-fno-omit-frame-pointer']
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


def run_windows(args):
    source = preprocess(args.source_root / "src/TimeSeries.C", args.cxx, args.raja)
    setters = []
    for match in re.finditer(r"void TimeSeries::set_window\s*\([^)]*\)\s*\{", source):
        setters.append(source[match.start():end(source, match.start())])
    if len(setters) == 1:
        header = preprocess(args.source_root / "src/TimeSeries.h", args.cxx, args.raja)
        for match in re.finditer(r"void set_window\s*\([^)]*\)\s*\{", header):
            method = header[match.start():end(header, match.start())]
            setters.append(method.replace("void set_window", "void TimeSeries::set_window", 1))
    assert len(setters) == 2, "two- and four-value window setters are required"
    defaults = []
    for field in ("m_winL", "m_winR", "m_winL2", "m_winR2"):
        match = re.search(r"\b" + field + r"\(([^)]+)\)", source)
        assert match, "window state must be initialized: " + field
        defaults.append(field + "=" + match.group(1))
    start = source.index("TimeSeries::misfit(")
    match = re.search(r"float_sw4\s+wghx\s*,\s*wghy\s*,\s*wghz\s*;", source[start:])
    assert match
    start += match.start()
    match = re.search(r"if\s*\(\s*m_use_win\s*\)", source[start:])
    assert match
    start += match.start()
    weight = source[start:end(source, start)]
    helper = ""
    if "window_weight(t)" in weight:
        match = re.search(r"float_sw4 TimeSeries::window_weight\([^)]*\) const\s*\{", source)
        assert match
        helper = source[match.start():end(source, match.start())]
    if "receiver_window_weight(" in weight:
        match = re.search(r"static float_sw4 receiver_window_weight\([^)]*\)\s*\{", source)
        assert match
        helper = source[match.start():end(source, match.start())]
    product = ""
    match = re.search(r"^float_sw4 TimeSeries::product_wgh\(", source, re.M)
    if match:
        product = source[match.start():end(source, match.start())]
    stub = r'''
#include <cmath>
#include <iostream>
#include <stdexcept>
#define float_sw4 double
using namespace std;
static double receiver_window_weight(double,double,double,double,double,double);
class TimeSeries {
 public:
  bool m_use_win=false;
  double m_dt=0.01, m_t0=3.9;
  bool m_use_x=true, m_use_y=true, m_use_z=true;
  int mLastTimeStep=5;
  double data[3][6], *mRecordedSol[3];
  TimeSeries() { for(int q=0;q<3;++q) { mRecordedSol[q]=data[q]; for(int i=0;i<6;++i) data[q][i]=1; } }
  double product_wgh(TimeSeries&) const;
  double DEFAULTS;
  void set_window(double, double);
  void set_window(double, double, double, double);
  double window_weight(double) const;
  double weight(double t) const {
    double itau=1/(5*m_dt), wghv=0.42, wghx=wghv, wghy=wghv, wghz=wghv;
    WEIGHT
    return wghx;
  }
};
'''.replace("DEFAULTS", ", ".join(defaults)).replace("WEIGHT", weight)
    main = r'''
static void close(double actual, double expected) {
  if (std::abs(actual-expected)>1e-12)
    throw std::runtime_error("window weight differs: actual="+std::to_string(actual)+
                             " expected="+std::to_string(expected));
}
static double interval(double t, double left, double right) {
  return 0.5*(std::tanh((t-left)/0.05)-std::tanh((t-right)/0.05));
}
int main() {
  TimeSeries single;
  close(single.weight(1.5),0.42);
  single.set_window(1,2);
  for(double t : {0.0,1.0,1.5,2.0,3.0,5.0}) close(single.weight(t),interval(t,1,2));
  single.set_window(1,2,4,5);
  for(double t : {0.0,1.0,1.5,2.0,3.0,4.0,4.5,5.0,6.0})
    close(single.weight(t),std::fmax(interval(t,1,2),interval(t,4,5)));
  single.set_window(1,3,2,4);
  for(double t : {0.0,1.0,2.0,2.5,3.0,4.0,5.0})
    close(single.weight(t),std::fmax(interval(t,1,3),interval(t,2,4)));
  single.set_window(1,2);
  close(single.weight(3.5),interval(3.5,1,2));
#ifdef HAVE_PRODUCT
  single.set_window(1,2,4,5);
  double expected=0;
  for(int i=0;i<=single.mLastTimeStep;++i) {
    double t=single.m_t0+i*single.m_dt;
    expected+=std::fmax(interval(t,1,2),interval(t,4,5));
  }
  close(single.product_wgh(single),3*expected);
  single.m_use_x=false;
  close(single.product_wgh(single),2*expected);
#endif
  std::cout<<"PASS: native single/two-window weights, overlap, and mode reset\n";
}
'''
    with tempfile.TemporaryDirectory(prefix="sw4-window-unit-") as name:
        root = Path(name)
        unit = root / "windows.cpp"
        unit.write_text(("#define HAVE_PRODUCT\n" if product else "") + stub + "\n".join(setters) + helper + product + main)
        subprocess.run([args.cxx, "-std=c++11", str(unit), "-o", str(root/"windows")], check=True)
        subprocess.run([str(root/"windows")], check=True)



def emit_native_main(root, destination):
    source = (root / "src/main.C").read_text()
    header = Path(__file__).resolve().with_name("receiver_objective_native.h")
    source = source.replace("using namespace std;",
                            "using namespace std;\n#include " + json.dumps(str(header)), 1)
    source, count = re.subn(r"simulation\.solve\s*\([^;]+;", "receiver_objective_checks(simulation);",
                            source, flags=re.S)
    assert count in (1, 2), count
    destination.write_text(source)


def execute(command, cwd):
    print(shlex.join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def remap(tokens, mapping):
    return [mapping.get(token, token) for token in tokens]


def build_native(source, dest, build_dir=None, make_optdir=None, revision=None):
    source = source.resolve()
    dest = dest.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    cuda = bool(make_optdir) or "CUDA_DEFINES" in (build_dir / "CMakeFiles/sw4.dir/flags.make").read_text()
    main = dest / ("native_main.cu" if cuda else "native_main.C")
    obj = dest / "native_main.o"
    link_obj = dest / "native_link.o"
    binary = dest / "objectives-native"
    emit_native_main(source, main)

    if make_optdir:
        # Standalone RAJA records its compile and link commands in the Make dry run.
        result = subprocess.run(
            ["make", "-n", "-W", "src/main.C", "sw4", "optdir=" + make_optdir,
             "GIT_VERSION=" + revision], cwd=source, text=True, capture_output=True, check=True)
        (dest / "make-commands.txt").write_text(result.stdout)
        (dest / "make-query-stderr.log").write_text(result.stderr)
        lines = result.stdout.splitlines()
        compile_line = next(line for line in lines if " -c ../src/main.C" in line)
        device_line = next(line for line in lines if " -dlink " in line)
        link_line = next(line for line in lines if " -o sw4 main.o file_link.o " in line)
        cwd = source / make_optdir
        command = shlex.split(compile_line.split(";", 1)[1])
        execute(remap(command, {"../src/main.C": str(main)}) + ["-o", str(obj)], cwd)
        mapping = {"main.o": str(obj), "file_link.o": str(link_obj), "sw4": str(binary)}
        execute(remap(shlex.split(device_line.split(";", 1)[1]), mapping), cwd)
        execute(remap(shlex.split(link_line.split(";", 1)[1]), mapping), cwd)
        return binary

    cwd = build_dir.resolve()
    target = cwd / "CMakeFiles/sw4.dir"
    flags = {}
    for line in (target / "flags.make").read_text().splitlines():
        if " = " in line and not line.startswith("#"):
            key, value = line.split(" = ", 1)
            flags[key] = value
    original_obj = "CMakeFiles/sw4.dir/" + ("gpu_sources/main.cu.o" if cuda else "src/main.C.o")
    rules = (target / "build.make").read_text().splitlines()
    compile_line = next(line.strip() for line in rules
                        if line.startswith("\t") and " -c " in line and " -o " + original_obj in line)
    for key, value in flags.items():
        compile_line = compile_line.replace("$(" + key + ")", value)
    assert "$(" not in compile_line, compile_line
    command = shlex.split(compile_line)
    command[command.index("-c") + 1] = str(main)
    execute(remap(command, {original_obj: str(obj), original_obj + ".d": str(dest / "native_main.d")}), cwd)

    mapping = {original_obj: str(obj), "CMakeFiles/sw4.dir/cmake_device_link.o": str(link_obj),
               "bin/sw4": str(binary)}
    if cuda:
        response = dest / "deviceObjects.rsp"
        response.write_text((target / "deviceObjects1.rsp").read_text().replace(original_obj, str(obj)))
        mapping["CMakeFiles/sw4.dir/deviceObjects1.rsp"] = str(response)
        execute(remap(shlex.split((target / "dlink.txt").read_text()), mapping), cwd)
    execute(remap(shlex.split((target / "link.txt").read_text()), mapping), cwd)
    return binary


def run_numerics(args):
    rows=check_numerics(args.source_root.resolve(),args.raja,args.work_dir/'candidate',args.cxx,args.asan)
    if args.cpu_reference:
        reference=check_numerics(args.cpu_reference.resolve(),False,args.work_dir/'cpu-reference',args.cxx,args.asan)
        for a,b in zip(rows,reference):
            assert a[:4]==b[:4]
            for i in [4,5,7,9,11,12,14,16,17]:
                assert abs(a[i]-b[i])<=2e-5*max(.001,abs(a[i]),abs(b[i])),('CPU parity',a[:4],i,a[i],b[i])
        print('PASS CPU reference parity')


def run_native_build(args):
    revision = args.revision
    if args.make_optdir and not revision:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                           cwd=args.source_root, text=True).strip()
    print(build_native(args.source_root, args.work_dir, args.build_dir,
                       args.make_optdir, revision))


def run_taylor(args):
    """Exercise the actual SW4 forward solve, material adjoint and objective."""
    import os
    import struct
    from receiver_template_startup import launch

    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    sw4 = args.sw4.resolve(strict=True)
    mopt = args.sw4mopt.resolve(strict=True)
    deck = """grid h=200 x=6000 y=6000 z=4000
time t=2.8 utcstart=01/01/2100:00:00:00.000
supergrid gp=6
developer cfl=0.7
block vp=6000 vs=3464 rho=2700
source x=2600 y=2800 z=800 Mzz=1 m0=1e16 t0=0.75 type=Ricker freq=2
"""
    station = "x=4000 y=3400 z=0"

    def execute_case(name, binary, inputs):
        root = work / name
        root.mkdir()
        input_path = root / "case.in"
        input_path.write_text(inputs)
        command = [args.launcher, "--exact", "--nodes=1", f"--ntasks={args.tasks}",
                   "--cpus-per-task=1", "--cpu-bind=cores", str(binary), str(input_path)]
        code, _ = launch(command, root, dict(os.environ, OMP_NUM_THREADS="1"), timeout=180)
        log = (root / "sw4.log").read_text()
        assert code == 0, log[-8000:]
        return root, log

    forward, _ = execute_case(
        "observed", sw4, deck + "fileio path=.\n" +
        f"rec {station} file=arrival usgsformat=1 sacformat=0 nsew=0\n")
    data = (forward / "arrival.txt").read_text().splitlines()
    assert len(data) > 20
    (forward / "zero.txt").write_text("\n".join(
        line if line.startswith("#") or not line.strip()
        else line.split()[0] + " 0 0 0" for line in data) + "\n")
    base = deck + f"""fileio path=run obspath={forward} temppath=tmp pfs=0
mparcart nx=1 ny=1 nz=1 type=vponly init=0
mscalefactors vp=6000
observation {station} file=arrival shift=0.08
"""
    invalid = "observation x=3200 y=4200 z=0 file=zero\n"

    def measure(label, receivers, excluded):
        run = base + receivers
        grad, grad_log = execute_case(
            label + "-gradient", mopt, run +
            "mrun task=gradtest misfit=crosscorr tsoutput=off savesteps=250\n")
        # The shipped gradient test independently computes a centered difference.
        text = (grad / "run/GradientTest.txt").read_text()
        rows = [[float(x) for x in line.split()] for line in text.splitlines() if line.strip()]
        assert len(rows) == 1 and len(rows[0]) >= 5, text
        index, _, _, adjoint, centered = rows[0][:5]
        assert index == 0 and math.isfinite(adjoint) and abs(adjoint) > 1e-10, rows
        assert abs(adjoint - centered) <= 2e-3 * abs(adjoint), rows

        assert grad_log.count(f"Travel-time window selection: 1 accepted, {excluded} excluded") == 1
        assert "Model lost " not in grad_log
        base_score = re.search(r"Unperturbed objective function f=(\S+)", grad_log)
        assert base_score
        return float(base_score[1]), rows[0]

    good = measure("resolved-only", "", 0)
    mixed = measure("mixed", invalid, 1)
    assert good == mixed, (good, mixed)
    base_score, (_, _, _, adjoint, centered, *_) = mixed

    # mfsurf's legacy var=rho maps to index 0, the sole Vp coefficient here.
    curve, curve_log = execute_case(
        "mixed-curve", mopt, base + invalid +
        "mrun task=func1d misfit=crosscorr tsoutput=off\n"
        "mfsurf i=1 j=1 k=1 var=rho pmin=0 pmax=40 npts=5\n")
    assert curve_log.count("Travel-time window selection: 1 accepted, 1 excluded") == 1
    assert "Model lost " not in curve_log
    payload = (curve / "run/fsurf.bin").read_bytes()
    ndim, count, lo, hi = struct.unpack_from("=iidd", payload)
    assert (ndim, count, lo, hi) == (1, 5, 0, 40)
    scores = struct.unpack_from(f"={count}d", payload, struct.calcsize("=iidd"))
    assert all(math.isfinite(value) for value in scores) and scores[0] > 0
    assert math.isclose(scores[0], base_score, rel_tol=1e-6)
    steps = (40, 20, 10)
    remainder = [abs(scores[h // 10] - scores[0] - h * adjoint) for h in steps]
    orders = [math.log(a / b, 2) for a, b in zip(remainder, remainder[1:])]
    assert all(1.7 <= order <= 2.3 for order in orders), (steps, remainder, orders)
    mixed = dict(adjoint=adjoint, centered=centered, scores=scores, steps=steps,
                 remainder=remainder, orders=orders)
    (work / "taylor-results.json").write_text(json.dumps(mixed, indent=2) + "\n")
    print("PASS model-level Taylor test and zero contribution from unresolved window")
    print("Taylor orders:", *mixed["orders"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    source = argparse.ArgumentParser(add_help=False)
    source.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    extraction = argparse.ArgumentParser(add_help=False, parents=[source])
    extraction.add_argument("--cxx", default="c++")
    extraction.add_argument("--raja", action="store_true", help="Select unified GPU implementation")

    windows = modes.add_parser("windows", parents=[extraction],
                               help="Check single and two-window weights and mode reset")
    windows.set_defaults(check=run_windows)

    numerics = modes.add_parser("numerics", parents=[extraction],
                                help="Check objectives, derivatives, and receiver data adjoints")
    numerics.add_argument("--work-dir", type=Path, required=True)
    numerics.add_argument("--cpu-reference", type=Path)
    numerics.add_argument("--asan", action="store_true",
                          help="Enable address and float-to-int overflow sanitizers")
    numerics.set_defaults(check=run_numerics)

    native = modes.add_parser("native-build", parents=[source],
                              help="Link a driver against existing SW4 objects")
    native.add_argument("--work-dir", type=Path, required=True)
    build = native.add_mutually_exclusive_group(required=True)
    build.add_argument("--build-dir", type=Path, help="Existing CMake Makefiles build directory")
    build.add_argument("--make-optdir", help="Existing standalone RAJA Make output directory")
    native.add_argument("--revision", help="Version for an archived RAJA source tree")
    native.set_defaults(check=run_native_build)


    taylor = modes.add_parser("taylor", help="Run a native CPU sw4mopt material-gradient Taylor test")
    taylor.add_argument("--work-dir", type=Path, required=True)
    taylor.add_argument("--sw4", type=Path, required=True, help="HDF5-enabled native CPU forward solver")
    taylor.add_argument("--sw4mopt", type=Path, required=True, help="Matching native CPU inversion solver")
    taylor.add_argument("--launcher", default="srun")
    taylor.add_argument("--tasks", type=int, default=2)
    taylor.set_defaults(check=run_taylor)

    args = parser.parse_args()
    args.check(args)


if __name__ == "__main__":
    main()
