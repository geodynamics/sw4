// Test-only native receiver-objective checks; included by a generated SW4 main.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

static double objective_wave(double t, int q) {
  double x = (t - 3.1) / .45;
  return (1 + .2 * q) * std::exp(-x*x) *
         (std::cos((4+q)*x) + .13*std::sin(2*x));
}
static void objective_fill(TimeSeries& ts, int n, double start, double dt,
                           double perturb, bool observed, bool zero = false) {
  ts.allocateRecordingArrays(n-1, start, dt);
  std::vector<float_sw4> sample(3);
  for (int i=0; i<n; ++i) {
    for (int q=0; q<3; ++q) {
      double value = observed ? objective_wave(start+i*dt,q) :
          1.13*objective_wave(start+i*dt-.137,q) + .01*std::sin(.1*i+q);
      sample[q] = zero ? 0 : value + perturb*std::sin(.13*i+.7*q);
    }
    ts.recordData(sample);
  }
}
static void receiver_objective_checks(EW& ew) {
  int errors=0;
  for (int storage=0; storage<2; ++storage) {
    for (int mode=0; mode<10; ++mode) {
      const int n=301, no=801;
      double start = mode==6 ? .005 : mode==7 ? 5.505 : 1.5;
      TimeSeries u(&ew,"test-u","test-u",TimeSeries::Velocity,false,!storage,
                   storage,"test.h5",1800,1800,0,false,0,1,true,0);
      TimeSeries o(&ew,"test-o","test-o",TimeSeries::Velocity,false,!storage,
                   storage,"test.h5",1800,1800,0,false,0,1,true,0);
      TimeSeries d(&ew,"test-d","test-d",TimeSeries::Velocity,false,!storage,
                   storage,"test.h5",1800,1800,0,false,0,1,true,0);
      if (mode==1 || mode==4) u.set_window(2.2,3.9);
      if (mode==2 || mode==5 || mode==9) u.set_window(2.1,2.9,3.2,4.1);
      if (mode==3) u.set_window(2,3.5,3,4.2);
      if (mode==4 || mode==5) u.exclude_component(false,true,true);
      if (mode==6 || mode==7) u.set_window(-1,10);
      objective_fill(u,n,start,.01,0,false);
      objective_fill(o,no,0,mode==9?.0125:.01,0,true,mode==8);
      objective_fill(d,n,start,.01,0,true,true);
      if (!u.myPoint()) continue;
      int before=errors;
      auto close=[&](double a,double b,double tolerance,const char* name) {
        if (!std::isfinite(a) || !std::isfinite(b) ||
            std::abs(a-b)>tolerance*std::max(.001,std::max(std::abs(a),std::abs(b)))) {
          std::cout<<"NUMERIC_FAIL "<<storage<<" "<<mode<<" "<<name
                   <<" actual="<<a<<" expected="<<b<<std::endl; ++errors;
        }
      };
      double g,h,h1,ag,ah,ah1,bg,bh,bh1;
      double J=u.misfit(o,&d,g,h,h1),old=o.get_shift(),e=1e-5;
      o.set_shift(old+e);double jp=u.misfit(o,nullptr,ag,ah,ah1);
      o.set_shift(old-e);double jm=u.misfit(o,nullptr,bg,bh,bh1);
      o.set_shift(old);
      close(g,(jp-jm)/(2*e),1e-5,"L2 shift derivative");
      close(h,(ag-bg)/(2*e),2e-5,"L2 shift curvature");
      double dot=0;
      for(int q=0;q<3;++q)for(int i=0;i<n;++i)
        dot+=d.getRecordingArray()[q][i]*std::sin(.13*i+.7*q);
      double step=storage?1e-3:1e-5;
      objective_fill(u,n,start,.01,step,false);jp=u.misfit(o,nullptr,ag,ah,ah1);
      objective_fill(u,n,start,.01,-step,false);jm=u.misfit(o,nullptr,bg,bh,bh1);
      objective_fill(u,n,start,.01,0,false);
      close(dot,(jp-jm)/(2*step),storage?5e-3:2e-5,"L2 data-adjoint dot product");
      double tt=0,ttdot=0;
      // Travel-time derivative checks require an interior, isolated maximum.
      if(mode!=6 && mode!=7) {
        tt=u.misfit2(o,&d);
        for(int q=0;q<3;++q)for(int i=0;i<n;++i)
          ttdot+=d.getRecordingArray()[q][i]*std::sin(.13*i+.7*q);
        objective_fill(u,n,start,.01,step,false);jp=u.misfit2(o,nullptr);
        objective_fill(u,n,start,.01,-step,false);jm=u.misfit2(o,nullptr);
        objective_fill(u,n,start,.01,0,false);
        close(ttdot,(jp-jm)/(2*step),storage?5e-3:2e-5,"travel-time data-adjoint dot product");
      }
      std::cout<<std::setprecision(17)<<"NATIVE_OBJECTIVE "<<storage<<" "<<mode
               <<" "<<J<<" "<<g<<" "<<h<<" "<<dot<<" "<<tt<<" "<<ttdot
               <<" "<<(errors-before)<<std::endl;
    }
  }
  int total=0;
  MPI_Allreduce(&errors,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(total) MPI_Abort(MPI_COMM_WORLD,91);
  int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0)std::cout<<"PASS native receiver objectives and data-adjoint checks"<<std::endl;
}
