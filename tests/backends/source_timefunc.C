// On-demand checks of the real source constructor, preparation, copying and forcing.
#include "EW.h"
#include "Source.h"
#include "GridPointSource.h"
#include "Filter.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

static std::vector<double> forcing(EW& ew, Source& source, double time, bool moment)
{
   std::vector<GridPointSource*> points;
   source.set_grid_point_sources4(&ew, points);
   double local[9]={}, total[9]={};
   for (auto point: points)
   {
      double f[3]; point->getFxyz(time,f);
      const double h=ew.mGridSize[point->m_grid];
      const double offset[3]={(point->m_i0-1)*h-source.getX0(),
                             (point->m_j0-1)*h-source.getY0(),
                             (point->m_k0-1)*h+ew.m_zmin[point->m_grid]-source.getZ0()};
      for (int d=0;d<(moment ? 3 : 1);++d)
         for (int c=0;c<3;++c)
            local[3*d+c]+=f[c]*h*h*h*(moment ? offset[d] : 1.0);
      delete point;
   }
   MPI_Allreduce(local,total,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   return std::vector<double>(total,total+(moment ? 9 : 3));
}

static void check(EW& ew, int channels, bool filtered)
{
   const int n=64;
   const double dt=.01, start=.25;
   const timeDep kind=channels==6 ? iDiscrete6moments : channels==3 ? iDiscrete3forces : iDiscrete;
   std::vector<double> raw(channels*(n+1)), expected;
   for(int c=0;c<channels;++c)
   {
      raw[c*(n+1)]=start;
      for(int i=0;i<n;++i)
         raw[c*(n+1)+i+1]=(c+1)*std::pow(std::sin(std::acos(-1.)*i/(n-1)),2)
                           +.1*c*std::sin(2*std::acos(-1.)*i/(n-1));
   }
   int samples=n;
   std::unique_ptr<Source> source;
   if(channels==6)
      source.reset(new Source(&ew,1/dt,start,1800,1800,1000,1,0,0,1,0,1,
                              kind,"moments",false,1,raw.data(),raw.size(),&samples,1));
   else
      source.reset(new Source(&ew,1/dt,start,1800,1800,1000,1,2,-.5,
                              kind,"force",false,1,raw.data(),raw.size(),&samples,1));
   // Copies of raw histories and prepared histories must preserve their representation.
   std::unique_ptr<Source> before(source->copy("before"));
   Filter filter(lowPass,4,2,0,12);
   filter.computeSOS(dt);
   const int padding=filtered ? static_cast<int>(std::ceil(filter.estimatePrecursor()/dt)) : 0;
   const int extent=n+2*padding;
   for(int c=0;c<channels;++c)
   {
      std::vector<double> signal(extent,0);
      std::copy(raw.begin()+c*(n+1)+1,raw.begin()+(c+1)*(n+1),signal.begin()+padding);
      if(filtered) filter.evaluate(extent,signal.data(),signal.data());
      expected.insert(expected.end(),signal.begin(),signal.end());
   }
   source->prepareTimeFunc(filtered,ew.getTimeStep(),ew.getNumberOfSteps(),&filter);
   before->prepareTimeFunc(filtered,ew.getTimeStep(),ew.getNumberOfSteps(),&filter);
   std::unique_ptr<Source> after(source->copy("after"));
   source->prepareTimeFunc(filtered,ew.getTimeStep(),ew.getNumberOfSteps(),&filter);
   after->prepareTimeFunc(filtered,ew.getTimeStep(),ew.getNumberOfSteps(),&filter);
   double maxerror=0;
   for(int i=0;i<extent;++i)
   {
      const double time=start-padding*dt+i*dt;
      for(Source* item: {source.get(),before.get(),after.get()})
      {
         const auto values=forcing(ew,*item,time,channels==6);
         const int history[]={0,1,2,1,3,4,2,4,5};
         for(size_t c=0;c<values.size();++c)
         {
            const int component=channels==6 ? history[c] : c;
            const double oracle=channels==1 ? expected[i]*(c==0 ? 1 : c==1 ? 2 : -.5) : expected[component*extent+i];
            maxerror=std::max(maxerror,std::abs(values[c]-oracle));
         }
      }
   }
   if(maxerror>2e-10) throw std::runtime_error("Source history error "+std::to_string(maxerror));
   if(!ew.getRank()) std::cout<<"PASS: "<<channels<<" histories, filtered="<<filtered<<" max_error="<<maxerror<<"\n";
}

int main(int argc,char**argv)
{
   MPI_Init(&argc,&argv);
   try {
      if(argc!=2) throw std::runtime_error("Usage: sw4_source_check case.in");
      std::vector<std::vector<Source*>> sources;
      std::vector<std::vector<TimeSeries*>> receivers;
      EW ew(argv[1],sources,receivers);
      ew.setupRun(sources);
      if(!ew.isInitialized()) throw std::runtime_error("Setup failed");
      for(int channels: {1,3,6}) for(bool filtered: {false,true}) check(ew,channels,filtered);
   } catch(const std::exception& e) {
      std::cerr<<e.what()<<"\n"; MPI_Abort(MPI_COMM_WORLD,1);
   }
   MPI_Finalize();
}
