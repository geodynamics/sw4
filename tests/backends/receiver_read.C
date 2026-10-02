// On-demand timing round trip through the real SAC, USGS and rechdf5 readers.
#include "EW.h"
#include "TimeSeries.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
   MPI_Init(&argc,&argv);
   try {
      if(argc!=5) throw std::runtime_error("Usage: sw4_receiver_check case.in output-directory displacement|velocity nsew");
      std::vector<std::vector<Source*>> sources;
      std::vector<std::vector<TimeSeries*>> receivers;
      EW ew(argv[1],sources,receivers);
      ew.setupRun(sources);
      if(!ew.isInitialized()) throw std::runtime_error("Setup failed");
      const bool velocity=std::string(argv[3])=="velocity";
      const bool geographic=std::string(argv[4])=="1";
      const auto mode=velocity ? TimeSeries::Velocity : TimeSeries::Displacement;
      const std::string suffix=velocity ? "v" : "";
      const std::string prefix=std::string(argv[2])+"/ascii";
      TimeSeries sac(&ew,"sac","station",mode,true,true,false,"",
                     1800,1800,200,false,0,1,!geographic);
      TimeSeries text(&ew,prefix,"station",mode,false,true,false,"",
                      1800,1800,200,false,0,1,!geographic);
      TimeSeries hdf(&ew,"hdf","station",mode,false,true,true,"",
                     1800,1800,200,false,0,1,!geographic);
      sac.readSACfiles(&ew,(prefix+(geographic ? ".e" : ".x")+suffix).c_str(),
                      (prefix+(geographic ? ".n" : ".y")+suffix).c_str(),
                      (prefix+(geographic ? ".u" : ".z")+suffix).c_str(),false);
      text.readFile(&ew,false);
#ifdef USE_HDF5
      if(!hdf.readSACHDF5(&ew,std::string(argv[2])+"/receivers.h5",false))
         throw std::runtime_error("Receiver HDF5 read failed");
#else
      throw std::runtime_error("This check requires USE_HDF5");
#endif
      double local_error=0, global_error=0;
      if(text.myPoint()) {
         if(text.getLastTimeStep()<1) throw std::runtime_error("Reader did not load samples");
         for(TimeSeries* series: {&sac,&text,&hdf}) {
            if(series->getLastTimeStep()!=text.getLastTimeStep())
               throw std::runtime_error("Reader sample counts differ");
            const double start=series->getStartTime()+series->getTimeShift();
            local_error=std::max(local_error,std::abs(start));
            local_error=std::max(local_error,std::abs(series->getDt()-text.getDt()));
         }
         for(int component=0;component<3;++component) {
            double peak=0;
            for(int i=0;i<=text.getLastTimeStep();++i)
               peak=std::max(peak,std::abs(text.getRecordingArray()[component][i]));
            const double tolerance=1e-12+2e-6*peak;
            for(int i=0;i<=text.getLastTimeStep();++i) {
               const double reference=text.getRecordingArray()[component][i];
               if(!std::isfinite(reference) ||
                  !std::isfinite(sac.getRecordingArray()[component][i]) ||
                  !std::isfinite(hdf.getRecordingArray()[component][i]) ||
                  std::abs(sac.getRecordingArray()[component][i]-reference)>tolerance ||
                  std::abs(hdf.getRecordingArray()[component][i]-reference)>tolerance)
                  throw std::runtime_error("Reader component "+std::to_string(component)+
                     " differs: SAC="+std::to_string(std::abs(sac.getRecordingArray()[component][i]-reference))+
                     " HDF5="+std::to_string(std::abs(hdf.getRecordingArray()[component][i]-reference))+
                     " reference="+std::to_string(reference));
            }
         }
      }
      MPI_Allreduce(&local_error,&global_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      if(global_error>1e-8) throw std::runtime_error("Reader timing differs: "+std::to_string(global_error));
      if(!ew.getRank()) std::cout<<"PASS: SAC/USGS/rechdf5 timing round trip, max_error="<<global_error<<"\n";
   } catch(const std::exception& e) {
      std::cerr<<e.what()<<"\n"; MPI_Abort(MPI_COMM_WORLD,1);
   }
   MPI_Finalize();
}
