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
      if(argc!=3) throw std::runtime_error("Usage: sw4_receiver_check case.in output-directory");
      std::vector<std::vector<Source*>> sources;
      std::vector<std::vector<TimeSeries*>> receivers;
      EW ew(argv[1],sources,receivers);
      ew.setupRun(sources);
      if(!ew.isInitialized()) throw std::runtime_error("Setup failed");
      const std::string prefix=std::string(argv[2])+"/ascii";
      TimeSeries sac(&ew,"sac","station",TimeSeries::Displacement,true,true,false,"",
                     1800,1800,200,false,0,1,true);
      TimeSeries text(&ew,prefix,"station",TimeSeries::Displacement,false,true,false,"",
                      1800,1800,200,false,0,1,true);
      TimeSeries hdf(&ew,"hdf","station",TimeSeries::Displacement,false,true,true,"",
                     1800,1800,200,false,0,1,true);
      sac.readSACfiles(&ew,(prefix+".x").c_str(),(prefix+".y").c_str(),(prefix+".z").c_str(),false);
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
         // The vertical component avoids any unrelated coordinate-basis rotation.
         for(int i=0;i<=text.getLastTimeStep();++i) {
            const double reference=text.getRecordingArray()[2][i];
            const double tolerance=1e-12+2e-6*std::abs(reference);
            if(std::abs(sac.getRecordingArray()[2][i]-reference)>tolerance ||
               std::abs(hdf.getRecordingArray()[2][i]-reference)>tolerance)
               throw std::runtime_error("Reader vertical samples differ");
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
