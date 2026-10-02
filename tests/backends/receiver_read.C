// On-demand timing round trip through the real SAC, USGS and rechdf5 readers.
#include "EW.h"
#include "TimeSeries.h"
#include "ReceiverComponents.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
   MPI_Init(&argc,&argv);
   try {
      if(argc!=5 && argc!=6) throw std::runtime_error("Usage: sw4_receiver_check case.in output-directory quantity nsew [downsample]");
      std::vector<std::vector<Source*>> sources;
      std::vector<std::vector<TimeSeries*>> receivers;
      EW ew(argv[1],sources,receivers);
      ew.setupRun(sources);
      if(!ew.isInitialized()) throw std::runtime_error("Setup failed");
      const std::string quantity=argv[3];
      const std::vector<std::pair<std::string,TimeSeries::receiverMode>> modes={
         {"displacement",TimeSeries::Displacement},{"velocity",TimeSeries::Velocity},
         {"div",TimeSeries::Div},{"curl",TimeSeries::Curl},{"strains",TimeSeries::Strains},
         {"displacementgradient",TimeSeries::DisplacementGradient}};
      auto found=std::find_if(modes.begin(),modes.end(),[&](const std::pair<std::string,TimeSeries::receiverMode>& item){return item.first==quantity;});
      if(found==modes.end()) throw std::runtime_error("Unknown quantity");
      const auto mode=found->second;
      const bool geographic=std::string(argv[4])=="1";
      const int downsample=argc==6 ? std::stoi(argv[5]):1;
      const auto components=sw4::receiver_components(mode,!geographic);
      const std::string prefix=std::string(argv[2])+"/ascii";
      TimeSeries sac(&ew,"sac","station",mode,true,true,false,"",
                     1800,1800,200,false,0,1,!geographic);
      TimeSeries text(&ew,prefix,"station",mode,false,true,false,"",
                      1800,1800,200,false,0,1,!geographic);
      TimeSeries hdf(&ew,"hdf","station",mode,false,true,true,"",
                     1800,1800,200,false,0,downsample,!geographic);
      std::vector<std::string> files;
      for(const auto& component:components) files.push_back(prefix+"."+component.suffix);
      if(!sac.readSACcomponents(&ew,files,false)) throw std::runtime_error("SAC receiver read failed");
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
            if(series->getLastTimeStep()!=(series==&hdf ? (text.getLastTimeStep()/downsample)*downsample:text.getLastTimeStep()))
               throw std::runtime_error("Reader sample counts differ");
            const double start=series->getStartTime()+series->getTimeShift();
            local_error=std::max(local_error,std::abs(start));
            local_error=std::max(local_error,std::abs(series->getDt()-text.getDt()));
         }
         for(size_t component=0;component<components.size();++component) {
            double peak=0;
            for(int i=0;i<=text.getLastTimeStep();++i)
               peak=std::max(peak,std::abs(text.getRecordingArray()[component][i]));
            const double tolerance=1e-12+2e-6*peak;
            for(int i=0;i<=text.getLastTimeStep();++i) {
               const double reference=text.getRecordingArray()[component][i];
               if(!std::isfinite(reference) ||
                  !std::isfinite(sac.getRecordingArray()[component][i]) ||
                  (i<=hdf.getLastTimeStep() && i%downsample==0 && !std::isfinite(hdf.getRecordingArray()[component][i])) ||
                  std::abs(sac.getRecordingArray()[component][i]-reference)>tolerance ||
                  (i<=hdf.getLastTimeStep() && i%downsample==0 && std::abs(hdf.getRecordingArray()[component][i]-reference)>tolerance))
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
