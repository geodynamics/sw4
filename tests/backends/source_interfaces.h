// Independent integral/first-moment oracle for actual Cartesian refinement grids.
#include "Source.h"
#include "GridPointSource.h"
#include "Filter.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

static int sw4_check_source_interfaces(EW& ew)
{
   const double norm[]={17./48,59./48,43./48,49./48};
   double max_error=0;
   int checks=0;
   for(int fine=1;fine<ew.mNumberOfCartesianGrids;++fine) {
      int coarse=fine-1;
      const double interface=ew.m_zmin[coarse];
      const double fine_h=ew.mGridSize[fine],coarse_h=ew.mGridSize[coarse];
      const int fine_nz=ew.m_global_nz[fine];
      const double source_x=(ew.m_global_nx[fine]-1)*fine_h/2-3;
      const double source_y=(ew.m_global_ny[fine]-1)*fine_h/2+3;
      REQUIRE2(std::abs(interface-(ew.m_zmin[fine]+(fine_nz-1)*fine_h))<1e-9,
               "Interface oracle requires adjacent Cartesian grids");
      // Generate depths from the grids, including all six special stencils,
      // remote controls, and exact/epsilon-sided interface positions.
      std::vector<double> depths;
      for(int kc: {fine_nz-3,fine_nz-2,fine_nz-1,fine_nz-6})
         depths.push_back(ew.m_zmin[fine]+(kc-1+.46)*fine_h);
      for(int kc: {1,2,3,6}) depths.push_back(interface+(kc-1+.47)*coarse_h);
      for(double epsilon: {-1e-7,0.,1e-7}) depths.push_back(interface+epsilon);
      for(size_t position=0;position<depths.size();++position) {
         const double z=depths[position];
         int i,j,k,g;
         int local=-1;
         if(ew.computeNearestGridPoint2(i,j,k,g,source_x,source_y,z)) local=g;
         MPI_Allreduce(&local,&g,1,MPI_INT,MPI_MAX,ew.m_1d_communicator);
         REQUIRE2(g>=0,"Source grid not found");
         int kc=std::floor((z-ew.m_zmin[g])/ew.mGridSize[g]+1);
         kc=std::max(1,std::min(kc,ew.m_global_nz[g]-1));
         if(position<4) REQUIRE2(g==fine && kc==fine_nz-(position<3 ? 3-position:6),
                                "Missed fine interface source stencil");
         if(position>=4 && position<8) REQUIRE2(g==coarse && kc==(position<7 ? position-3:6),
                                               "Missed coarse interface source stencil");
         for(bool moment: {false,true}) {
            const int channels=moment ? 6:3,n=64;
            const double dt=.01,start=.25;
            std::vector<double> raw(channels*(n+1));
            for(int c=0;c<channels;++c) {
               raw[c*(n+1)]=start;
               for(int t=0;t<n;++t) raw[c*(n+1)+1+t]=c+1;
            }
            int samples=n;
            std::unique_ptr<Source> source;
            if(moment) source.reset(new Source(&ew,1/dt,start,source_x,source_y,z,1,0,0,1,0,1,
                 iDiscrete6moments,"interface",false,1,raw.data(),raw.size(),&samples,1));
            else source.reset(new Source(&ew,1/dt,start,source_x,source_y,z,1,1,1,
                 iDiscrete3forces,"interface",false,1,raw.data(),raw.size(),&samples,1));
            Filter filter(lowPass,4,2,0,12); filter.computeSOS(dt);
            source->prepareTimeFunc(false,ew.getTimeStep(),ew.getNumberOfSteps(),&filter);
            std::vector<GridPointSource*> points;
            source->set_grid_point_sources4(&ew,points);
            int ranks, local_owner=points.empty() ? 0:1, owners=0;
            MPI_Comm_size(ew.m_1d_communicator,&ranks);
            MPI_Allreduce(&local_owner,&owners,1,MPI_INT,MPI_SUM,ew.m_1d_communicator);
            REQUIRE2((ranks!=2 && ranks!=4) || owners>1,"Source oracle did not cross MPI ownership");
            double local_integrals[12]={},total[12]={};
            for(auto point:points) {
               const int pg=point->m_grid,pk=point->m_k0,nz=ew.m_global_nz[pg];
               const double h=ew.mGridSize[pg];
               double weight=h*h*h;
               if(pk<=4) weight*=norm[pk-1];
               else if(pk>=nz-3) weight*=norm[nz-pk];
               double force[3]; point->getFxyz(start+.2,force);
               const double offset[]={(point->m_i0-1)*h-source_x,(point->m_j0-1)*h-source_y,
                  (pk-1)*h+ew.m_zmin[pg]-z};
               for(int c=0;c<3;++c) {
                  local_integrals[c]+=force[c]*weight;
                  for(int d=0;d<3;++d) local_integrals[3+3*d+c]+=force[c]*weight*offset[d];
               }
               delete point;
            }
            MPI_Allreduce(local_integrals,total,12,MPI_DOUBLE,MPI_SUM,ew.m_1d_communicator);
            const int tensor[]={0,1,2,1,3,4,2,4,5};
            for(int c=0;c<12;++c) {
               const double expected=c<3 ? (moment ? 0.:c+1):
                  (moment ? tensor[c-3]+1.:0.);
               max_error=std::max(max_error,std::abs(total[c]-expected));
            }
            ++checks;
         }
         if(!ew.getRank()) std::cout<<"INTERFACE_CASE grid="<<g<<" kc="<<kc
            <<" Nz="<<ew.m_global_nz[g]<<" h="<<ew.mGridSize[g]<<" z="<<z<<std::endl;
      }
   }
   REQUIRE2(checks>0 && max_error<2e-9,"Interface integral/moment error "<<max_error);
   if(!ew.getRank()) std::cout<<"PASS: interface forces and nine moments, cases="<<checks
      <<", max_error="<<max_error<<std::endl;
   return 0;
}
