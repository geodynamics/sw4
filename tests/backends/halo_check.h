// Explicit test executable only: exercise the production halo transport.
#include <array>
#include <iostream>
#include <vector>

inline double sw4_halo_value(int phase, int rank, int c, int i, int j, int k)
{
   return phase*1e11 + rank*1e8 + c*1e6 + (k+100)*10000 + (j+100)*100 + i+100;
}

static int sw4_check_halos(EW& ew)
{
   int rank, size;
   MPI_Comm_rank(ew.m_cartesian_communicator, &rank);
   MPI_Comm_size(ew.m_cartesian_communicator, &size);
   int dimensions[2],periods[2],coordinates[2];
   MPI_Cart_get(ew.m_cartesian_communicator,2,dimensions,periods,coordinates);
   const int padding=ew.getNumberOfParallelPaddingPoints();
   long long errors=0;
   // GPU communication supports both component-major and interleaved layouts.
   // Native C kernels use component-major arrays exclusively.
#if defined(SW4_USE_RAJA)
   const std::vector<bool> layouts={true,false};
#else
   const std::vector<bool> layouts={true};
#endif
   for(bool component_major: layouts) {
#if defined(SW4_USE_RAJA)
      if(ew.m_croutines!=component_major) {
         ew.m_croutines=component_major;
         Sarray::m_corder=component_major;
         ew.setupMPICommunications();
      }
#endif
      for(int g=0;g<ew.mNumberOfGrids;++g) {
         const int ib=ew.m_iStart[g],ie=ew.m_iEnd[g];
         const int jb=ew.m_jStart[g],je=ew.m_jEnd[g];
         const int kb=ew.m_kStart[g],ke=ew.m_kEnd[g];
         const std::array<int,4> owned={{ib+(ew.m_neighbor[0]!=MPI_PROC_NULL ? padding:0),
            ie-(ew.m_neighbor[1]!=MPI_PROC_NULL ? padding:0),
            jb+(ew.m_neighbor[2]!=MPI_PROC_NULL ? padding:0),
            je-(ew.m_neighbor[3]!=MPI_PROC_NULL ? padding:0)}};
         std::vector<int> bounds(4*size);
         MPI_Allgather(owned.data(),4,MPI_INT,bounds.data(),4,MPI_INT,ew.m_cartesian_communicator);
         for(int nc: {1,3,4,21}) {
            Sarray field(nc,ib,ie,jb,je,kb,ke);
            for(int phase=0;phase<3;++phase) {
#if defined(SW4_USE_RAJA)
               // A pending device producer must complete before MPI packs it.
               auto data=field.c_ptr();
               const int ni=ie-ib+1,nj=je-jb+1,nk=ke-kb+1;
               RAJA::RangeSegment points(0,field.m_npts),dummy(0,1);
               RAJA::kernel<BUFFER_POL>(RAJA::make_tuple(points,dummy),
                  [=] RAJA_DEVICE(int index,int unused) {
                  int cell=component_major ? index%(ni*nj*nk):index/nc;
                  int c=component_major ? index/(ni*nj*nk)+1:index%nc+1;
                  int i=ib+cell%ni,j=jb+(cell/ni)%nj,k=kb+cell/(ni*nj);
                  data[index]=phase*1e11+rank*1e8+c*1e6+(k+100)*10000+(j+100)*100+i+100;
               });
#else
               for(int k=kb;k<=ke;++k) for(int j=jb;j<=je;++j)
                  for(int i=ib;i<=ie;++i) for(int c=1;c<=nc;++c)
                     field(c,i,j,k)=sw4_halo_value(phase,rank,c,i,j,k);
#endif
               ew.communicate_array(field,g);
               for(int j=jb;j<=je;++j) for(int i=ib;i<=ie;++i) {
                  int owner=-1;
                  for(int r=0;r<size;++r)
                     if(i>=bounds[4*r] && i<=bounds[4*r+1] &&
                        j>=bounds[4*r+2] && j<=bounds[4*r+3]) {
                        if(owner!=-1) ++errors;
                        owner=r;
                     }
                  if(owner<0) { ++errors; continue; }
                  for(int k=kb;k<=ke;++k) for(int c=1;c<=nc;++c)
                     if(field(c,i,j,k)!=sw4_halo_value(phase,owner,c,i,j,k)) ++errors;
               }
            }
         }
      }
   }
   long long total;
   MPI_Allreduce(&errors,&total,1,MPI_LONG_LONG,MPI_SUM,ew.m_cartesian_communicator);
   if(!rank) std::cout<<(total ? "FAIL" : "PASS")
      <<": 1/3/4/21-component halo sentinels, layouts="<<layouts.size()
      <<", ranks="<<size<<", dimensions="<<dimensions[0]<<"x"<<dimensions[1]<<", errors="<<total<<std::endl;
   return total ? 1:0;
}
