// Polynomial oracle for all Cartesian anisotropic rows, including both closures.
#include "EW.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

void innerloopanisgstrvc_ci(int,int,int,int,int,int,int,float_sw4*,float_sw4*,
                          float_sw4*,int*,float_sw4*,float_sw4*,float_sw4*,
                          float_sw4,float_sw4*,float_sw4*,float_sw4*);

int main(int argc,char**argv)
{
   MPI_Init(&argc,&argv);
   try {
      if(argc!=2) throw std::runtime_error("Usage: sw4_anisotropic_check case.in");
      std::vector<std::vector<Source*>> sources;
      std::vector<std::vector<TimeSeries*>> receivers;
      EW ew(argv[1],sources,receivers);
      float_sw4 acof[384],ghcof[6],bop[24],bope[48],sbop[6];
      ew.GetStencilCoefficients(acof,ghcof,bop,bope,sbop);
      constexpr int first=-2,last=26,n=last-first+1,nk=24,cells=n*n*n;
      constexpr double h=.5,lambda=3,mu=2;
      auto index=[](int m,int i,int j,int k) {
         return (m-1)*cells+(i-first)+n*(j-first)+n*n*(k-first);
      };
      std::vector<float_sw4> u(3*cells),lu(3*cells),c(21*cells),stretch(n,1);
      // Voigt order is xx, yy, zz, yz, xz, xy; SW4 packs these pairs.
      const int packed[21][2]={{0,0},{0,5},{0,4},{0,1},{0,3},{0,2},
         {5,5},{5,4},{5,1},{5,3},{5,2},{4,4},{4,1},{4,3},{4,2},
         {1,1},{1,3},{1,2},{3,3},{3,2},{2,2}};
      const int voigt[3][3]={{0,5,4},{5,1,3},{4,3,2}};
      int sides[6]={0,0,0,0,1,1};
      double error=0;
      for(int material=0;material<2;++material) {
         double tensor[6][6];
         for(int a=0;a<6;++a) for(int b=0;b<6;++b) {
            if(material==0)
               tensor[a][b]=a<3 && b<3 ? lambda+(a==b ? 2*mu : 0) : (a==b ? mu : 0);
            else {
               // Symmetric, strictly diagonally dominant anisotropic stiffness.
               const double diagonal[6]={10,12,14,3,4,5};
               tensor[a][b]=a==b ? diagonal[a] : .02*(a+b+1);
            }
         }
         for(int m=1;m<=21;++m)
            std::fill(c.begin()+(m-1)*cells,c.begin()+m*cells,
                      tensor[packed[m-1][0]][packed[m-1][1]]);
         for(int component=1;component<=3;++component)
            for(int axis=0;axis<3;++axis) for(int second=axis;second<3;++second) {
               std::fill(u.begin(),u.end(),0);std::fill(lu.begin(),lu.end(),0);
               for(int k=first;k<=last;++k) for(int j=first;j<=last;++j) for(int i=first;i<=last;++i) {
                  const double xyz[3]={(i-1)*h,(j-1)*h,(k-1)*h};
                  u[index(component,i,j,k)]=xyz[axis]*xyz[second];
               }
               innerloopanisgstrvc_ci(first,last,first,last,first,last,nk,u.data(),lu.data(),
                                     c.data(),sides,acof,bope,ghcof,h,stretch.data(),stretch.data(),stretch.data());
               // L_m = C_(m,a,l,b) d_a d_b u_l, derived independently of the stencil.
               for(int k=1;k<=nk;++k) for(int m=1;m<=3;++m) {
                  const double expected=tensor[voigt[m-1][axis]][voigt[component-1][second]]
                                       +tensor[voigt[m-1][second]][voigt[component-1][axis]];
                  const double actual=lu[index(m,12,13,k)];
                  if(!std::isfinite(actual)) throw std::runtime_error("Nonfinite anisotropic residual");
                  error=std::max(error,std::abs(actual-expected));
               }
            }
      }
      if(error>2e-10) throw std::runtime_error("Anisotropic polynomial residual error "+std::to_string(error));
      if(!ew.getRank()) std::cout<<"PASS: all anisotropic interior/top/bottom rows, max_error="<<error<<"\n";
   } catch(const std::exception& e) {
      std::cerr<<e.what()<<"\n";MPI_Abort(MPI_COMM_WORLD,1);
   }
   MPI_Finalize();
}
