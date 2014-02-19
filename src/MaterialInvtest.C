#include "EW.h"
#include "F77_FUNC.h"
#include "MaterialInvtest.h"

extern "C" {
   void F77_FUNC(invtestmtrl,INVTESTMTRL)(  int*, int*, int*, int*, int*, int*, double*, double*,
				       double*, double*, double*, int* );
   void F77_FUNC(invtestmtrlc,INVTESTMTRLC)( int*, int*, int*, int*, int*, int*, double*, double*,
					double*, double*, double*, double*, int* );
}

//-----------------------------------------------------------------------
MaterialInvtest::MaterialInvtest( EW* a_ew, int nr ) :
   mEW(a_ew),
   m_nr(nr)
{
   mCoversAllPoints = true;
}

//-----------------------------------------------------------------------
void MaterialInvtest::set_material_properties( std::vector<Sarray> & rho, std::vector<Sarray> & cs,
					       std::vector<Sarray> & cp, std::vector<Sarray>& xis,
					       std::vector<Sarray>& xip)
{
   int ifirst, ilast, jfirst, jlast, kfirst, klast;
   double* rho_ptr, *cs_ptr, *cp_ptr;

   for( int g = 0 ; g < mEW->mNumberOfCartesianGrids; g++) 
   {
      rho_ptr   = rho[g].c_ptr();
      cs_ptr    = cs[g].c_ptr();
      cp_ptr    = cp[g].c_ptr();
      ifirst = mEW->m_iStart[g];
      ilast  = mEW->m_iEnd[g];
      jfirst = mEW->m_jStart[g];
      jlast  = mEW->m_jEnd[g];
      kfirst = mEW->m_kStart[g];
      klast  = mEW->m_kEnd[g];
      double h = mEW->mGridSize[g]; 
      double zmin = mEW->m_zmin[g];
      F77_FUNC(invtestmtrl,INVTESTMTRL)( &ifirst, &ilast, &jfirst, &jlast,
					 &kfirst, &klast, rho_ptr, cs_ptr,
					 cp_ptr, &h, &zmin, &m_nr);
   }
   if( mEW->topographyExists() )
   {
      int g = mEW->mNumberOfGrids-1;
      rho_ptr   = rho[g].c_ptr();
      cs_ptr    = cs[g].c_ptr();
      cp_ptr    = cp[g].c_ptr();
      ifirst = mEW->m_iStart[g];
      ilast  = mEW->m_iEnd[g];
      jfirst = mEW->m_jStart[g];
      jlast  = mEW->m_jEnd[g];
      kfirst = mEW->m_kStart[g];
      klast  = mEW->m_kEnd[g];
      double* x_ptr = mEW->mX.c_ptr();
      double* y_ptr = mEW->mY.c_ptr();
      double* z_ptr = mEW->mZ.c_ptr();
      F77_FUNC(invtestmtrlc,INVTESTMTRLC)( &ifirst, &ilast, &jfirst, &jlast,
					   &kfirst, &klast, rho_ptr, cs_ptr,
					   cp_ptr, x_ptr, y_ptr, z_ptr, &m_nr);
   }
}
