
#include "EW.h"
#include "MaterialVolimagefile.h"

//-----------------------------------------------------------------------
MaterialVolimagefile::MaterialVolimagefile( EW* a_ew, bool rhomula, std::string path, std::string rhofile, 
	std::string mufile, std::string lambdafile, std::string qpfile, std::string qsfile ) :
   mEW(a_ew)
{
   mCoversAllPoints = true;
   m_rhomula = rhomula;
   m_path = path;
   m_rho = rhofile;
   m_mu = mufile;
   m_lambda = lambdafile;
   m_qp = qpfile;
   m_qs = qsfile;
}

//-----------------------------------------------------------------------
void MaterialVolimagefile::set_material_properties( std::vector<Sarray> & rho, std::vector<Sarray> & cs,
				std::vector<Sarray> & cp, std::vector<Sarray>& xis, std::vector<Sarray>& xip)
{
   mEW->read_volimage( m_path, m_rho, rho );
   mEW->read_volimage( m_path, m_mu, cs );
   mEW->read_volimage( m_path, m_lambda, cp );
   //   mEW->check_min_max_int( rho );

   mEW->communicate_arrays( rho );
   mEW->communicate_arrays( cs );
   mEW->communicate_arrays( cp );
   mEW->update_curvilinear_cartesian_interface( rho );
   mEW->update_curvilinear_cartesian_interface( cs );
   mEW->update_curvilinear_cartesian_interface( cp );

   if( xis[0].is_defined() )
   {
      mEW->read_volimage( m_path, m_qs, xis );
      mEW->communicate_arrays( xis );
      mEW->update_curvilinear_cartesian_interface( xis );
   }
   if( xip[0].is_defined() )
   {
      mEW->read_volimage( m_path, m_qp, xip );
      mEW->communicate_arrays( xip );
      mEW->update_curvilinear_cartesian_interface( xip );
   }


   if( m_rhomula )
   {
      for( int g = 0 ; g < mEW->mNumberOfGrids; g++) 
      {
	 // Convert to cp,cs,rho, here cs is mu, cp is lambda
	 double* rhoptr = rho[g].c_ptr();
	 double* csptr  = cs[g].c_ptr();
	 double* cpptr  = cp[g].c_ptr();
	 for( size_t i = 0 ; i < rho[g].npts() ; i++ )
	 {
            if( rhoptr[i] != -1 )
	    {
	       csptr[i] = csptr[i]/rhoptr[i];
	       cpptr[i] = 2*csptr[i] + cpptr[i]/rhoptr[i];
	       csptr[i] = sqrt(csptr[i]);
	       cpptr[i] = sqrt(cpptr[i]);
	    }
	 }
      }
   }
}
