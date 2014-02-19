#include "ConvParOutput.h"

//-----------------------------------------------------------------------
ConvParOutput::ConvParOutput( EW& ew, int n, int nvar, int myrank, int verbose,
			      bool cg ):
   m_n(n),
   m_nvar(nvar),
   m_myrank(myrank),
   m_verbose(verbose),
   m_ind(0),
   m_cg(cg)
{
   m_convfile = ew.getOutputPath() + "convergence.log";
   m_parafile = ew.getOutputPath() + "parameters.log";

   m_stepfileio = true;

   if( m_stepfileio && m_myrank == 0 )
   {
      m_fd  = fopen(m_convfile.c_str(),"w");
      m_fdx = fopen(m_parafile.c_str(),"w");
   }
   if( m_n > 11 )
      m_nsize = m_nvar + m_n - 11;
   else
      m_nsize = m_nvar;
}

//-----------------------------------------------------------------------
void ConvParOutput::print_dfmsg( double f, double* df, double* sf  )
{
   if( m_myrank == 0 )
   {
      cout << "-----------------------------------------------------------------------" << endl;
      cout << "Misfit= "  << f << endl;
      cout << " Scaled gradient = " ;
      for( int i=0 ; i < m_nvar ; i++ )
      {
	 cout << df[i]*sf[i] << " ";
	 if( i == 5 )
	    cout << endl << "      " ;
      }
	 // Observation shifts
      for( int i=11 ; i < m_n ; i++ )
      {
	 cout << "   " << df[i]*sf[i] << " " ;
	 if( i % 5 == 0 )
	    cout << endl <<  "   " ;
      }
      cout << endl;
   }
}

//-----------------------------------------------------------------------
void ConvParOutput::print_xmsg( double* x, double rnorm, double dxnorm,
				int j, int k )
{
   if( m_myrank == 0 )
   {
      cout << "-----------------------------------------------------------------------" << endl;
      if( m_cg )
	 cout << " it=" << j << "," << k << " dfnorm= " << rnorm << " dxnorm= " << dxnorm << endl;
      else
	 cout << " it=" << j << " dfnorm= " << rnorm << " dxnorm= " << dxnorm << endl;
      cout << "  x = " ;
      for( int i=0 ; i < m_nvar ; i++ )
      {
	 cout << x[i] << " ";
	 if( i==5 )
	    cout << endl << "      " ;
      }
      for( int i=11 ; i < m_n ; i++ )
      {
	 cout << x[i] << " ";
	 if( i % 5 == 0 )
	    cout << endl << "      " ;
      }
      cout << endl;
   }
}

//-----------------------------------------------------------------------
void ConvParOutput::print_vector( double* v, const char* msg, int m  )
{
   if( m_myrank == 0 && m_verbose > m )
   {
      cout << msg << endl;
      for( int i=0 ; i < m_nvar ; i++ )
      {
	 cout << v[i] << " " ;
         if( i % 5 == 0 )
	    cout << endl;
      }
      for( int i=11 ; i < m_n ; i++ )
      {
	 cout << v[i] << " " ;
         if( i % 5 == 0 )
	    cout << endl;
      }
      cout << endl;
   }
}

//-----------------------------------------------------------------------
void ConvParOutput::print_scalar( double s, const char* msg, int m  )
{
   if( m_myrank == 0 && m_verbose > m )
      cout << msg << " " << s << endl;
}

//-----------------------------------------------------------------------
void ConvParOutput::save_step( double f, double* df, double* sf,
			       double* x, double rnorm, double dxnorm,
			       int nreductions, int j, int k )
{
   if( m_myrank == 0 )
   {
      if( m_stepfileio )
      {
	 if( m_cg )
	 {
	    fprintf(m_fd, "%i %i %15.7g %15.7g %15.7g %i\n", j,k, rnorm, dxnorm, f, nreductions );
	    fprintf(m_fdx,"%i %i ",j,k);
	 }
	 else
	 {
	    fprintf(m_fd, "%i %15.7g %15.7g %15.7g %i\n", j,rnorm, dxnorm, f, nreductions );
	    fprintf(m_fdx,"%i ",j);
	 }
	 for( int i=0; i < m_nvar ;i++ )
	    fprintf(m_fdx," %15.7g ",x[i]);
	 for( int i=11; i < m_n ;i++ )
	    fprintf(m_fdx," %15.7g ",x[i]);
	 fprintf(m_fdx,"\n");
         fflush(m_fd);
	 fflush(m_fdx);
      }
      else
      {
         m_iconvdata.push_back(j);
         m_iconvdata.push_back(k);
         m_iconvdata.push_back(nreductions);
	 //	 m_iconvdata[3*m_ind]  = j;
	 //	 m_iconvdata[3*m_ind+1]= k;
	 //	 m_iconvdata[3*m_ind+2]= nreductions;
         m_convdata.push_back(rnorm);
         m_convdata.push_back(dxnorm);
         m_convdata.push_back(f);
	 //	 m_convdata[3*m_ind]   = rnorm;
	 //	 m_convdata[3*m_ind+1] = dxnorm;
	 //	 m_convdata[3*m_ind+2] = f;
	 for( int i = 0 ; i < m_nvar ; i++ )
            m_paradata.push_back(x[i]);
	    //	    m_paradata[m_nsize*m_ind+i] = x[i];
	 for( int i = 11 ; i < m_n ; i++ )
            m_paradata.push_back(x[i]);
	 //	    m_paradata[m_nsize*m_ind+i] = x[i];
	 m_ind++;
      }
   }
}

//-----------------------------------------------------------------------
void ConvParOutput::finish()
{
   if( m_myrank == 0 )
   {
      if( m_stepfileio )
      {
	 fclose(m_fd);
	 fclose(m_fdx);
      }
      else
      {
	 m_fd = fopen(m_convfile.c_str(),"w");
	 m_fdx= fopen(m_parafile.c_str(),"w");
         for( int i= 0 ; i < m_ind ; i ++ )
	 {
            if( m_cg )
	    {
	       fprintf(m_fd, "%i %i %15.7g %15.7g %15.7g %i \n", m_iconvdata[3*i],m_iconvdata[3*i+1], m_convdata[3*i],
		       m_convdata[3*i+1], m_convdata[3*i+2], m_iconvdata[3*i+2] );
	       fprintf(m_fdx,"%i %i ",m_iconvdata[3*i],m_iconvdata[3*i+1]);
	    }
	    else
	    {
	       fprintf(m_fd, "%i %15.7g %15.7g %15.7g %i \n", m_iconvdata[3*i], m_convdata[3*i],
		       m_convdata[3*i+1], m_convdata[3*i+2], m_iconvdata[3*i+2] );
	       fprintf(m_fdx,"%i ",m_iconvdata[3*i]);
	    }
	    for( int j= 0; j < m_nsize ;j++ )
	       fprintf(m_fdx," %15.7g ",m_paradata[m_nsize*i+j]);
	    fprintf(m_fdx,"\n");
	 }
	 fclose(m_fd);
	 fclose(m_fdx);
         m_ind = 0;
      }
   }
}

//-----------------------------------------------------------------------
ConvParOutput::~ConvParOutput()
{

}
