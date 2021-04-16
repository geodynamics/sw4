#include <cstring>
#include <unistd.h>
#include <fcntl.h>

#include "MaterialParCart.h"
#include "MaterialParCartesian.h"
#include "MaterialParCartesianVels.h"
#include "MaterialParCartesianVp.h"
#include "MaterialParCartesianVsVp.h"
#include "MaterialParAllpts.h"

#include "EW.h"

#include "SfileOutput.h"
#include "Mopt.h"

//-----------------------------------------------------------------------
Mopt::Mopt( EW* a_ew )
{
   m_ew = a_ew;
   m_opttest = 1;

   m_typrhoset = false;
   m_typmuset = false;
   m_typlambdaset = false;
   m_rhosffactor = 1;
   m_musffactor = 1;
   m_lambdasffactor = 1;
   m_nspar = 0;

   m_scales_file_given = false;
   m_rhoscale    = 1;
   m_muscale     = 1;
   m_lambdascale = 1;
   m_misfitscale = 1;
   m_vsscale     = 1;
   m_vpscale     = 1;
   MPI_Comm_rank( m_ew->m_1d_communicator, &m_myrank );
   m_optmethod = 1;
   m_nbfgs_vectors = 10;
   m_ihess_guess = 2;
   m_maxit = 10;
   m_maxsubit = 0;
   m_dolinesearch = true;
   m_fletcher_reeves = false;
   m_wolfe = false;
   m_mcheck = false;
   m_output_ts = false;
   m_test_regularizer = false;
   m_do_profiling = false;
   m_use_pseudohessian = false;
   m_tolerance = 1e-12;
   m_var    = 0;
   m_var2   = 0;
   m_itest  = 1;
   m_jtest  = 1;
   m_ktest  = 1;
   m_itest2 = 2;
   m_jtest2 = 1;
   m_ktest2 = 1;
   m_pmin = -300;
   m_pmax =  300;
   m_pmin2= -300;
   m_pmax2=  300;
   m_nsurfpts  = 10;
   m_nsurfpts2 = 10;
   m_misfit1d_images = false;
   m_path = "./";
   m_reg_coeff = 0.0;
   m_nstot = 0;
   m_sfs = NULL;
   m_sfm = NULL;
   m_xs0 = NULL;
   m_misfit = L2;
   m_mpcart0 = NULL;
   m_mp = NULL;
   m_nsteps_in_memory = 10;
   m_write_dfm = false;
}  

//-----------------------------------------------------------------------
bool Mopt::parseInputFileOpt( std::string filename )
{
   char buffer[256];
   ifstream inputFile;
   MPI_Barrier(m_ew->m_1d_communicator);

   inputFile.open(filename.c_str());
   if (!inputFile.is_open())
   {
      if (m_myrank == 0)
	 cerr << endl << "ERROR OPENING INPUT FILE: " << filename << endl << endl;
      return false;
   }
   while (!inputFile.eof())
   {    
      inputFile.getline(buffer, 256);
      if( strlen(buffer) > 0 && !(startswith("#",buffer)||startswith("\n",buffer)||startswith("\r",buffer) ) )
      {
	 if( startswith("mparcart", buffer) )
	    processMaterialParCart( buffer );
         else if( startswith("mpallpts",buffer) )
	    processMaterialAllpts( buffer );
         else if( startswith("mrun",buffer) )
	    processMrun(buffer);
         else if( startswith("mscalefactors",buffer) )
	    processMscalefactors(buffer);
         else if( startswith("lbfgs",buffer) )
	    processLBFGS(buffer);
         else if( startswith("nlcg",buffer) )
	    processNLCG(buffer);
         else if( startswith("mfsurf",buffer) )
	    processMfsurf( buffer );
         else if( startswith("mimagehdf5",buffer) )
	    processMimage( buffer, true);
         else if( startswith("mimage",buffer) )
	    processMimage( buffer, false);
	 else if (startswith("m3dimage", buffer))
	   processM3Dimage(buffer);
	 else if (startswith("sfileoutput", buffer))
	   processSfileoutput(buffer);
         else if( startswith("mtypx",buffer) )
	    processMtypx( buffer );
         else if( startswith("fileio",buffer) )
	    processMfileio( buffer );
         else if( startswith("regularize",buffer) )
	    processMregularize( buffer );
         //         else if( startswith("refinement",buffer) )
         //	    CHECK_INPUT(false,"ERROR: sw4mopt does not support mesh refinement");
      }
   }
   inputFile.close();
   MPI_Barrier(m_ew->m_1d_communicator);
   m_ew->create_directory(m_path);
   m_mp->set_path(m_path);

// wait until all processes have read the input file
   if( m_ew->getVerbosity() >=3 && m_myrank == 0 )
      cout << "********parseInputFileOpt: Done reading the input file*********" << endl;
   return true;
}

//-----------------------------------------------------------------------
void Mopt::badOption(string name, char* option) const
{
   if (m_myrank == 0)
      cout << "\tWarning: ignoring " << name << " line option '" << option << "'" << endl;
}

//-----------------------------------------------------------------------
bool Mopt::startswith(const char begin[], char *line)
{
  int lenb = strlen(begin);

  // We ignore any preceeding whitespace
  while (strncmp(line, " ", 1) == 0 || strncmp(line, "\t", 1) == 0)
    line++;
    
  if (strncmp(begin, line, lenb) == 0)
     return true;
  else 
     return false;
}

//-----------------------------------------------------------------------
void Mopt::processMaterialAllpts( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mpallpts", token) == 0,
	       "ERROR: not an mpallpts line: " << token);
   token = strtok(NULL, " \t");
   int variables = 0;

   char file[256]= " \0"; //shut up memory checker
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("type=",token) )
      {
	 token += 5;
	 int len=strlen(token);
	 if( strncmp("velocity",token,len)==0 )
	    variables = 1;
	 else if( strncmp("vsvp",token,len)==0 )
	    variables = 2;
	 else if( strncmp("vponly",token,len)==0 )
	    variables = 3;
	 else
	    variables = 0;
      }
      else
      {
	 badOption("mpallpts", token);
      }
      token = strtok(NULL, " \t");
   }
   if( m_mp != NULL )
      cout << "Error: more than one material parameterization command" << endl;

   m_mp = new MaterialParAllpts( m_ew, file, variables );
}

//-----------------------------------------------------------------------
void Mopt::processMaterialParCart( char* buffer )
{
   int nr=1;
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mparcart", token) == 0,
	       "ERROR: not an mparcart line: " << token);
   token = strtok(NULL, " \t");

   bool vel = false, vponly=false, vsvp=false, mparcartfile =false;
   bool shared=false;
   int nx=3, ny=3, nz=3, init=0;
   double ratio=1.732, gamma=1, amp=0.1, omega=2*M_PI;
   char file[256]= " \0"; //shut up memory checker

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("type=",token) )
      {
	 token += 5;
	 int len=strlen(token);
	 vel    =strncmp("velocity",token,len)==0;
	 vponly =strncmp("vponly",token,len)==0;
	 vsvp   =strncmp("vsvp",token,len)==0;
      }
      else if( startswith("nx=",token) )
      {
         token += 3;
         nx = atoi(token);
      }
      else if( startswith("ny=",token) )
      {
         token += 3;
         ny = atoi(token);
      }
      else if( startswith("nz=",token) )
      {
         token += 3;
         nz = atoi(token);
      }
      else if( startswith("init=",token) )
      {
         token += 5;
	 if( strcmp(token,"0")==0 )
            init = 0;
	 else if( strcmp(token,"sineperturbation")==0 )
	    init = 1;
	 else if( strcmp(token,"interpolate")== 0 )
	    init = 3;
	 else
	 {
            strncpy(file,token,256);
	    init = 2;
	 }
      }
      else if( startswith("filetype=",token) )
      {
	 token += 9;
	 CHECK_INPUT(strcmp(token,"mpar")==0 || strcmp(token,"mpc")==0, 
                   "ERROR: processing mparcart, file type " << token << "not recognized" );
	 if( strcmp(token,"mpc")== 0 )
	 {
	    mparcartfile = true;
	 }
      }
      else if( startswith("amplitude=",token) )
      {
         token += 10;
         amp = atof(token);
      }
      else if( startswith("omega=",token) )
      {
         token += 6;
         omega = atof(token);
      }
      else if( startswith("shared=",token) )
      {
         token += 7;
         shared = strcmp(token,"1") == 0   ||
                  strcmp(token,"yes") == 0 || 
                  strcmp(token,"on") == 0;
      }
      else
      {
	 badOption("mparcart", token);
      }
      token = strtok(NULL, " \t");
   }
   // if mparcartfiletype and file name given, then ok to read mparcart file.
   // Ignore filetype if init does not provide a file name.
   if( mparcartfile && init == 2 )
      init = 4;

   if( m_mp != NULL )
      cout << "Error: more than one material parameterization command" << endl;

   // Make sure the material grid is coarser than the global grid
   VERIFY2( nx <= m_ew->m_global_nx[0] && ny <= m_ew->m_global_ny[0] && nz <= m_ew->m_global_nz[0], "MaterialParCart: The material grid must be coarser than the global grid")

   int varcase=1;
   if( vel )
      varcase=2;
   else if( vsvp )
      varcase=3;
   else if( vponly )
      varcase=4;
   if( !shared )
      m_mp = new MaterialParCart( m_ew, nx, ny, nz, init, varcase, file, amp, omega );
   else
   {
      if (vel)
         m_mp = new MaterialParCartesianVels( m_ew, nx, ny, nz, init, file );
      else if( vponly )
         m_mp = new MaterialParCartesianVp( m_ew, nx, ny, nz, init, file, ratio, gamma, true );
      else if( vsvp )
         m_mp = new MaterialParCartesianVsVp( m_ew, nx, ny, nz, init, file );
      else
         m_mp = new MaterialParCartesian( m_ew, nx, ny, nz, init, file );
   }
   // use for material projection only:
   m_mpcart0 = new MaterialParCartesian( m_ew, nx, ny, nz, 0, "");
} // end processMaterialParCart

//-----------------------------------------------------------------------
void Mopt::processMrun( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mrun", token) == 0,
	       "ERROR: not an mrun line: " << token);
   token = strtok(NULL, " \t");
   bool quiet = true;
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("task=",token) )
      {
         token += 5;
	 if( strcmp(token,"minvert") == 0 )
	    m_opttest = 1;
	 else if( strcmp(token,"gradtest") == 0 )
	    m_opttest = 2;
	 else if( strcmp(token,"hesstest") == 0 )
	    m_opttest = 3;
	 else if( strcmp(token,"func1d") == 0 )
	    m_opttest = 4;
	 else if( strcmp(token,"func2d") == 0 )
	    m_opttest = 5;
	 else if( strcmp(token,"forward") == 0 )
	    m_opttest = 6;
	 else if( strcmp(token,"projectmtrl") == 0 )
	    m_opttest = 7;
	 else if( strcmp(token,"testmtrl") == 0 )
	    m_opttest = 8;
	 else if( strcmp(token,"computegrad") == 0 )
	    m_opttest = 9;
	 else if( strcmp(token,"minvert+src11") == 0 )
	 {
	    m_opttest = 1;
            m_nspar = 11;
	 }
	 else if( strcmp(token,"minvert+src10") == 0 )
	 {
	    m_opttest = 1;
            m_nspar = 10;
	 }
	 else if( strcmp(token,"minvert+src9") == 0 )
	 {
	    m_opttest = 1;
            m_nspar = 9;
	 }
	 else if( strcmp(token,"minvert+src6") == 0 )
	 {
	    m_opttest = 1;
            m_nspar = 6;
	 }
	 else
	    cout << "ERROR: mrun task=" << token << " not recognized " << endl;
      }
      else if( startswith("mcheck=",token) )
      {
         token += 7;
         CHECK_INPUT(strcmp("on",token)== 0||strcmp("off",token)== 0,"ERROR, mrun mcheck= " << token <<
		     " not understood" << endl);
	 m_mcheck = strcmp("on",token)== 0;
      }
      else if( startswith("tsoutput=",token) )
      {
         token += 9;
         CHECK_INPUT(strcmp("on",token)== 0||strcmp("off",token)== 0,"ERROR, mrun tsoutput= " << token <<
		     " not understood" << endl);
	 m_output_ts = strcmp("on",token)== 0;
      }
      else if( startswith("quiet=",token) )
      {
	 token += 6;
         CHECK_INPUT(strcmp("yes",token)== 0||strcmp("no",token)== 0,"ERROR, mrun quiet= " << token <<
		     " not understood" << endl);
	 m_ew->setQuiet(strcmp("yes",token)==0);
      }
      else if( startswith("misfit=",token) )
      {
	 token += 7;
	 if( strcmp(token,"l2")==0 )
	    m_misfit = L2;
	 else if( strcmp(token,"crosscorr")==0 )
	    m_misfit = CROSSCORR;
	 else
	    CHECK_INPUT(false,"ERROR, mrun misfit= " << token << " not understood" << endl);
      }
      else if( startswith("savesteps=",token) )
      {
 	 token += 10;
	 m_nsteps_in_memory = atoi(token);
      }
      else if( startswith("profiling=",token) )
      {
	 token += 10;
	 int n = strlen(token);
	 if( strncmp("yes",token,n)== 0 || strncmp("on",token,n)==0 || strncmp("1",token,n)== 0)
	    m_do_profiling = true;
      }
      else if( startswith("pseudohessian=",token) )
      {
         token += 14;
	 int n = strlen(token);
	 if( strncmp("yes",token,n)== 0 || strncmp("on",token,n)==0 || strncmp("1",token,n)== 0)
	    m_use_pseudohessian = true;
      }
      else if( startswith("zerosrc=",token) )
      {
         token += 8;
	 int n = strlen(token);
	 if( strncmp("yes",token,n)== 0 || strncmp("on",token,n)==0 || strncmp("1",token,n)== 0 )
            m_ew->set_zerograd();
         //	    m_zerograd_at_src = true;
      }
      else if( startswith("zerosrcpad=",token) )
      {
         token += 11;
	 int n = strlen(token);
         m_ew->set_zerograd_pad(atoi(token));
         //         m_zerograd_pad = atoi(token);
      }
      else if( startswith("writedfm=",token) )
      {
         token += 9;
	 int n = strlen(token);
	 if( strncmp("yes",token,n)== 0 || strncmp("1",token,n)== 0 || strncmp("on",token,n)==0 )
	    m_write_dfm = true;
      }
      else
         badOption("mrun",token);
      token = strtok(NULL," \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::processMscalefactors( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mscalefactors", token) == 0,
	       "ERROR: not an mscalefactors line: " << token);
   token = strtok(NULL, " \t");
   double misfitscale=1;
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("rho=",token) )
      {
         token += 4;
	 m_rhoscale = atof(token);
      }
      else if( startswith("mu=",token) )
      {
         token += 3;
	 m_muscale = atof(token);
      }
      else if( startswith("lambda=",token) )
      {
         token += 7;
	 m_lambdascale = atof(token);
      }
      else if( startswith("misfit=",token) )
      {
         token += 7;
	 m_misfitscale = atof(token);
      }
      else if( startswith("vp=",token) )
      {
	 token += 3;
	 m_vpscale = atof(token);
      }
      else if( startswith("vs=",token) )
      {
	 token += 3;
	 m_vsscale = atof(token);
      }
      else if( startswith("file=",token) )
      {
         token += 5;
	 m_scales_fname = token;
         m_scales_file_given = true;
      }
      else
         badOption("mscalefactors",token);
      token = strtok(NULL," \t");
   }
   double imf = 1/sqrt(m_misfitscale);
   m_rhoscale    *= imf;
   m_muscale     *= imf;
   m_lambdascale *= imf;
   m_vsscale     *= imf;
   m_vpscale     *= imf;
}

//-----------------------------------------------------------------------
void Mopt::processMfileio( char* buffer )
{
   char* path = 0;
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("fileio", token) == 0, "ERROR: not a fileio line...: " << token);
   token = strtok(NULL, " \t");

   string err = "FileIO Error: ";

   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
	 break;
      if(startswith("path=", token)) 
      {
	 token += 5;
	 m_path = token;
	 m_path += '/';
      }
      token = strtok(NULL, " \t");
   }   
}

//-----------------------------------------------------------------------
void Mopt::processMregularize( char* buffer )
{
   char* path = 0;
   char* token = strtok(buffer, " \t");
   int n;
   CHECK_INPUT(strcmp("regularize", token) == 0, "ERROR: not a regularize line...: " << token);
   token = strtok(NULL, " \t");
   //   double scale_coeff = 1;

   string err = "Regularize Error: ";

   while (token != NULL)
   {
      if (startswith("#", token) || startswith(" ", buffer))
	 break;
      if(startswith("coeff=", token)) 
      {
	 token += 6;
         //	 scale_coeff = atof(token);
         m_reg_coeff = atof(token);
      }
      else if( startswith("testmode=",token))
      {
	 token += 9;
	 n = strlen(token);
	 if( strncmp(token,"yes",n)==0 || strncmp(token,"on",n)==0 || strncmp(token,"1",n)==0 )
	    m_test_regularizer = true;
	 else
	    m_test_regularizer = false;
      }
      else
         badOption("regularize",token);
      token = strtok(NULL, " \t");
   }
   //   m_reg_coeff = scale_coeff;
   
}// end processMregularize

//-----------------------------------------------------------------------
void Mopt::get_scalefactors( double& rhoscale, double& muscale,
			     double& lambdascale )
{
   rhoscale = m_rhoscale;
   muscale = m_muscale;
   lambdascale = m_lambdascale;
}

//-----------------------------------------------------------------------
void Mopt::set_sscalefactors( )
{
   int my_nmpars, nmpard, nmpard_global;
   m_mp->get_nr_of_parameters( my_nmpars, nmpard, nmpard_global );
   m_nstot = m_nspar + my_nmpars;
   if( m_nstot > 0 )
   {
      m_sfs = new double[m_nstot];
      double *sfs = &m_sfs[m_nspar];

// first set the source scale factors (if any)
      set_sourcescalefactors( m_nspar, m_sfs );

      if( !m_scales_file_given )
	 m_mp->set_scalefactors( my_nmpars, sfs, m_rhoscale, m_muscale, m_lambdascale, m_vsscale, m_vpscale );
      else
      {
	 int errflag = 0;
	 if( m_myrank == 0 )
	 {
	    cout << "Reading scale factor file " << " my_nmpars = " << my_nmpars << endl;
	    string fname = m_path + m_scales_fname;
	 //	 int fd=open(m_scales_fname.c_str(),O_RDONLY );
	    int fd=open(fname.c_str(),O_RDONLY );
	    VERIFY2( fd != -1, "ERROR reading scale factors. File " << m_scales_fname << " could not be opened"<<endl);
	    int nmpars_read;
	    size_t nr = read(fd,&nmpars_read,sizeof(int));
	    if( nmpars_read == my_nmpars && nr == sizeof(int) )
	    {
	       nr = read(fd,sfs,my_nmpars*sizeof(double));
	       if( nr != my_nmpars*sizeof(double) )
		  errflag = 2;
	    }
	    else if( nmpars_read != my_nmpars )
	       errflag = 3;
	    else
	       errflag = 1;
	    close(fd);

	    double imf = 1/sqrt(m_misfitscale);
	    for( int i=0 ; i < my_nmpars ; i++ )
	       sfs[i] *= imf;
	 }
	 MPI_Bcast( sfs, my_nmpars, MPI_DOUBLE, 0, m_ew->m_1d_communicator );
	 VERIFY2( errflag == 0, "Error no " << errflag << " in Mopt::set_sscalefactors");
      }
   }
}

//-----------------------------------------------------------------------
void Mopt::set_baseMat(double* xs, double* xm )
{
  int my_nmpars, nmpard, nmpard_global;
  m_mp->get_nr_of_parameters( my_nmpars, nmpard, nmpard_global );
  m_nstot = m_nspar + my_nmpars;
  m_xs0 = new double[m_nstot];
  m_xm0 = new double[nmpard];

  for( int i=0 ; i < my_nmpars ; i++)
     m_xs0[i]  = xs[i];
  for( int i=0 ; i < nmpard ; i++)
     m_xm0[i] = xm[i];
}

//-----------------------------------------------------------------------
void Mopt::set_dscalefactors()
{
   int nmpars, nmpard, nmpard_global;
   m_mp->get_nr_of_parameters( nmpars, nmpard, nmpard_global );
   if( nmpard > 0 )
   {
      m_sfm = new double[nmpard];
      if( !m_scales_file_given )
	 m_mp->set_scalefactors( nmpard, m_sfm, m_rhoscale, m_muscale, m_lambdascale, m_vsscale, m_vpscale );
      else
      {
   // Not yet implemented
	 for( int i=0 ; i<nmpard ;i++)
	    m_sfm[i]=1;
      }
   }
}

//-----------------------------------------------------------------------
void Mopt::set_sourcescalefactors( int nspar, double* sfs )
{
   double sfall[11];
   m_ew->get_scalefactors( sfall );
   if( nspar == 11 )
      for( int i=0; i<11 ;i++ )
	 sfs[i] = sfall[i];
   else if( nspar == 10 )
      for( int i=0; i<10 ;i++ )
	 sfs[i] = sfall[i];
   else if( nspar == 9 )
      for( int i=0; i<9 ;i++ )
	 sfs[i] = sfall[i];
   else if( nspar == 6 )
      for( int i=0 ; i< 6 ;i++)
	 sfs[i] = sfall[i+3];
   else if( nspar !=0 )
      cout << "Error in set_sourcescalefactors, nspar = " << nspar
	   << " undefined case"<< endl;
}

//-----------------------------------------------------------------------
void Mopt::processLBFGS( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("lbfgs", token) == 0,
	       "ERROR: not an lbfgs line: " << token);
   token = strtok(NULL, " \t");
   m_optmethod = 1;
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("nvectors=",token) )
      {
         token += 9;
	 m_nbfgs_vectors = atoi(token);
	 CHECK_INPUT( m_nbfgs_vectors > 0, 
           "ERROR: lbfgs, nvectors must be positive, not " << m_nbfgs_vectors << endl );
      }
      else if( startswith("ihess0=",token) )
      {
         token += 7;
         if( strcmp(token,"scale-factors")== 0 )
	    m_ihess_guess = 1;
	 else if( strcmp(token,"gamma")==0 )
	    m_ihess_guess = 2;
         else
	    CHECK_INPUT(false,"ERROR: lbfgs ihess0="  << token << " not understood " << endl);
      }
      else if( startswith("maxit=",token) )
      {
         token += 6;
	 m_maxit = atoi(token);
      }
      else if( startswith("tolerance=",token) )
      {
         token += 10;
	 m_tolerance = atof(token);
      }
      else if( startswith("linesearch=",token) )
      {
         token += 11;
         CHECK_INPUT( strcmp(token,"on")==0 ||  strcmp(token,"off")==0 || strcmp(token,"wolfe")==0,
		      "lbfgs: ERROR: linesearch=" << token << "not understood");
         if( strcmp(token,"wolfe")== 0 )
	 {
            m_dolinesearch = true;
	    m_wolfe = true;
 	 }
         else 
	 {
	    m_dolinesearch = strcmp(token,"on")==0;
            m_wolfe = false;
	 }
      }
      else
         badOption("lbfgs",token);
      token = strtok(NULL," \t");
   }
}


//-----------------------------------------------------------------------
void Mopt::processNLCG( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("nlcg", token) == 0,
	       "ERROR: not an nlcg line: " << token);
   token = strtok(NULL, " \t");
   m_optmethod = 2;
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("maxsubit=",token) )
      {
	 token += 9;
	 m_maxsubit = atoi(token);
      }
      else if( startswith("maxit=",token) )
      {
         token += 6;
	 m_maxit = atoi(token);
      }
      else if( startswith("tolerance=",token) )
      {
         token += 10;
	 m_tolerance = atof(token);
      }
      else if( startswith("linesearch=",token) )
      {
         token += 11;
         CHECK_INPUT( strcmp(token,"on")==0 ||  strcmp(token,"off")==0,
		      "nlcg: ERROR: linesearch=" << token << "not understood");
         m_dolinesearch = strcmp(token,"on")==0;
      }
      else if( startswith("subtype=",token) )
      {
         token += 8;
         m_fletcher_reeves = strcmp(token,"fletcher-reeves")==0
	 || strcmp(token,"Fletcher-Reeves")==0;
         if( !m_fletcher_reeves )
	    CHECK_INPUT( strcmp(token,"Polak-Ribiere")==0 || strcmp(token,"polak-ribiere")== 0,
			 "ERROR: nclg subtype " << token << " not understood" << endl);
      }
      else
         badOption("nlcg",token);
      token = strtok(NULL," \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::processMfsurf( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mfsurf", token) == 0,
	       "ERROR: not an mfsurf line: " << token);
   token = strtok(NULL, " \t");

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("var=",token) )
      {
         token += 4;
         CHECK_INPUT( strcmp(token,"rho")==0 || strcmp(token,"mu")==0 || strcmp(token,"lambda")==0
		      || strcmp(token,"vp")==0 || strcmp(token,"vs")==0,
		      "ERROR: mfsurf, var= " << token << " not recognized"<<endl);
	 if( strcmp(token,"rho")==0 )
	    m_var = 0;
	 else if( strcmp(token,"mu")==0 || strcmp(token,"vs")==0 )
	    m_var = 1;
	 else
	    m_var = 2;
      }
      else if( startswith("var2=",token) )
      {
         token += 5;
         CHECK_INPUT( strcmp(token,"rho")==0 || strcmp(token,"mu")==0 || strcmp(token,"lambda")==0,
		      "ERROR: mfsurf, var2= " << token << " not recognized"<<endl);
	 if( strcmp(token,"rho")==0 )
	    m_var2 = 0;
	 else if( strcmp(token,"mu")==0 )
	    m_var2 = 1;
	 else
	    m_var2 = 2;
      }
      else if( startswith("i=",token) )
      {
         token += 2;
         m_itest = atoi(token);
      }
      else if( startswith("j=",token) )
      {
         token += 2;
         m_jtest = atoi(token);
      }
      else if( startswith("k=",token) )
      {
         token += 2;
         m_ktest = atoi(token);
      }
      else if( startswith("i2=",token) )
      {
         token += 3;
         m_itest2 = atoi(token);
      }
      else if( startswith("j2=",token) )
      {
         token += 3;
         m_jtest2 = atoi(token);
      }
      else if( startswith("k2=",token) )
      {
         token += 3;
         m_ktest2 = atoi(token);
      }
      else if( startswith("pmin=",token))
      {
         token += 5;
	 m_pmin = atof(token);
      }
      else if( startswith("pmax=",token))
      {
         token += 5;
	 m_pmax = atof(token);
      }
      else if( startswith("pmin2=",token))
      {
         token += 6;
	 m_pmin2 = atof(token);
      }
      else if( startswith("pmax2=",token))
      {
         token += 6;
	 m_pmax2 = atof(token);
      }
      else if( startswith("npts=",token))
      {
         token += 5;
	 m_nsurfpts = atoi(token);
      }
      else if( startswith("npts2=",token))
      {
         token += 6;
	 m_nsurfpts2 = atoi(token);
      }
      else if( startswith("imagesave=",token) )
      {
	 token += 10;
	 m_misfit1d_images = (strcmp(token,"1")==0)||(strcmp(token,"yes")==0);
      }
      else
         badOption("mfsurf",token);
      token = strtok(NULL," \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::processMimage( char* buffer, bool use_hdf5)
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mimage", token) == 0,
	       "ERROR: not an mimage line: " << token);

   int iter=1, iterInterval=0;
   Image::ImageMode mode=Image::RHO;
   Image::ImageOrientation locationType=Image::UNDEFINED;
   string filePrefix="mimage";
   double coordValue;
   int gridPointValue;
   bool coordWasSet = false;
   bool use_double  = false;

   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("iter=",token) )
      {
	 token += 5; // skip iter=
	 iter = atoi(token);
      }
      else if (startswith("iterInterval=", token) )
      {
	 token += 13; // skip iterInterval=
	 iterInterval = atoi(token);
      }
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if (startswith("mode=", token))
      {
	 token += 5; // skip mode=
	 if (strcmp(token, "rho") == 0)   mode = Image::RHO;
	 else if (strcmp(token, "lambda") == 0)   mode = Image::LAMBDA;
	 else if (strcmp(token, "mu") == 0)   mode = Image::MU;
	 else if (strcmp(token, "p") == 0)   mode = Image::P;
	 else if (strcmp(token, "s") == 0)   mode = Image::S;
	 else if (strcmp(token, "gradrho") == 0)   mode = Image::GRADRHO;
	 else if (strcmp(token, "gradmu") == 0)   mode = Image::GRADMU;
	 else if (strcmp(token, "gradlambda") == 0)   mode = Image::GRADLAMBDA;
	 else if (strcmp(token, "gradp") == 0)   mode = Image::GRADP;
	 else if (strcmp(token, "grads") == 0)   mode = Image::GRADS;
         else
	    CHECK_INPUT( false, "Processing mimage command: mode must be one of the following: " << endl
			 << "rho|mu|lambda|p|s|gradrho|gradmu|gradlambda|gradp|grads" << endl
			 << "not " << token << endl );
      }
      else if( startswith("precision=",token) )
      {
	 token += 10;
         CHECK_INPUT( strcmp(token,"double")==0 || strcmp(token,"float")==0, "Processing mimage command: "
		      << " precision must be float or double, not " << token << endl );
	 use_double =  strcmp(token,"double")==0;
      }
      else if (startswith("x=", token))
      {
	 token += 2; // skip x=
	 if ( !coordWasSet )
	 {
	    coordWasSet  = true;
	    locationType = Image::X;
	    coordValue   = atof(token);
	 }
      }
      else if (startswith("y=", token))
      {
	 token += 2; // skip y=
	 if ( !coordWasSet )
	 {
	    coordWasSet  = true;
	    locationType = Image::Y;
	    coordValue   = atof(token);
	 }
      }
      else if (startswith("z=", token))
      {
	 token += 2; // skip z=
	 if ( !coordWasSet )
	 {
	    coordWasSet  = true;
	    locationType = Image::Z;
	    coordValue   = atof(token);
	 }
      }
      else
	 badOption("mimage", token);
      token = strtok(NULL, " \t");
   }
   CHECK_INPUT( coordWasSet, "ERROR: Processing image command: one of the coordinate (x,y,z) option must be set " <<
		" to determine the image's 2D plane" << endl);
   double time=-1, timeInterval=-1;
   Image* i = new Image(m_ew, time, timeInterval, iter, iterInterval, 
		 filePrefix, mode, locationType, coordValue, use_double, use_hdf5);
   i->computeGridPtIndex();
   i->allocatePlane();
   m_image_files.push_back(i);
}

//-----------------------------------------------------------------------
void Mopt::processSfileoutput( char* buffer )
{
   int cycle=0, cycleInterval=0;
   int sampleFactorV = 1;
   int sampleFactorH = 1;
   float_sw4 time=0.0, timeInterval=0.0;
   bool timingSet = false;
   float_sw4 tStart = -999.99;
   string filePrefix="sfileoutput";
   bool use_double = false;
  
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("sfileoutput", token) == 0, "ERROR: Not a sfileoutput line...: " << token );

   token = strtok(NULL, " \t");
   string err = "sfileoutput Error: ";
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	 // Ignore commented lines and lines with just a space.
	 break;
      /* if (startswith("time=", token) ) */
      /* { */
	 /* token += 5; // skip time= */
	 /* CHECK_INPUT( atof(token) >= 0.,"Processing sfileoutput command: time must be a non-negative number, not: " << token); */
	 /* time = atof(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("timeInterval=", token) ) */
      /* { */
	 /* token += 13; // skip timeInterval= */
	 /* CHECK_INPUT( atof(token) >= 0.,"Processing sfileoutput command: timeInterval must be a non-negative number, not: " << token); */
	 /* timeInterval = atof(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("startTime=", token) ) */
      /* { */
	 /* token += 10; // skip startTime= */
	 /* tStart = atof(token); */
      /* } */
      else if (startswith("sampleFactorH=", token) )
      {
	 token += 14; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactorH must be a positive integer, not: " << token);
	 sampleFactorH = atoi(token);
      }
      else if (startswith("sampleFactorV=", token) )
      {
	 token += 14; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactorV must be a positive integer, not: " << token);
	 sampleFactorV= atoi(token);
      }
      else if (startswith("sampleFactor=", token) )
      {
	 token += 13; 
	 CHECK_INPUT( atoi(token) >= 1,"Processing sfileoutput command: sampleFactor must be a positive integer, not: " << token);
	 sampleFactorH = sampleFactorV = atoi(token);
      }
      /* else if (startswith("cycle=", token) ) */
      /* { */
	 /* token += 6; // skip cycle= */
	 /* CHECK_INPUT( atoi(token) >= 0.,"Processing sfileoutput command: cycle must be a non-negative integer, not: " << token); */
	 /* cycle = atoi(token); */
	 /* timingSet = true; */
      /* } */
      /* else if (startswith("cycleInterval=", token) ) */
      /* { */
	 /* token += 14; // skip cycleInterval= */
	 /* CHECK_INPUT( atoi(token) >= 0.,"Processing sfileoutput command: cycleInterval must be a non-negative integer, not: " << token); */
	 /* cycleInterval = atoi(token); */
	 /* timingSet = true; */
      /* } */
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if( startswith("precision=",token) )
      {
	 token += 10;
	 CHECK_INPUT( startswith("double",token) || startswith("float",token),
		      "Processing sfileoutput command: precision must be float or double, not '" << token );
	 use_double =  startswith("double",token);
      }
      else
      {
	 badOption("sfileoutput", token);
      }
      token = strtok(NULL, " \t");
   }

   SfileOutput* sfile = new SfileOutput( m_ew, time, timeInterval, cycle, cycleInterval, 
     		       tStart, filePrefix, sampleFactorH, sampleFactorV, use_double);
   sfile->setup_images( );
   m_sfiles.push_back(sfile);
}



//-----------------------------------------------------------------------
void Mopt::processM3Dimage( char* buffer )
{
   int iter=-1, iterInterval=0;
   Image3D::Image3DMode mode=Image3D::RHO;
   double time=0.0, timeInterval=0.0;
   bool timingSet = false;
   double tStart = -999.99;
   string filePrefix="m3dimage";
   bool use_double = false;
  
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("m3dimage", token) == 0, "ERROR: Not a m3dimage line...: " << token );

   token = strtok(NULL, " \t");
   string err = "m3dimage Error: ";
   while (token != NULL)
   {
     // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
      {
	 // Ignore commented lines and lines with just a space.
	 break;
      }
      if (startswith("iter=", token) )
      {
	 token += 5; // skip iter=
	 CHECK_INPUT( atoi(token) >= 0.,"Processing m3dimage command: iter must be a non-negative integer, not: " << token);
	 iter = atoi(token);
	 timingSet = true;
      }
      else if (startswith("iterInterval=", token) )
      {
	 token += 13; // skip iterInterval=
	 CHECK_INPUT( atoi(token) >= 0.,"Processing m3dimage command: iterInterval must be a non-negative integer, not: " << token);
	 iterInterval = atoi(token);
	 timingSet = true;
      }
      else if (startswith("file=", token))
      {
	 token += 5; // skip file=
	 filePrefix = token;
      }
      else if (startswith("mode=", token))
      {
	 token += 5; // skip mode=
	 if (strcmp(token, "rho") == 0)   mode = Image3D::RHO;
	 else if (strcmp(token, "p") == 0)   mode = Image3D::P;
	 else if (strcmp(token, "s") == 0)   mode = Image3D::S;
	 else if (strcmp(token, "mu") == 0)   mode = Image3D::MU;
	 else if (strcmp(token, "lambda") == 0)   mode = Image3D::LAMBDA;
	 else if (strcmp(token, "gradrho") == 0)   mode = Image3D::GRADRHO;
	 else if (strcmp(token, "gradp") == 0)   mode = Image3D::GRADP;
	 else if (strcmp(token, "grads") == 0)   mode = Image3D::GRADS;
	 else if (strcmp(token, "gradmu") == 0)   mode = Image3D::GRADMU;
	 else if (strcmp(token, "gradlambda") == 0)   mode = Image3D::GRADLAMBDA;
	 else
	 {
	    //	    mode = static_cast<Image3D::Image3DMode>(atoi(token));
	    CHECK_INPUT(0,"Processing m3dimage command: " << "mode must be one of the following: " << endl <<
			"\t rho|p|s|mu|lambda|gradrho|gradp|grads|gradmu|gradlambda *not: "<< token );
	 }
      }
      else if( startswith("precision=",token) )
      {
	 token += 10;
	 CHECK_INPUT( startswith("double",token) || startswith("float",token),
		      "Processing m3dimage command: precision must be float or double, not '" << token );
	 use_double =  startswith("double",token);
      }
      else
      {
	 badOption("m3dimage", token);
      }
      token = strtok(NULL, " \t");
   }

   CHECK_INPUT( timingSet, "Processing m3dimage command: " << 
		"at least one timing mechanism must be set: iter, iterInterval"  << endl );
   Image3D* im3 = new Image3D( m_ew, time, timeInterval, iter, iterInterval, 
 			       tStart, filePrefix, mode, use_double );
   im3->setup_images( );
   m_3dimage_files.push_back(im3);
}

//-----------------------------------------------------------------------
void Mopt::processMtypx( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mtypx", token) == 0,
	       "ERROR: not an mtypx line: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("rho=",token) )
      {
         token += 4;
	 m_typrhoset = true;
	 m_typrho = atof(token);
      }
      else if( startswith("mu=",token) )
      {
         token += 3;
	 m_typmuset = true;
	 m_typmu = atof(token);
      }
      else if( startswith("lambda=",token) )
      {
         token += 7;
	 m_typlambdaset = true;
	 m_typlambda = atof(token);
      }
      else if( startswith("rhosffactor=",token) )	 
      {
         token += 12;
	 m_rhosffactor = atof(token);
      }
      else if( startswith("musffactor=",token) )	 
      {
         token += 11;
	 m_musffactor = atof(token);
      }
      else if( startswith("lambdasffactor=",token) )	 
      {
         token += 15;
	 m_lambdasffactor = atof(token);
      }
      else
	 badOption("mtypx", token);
      token = strtok(NULL, " \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::set_typx( int nmpar, double* sf, double* typx )
{
   for( int i=0 ; i < nmpar ; i += 3 )
   {
      if( m_typrhoset )
	 typx[i] = m_typrho;
      else
	 typx[i] = m_rhosffactor*sf[i];
      
      if( m_typmuset )
	 typx[i+1] = m_typmu;
      else
	 typx[i+1] = m_musffactor*sf[i+1];

      if( m_typlambdaset )
	 typx[i+2] = m_typlambda;
      else
	 typx[i+2] = m_lambdasffactor*sf[i+2];
   }
}

//-----------------------------------------------------------------------
void Mopt::init_pseudohessian( vector<Sarray>& ph )
{
   if( m_use_pseudohessian )
   {
      for( int g=0 ; g < m_ew->mNumberOfCartesianGrids ; g++ )
      {
         ph[g].define(3,m_ew->m_iStart[g],m_ew->m_iEnd[g],
                        m_ew->m_jStart[g],m_ew->m_jEnd[g],
                        m_ew->m_kStart[g],m_ew->m_kEnd[g] );
         ph[g].set_to_zero();
      }
   }
}

//-----------------------------------------------------------------------
int Mopt::get_pseudo_hessian_case( )
{
   if( m_use_pseudohessian )
      return m_mp->get_varcase();
   else
      return 0;
}
