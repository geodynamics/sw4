#include "EW.h"
#include "DataPatches.h"
#include <cstring>
#include "version.h"
//#include "MaterialParameterization.h"
#include "MaterialParAllpts.h"
#include "MaterialParCartesian.h"
#include "Mopt.h"

#include <fcntl.h>
#include <unistd.h>

void lbfgs( EW& simulation, int nspar, int nmpars, double* xs, double* sf, 
	    int nmpard, double* xm, double* sfm,
	    vector<Source*>& GlobalSources,
	    vector<TimeSeries*>& GlobalTimeSeries,
	    vector<TimeSeries*> & GlobalObservations,
	    int m, int myRank, MaterialParameterization* mp );


void usage(string thereason)
{
  cout << endl
       << "sbp4mopt - Summation by parts 4th order inverse seismic wave solver"  << endl << endl
       << "Usage: sbp4mopt [-v] file.in" << endl
       << "\t -v:      prints out the version info" << endl
       << "\t file.in: an input file" << endl << endl
       << "Reason for message: " << thereason << endl;
}



//-----------------------------------------------------------------------
void set_source_pars( int nspar, double srcpars[11], double* xs )
{
   // nspar =11, all parameters=(x0,y0,z0,m_{xx}-m_{zz},t0,freq)
   // nspar =10, no freq., parameters=(x0,y0,z0,m_{xx}-m_{zz},t0)
   // nspar = 9 , no freq or t0, parameters=(x0,y0,z0,m_{xx}-m_{zz})
   // nspar = 0, source not part of inversion.
   if( nspar == 11 )
      for( int i=0; i<11 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar == 10 )
      for( int i=0; i<10 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar == 9 )
      for( int i=0; i<9 ;i++ )
	 srcpars[i] = xs[i];
   else if( nspar !=0 )
      cout << "Error in set_source_pars, nspar = " << nspar
	   << " undefined case"<< endl;
}

//-----------------------------------------------------------------------
void get_source_pars( int nspar, double srcpars[11], double* xs )
{
   // nspar =11, all parameters=(x0,y0,z0,m_{xx}-m_{zz},t0,freq)
   // nspar =10, no freq., parameters=(x0,y0,z0,m_{xx}-m_{zz},t0)
   // nspar = 9 , no freq or t0, parameters=(x0,y0,z0,m_{xx}-m_{zz})
   // nspar = 0, source not part of inversion.
   if( nspar == 11 )
      for( int i=0; i<11 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar == 10 )
      for( int i=0; i<10 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar == 9 )
      for( int i=0; i<9 ;i++ )
	 xs[i] = srcpars[i];
   else if( nspar !=0 )
      cout << "Error in get_source_pars, nspar = " << nspar
	   << " undefined case"<< endl;
}

//-----------------------------------------------------------------------
void compute_f( EW& simulation, int nspar, int nmpars, double* xs,
		int nmpard, double* xm,
		vector<Source*>& GlobalSources,
		vector<TimeSeries*>& GlobalTimeSeries,
		vector<TimeSeries*>& GlobalObservations,
		double& mf, MaterialParameterization* mp )
//-----------------------------------------------------------------------
// Compute misfit.
//
// Input: simulation - Simulation object
//        nspar  - Number of source parameters
//        nmpars - Number of shared material parameters
//        xs     - Vector of shared unknown, xs[0..nspar-1] are the source
//                 parameters, xs[nspar..nmpars+nspar-1] are the material parameters.
//        nmpard - Number of distributed material parameters
//        xm     - Vector of distributed unknown, of size nmpard.
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        mp - Pointer to object describing the parameterization of the material.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         mf               - The misfit.
//-----------------------------------------------------------------------
{
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy(" ");

   // fetch all 11 source parameters, and merge in the unknowns
   double srcpars[11];
   src[0]->get_parameters( srcpars );
   set_source_pars( nspar, srcpars, xs );
   src[0]->set_parameters( srcpars );

// Translate one-dimensional parameter vector xm to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);

//New
   mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda );
// Old
 //   simulation.parameters_to_material( nm, xm, rho, mu, lambda );

// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng);
   simulation.solve_allpars( src, rho, mu, lambda, GlobalTimeSeries, U, Um, upred_saved, ucorr_saved );

// Compute misfit
   mf = 0;
   double dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      mf += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], NULL, dshift, ddshift, dd1shift );

   double mftmp = mf;
   MPI_Allreduce(&mftmp,&mf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

// Give back memory
   for( unsigned int g=0 ; g < ng ; g++ )
   {
      delete upred_saved[g];
      delete ucorr_saved[g];
   }
   delete src[0];
}

//-----------------------------------------------------------------------
void compute_f_and_df( EW& simulation, int nspar, int nmpars, double* xs,
		       int nmpard, double* xm,
		       vector<Source*>& GlobalSources,
		       vector<TimeSeries*>& GlobalTimeSeries,
		       vector<TimeSeries*>& GlobalObservations, 
		       double& f, double* dfs, double* dfm, int myrank,
		       MaterialParameterization* mp )

//-----------------------------------------------------------------------
// Compute misfit and its gradient.
//
// Input: simulation - Simulation object
//        nspar  - Number of source parameters (shared among processors).
//        nmpars - Number of shared material parameters.
//        xs     - Vector of shared unknown, xs[0..nspar-1] are the source
//                 parameters, xs[nspar..nmpars+nspar-1] are the material parameters.
//        nmpard - Number of distributed material parameters
//        xm     - Vector of distributed unknown, of size nmpard.
//        GlobalSources - The single source object
//        GlobalTimeSeries   - TimeSeries objects, number of objects and
//                    locations should agree with the GlobalObservations vector.
//        GlobalObservations - The observed data at receivers.
//        mp - Pointer to object describing the parameterization of the material.
//
// Output: GlobalTimeSeries - The solution of the forward problem at the stations.
//         f                - The misfit.
//         dfs              - Gradient wrt to the source of misfit.
//         dfm              - Gradient wrt to the material of misfit.
//-----------------------------------------------------------------------
{
   int verbose = 0;
   vector<Source*> src(1);
   src[0] = GlobalSources[0]->copy(" ");

   // fetch all 11 source parameters, and merge in the unknowns
   double srcpars[11];
   src[0]->get_parameters( srcpars );
   set_source_pars( nspar, srcpars, xs );
   src[0]->set_parameters( srcpars );

// Translate one-dimensional parameter vector xm to material data (rho,mu,lambda)
   int ng = simulation.mNumberOfGrids;
   vector<Sarray> rho(ng), mu(ng), lambda(ng);
//New
   mp->get_material( nmpard, xm, nmpars, &xs[nspar], rho, mu, lambda );
// Old   
//   simulation.parameters_to_material( nmpar, xm, rho, mu, lambda );


// Run forward problem with guessed source, upred_saved,ucorr_saved are allocated
// inside solve_allpars. U and Um are final time solutions, to be used as 'initial' data
// when reconstructing U backwards.
   vector<DataPatches*> upred_saved(ng), ucorr_saved(ng);
   vector<Sarray> U(ng), Um(ng);
   simulation.solve_allpars( src, rho, mu, lambda, GlobalTimeSeries, U, Um, upred_saved, ucorr_saved );

// Compute misfit, 'diffs' will hold the source for the adjoint problem
   vector<TimeSeries*> diffs;
   for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
   {
      //      GlobalTimeSeries[m]->writeFile();
      TimeSeries *elem = GlobalTimeSeries[m]->copy( &simulation, "diffsrc" );
      diffs.push_back(elem);
   }
   f = 0;
   double dshift, ddshift, dd1shift;
   for( int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      f += GlobalTimeSeries[m]->misfit( *GlobalObservations[m], diffs[m], dshift, ddshift, dd1shift );

   double mftmp = f;
   MPI_Allreduce(&mftmp,&f,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   if( myrank == 0 && verbose >= 1 )
   {
      cout.precision(16);  
      cout << " Misfit is = " << f << endl;
   }
// Get gradient by solving the adjoint problem:
// Old
//   simulation.solve_backward_allpars( src, rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfs, nmpar, dfm );

// New
   double dfsrc[11];
   vector<Sarray> gRho(ng), gMu(ng), gLambda(ng);
   simulation.solve_backward_allpars( src, rho, mu, lambda,  diffs, U, Um, upred_saved, ucorr_saved, dfsrc, gRho, gMu, gLambda );
   get_source_pars( nspar, dfsrc, dfs );   
   mp->get_gradient( nmpard, xm, nmpars, &xs[nspar], &dfs[nspar], dfm, gRho, gMu, gLambda );
   
// diffs no longer needed, give back memory
   for( unsigned int m = 0 ; m < GlobalTimeSeries.size() ; m++ )
      delete diffs[m];
   diffs.clear();

// Give back memory
   for( unsigned int g=0 ; g < ng ; g++ )
   {
      delete upred_saved[g];
      delete ucorr_saved[g];
   }
   delete src[0];
}


//-----------------------------------------------------------------------
void restrict( int active[6], int wind[6], double* xm, double* xmi )
{
   int ni = (wind[1]-wind[0]+1);
   int nj = (wind[3]-wind[2]+1);
   int nia= (active[1]-active[0]+1);
   int nja= (active[3]-active[2]+1);

   for( int k=wind[4] ; k<=wind[5] ; k++ )
      for( int j=wind[2] ; j<=wind[3] ; j++ )
	 for( int i=wind[0] ; i<=wind[1] ; i++ )
	 {
            size_t indi = (i-wind[0]) + ni*(j-wind[2]) + ni*nj*(k-wind[4]);
	    size_t ind  = (i-active[0]) + nia*(j-active[2]) + nia*nja*(k-active[4]);
	    xmi[3*indi]   = xm[3*ind];
	    xmi[3*indi+1] = xm[3*ind+1];
	    xmi[3*indi+2] = xm[3*ind+2];
	 }
}

//-----------------------------------------------------------------------
void gradient_test( EW& simulation, vector<Source*>& GlobalSources, vector<TimeSeries*>& GlobalTimeSeries,
		   vector<TimeSeries*>& GlobalObservations, int nspar, int nmpars, double* xs, int nmpard, double* xm,
		   int myRank, MaterialParameterization* mp )
{
   int ns = nspar+nmpars;
   double* dfs;
   if( ns > 0 )
      dfs = new double[ns];
   double* dfm;
   if( nmpard > 0 )
      dfm = new double[nmpard];

   int sharedpars = 1;

   double f, fp;
      
   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
	   		        GlobalObservations, f, dfs, dfm, myRank, mp );
   if( myRank == 0 )
      cout << "Initial f = " << f << " " << endl;

   double h=1e-6;

   if( ns>0 )
   {
      if( myRank == 0 )
	 cout << "Gradient testing shared parameters :" << endl;
      for( int ind=0 ; ind < ns ; ind++ )
      {
	 xs[ind] += h;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fp, mp );
	 double dfnum = (fp-f)/h;
	 double dfan  = dfs[ind];
         double relerr = fabs(dfan-dfnum)/(fabs(dfan)+1e-10);
	 if( myRank == 0 && relerr > 1e-6 )
	 {
            cout << " ind = " << ind << "f = " << fp << " h= " << h << " dfan = " << dfan
		 << " dfnum = " << dfnum << " err = " << dfan-dfnum << endl;
	 }
	 xs[ind] -= h;
	 if( (ind % 100) && myRank == 0 )
	    cout << "Done ind = " << ind << endl;
      }
   }
   int nms, nmd, nmpard_global;
   mp->get_nr_of_parameters( nms, nmd, nmpard_global ) ;
   if( nmpard_global > 0 )
   {
      if( myRank == 0 )
	 cout << "Gradient testing distributed parameters :" << endl;
      for( size_t indg = 0 ; indg < nmpard_global ; indg++ )
      {
         ssize_t ind = mp->local_index(indg);
	 if( ind >=0 )
	    xm[ind] += h;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, fp, mp );
	 double dfnum = (fp-f)/h;
	 double dfan;
	 if( ind >=0 )
	    dfan = dfm[ind];

	 if( myRank == 0 )
	    cout << "f = " << fp << " h= " << h << " dfan = " << dfan << " dfnum = " << dfnum << " err = " << dfan-dfnum << endl;
         if( ind >= 0 )
	    xm[ind] -= h;
      }
   }
   if( ns > 0 )
      delete[] dfs;
   if( nmpard > 0 )
      delete[] dfm;
}

//-----------------------------------------------------------------------
void hessian_test( EW& simulation, vector<Source*>& GlobalSources, vector<TimeSeries*>& GlobalTimeSeries,
		   vector<TimeSeries*>& GlobalObservations, int nspar, int nmpars, double* xs, int nmpard, double* xm,
		   int myRank, MaterialParameterization* mp )
{
	      // Hessian test
   int ns = nspar+nmpars;
   double f;
   double* dfs;
   if( ns > 0 )
      dfs = new double[ns];
   double* dfm;
   if( nmpard > 0 )
      dfm = new double[nmpard];

   compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, f, dfs, dfm, myRank, mp );

   double* dfsp;
   if( ns > 0 )
      dfsp= new double[ns]; 

   double* dfmp;
   if( nmpard > 0 )
      dfmp= new double[nmpard]; 

   int sharedpars = 1;
   double h =1e-6;
   int grid = 0;

   if( sharedpars == 1 )
   {
      double* hess = new double[ns*ns];
      double fp;
      if( myRank == 0 )
	 cout << "Hessian computation of shared parameters :" << endl;
      
      for( int ind=0 ; ind < ns ; ind++ )
      {
	 xs[ind] += h;
	 compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			   GlobalObservations, fp, dfsp, dfmp, myRank, mp );
	 for( int p= 0 ; p < ns ; p++ )
	    hess[p+ns*ind] = (dfsp[p]-dfs[p])/h;
	 xs[ind] -= h;
         if( myRank == 0 )
	    cout << " done "  << ind << endl;
      }
      if( myRank == 0 )
      {
	 int fid = open( "hessdir/hessians.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	 if (fid == -1 )
	 {
	    VERIFY2(0, "ERROR: error opening file hessians.bin for writing header");
	    exit(-1);
	 }
	 int prec = 8;
	 size_t nr=write(fid,&prec,sizeof(int));
	 int npatch=1;
	 nr=write(fid,&npatch,sizeof(int));
	 int nc = 1;
	 nr=write(fid,&nc,sizeof(int));
	 double gz=simulation.mGridSize[0];
	 nr=write(fid,&gz,sizeof(double));
	 int dims[6]={1,ns,1,ns,1,1};
	 nr=write(fid,dims,6*sizeof(int));
	 nr=write(fid,hess,ns*ns*sizeof(double));
	 close(fid);
      }
   }
   else
   {
      int active[6], activeg[6];
      activeg[0] = simulation.m_iStartActGlobal[grid];
      activeg[1] = simulation.m_iEndActGlobal[grid];
      activeg[2] = simulation.m_jStartActGlobal[grid];
      activeg[3] = simulation.m_jEndActGlobal[grid];
      activeg[4] = simulation.m_kStartActGlobal[grid];
      activeg[5] = simulation.m_kEndActGlobal[grid];
      active[0]  = simulation.m_iStartAct[grid];
      active[1]  = simulation.m_iEndAct[grid];
      active[2]  = simulation.m_jStartAct[grid];
      active[3]  = simulation.m_jEndAct[grid];
      active[4]  = simulation.m_kStartAct[grid];
      active[5]  = simulation.m_kEndAct[grid];
      int wind[6];
	      // wind is intersection of 'active' and 'interior of proc.'
      for( int i=0 ; i < 6  ; i++ )
	 wind[i] = active[i];
      if( active[1] >= active[0] )
      {
	 if( wind[0] < simulation.m_iStartInt[grid] )
	    wind[0] = simulation.m_iStartInt[grid];
	 if( wind[1] > simulation.m_iEndInt[grid] )
	    wind[1] = simulation.m_iEndInt[grid];
	 if( wind[0] > wind[1] )
	    wind[1] = wind[0]-1;
      }
      if( active[3] >= active[2] )
      {
	 if( wind[2] < simulation.m_jStartInt[grid] )
	    wind[2] = simulation.m_jStartInt[grid];
	 if( wind[3] > simulation.m_jEndInt[grid] )
	    wind[3] = simulation.m_jEndInt[grid];
	 if( wind[2] > wind[3] )
	    wind[3] = wind[2]-1;
      }
      if( active[5] >= active[4] )
      {
	 if( wind[4] < simulation.m_kStartInt[grid] )
	    wind[4] = simulation.m_kStartInt[grid];
	 if( wind[5] > simulation.m_kEndInt[grid] )
	    wind[5] = simulation.m_kEndInt[grid];
	 if( wind[4] > wind[5] )
	    wind[5] = wind[4]-1;
      }
      int start[3]   ={wind[0]-activeg[0],wind[2]-activeg[2],wind[4]-activeg[4]};
      int locsize[3] ={wind[1]-wind[0]+1,wind[3]-wind[2]+1,wind[5]-wind[4]+1};
      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
              
	      //              cout << myRank << " " << activeg[2] << " " << activeg[3] << " " <<  active[2] << " " << active[3] <<  " " << wind[2] << " " << wind[3] << endl;
     //              cout << myRank << " " << simulation .m_jStart[grid] << " " << simulation.m_jEnd[grid] << endl;
	      //              int start[3]   ={active[0]-activeg[0],active[2]-activeg[2],active[4]-activeg[4]};
	      //	      int locsize[3] ={active[1]-active[0]+1,active[3]-active[2]+1,active[5]-active[4]+1};
	      //	      int globsize[3]={activeg[1]-activeg[0]+1,activeg[3]-activeg[2]+1,activeg[5]-activeg[4]+1};
      int nptsbuf = 1000000;
      int iwrite = myRank == 0 || ( myRank % 8 == 0 );

      Parallel_IO* pio = new Parallel_IO(1,1,globsize,locsize,start,nptsbuf,0);
      int fid;
      if( myRank == 0 )
      {
		 // create file, write header
	 fid = open( "hessdir/hessian.bin", O_CREAT | O_TRUNC | O_WRONLY, 0660 ); 
	 if (fid == -1 )
	 {
	    VERIFY2(0, "ERROR: error opening file hessian.bin for writing header");
	    exit(-1);
	 }
	 int prec = 8;
	 size_t nr=write(fid,&prec,sizeof(int));
	 int npatch=1;
	 nr=write(fid,&npatch,sizeof(int));
	 int nc = 3;
	 nr=write(fid,&nc,sizeof(int));
	 double gz=simulation.mGridSize[0];
	 nr=write(fid,&gz,sizeof(double));
	 int dims[6]={activeg[0],globsize[0],activeg[2],globsize[1],activeg[4],globsize[2]};
	 nr=write(fid,dims,6*sizeof(int));
      }
      else
	 fid = open( "hessdir/hessian.bin", O_WRONLY );
      size_t offset = sizeof(double)+9*sizeof(int);
      size_t colsize = globsize[0]*((size_t)globsize[1])*globsize[2]*3;

	      // Verify IO
      int npts = (wind[1]-wind[0]+1)*(wind[3]-wind[2]+1)*(wind[5]-wind[4]+1);
      double* xmi;
      if( npts > 0 )
	 xmi = new double[3*npts];
	      //	          simulation.get_material_parameter( nmpar, xm );
	      //		  if( npts > 0 )
	      //		     restrict( active, wind, xm, xmi );
	      //			  pio->write_array( &fid, 3, xmi, offset, "double");
	      //			  close(fid);
	      //			  exit(0);
      for( int kper = activeg[4] ; kper <= activeg[5] ; kper++ )
	 for( int jper = activeg[2] ; jper <= activeg[3] ; jper++ )
	    for( int iper = activeg[0] ; iper <= activeg[1] ; iper++ )
	       for( int var = 0 ; var < 3 ; var++ )
	       {
	       // perturb material 
	       //	       simulation.perturb_mtrl(iper,jper,kper,h,grid,var);
	       //	       simulation.get_material_parameter( nmpard, xm );
		  ssize_t pind = mp->parameter_index(iper,jper,kper,grid,var);
		  if( pind >= 0 )
		     xm[pind] += h;
		  //	       mp->perturb_material(iper,jper,kper,grid,var,h,xs,xm);
		  compute_f_and_df( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
				 GlobalObservations, f, dfs, dfmp, myRank, mp );
		  for( int p= 0 ; p < nmpard ; p++ )
		     dfmp[p] = (dfmp[p]-dfm[p])/h;
			  // Save Hessian column
		  restrict( active, wind, dfmp, xmi );
		  pio->write_array( &fid, 3, xmi, offset, "double");
		  offset += colsize*sizeof(double);
		       // restore material 
	       //	       simulation.perturb_mtrl(iper,jper,kper,-h,grid,var);
	       //	       mp->perturb_material(iper,jper,kper,grid,var,-h,xs,xm);
		  if( pind >= 0 )
		     xm[pind] -= h;
		   
	       }
      close(fid);
      if( npts > 0 )
	 delete[] xmi;
   }
   if( nmpard > 0 )
   {
      delete[] dfm;
      delete[] dfmp;
   }
   if( ns > 0 )
   {
      delete[] dfs;
      delete[] dfsp;
   }
}


//-----------------------------------------------------------------------
void misfit_curve( int i, int j, int k, int var, double pmin, double pmax,
		   int npts, EW& simulation, MaterialParameterization* mp,
		   int nspar, int nmpars, double* xs, int nmpard, double* xm,
		   vector<Source*>& GlobalSources, vector<TimeSeries*>& GlobalTimeSeries,
		   vector<TimeSeries*>& GlobalObservations, int myRank )
{
   double* fcn = new double[npts];
   ssize_t ind=mp->parameter_index(i,j,k,0,var);
   double xoriginal = xs[ind];
   for( int m=0 ; m < npts ; m++ )
   {
      double p = pmin + static_cast<double>(m)/(npts-1)*(pmax-pmin);
      xs[ind] = xoriginal+p;
      double f;
      compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		 GlobalObservations, f, mp );
      fcn[m] = f;
   }
   if( myRank == 0 )
   {
      int fd=open("fsurf.bin", O_CREAT|O_TRUNC|O_WRONLY, 0660 );
      int dims=1;
      size_t nr = write(fd,&dims,sizeof(int));
      nr = write(fd,&npts,sizeof(int));
      nr = write(fd,&pmin,sizeof(double));
      nr = write(fd,&pmax,sizeof(double));
      nr = write(fd,fcn,npts*sizeof(double));
      close(fd);
   }
}

//-----------------------------------------------------------------------
void misfit_surface( int ix1, int jx1, int kx1, int ix2, int jx2, int kx2,
		     int varx1, int varx2, double pmin1, double pmax1,
		     double pmin2, double pmax2, int npts1, int npts2, EW& simulation,
		     MaterialParameterization* mp, int nspar, int nmpars,
		     double* xs, int nmpard, double* xm,
		     vector<Source*>& GlobalSources, vector<TimeSeries*>& GlobalTimeSeries,
		     vector<TimeSeries*>& GlobalObservations, int myRank )
{
   double* fcn = new double[npts1*npts2];
   ssize_t ind1=mp->parameter_index(ix1,jx1,kx1,0,varx1);
   ssize_t ind2=mp->parameter_index(ix2,jx2,kx2,0,varx2);
   double xoriginal1 = xs[ind1];
   double xoriginal2 = xs[ind2];
   for( int m1=0 ; m1 < npts1 ; m1++ )
   {
      for( int m2=0 ; m2 < npts2 ; m2++ )
      {
	 double p1 = pmin1 + static_cast<double>(m1)/(npts1-1)*(pmax1-pmin1);
	 xs[ind1] = xoriginal1+p1;
	 double p2 = pmin2 + static_cast<double>(m2)/(npts2-1)*(pmax2-pmin2);
	 xs[ind2] = xoriginal2+p2;
	 double f;
	 compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
		    GlobalObservations, f, mp );
	 fcn[m1+npts1*m2] = f;
      }
      if( myRank == 0 )
	 cout << "Done m = " << m1 << endl;
   }
   if( myRank == 0 )
   {
      int fd=open("fsurf.bin", O_CREAT|O_TRUNC|O_WRONLY, 0660 );
      int dims=2;
      size_t nr = write(fd,&dims,sizeof(int));
      nr = write(fd,&npts1,sizeof(int));
      nr = write(fd,&npts2,sizeof(int));
      nr = write(fd,&pmin1,sizeof(double));
      nr = write(fd,&pmax1,sizeof(double));
      nr = write(fd,&pmin2,sizeof(double));
      nr = write(fd,&pmax2,sizeof(double));
      nr = write(fd,fcn,npts1*npts2*sizeof(double));
      close(fd);
   }
}

//-----------------------------------------------------------------------
int start_minv( int argc, char **argv, string& input_file, int& myRank,
		int& nProcs )
{
  stringstream reason;

  // Initialize MPI...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // mpi2 adds on four more args...  [-p4pg, dir, -p4wd, dir]
  int mpi2args = 4;
  if (argc != 2 && argc != 3 && argc != (2+mpi2args) && argc != (3+mpi2args) )
    {
      reason << "Wrong number of args (1-2), not: " << argc-1 << endl; 
      
      if (myRank == 0)
      {
	for (int i = 0; i < argc; ++i)
	  cout << "Argv[" << i << "] = " << argv[i] << endl;

	usage(reason.str());
      }
// Stop MPI
      MPI_Finalize();
      return 1;
    }
  else if (strcmp(argv[1],"-v") == 0 )
  {
     if (myRank == 0)
        cout << ewversion::getVersionInfo() << endl;
// Stop MPI
     MPI_Finalize();
     return 2;
  }
  else
     if (argc == 1) 
     {
        reason  << "ERROR: ****No input file specified!" << endl;
        for (int i = 0; i < argc; ++i)
           reason << "Argv[" << i << "] = " << argv[i] << endl;
        
        if (myRank == 0) usage(reason.str());
// Stop MPI
	MPI_Finalize();
        return 1;
     }

  else
    input_file = argv[1];

  if (myRank == 0) 
  {
    cout << ewversion::getVersionInfo() << endl;
    cout << "Input file: " << input_file << endl;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  return 0;
}

//-----------------------------------------------------------------------
int main(int argc, char **argv)
{
  string fileName;
  int myRank, nProcs;
  int status = start_minv( argc, argv, fileName, myRank, nProcs );
  
  if( status == 0 )
  {
// Save the source description here
     vector<Source*> GlobalSources; 
// Save the time series here
     vector<TimeSeries*> GlobalTimeSeries;
     vector<TimeSeries*> GlobalObservations;

// make a new simulation object by reading the input file 'fileName'
     EW simulation(fileName, GlobalSources, GlobalObservations, true );

     if (!simulation.wasParsingSuccessful())
     {
	if (myRank == 0)
	   cout << "Error: there were problems parsing the input file" << endl;
	status = 1;
     }
     else
     {
// get the simulation object ready for time-stepping
	simulation.setupRun( GlobalSources );
	if (!simulation.isInitialized())
	{ 
	   if (myRank == 0)
	      cout << "Error: simulation object not ready for time stepping" << endl;
	   status=1;
	}
	else if( GlobalSources.size() != 1 )
	{
	   if (myRank == 0)
	      cout << "Source optmization only implemented for a single source" << endl;
	}
	else
	{
// Successful initialization

	   simulation.setQuiet(true);
	   //	   simulation.setQuiet(false);
// Make observations aware of the utc reference time, if set.
// Filter observed data if required
	   for( int m = 0; m < GlobalObservations.size(); m++ )
	   {
	      //	      simulation.set_utcref( *GlobalObservations[m] );
	      if( simulation.m_prefilter_sources && simulation.m_filter_observations )
	      {
		 GlobalObservations[m]->filter_data( simulation.m_filterobs_ptr );
		 GlobalObservations[m]->writeFile( "_fi" );
	      }
	   }

//  First copy observations to GlobalTimeSeries, later, solve will insert 
//  the simulation time step and start time into GlobalTimeSeries.
	   for( int m = 0; m < GlobalObservations.size(); m++ )
	   {
	      string newname = "_out";
	      TimeSeries *elem = GlobalObservations[m]->copy( &simulation, newname, true );
	      GlobalTimeSeries.push_back(elem);
	   }

	   if (myRank == 0)
	   {
	      cout << "Running sw4mopt on " <<  nProcs << " processors..." << endl
		   << "Writing output to directory: " 
		   << simulation.getOutputPath() << endl;
	   }

           int mtrlpar = 1;
// Configure optimizer 
           Mopt* mopt = new Mopt( &simulation );
	   mopt->parseInputFileOpt( fileName );

// Select material parameterization
           MaterialParameterization* mp = mopt->m_mp;

           int nmpars, nmpard, nmpard_global;
	   mp->get_nr_of_parameters( nmpars, nmpard, nmpard_global );

	   double* xm=NULL;
           if( nmpard > 0 )
	      xm = new double[nmpard];
	   
// Perturb material if gradient testing
//	   simulation.perturb_mtrl();

// Default initial guess material = the material given in the input file
//           int nmpar;


	   //	   simulation.get_nr_of_material_parameters( nmpar );
	   //           cout << "nmpar " << nmpar << endl;
	   //	   double* xm=NULL;
	   //           if( nmpar > 0 )
	   //	   {
	   //	      xm = new double[nmpar];
	   //	      simulation.get_material_parameter( nmpar, xm );
	   //	   }


// figure out how many parameters we need
	   int maxit, maxrestart, varcase=0, stepselection=0;
	   bool dolinesearch, fletcher_reeves=true, testing;
	   double tolerance;

	   simulation.get_cgparameters( maxit, maxrestart, tolerance, fletcher_reeves, stepselection,
					dolinesearch, varcase, testing );
           varcase = 4;
	   int nspar=11;
           if( varcase == 0 )
	      nspar = 11;
	   else if( varcase == 1 )
	      nspar = 10;
	   else if( varcase == 2 )
	      nspar = 9;
           else if( varcase == 4 )
	      nspar = 0;
	   else
	      if( myRank == 0 )
		 cout << "Error: varcase = " << varcase << " not defined for sw4mopt" << endl;

           int ns = nmpars + nspar;
	   double *xs = new double[ns];

// Default initial guess, the input source, stored in GlobalSources[0]
           double xspar[11];
	   GlobalSources[0]->get_parameters(xspar);
	   set_source_pars( nspar, xspar, xs );

// Initialize the material parameters
           //           mp->read_parameters("mtrl5x5x5.bin",nmpars,xs);
	   //           mp->parameters_from_basematerial( nmpard, xm, nmpars, xs );

           mp->get_parameters(nmpard,xm,nmpars,xs,simulation.mRho,simulation.mMu,simulation.mLambda );
	   mp->write_parameters("mtrlpar-init.bin",nmpars,xs);

// Scaling factors, for the first tests set equal to one.
	   double* sf  = NULL;
           double* sfm = NULL;

	   // Values for test problem
	   double rhoscale=1, muscale=3.24, lascale=66.95, fscale=400;
           mopt->get_scalefactors( rhoscale, muscale, lascale, fscale );
	   //           cout << "sfactors " << rhoscale << " " << muscale << " " << lascale << " " << fscale << endl;
           if( ns > 0 )
	   {
	      sf  = new double[ns];
              for( int i=0 ; i < ns ; i++ )
		 sf[i] = 1;
	      double isfscale=1/sqrt(fscale);
	      for( int i=0 ; i < nmpars ; i += 3 )
	      {
		 sf[nspar+i]   = isfscale*rhoscale;
		 sf[nspar+i+1] = isfscale*muscale;
		 sf[nspar+i+2] = isfscale*lascale;
	      }
	   }
           if( nmpard > 0 )
	   {
              sfm = new double[nmpard];
	      for( int i=0 ; i<nmpard ;i++)
		 sfm[i]=1;
	   //	   simulation.get_scale_factors( nmpar, sfm );
	   }

// Output source initial guess
	   if( myRank == 0 && nspar > 0 )
	   {
	      cout << "Initial source guess : \n";
              char* names[11] = {" x0 = "," y0 = ", " z0 = "," mxx = ", " mxy = ", " mxz = ",
			   " myy = ", " myz = ", " mzz = ", " t0 = ", " freq = " };
              for( int i=0 ; i < nspar ; i++ )
	      {
		 cout << names[i] << xs[i] ;
                 if ( (i+1) % 3 == 0 || i == nspar-1 )
		    cout << endl;
	      }
	   }

           if( mopt->m_opttest == 2 )
	   {
// Gradient_test compares the computed gradient with the gradient obtained by numerical differentiation.
              gradient_test( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nspar, nmpars, xs,
			    nmpard, xm, myRank, mp );
	   }
	   else if( mopt->m_opttest == 3 )
	   {
// Hessian_test outputs the Hessian computed by numerical differentiation.
              hessian_test( simulation, GlobalSources, GlobalTimeSeries, GlobalObservations, nspar, nmpars, xs,
			    nmpard, xm, myRank, mp );
	   }
	   else if( mopt->m_opttest == 4 )
	   {
// Compute and save a one dimensional cut through the objective function
              int npts = 100;
	      int ix=2, jx=2, kx=2, varx=2;
	      double pmin=-2,pmax=2;
              misfit_curve( ix, jx, kx, varx, pmin, pmax, npts, simulation, mp, nspar, nmpars, xs,
			     nmpard, xm, GlobalSources, GlobalTimeSeries, GlobalObservations, myRank );
	   }
	   else if( mopt->m_opttest == 5 )
	   {
// Compute and save a two dimensional cut through the objective function
              int ix1=2, jx1=2, kx1=2, ix2=2, jx2=2, kx2=2, varx1=1, varx2=2;
	      double pmin1=-0.5,pmax1=0.5,pmin2=-1,pmax2=1;
	      int npts1=30,npts2=30;
              misfit_surface( ix1, jx1, kx1, ix2, jx2, kx2, varx1, varx2, pmin1, pmax1,
			      pmin2, pmax2, npts1, npts2, simulation, mp, nspar, nmpars,
			      xs, nmpard, xm, GlobalSources, GlobalTimeSeries, GlobalObservations,
			      myRank );
	   }
           else if( mopt->m_opttest == 6 )
	   {
// Solve forward problem to generate synthetic data
              double f;
	      compute_f( simulation, nspar, nmpars, xs, nmpard, xm, GlobalSources, GlobalTimeSeries,
			 GlobalObservations, f, mp );
              for( int m=0 ; m < GlobalTimeSeries.size() ; m++ )
		 GlobalTimeSeries[m]->writeFile( );
	   }
           else if( mopt->m_opttest == 1 )
	   {
// Run optimizer (default)
	      int method, bfgs_m;
	      simulation.get_optmethod( method, bfgs_m );
	      lbfgs( simulation, nspar, nmpars, xs, sf, nmpard, xm, sfm, GlobalSources, GlobalTimeSeries,
		     GlobalObservations, bfgs_m, myRank, mp );
	   }
           else
	      if( myRank == 0 )
		 cout << "ERROR: m_opttest = " << mopt->m_opttest << " is not a valid choice" << endl;

	   if( myRank == 0 )
	   {
	      cout << "============================================================" << endl
		   << " sw4opt ( Material/Source estimation solver) finished! " << endl
		   << "============================================================" << endl;
	   }
	}
     }
  }
  else if( status == 1 )
     cout  << "============================================================" << endl
	   << "The execution on proc " << myRank << " was unsuccessful." << endl
	   << "============================================================" << endl;
  if( status == 2 )
     status = 0;
// Stop MPI
  MPI_Finalize();
  return 0;
  // Note: Always return 0, to avoid having one error message per process from LC slurmd.
  //  return status;
} 


   
