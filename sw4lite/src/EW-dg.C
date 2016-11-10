#include "sw4.h"

#include "EW.h"

#include <sstream>
#include <fstream>

#ifdef SW4_OPENMP
#include <omp.h>
#endif

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "Source.h"
#include "GridPointSource.h"
#include "CheckPoint.h"

#include "DG_FUNC.h"

void EW::set_dg_orders( int qu, int qv)
{
    m_qu = qu;
    m_qv = qv;
    m_nint = (qu+1);
}

//-----------------------------------------------------------------------
void EW::timeStepLoopdGalerkin()
{
        // DG timestepping routine
   float_sw4 time_begin_solve = MPI_Wtime();
#ifdef SW4_CROUTINES
   std::cout << "ERROR, DG requires compilation with fortran" << endl;
   exit(2);
#else   
#ifdef SW4_OPENMP
#pragma omp parallel
   {
      if( omp_get_thread_num() == 0 &&  m_myrank == 0  )
      {
	cout << endl << "***** Number of MPI-tasks: " << m_nprocs << " *******" << endl;
	 int nth=omp_get_num_threads();
	 cout <<  "****** Using OpenMP with " << nth << " thread";
	 if( nth > 1 )
	    cout << "s";
	 cout << " per MPI task *******" << endl<< endl;
      }
   }
#endif

    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];

        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    int ni = ilast-ifirst+1;
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    int nj = jlast-jfirst+1;
    if (m_global_nz[g] == klast){klast = klast-1;}
    int nk = klast-kfirst+1;
        // What order are we running with? 
    int q_u = m_qu;
    int q_v = m_qv;
    int n_int = m_nint;
    int n_int_id = 3*(q_u+1); // Initial data uses higher order quadrature.
    
    size_t mat_size;
    mat_size = ((size_t) 3)*(m_qu+1)*(m_qu+1)*(m_qu+1)*(ni*nj*nk);
    double * udg = new double[mat_size];
    double * vdg = new double[mat_size];
    double * updg = new double[mat_size];
    double * vpdg = new double[mat_size];
    double * utdg = new double[mat_size];
    double * vtdg = new double[mat_size];
    double * force_v = new double[mat_size];
    double * force_u = new double[mat_size];
    mat_size = ((size_t) 3)*(q_v+1)*(q_v+1)*(q_v+1)*(q_v+1)*(q_v+1)*(q_v+1)*3;
    double * MV = new double[mat_size];
    double * MU = new double[mat_size];
    double * SV = new double[mat_size];
    double * SU = new double[mat_size];
    mat_size = ((size_t) (n_int+1))*(n_int+1)*pow(q_v+1,3)*3*6;
    double * LV = new double[mat_size];
    mat_size = ((size_t) (n_int+1))*(n_int+1)*pow(q_u+1,3)*3*3*6;
    double * LU = new double[mat_size];
    mat_size = ((size_t) n_int+1)*(n_int+1)*3*6*ni*nj*nk;
    double * w_in_all_faces = new double[mat_size];
    double * w_out_all_faces = new double[mat_size];
    double * w_star_all_faces = new double[mat_size];
    double * v_in_all_faces = new double[mat_size];
    double * v_out_all_faces = new double[mat_size];
    double * v_star_all_faces = new double[mat_size];

        // Arrays for communication
    size_t xside_size, yside_size;
        // X-sides 2 sides x 2 fields
    xside_size = ((size_t) n_int+1)*(n_int+1)*3*6*nj*nk*2;
    double * x_out_b = new double[xside_size];
    double * x_out_e = new double[xside_size];
    double * x_in_b = new double[xside_size];
    double * x_in_e = new double[xside_size];
        // Y-sides 2 sides x 2 fields
    yside_size = ((size_t) n_int+1)*(n_int+1)*3*6*ni*nk*2;
    double * y_out_b = new double[yside_size];
    double * y_out_e = new double[yside_size];
    double * y_in_b = new double[yside_size];
    double * y_in_e = new double[yside_size];
    MPI_Status status;
    int xtag1 = 345;
    int xtag2 = 346;
    int ytag1 = 347;
    int ytag2 = 348;
    
        // On exit MU and MV has been replaced by inv(MU) and inv(MV)
    assemble(MU,MV,SU,SV,LU,LV);
        // Initial data
    initialData(udg,vdg);
    
    double t = 0.0;
        // Compute error in initial data

    if( m_myrank == 0 ){
        cout << endl << "****** Computing the error in initial data ******" << endl;
    }
    computeError(udg,vdg,t);
    double mu, lambda,rho;
        // Assume constant material, sample it in middle of domain
    int imid = (ifirst+ilast)/2;
    int jmid = (jfirst+jlast)/2; 
    int kmid = (kfirst+klast)/2; 
    mu = mMu[g](imid,jmid,kmid);
    lambda = mLambda[g](imid,jmid,kmid);
    rho = mRho[g](imid,jmid,kmid);
        // Compute timesteps
    double dt = mCFL/(sqrt((2.*mu+lambda)/rho)*(q_u+1.5))*h;
    int nsteps = ((int) mTmax/dt)+1;
    dt = mTmax/(1.0*nsteps);
    mDt = dt;
    if( m_myrank == 0 ){cout << endl << "****** CFL= " << mCFL << " time step=" << mDt << endl;}

    double rec_coord[3];
    rec_coord[0] = 5.0+0.5*h;
    rec_coord[1] = 5.0+0.5*h;
    rec_coord[2] = 6.0+0.5*h;
        // Look for a recorder
        // Is this recorder inside this processor?
        // We always push into + quadrant if we are on a face/edge/vertex
    double tol = 1e-10;
    int inx,iny,inz;
    double r1,r2,r3;
    double rtmp;
    inx = rec_coord[0]/h;
    rtmp = rec_coord[0];
    if (fabs(rtmp - inx*h) < tol) {
        rtmp = rtmp + 0.1*h;
    }
    inx = rtmp/h+1;
    iny = rec_coord[1]/h;
    rtmp = rec_coord[1];
    if (fabs(rtmp - iny*h) < tol) {
        rtmp = rtmp + 0.1*h;
    }
    iny = rtmp/h+1;
    inz = rec_coord[2]/h;
    rtmp = rec_coord[2];
    if (fabs(rtmp - inz*h) < tol) {
        rtmp = rtmp + 0.1*h;
    }
    inz = rtmp/h+1;
    int got_rec = 0;
    
    double * px_rec = new double[q_u+1];
    double * py_rec = new double[q_u+1];
    double * pz_rec = new double[q_u+1];
    double * prx_rec = new double[q_u+1];
    double * pry_rec = new double[q_u+1];
    double * prz_rec = new double[q_u+1];
    int ncomponents = 3;
    double urec[ncomponents];
    double vrec[ncomponents];
    
    ofstream myfile;    
    if ( inx >= ifirst && inx <= ilast &&
         iny >= jfirst && iny <= jlast){
        got_rec = 1;
        char fn [100];
	snprintf (fn, sizeof fn, "error_%05d.txt", klast);
        myfile.open (fn);

        r1 = 2.0*((rec_coord[0] - (inx-1.0)*h)/h)-1.0;
        r2 = 2.0*((rec_coord[1] - (iny-1.0)*h)/h)-1.0;
        r3 = 2.0*((rec_coord[2] - (inz-1.0)*h)/h)-1.0;
        cout << endl << "****** Recorder found in proc " << m_myrank << "  ******" << endl;
        cout.precision(7);
	cout << fixed << "x,y,z computed location : " << ((inx-0.5)+0.5*r1)*h << " "
	     << ((iny-0.5)+0.5*r2)*h << " "<<  ((inz-0.5)+0.5*r3)*h << endl << endl ;
        eval_legendre(px_rec,prx_rec,&q_u,&r1);
        eval_legendre(py_rec,pry_rec,&q_u,&r2);
        eval_legendre(pz_rec,prz_rec,&q_u,&r3);
        get_recorder(udg,vdg,urec,vrec,
                     &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u,
                     px_rec,prx_rec,py_rec,pry_rec,pz_rec,prz_rec,
                     &inx,&iny,&inz,&ncomponents);
    }

    int got_src = 0;    
    double * px_src;
    double * py_src;
    double * pz_src;
    double * prx_src;
    double * pry_src;
    double * prz_src;
    double * point_src;
    double src_coord[3], f_amp[3];
    int inxs,inys,inzs;

    if(m_globalUniqueSources.size() > 0){
        vector<Source*>& a_sources = m_globalUniqueSources;
        Source source = *a_sources[0];
        src_coord[0] = source.getX0();
        src_coord[1] = source.getY0();
        src_coord[2] = source.getZ0();
        source.getForces(f_amp[0],f_amp[1],f_amp[2]);
            // Look for a src
            // Is this src inside this processor?
            // We always push into + quadrant if we are on a face/edge/vertex
        inxs = src_coord[0]/h;
        rtmp = src_coord[0];
        if (fabs(rtmp - inxs*h) < tol) {
            rtmp = rtmp + 0.1*h;
        }
        inxs = rtmp/h+1;
        inys = src_coord[1]/h;
        rtmp = src_coord[1];
        if (fabs(rtmp - inys*h) < tol) {
            rtmp = rtmp + 0.1*h;
        }
        inys = rtmp/h+1;
        inzs = src_coord[2]/h;
        rtmp = src_coord[2];
        if (fabs(rtmp - inzs*h) < tol) {
            rtmp = rtmp + 0.1*h;
        }
        inzs = rtmp/h+1;
        px_src = new double[q_u+1];
        py_src = new double[q_u+1];
        pz_src = new double[q_u+1];
        prx_src = new double[q_u+1];
        pry_src = new double[q_u+1];
        prz_src = new double[q_u+1];
        point_src = new double[(q_v+1)*(q_v+1)*(q_v+1)];
        if ( inxs >= ifirst && inxs <= ilast &&
             inys >= jfirst && inys <= jlast){
            got_src = 1;
            cout << endl << "****** Source found in proc " << m_myrank << "  ******" << endl;
            cout << "Source coordinates x0: " << src_coord[0]
        	 << " y0: " << src_coord[1] << " z0: "
        	 << src_coord[2] << endl;
            r1 = 2.0*((src_coord[0] - (inxs-1.)*h)/h)-1.0;
            r2 = 2.0*((src_coord[1] - (inys-1.)*h)/h)-1.0;
            r3 = 2.0*((src_coord[2] - (inzs-1.)*h)/h)-1.0;
            cout << fixed << "r1, r1 and r3 local coord : " << r1 << " "<< r2 << " "<< r3 << endl << endl ;
            eval_legendre(px_src,prx_src,&q_u,&r1);
            eval_legendre(py_src,pry_src,&q_u,&r2);
            eval_legendre(pz_src,prz_src,&q_u,&r3);
            get_dirac_source(point_src,px_src,py_src,pz_src,&q_u);
        }
    }

    double df;
        //int ntay = 4*(((int) (q_u-1)/4)+1);
    int ntay = max(4,q_u+1);
    if( m_myrank == 0 ){
        cout << "Starting DG-solver Using " << ntay << " time derivatives and degree "
             << q_u << " polynomials."<<  endl;
    }

    int nsrc;
    double * tg;
    double * ct;
    double * fsrc;
    double * source_tay_coeff;
    
    if(got_src == 1){
            //DEAA What can we get away with here?
        nsrc = ntay+3;
        tg = new double[nsrc+1];
        ct = new double[(nsrc+1)*(ntay+1)];
        fsrc = new double[nsrc+1];
        source_tay_coeff = new double[nsrc+1];
        
        set_tay_weights(tg,ct,&nsrc,&ntay);
        get_source_tay_coeff(source_tay_coeff,tg,ct,&nsrc,&ntay,&t,&dt);
    }

        // Do we need to set boundary conditions?
    int sbx_b,sbx_e,sby_b,sby_e;
    sbx_b = (ifirst == 1) ? 1 : 0;                
    sbx_e = (ilast == (m_global_nx[g]-1)) ? 1 : 0;
    sby_b = (jfirst == 1) ? 1 : 0;                
    sby_e = (jlast == (m_global_ny[g]-1)) ? 1 : 0;
    
    for (int it = 1; it <= nsteps; it++){
        t = (it-1)*dt;
        if (got_src == 1){
             get_source_tay_coeff(source_tay_coeff,tg,ct,&nsrc,&ntay,&t,&dt);
         }
        df = dt;
            // Start the update to u(t+dt) = updg = dt*utdg+dt^2/2* (d\dt)(utdg) ...
        start_next_step(updg,vpdg,udg,vdg,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u);
        
        for (int itay = 1; itay <= ntay; itay++){

            build_w_and_v(udg,vdg,w_in_all_faces,v_in_all_faces);

            pass_outside_fluxes(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,
                                &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int);
                // Copy the data that is to be sent
            get_comm_sides(x_in_b,x_in_e,y_in_b,y_in_e,v_in_all_faces,w_in_all_faces,
                           &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int);
                // X-direction communication
            MPI_Sendrecv(x_in_e,xside_size,MPI_DOUBLE, m_neighbor[1], xtag1,
                         x_out_b,xside_size, MPI_DOUBLE, m_neighbor[0], xtag1,
                         m_cartesian_communicator, &status);
            MPI_Sendrecv(x_in_b,xside_size,MPI_DOUBLE, m_neighbor[0], xtag2,
                         x_out_e,xside_size, MPI_DOUBLE, m_neighbor[1], xtag2,
                         m_cartesian_communicator, &status);
                // Y-direction communication
            MPI_Sendrecv(y_in_e,yside_size,MPI_DOUBLE, m_neighbor[3], ytag1,
                         y_out_b,yside_size, MPI_DOUBLE, m_neighbor[2], ytag1,
                         m_cartesian_communicator, &status);
            MPI_Sendrecv(y_in_b,yside_size,MPI_DOUBLE, m_neighbor[2], ytag2,
                         y_out_e,yside_size, MPI_DOUBLE, m_neighbor[3], ytag2,
                         m_cartesian_communicator, &status);
            put_comm_sides(x_out_b,x_out_e,y_out_b,y_out_e,v_out_all_faces,w_out_all_faces,
                           &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int);
                // Overwrite with boundary conditions where needed
            set_boundary_conditions(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,
                                    &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int,
                                    &sbx_b,&sbx_e,&sby_b,&sby_e,&m_dg_bc_type);
                // DEAA: This needs to take the flux parameters as input.
            numerical_fluxes(w_star_all_faces,v_star_all_faces,w_in_all_faces,v_in_all_faces,
                             w_out_all_faces,v_out_all_faces);
            compute_surface_integrals(v_in_all_faces,v_star_all_faces,w_star_all_faces,
                                      force_u,force_v,LU,LV,&h,&q_u,&q_v,&n_int,
                                      &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast);
             if (got_src == 1){
                     // Add forcing
                 double stc = source_tay_coeff[itay-1];
                 add_dirac_source(force_v,point_src,f_amp,&stc,
                                  &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,
                                  &inxs,&inys,&inzs);
             }
            compute_time_derivatives(utdg,vtdg,udg,vdg,force_u,force_v,
                                     MU,MV,SU,SV,
                                     &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,
                                     &q_v,&q_u);
            taylor_swap(utdg,vtdg,updg,vpdg,udg,vdg,&df,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u);
            df = df*dt/(1+itay);
        }
            // Swap the new into the old
        start_next_step(udg,vdg,updg,vpdg,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u);
        t = t+dt;
        if(got_rec == 1){
            get_recorder(udg,vdg,urec,vrec,
                         &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u,
                         px_rec,prx_rec,py_rec,pry_rec,pz_rec,prz_rec,
                         &inx,&iny,&inz,&ncomponents);
            cout.precision(16);
            double ups[3];
            double x_err = rec_coord[0];
            double y_err = rec_coord[1];
            double z_err = rec_coord[2];
            get_exact_point_source_dG(ups,t,x_err,y_err,z_err);
            myfile << fixed << t << " " << ups[2] << " " << urec[2] << " " << ups[2]-urec[2] << endl;
            char fn [100];
            snprintf (fn, sizeof fn, "ul_%05d_%05d.txt",klast, it);

            ofstream linefile;    
            linefile.open (fn);
            linefile.precision(16);
                // record at the center of each element along the source.
            for (int iline = kfirst; iline <= klast; iline++){
                x_err = (inx-0.5)*h;
                y_err = (iny-0.5)*h;
                z_err = (iline-0.5)*h;
                get_exact_point_source_dG(ups,t,x_err,y_err,z_err);
                get_recorder(udg,vdg,urec,vrec,
                             &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&q_u,
                             px_rec,prx_rec,py_rec,pry_rec,pz_rec,prz_rec,
                             &inx,&iny,&iline,&ncomponents);
                linefile << fixed << " " << urec[0] << " " << urec[1] << " " << urec[2] <<
                    " " << ups[0] << " " << ups[1] << " " << ups[2] << endl;
            }
            linefile.close();
        }
        
    }

    if(m_myrank == 0){cout << "checking errors at final time: " << mTmax <<  endl;}
    computeError(udg,vdg,t);
    if( m_myrank == 0 ){
        float_sw4 time_end_solve = MPI_Wtime();
        cout << endl << "***** Solver execution time: " << time_end_solve - time_begin_solve <<
            " seconds ******" << endl;
    }

    delete [] udg,vdg,utdg,vtdg,updg,vpdg,MU,MV,SU,SV,LU,LV,force_u,force_v;
    delete [] x_out_b,x_out_e,y_out_b,y_out_e,x_in_b,x_in_e,y_in_b,y_in_e;
    delete [] w_in_all_faces,w_out_all_faces,w_star_all_faces,v_in_all_faces,
        v_out_all_faces,v_star_all_faces;

    delete [] px_rec, py_rec, pz_rec;
    delete [] prx_rec, pry_rec, prz_rec;
    if(got_rec == 1)
        myfile.close();


#endif
}


void EW::assemble(double* MU,double* MV,double* SU,double* SV,double* LU,double* LV)
{


    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];
        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    if (m_global_nz[g] == klast){klast = klast-1;}
    int n_int = m_nint;
    int q_u = m_qu;
    int q_v = m_qv;

        // Currently we only support constant coeff.
    double mu, lambda,rho;
        // Assume constant material, sample it in middle of domain
    int imid = (ifirst+ilast)/2;
    int jmid = (jfirst+jlast)/2; 
    int kmid = (kfirst+klast)/2; 
    mu = mMu[g](imid,jmid,kmid);
    lambda = mLambda[g](imid,jmid,kmid);
    rho = mRho[g](imid,jmid,kmid);
#ifdef SW4_CROUTINES
    std::cout << "ERROR, DG requires compilation with fortran" << endl;
    exit(2);
#else
    assemble_const_coeff(MU,MV,SU,SV,LU,LV,&q_u,&q_v,&n_int,&h,&lambda,&mu,&rho);
    factor(MU,MV,&q_u,&q_v);
#endif
}


//---------------------------------------------------------------------------
void EW::computeError(double* udg, double* vdg, double t)
{

    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];
        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    if (m_global_nz[g] == klast){klast = klast-1;}
    int n_int = m_nint;
    int q_u = m_qu;
    int q_v = m_qv;
    
    double * l2_err = new double[3];
    double * l2_err_tmp = new double[3];
    if(m_single_mode_problem == 1){
#ifdef SW4_CROUTINES
       std::cout << "ERROR, DG requires compilation with fortran" << endl;
       exit(2);
#else
       compute_single_mode_error(l2_err_tmp,udg,&ifirst,&ilast,&jfirst,&jlast,
                                  &kfirst,&klast,&h,&t,&q_u,&n_int);
#endif
        MPI_Allreduce( l2_err_tmp, l2_err, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if( m_myrank == 0 ){
            cout << "checking errors at time: " << t <<  endl;
            cout << "l2_err 1 " << sqrt(l2_err[0]) <<  endl;
            cout << "l2_err 2 " << sqrt(l2_err[1]) <<  endl;
            cout << "l2_err 3 " << sqrt(l2_err[2]) <<  endl;
        }
    }
    else{
            // Assume constant material, sample it in middle of domain
        int g = 0;
        int imid = (m_iStartInt[g]+m_iEndInt[g])/2;
        int jmid = (m_jStartInt[g]+m_jEndInt[g])/2;
        int kmid = (m_kStartInt[g]+m_kEndInt[g])/2;
        double rho   = mRho[g](imid,jmid,kmid);
        double beta  =  sqrt( mMu[g](imid,jmid,kmid)/rho);
        double alpha =  sqrt( (2*mMu[g](imid,jmid,kmid)+mLambda[g](imid,jmid,kmid))/rho);
        
        double x0 = 5.+h/2.;
        double y0 = 5.+h/2.;
        double z0 = 5.+h/2.;
        double fr = 1.0;
        double t0 = 0.0;
        double fx = 0.0;
        double fy = 0.0;
        double fz = 1.0;
        double parameters[11];
        parameters[0] = rho;
        parameters[1] = beta;
        parameters[2] = alpha;
        parameters[3] = x0;
        parameters[4] = y0;
        parameters[5] = z0;
        parameters[6] = fr;
        parameters[7] = t0;
        parameters[8] = fx;
        parameters[9] = fy;
        parameters[10] = fz;
#ifdef SW4_CROUTINES
	std::cout << "ERROR, DG requires compilation with fortran" << endl;
	exit(2);
#else
        compute_point_dirac_error(l2_err_tmp,udg,&ifirst,&ilast,&jfirst,&jlast,
                                  &kfirst,&klast,&h,&t,&q_u,&n_int,parameters);
#endif
        MPI_Allreduce( &l2_err_tmp[1], &l2_err[1], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( &l2_err_tmp[0], &l2_err[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if( m_myrank == 0 ){
            cout << "checking errors at time: " << t <<  endl;
            cout << "Max-err in u_z " << l2_err[0] <<  endl;
            cout << "l2_err in u_y " << sqrt(l2_err[1]) <<  endl;
            cout << "l2_err in u_z " << sqrt(l2_err[2]) <<  endl;
        }
    }

}



void EW::numerical_fluxes(double* w_star_all_faces, double* v_star_all_faces,
                          double* w_in_all_faces, double* v_in_all_faces,
                          double* w_out_all_faces, double* v_out_all_faces)
{

    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];
        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    if (m_global_nz[g] == klast){klast = klast-1;}
    int n_int = m_nint;
    int q_u = m_qu;
    int q_v = m_qv;

    double mu, lambda,rho;
        // Assume constant material, sample it in middle of domain
    int imid = (ifirst+ilast)/2;
    int jmid = (jfirst+jlast)/2; 
    int kmid = (kfirst+klast)/2; 
    mu = mMu[g](imid,jmid,kmid);
    lambda = mLambda[g](imid,jmid,kmid);
    rho = mRho[g](imid,jmid,kmid);
#ifdef SW4_CROUTINES
    std::cout << "ERROR, DG requires compilation with fortran" << endl;
    exit(2);
#else
    compute_numerical_fluxes(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,
                             v_star_all_faces,w_star_all_faces,
                             &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int,&lambda,&mu,&rho);
#endif
}


//---------------------------------------------------------------------------
void EW::initialData(double * udg, double * vdg)
{

    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];
        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    if (m_global_nz[g] == klast){klast = klast-1;}
    int n_int_id = 3*m_nint;
    int q_u = m_qu;
    int q_v = m_qv;

    int id_type = 0; // Default to homogenous data
    if(m_single_mode_problem == 1){
        id_type=1;
    }
#ifdef SW4_CROUTINES
    std::cout << "ERROR, DG requires compilation with fortran" << endl;
    exit(2);
#else
    get_initial_data(udg,vdg,&ifirst,&ilast,&jfirst,&jlast,
                     &kfirst,&klast,&h,&q_v,&q_u,&n_int_id,&id_type);
#endif
}

//---------------------------------------------------------------------------
void EW::build_w_and_v(double* udg, double* vdg, double* w_in_all_faces, double* v_in_all_faces)
{

    int g = 0;
    double h = mGridSize[g];
    int ifirst = m_iStartInt[g], ilast = m_iEndInt[g];
    int jfirst = m_jStartInt[g], jlast = m_jEndInt[g];
    int kfirst = m_kStartInt[g], klast = m_kEndInt[g];
        // Throw away the last gridpoint as number of elements are one less.
    if (m_global_nx[g] == ilast){ilast = ilast-1;}
    if (m_global_ny[g] == jlast){jlast = jlast-1;}
    if (m_global_nz[g] == klast){klast = klast-1;}
    int n_int = m_nint;
    int q_u = m_qu;
    int q_v = m_qv;

        // Currently we only support constant coeff.
    double mu, lambda,rho;
        // Assume constant material, sample it in middle of domain
    int imid = (ifirst+ilast)/2;
    int jmid = (jfirst+jlast)/2; 
    int kmid = (kfirst+klast)/2; 
    mu = mMu[g](imid,jmid,kmid);
    lambda = mLambda[g](imid,jmid,kmid);
    rho = mRho[g](imid,jmid,kmid);
#ifdef SW4_CROUTINES
    std::cout << "ERROR, DG requires compilation with fortran" << endl;
    exit(2);
#else
    build_my_v_const_coeff(vdg,v_in_all_faces,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&n_int);
    build_my_w_const_coeff(udg,w_in_all_faces,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&h,&q_u,&n_int,&lambda,&mu);
#endif
}

void EW::processdGalerkin(char* buffer)
{
    char* token = strtok(buffer, " \t");
    token = strtok(NULL, " \t");
    string err = "ERROR in processGalerkin: ";
    while (token != NULL)
    {
            // while there are still tokens in the string
        if (startswith("#", token) || startswith(" ", buffer))
                // Ignore commented lines and lines with just a space.
            break;
        if (startswith("order=", token))
        {
            token += 6; // skip order=
            CHECK_INPUT(atoi(token) >= 1, err << "order must be greater than 0: " << token);
            int qorder = atoi(token);
            set_dg_orders(qorder,qorder);
        }
        else if (startswith("singlemode=", token))
        {
            token += 11; // skip singlemode=
            CHECK_INPUT(atoi(token) >= 0, err << "Single mode should be 0 or 1: " << token);
            int sm = atoi(token);
            if(sm == 1){
                m_single_mode_problem = 1;
                m_dg_bc_type = 1;
                double rho = 1.0;
                double mu  = 1.0;
                double la  = 2.0;
                int g = 0;
                    mRho[g].set_value(rho);
                    mMu[g].set_value(mu);
                    mLambda[g].set_value(la);
            }
            else{
                m_single_mode_problem = 0;
                m_dg_bc_type = 0;
                double rho = 1.0;
                double mu  = 1.0;
                double la  = 1.0;
                int g = 0;
                mRho[g].set_value(rho);
                mMu[g].set_value(mu);
                mLambda[g].set_value(la);
            }
        }
        else
        {
            badOption("dgalerkin", token);
        }
        token = strtok(NULL, " \t");
    }
}

//-----------------------------------------------------------------------
void EW::get_exact_point_source_dG(double* u, double t, double x, double y, double z)
{
        // Assume constant material, sample it in middle of domain
    float_sw4 rho   = 1.0;
    float_sw4 beta  =  1.0;
    float_sw4 alpha =  sqrt(3.0);
    float_sw4 h   = mGridSize[0];
    
    float_sw4 x0 = 0.5*h+5.;
    float_sw4 y0 = 0.5*h+5.;
    float_sw4 z0 = 0.5*h+5.;
    float_sw4 fr = 1.0;
    float_sw4 time = t;
    float_sw4 fx = 0.0;
    float_sw4 fy = 0.0;
    float_sw4 fz = 1.0;
    float_sw4 eps = 1e-3*h;
    float_sw4 R = sqrt( (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) );
    float_sw4 A, B;
    if( R < eps )
        u[0] = u[1] = u[2] = 0.0;
    else
    {
            //else if (tD == iC6SmoothBump)
        A = ( 1/pow(alpha,2) * C6SmoothBump(time, fr*R, alpha) - 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) +
              3/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R*R*R)  ;
	
        B = ( 1/pow(beta,2) * C6SmoothBump(time, fr*R, beta) -
              1/pow(fr*R,2) * C6SmoothBump_x_T_Integral(time, fr*R, alpha, beta) ) / (4*M_PI*rho*R) ;
    }
    u[0] = ( (x - x0)*(x - x0)*fx + (x - x0)*(y - y0)*fy + (x - x0)*(z - z0)*fz )*A + fx*B;
    u[1] = ( (y - y0)*(x - x0)*fx + (y - y0)*(y - y0)*fy + (y - y0)*(z - z0)*fz )*A + fy*B;
    u[2] = ( (z - z0)*(x - x0)*fx + (z - z0)*(y - y0)*fy + (z - z0)*(z - z0)*fz )*A + fz*B;
}

