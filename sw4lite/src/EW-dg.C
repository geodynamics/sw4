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
    double df;
        //int ntay = 4*(((int) (q_u-1)/4)+1);
    int ntay = max(4,q_u+1);
    if( m_myrank == 0 ){
        cout << "Starting DG-solver Using " << ntay << " time derivatives and degree "
             << q_u << " polynomials."<<  endl;
    }
        // Do we need to set boundary conditions?
    int sbx_b,sbx_e,sby_b,sby_e;
    sbx_b = (ifirst == 1) ? 1 : 0;                
    sbx_e = (ilast == (m_global_nx[g]-1)) ? 1 : 0;
    sby_b = (jfirst == 1) ? 1 : 0;                
    sby_e = (jlast == (m_global_ny[g]-1)) ? 1 : 0;
    
    for (int it = 1; it <= nsteps; it++){
        t = (it-1)*dt;
        // if (got_src == 1){
        //     get_source_tay_coeff(source_tay_coeff,tg,ct,&nsrc,&ntay,&t,&dt);
        // }
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
                // if (got_src == 1){
                //         // Add forcing
                //     double stc = source_tay_coeff[itay-1];
                //     add_dirac_source(force_v,point_src,f_amp,&stc,
                //                      &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,
                //                      &inxs,&inys,&inzs);
                // }
            
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
    assemble_const_coeff(MU,MV,SU,SV,LU,LV,&q_u,&q_v,&n_int,&h,&lambda,&mu,&rho);
    factor(MU,MV,&q_u,&q_v);

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
        compute_single_mode_error(l2_err_tmp,udg,&ifirst,&ilast,&jfirst,&jlast,
                                  &kfirst,&klast,&h,&t,&q_u,&n_int);
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
        compute_point_dirac_error(l2_err_tmp,udg,&ifirst,&ilast,&jfirst,&jlast,
                                  &kfirst,&klast,&h,&t,&q_u,&n_int,parameters);
        MPI_Allreduce( &l2_err_tmp[1], &l2_err[1], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( &l2_err_tmp[0], &l2_err[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if( m_myrank == 0 ){
            cout << "checking errors at time: " << t <<  endl;
            cout << "l2_err 1 " << l2_err[0] <<  endl;
            cout << "l2_err 2 " << sqrt(l2_err[1]) <<  endl;
            cout << "l2_err 3 " << sqrt(l2_err[2]) <<  endl;
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
    compute_numerical_fluxes(v_out_all_faces,v_in_all_faces,w_out_all_faces,w_in_all_faces,
                             v_star_all_faces,w_star_all_faces,
                             &ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&n_int,&lambda,&mu,&rho);

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
    get_initial_data(udg,vdg,&ifirst,&ilast,&jfirst,&jlast,
                     &kfirst,&klast,&h,&q_v,&q_u,&n_int_id,&id_type);

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

    build_my_v_const_coeff(vdg,v_in_all_faces,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&q_v,&n_int);
    build_my_w_const_coeff(udg,w_in_all_faces,&ifirst,&ilast,&jfirst,&jlast,&kfirst,&klast,&h,&q_u,&n_int,&lambda,&mu);

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
                for( int g=0 ; g < mNumberOfGrids ; g++ )
                {
                    mRho[g].set_value(rho);
                    mMu[g].set_value(mu);
                    mLambda[g].set_value(la);
                }
            }
        }
        else
        {
            badOption("dgalerkin", token);
        }
        token = strtok(NULL, " \t");
    }
}


