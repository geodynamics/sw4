
extern "C" {
    void assemble(double* MU,double* MV,double* SU,double* SV,double* LU,double* LV,int* q_u,int* q_v,int* nint,double* h);
    void assemble_const_coeff(double* MU,double* MV,double* SU,double* SV,double* LU,double* LV,
                              int* q_u,int* q_v,int* nint,double* h, double* lambda_lame, double* mu_lame, double* rho);
    void factor(double* MU,double* MV,int* q_u,int* q_v);
    void compute_surface_integrals(double* v_in_all_faces,double* v_star_all_faces,double* w_star_all_faces,
                                   double* force_u,double* force_v,double* LU,double* LV,double* h,int* q_u,int* q_v,
                                   int* nint,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast);
    void get_initial_data(double* udg,double* vdg,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                          double* h,int* q_v,int* q_u,int* nint_id, int* id_type);
    void build_my_v(double* vdg,double* v_all_faces,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                    int* q_v,int* nint);
    void build_my_w(double* udg,double* w_all_faces,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                    double* h,int* q_u,int* nint);
    void build_my_v_const_coeff(double* vdg,double* v_all_faces,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                    int* q_v,int* nint);
    void build_my_w_const_coeff(double* udg,double* w_all_faces,int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                                double* h,int* q_u,int* nint, double* lambda_lame, double* mu_lame);
    void compute_single_mode_error(double* l2_err,double* udg,int* ifirst,int* ilast,int* jfirst,int* jlast,
                                   int* kfirst,int* klast,double* h,double* t,int* q_u,int* nint);
    void compute_point_dirac_error(double* l2_err,double* udg,int* ifirst,int* ilast,int* jfirst,int* jlast,
                                   int* kfirst,int* klast,double* h,double* t,int* q_u,int* nint, double* parameters);
    void start_next_step(double* updg,double* vpdg,double* udg,double* vdg,
                         int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* q_v,int* q_u);
    void taylor_swap(double* utdg,double* vtdg, double* updg,double* vpdg,double* udg,double* vdg,double* df, 
                         int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* q_v,int* q_u);
    void pass_outside_fluxes(double* v_out_all_faces,double* v_in_all_faces,double* w_out_all_faces,double* w_in_all_faces,
                             int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* nint);
    void set_boundary_conditions(double* v_out_all_faces,double* v_in_all_faces,double* w_out_all_faces,double* w_in_all_faces,
                                 int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* nint,
                                 int* sbx_b, int* sbx_e, int* sby_b, int* sby_e,int* bc_type);
    void compute_numerical_fluxes(double* v_out_all_faces,double* v_in_all_faces,double* w_out_all_faces,double* w_in_all_faces,
                                  double* v_star_all_faces,double* w_star_all_faces,
                                  int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* nint,
                                  double* lambda_lame, double* mu_lame, double* rho);
    void compute_time_derivatives(double* utdg,double* vtdg,double* udg,double* vdg,
                                  double* force_u,double* force_v,
                                  double* MU,double* MV,double* SU,double* SV,
                                  int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                                  int* q_v,int* q_u);
    void get_comm_sides(double* x_in_b,double* x_in_e,double* y_in_b,double* y_in_e,
                        double* v_in_all_faces,double* w_in_all_faces,
                        int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* nint);
    void put_comm_sides(double* x_out_b,double* x_out_e,double* y_out_b,double* y_out_e,
                        double* v_out_all_faces,double* w_out_all_faces,
                        int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* nint);
    void eval_legendre(double* P,double* Pr,int* qu,double* x);
    void get_recorder(double* udg,double* vdg,double* urec,double* vrec,
                      int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,
                      int* q_v,int* q_u,
                      double* Px,double* Prx,double* Py,double* Pry,double* Pz,double* Prz,int* ix,int* iy,int* iz,int* ncomponents);
    void set_tay_weights(double* tg,double* ct,int* nsrc,int* ntay);
    void get_source_tay_coeff(double* source_tay_coeff,double* tg,double* ct,int* nsrc,int* ntay,double* t,double* dt);
    void add_dirac_source(double* force_v,double* point_src,double* f_amp,double* source_tay_coeff,
                          int* ifirst,int* ilast,int* jfirst,int* jlast,int* kfirst,int* klast,int* q_v,
                          int* ix,int* iy,int* iz);
    void get_dirac_source(double* point_src,double* px,double* py,double* pz,int* q_v);
}

