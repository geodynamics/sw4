#ifndef MATERIALPFILE_FILE
#define MATERIALPFILE_FILE

#include <string>

#include "MaterialData.h"

class EW;

using namespace std;

class MaterialPfile : public MaterialData
{
 public:
   
   MaterialPfile( EW * a_ew,
	      const std::string file,
	      const std::string directory,
	      const int nstenc,
	      const double vpminppm,
	      const double vsminppm,
	      const double rhominppm,
	      const bool flatten,
	      const bool coords_geographic );

  ~MaterialPfile();

  void set_material_properties(std::vector<Sarray> & rho, 
			       std::vector<Sarray> & cs,
			       std::vector<Sarray> & cp, 
			       std::vector<Sarray> & xis, 
			       std::vector<Sarray> & xip);

  //  void getMinMaxBoundsZ(double& zmin, double& zmax);
   
 protected:
  inline bool inside( double lat, double lon, double z )
  {
    return m_latmin <= lat && lat <= m_latmax && m_lonmin <= lon && lon <= m_lonmax 
      && m_zmin <= z && z <= m_zmax;
       //      && m_elevmin <= elev && elev <= m_elevmax;
  }
  inline bool inside_cart( double x, double y, double z )
  {
    return m_xmin <= x && x <= m_xmax && m_ymin <= y && y <= m_ymax 
      && m_zmin <= z && z <= m_zmax;
    //      && m_elevmin <= elev && elev <= m_elevmax;
  }

   void read_pfile( );

   void sample_cart(double,double,double,double&,double&,double&,double&,double&,
	       double&,double&,bool&) ;

   void sample_latlon(double,double,double,double&,double&,double&,double&,double&,
	       double&,double&,bool&) ;

   EW* mEW;
   int m_nlat, m_nlon, m_nmaxdepth, m_nx, m_ny;
   int m_nstenc;
   double m_h;
   int     m_ksed, m_kmoho, m_k410, m_k660;
   double *m_lon, *m_lat, *m_x, *m_y;
   double *m_z, *m_vp, *m_vs, *m_rho, *m_qp, *m_qs;
   double  m_vpmin, m_vsmin, m_rhomin;
   double *m_st, *m_ct;
   string m_model_file, m_model_dir, m_model_name;
   bool m_qf;
   bool m_flatten;

   double m_latmin, m_latmax, m_lonmin, m_lonmax, m_elevmin, m_elevmax;
   double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
   bool m_coords_geographic, m_absoluteDepth;
};
#endif
