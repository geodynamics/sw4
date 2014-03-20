#ifndef MATERIALRFILE_FILE
#define MATERIALRFILE_FILE

#include <string>

#include "MaterialData.h"

class EW;

using namespace std;

class MaterialRfile : public MaterialData
{
 public:
   
   MaterialRfile( EW * a_ew,
		  const std::string file,
		  const std::string directory );

  ~MaterialRfile();

  void set_material_properties(std::vector<Sarray> & rho, 
			       std::vector<Sarray> & cs,
			       std::vector<Sarray> & cp, 
			       std::vector<Sarray> & xis, 
			       std::vector<Sarray> & xip);

  //  int get_material_pt( double x, double y, double z, double& rho, double& cs, double& cp,
  //		       double& qs, double& qp );

  //  void getMinMaxBoundsZ(double& zmin, double& zmax);
   
 protected:
  inline bool inside( double x, double y, double z )
  {
    return m_xminloc <= x && x <= m_xmaxloc && m_yminloc <= y && y <= m_ymaxloc 
      && m_zminloc <= z && z <= m_zmaxloc;
  }

   void read_rfile( );
   void fill_in_fluids( );
   int io_processor( );

   EW* mEW;

   std::string m_model_file, m_model_dir;

   bool m_use_attenuation;

   int m_npatches;
   // Index range of patch in file coordinate system
   vector<int> m_ifirst, m_ilast, m_jfirst, m_jlast, m_kfirst, m_klast, m_ni, m_nj, m_nk;

  // file coordinate system is x=(i-1)*m_hx[gr] + m_xmin[gr], in SW4 coordinates.
   vector<double> m_z0, m_hh, m_hv;
   double m_x0, m_y0;

   // xminloc, xmaxloc, etc. is the bounding box for the set of data patches in this processor.
   double m_xminloc, m_xmaxloc, m_yminloc, m_ymaxloc, m_zminloc, m_zmaxloc;

// 3-dimensional Sarrays
   vector<Sarray> mMaterial;

   //   int m_nlat, m_nlon, m_nmaxdepth, m_nx, m_ny;
   //   int m_nstenc;
   //   double m_h, m_dlon, m_dlat;
   //   int     m_ksed, m_kmoho, m_k410, m_k660;
   //   double *m_lon, *m_lat, *m_x, *m_y;
   //   double  m_vpmin, m_vsmin, m_rhomin;
   //   string m_model_file, m_model_dir, m_model_name;
   //   bool m_qf;
   //   double m_latmin, m_latmax, m_lonmin, m_lonmax, m_depthmin, m_depthmax;
   //   double m_xmin, m_xmax, m_ymin, m_ymax;
   //   bool m_coords_geographic;
};
#endif
