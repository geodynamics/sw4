//-*-c++-*-
#ifndef MATERIALBLOCK_H
#define MATERIALBLOCK_H

#include "MaterialData.h"

class EW;
class MaterialBlock : public MaterialData
{
 public:
   MaterialBlock( EW * a_ew,double rho, double vs, double vp, double xmin, double xmax, double ymin,
		  double ymax, double zmin, double zmax, double qs=-1, double qp=-1,
		  double freq=1 );

   virtual void set_material_properties( std::vector<Sarray> &rho, std::vector<Sarray> &cs,
					 std::vector<Sarray> &cp,
					 std::vector<Sarray>& xis, std::vector<Sarray>& xip);
   void set_gradients( double rhograd, double vsgrad, double vpgrad );
   void set_absoluteDepth( bool absDepth );

private:
  bool inside_block( double x, double y, double z );
  double m_rho, m_vp, m_vs, m_qp, m_qs, m_freq;
  double m_vpgrad, m_vsgrad, m_rhograd;
  double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax;
  double m_tol;
  bool m_absoluteDepth;
  EW *mEW; // where is this pointer needed?
};

#endif
