#include "MaterialProperty.h"

MaterialProperty:: MaterialProperty(int id, double vp0, double vp1, double vp2, double vs0, double vs1, double vs2, 
				    double rho0, double rho1, double rho2, double qp, double qs )
{
  m_materialID = id;
  m_vp0 = vp0;
  m_vp1 = vp1;
  m_vp2 = vp2;

  m_vs0 = vs0;
  m_vs1 = vs1;
  m_vs2 = vs2;

  m_rho0 = rho0;
  m_rho1 = rho1;
  m_rho2 = rho2;

  m_vp1o2 = 0.;
  m_vs1o2 = 0.;
  m_rho1o2 = 0.;

  m_qp = qp;
  m_qs = qs;
}

void MaterialProperty::setSqrtCoefficients( double vp1o2, double vs1o2, double rho1o2 )
{
  m_vp1o2 = vp1o2;
  m_vs1o2 = vs1o2;
  m_rho1o2 = rho1o2;
}



