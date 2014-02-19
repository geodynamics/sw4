#ifndef WPP_MATERIALPROPERTY_H
#define WPP_MATERIALPROPERTY_H

#include <string>

class MaterialProperty
{
public:
// -----------------------------------------------------------------
// class MaterialProperty
//
// Defines Vp, Vs, Rho and, optionally, Qp and Qs.
// ------------------------------------------------------------------
MaterialProperty(int id, double vp0, double vp1, double vp2, double vs0, double vs1, double vs2, 
		 double rho0, double rho1, double rho2, double qp, double qs );

void setSqrtCoefficients( double vp1o2, double vs1o2, double rho1o2 );

int m_materialID;
double m_vp0, m_vp1, m_vp2, m_vp1o2, m_vs0, m_vs1, m_vs2, m_vs1o2, m_rho0, m_rho1, m_rho2, m_rho1o2, m_qp, m_qs;

private:
MaterialProperty();
};
#endif
