#include <math.h>

#include "TestRayleighWave.h"
#include "F77_FUNC.h"

extern "C" {
double F77_FUNC(rvel,RVEL)(double *lambda, double *mu);
}

TestRayleighWave:: TestRayleighWave( double rho, double cs, double cp, double alpha_rad ) : 
  m_rho(rho), m_cs(cs), m_cp(cp), m_omega(2*M_PI/1000.0), m_alpha(alpha_rad) 
// m_omega must be assigned before the exact solution can be evaluated
{
  double xi;
  
  m_mu = m_cs*m_cs*m_rho;
  m_lambda = m_cp*m_cp*m_rho-2*m_mu;
  xi = F77_FUNC(rvel,RVEL)( &m_lambda, &m_mu );
  m_cr = xi*m_cs;
}
