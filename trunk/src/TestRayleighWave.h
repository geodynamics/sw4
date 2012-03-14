//-*-c++-*-
#ifndef TEST_RAYLEIGH_WAVE_H
#define TEST_RAYLEIGH_WAVE_H

#include <iostream>

class TestRayleighWave
{
public:

TestRayleighWave( double rho, double cs, double cp, 
		  double kx, double ky ) : 
  m_rho(rho),m_cs(cs),m_cp(cp),m_omega(2*M_PI/1000.0),m_kx(kx),m_ky(ky) 
// m_omega must be assigned before the exact solution can be evaluated
  {
    m_mu = m_cs*m_cs*m_rho;
    m_lambda = m_cp*m_cp*m_rho-2*m_mu;
    double muolap2mu = m_mu/(2*m_mu+m_lambda);
    double xi = 0.8, er=1;
    int it = 0;
    double xip, rat;
    while( it < 100 && er > 1e-14 )
    {
      rat=(1-muolap2mu*xi)/(1-xi);
      xip = xi - (sqrt((1-xi)*(1-muolap2mu*xi))-(1-0.5*xi)*(1-0.5*xi))/( 
        ( -0.5*sqrt( rat )-muolap2mu*sqrt(1/rat) +1-0.5*xi) );
      er = fabs(xip-xi);
      xi = xip;
      it++;
    }
    if( er > 1e-14 )
    {
      std::cout << "TestRayleighWave: Error could not compute cr " << std::endl;
      std::cout << "Error = " << er << " after " << it << " newton iterations " << std::endl;
    }
    else
    {
      m_cr = xi*m_cs;
    }
  }

double m_rho, m_cp, m_cs, m_cr, m_lambda, m_mu, m_omega, m_kx, m_ky;

private:
TestRayleighWave(const TestRayleighWave&);
TestRayleighWave& operator=(const TestRayleighWave&);

};

#endif
