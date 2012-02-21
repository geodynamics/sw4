//-*-c++-*-
#ifndef TEST_POINT_SOURCE_H
#define TEST_POINT_SOURCE_H

class TestPointSource
{
public:

   TestPointSource( double rho, double cs, double cp ) : 
                            m_rho(rho),m_cs(cs),m_cp(cp)
{
   m_mu = m_cs*m_cs*m_rho;
   m_lambda = m_cp*m_cp*m_rho-2*m_mu;
}

double m_rho, m_cp, m_cs, m_lambda, m_mu;

private:
TestPointSource(const TestPointSource&);
TestPointSource& operator=(const TestPointSource&);

};

#endif
