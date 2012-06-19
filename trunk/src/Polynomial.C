// -*-c++-*-
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "Polynomial.h"

using namespace std;

Polynomial::Polynomial()
{
  for (int q=0; q<3; q++)
  {
    m_c[q] = 0.;
  }
} // end default constructor

Polynomial::Polynomial(double c[3])
{
  for (int q=0; q<3; q++)
  {
    m_c[q] = c[q];
  }
} // end constructor

double Polynomial::coeff(unsigned int q)
{
  return m_c[q];
}

// output all coefficients
ostream& operator<<( ostream& output, const Polynomial& s )
{
  output << "s^0: " << s.m_c[0] << ", s^1: " << s.m_c[1] << ", s^2: " << s.m_c[2];
}
