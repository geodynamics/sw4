// -*-c++-*-
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "SecondOrderSection.h"

using namespace std;

SecondOrderSection::SecondOrderSection()
{
  for (int q=0; q<3; q++)
  {
    m_n.m_c[q] = 0;
    m_d.m_c[q] = 0;
  }
  
} // end default constructor

SecondOrderSection::SecondOrderSection(Polynomial &nom, Polynomial &denom)
{
  for (int q=0; q<3; q++)
  {
    m_n.m_c[q] = nom.m_c[q];
    m_d.m_c[q] = denom.m_c[q];
  }
  
} // end constructor

double SecondOrderSection::numer(unsigned int q)
{
  return m_n.m_c[q];
}

double SecondOrderSection::denom(unsigned int q)
{
  return m_d.m_c[q];
}

