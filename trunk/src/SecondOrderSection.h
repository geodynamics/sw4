// -*-c++-*-
#ifndef SECOND_ORDER_SECTION_H
#define SECOND_ORDER_SECTION_H

using namespace std;
#include "Polynomial.h"

class SecondOrderSection{
public:
SecondOrderSection(Polynomial &n, Polynomial &d);
double numer(unsigned int q);
double denom(unsigned int q);

// for efficiency and simplicity reasons, we make the coefficients public
Polynomial m_n, m_d;

private:   
SecondOrderSection();

};

#endif
