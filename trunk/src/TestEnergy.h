//-*-c++-*-
#ifndef TEST_ENERGY_H
#define TEST_ENERGY_H

#include <stdlib.h>

class TestEnergy
{
public:

   TestEnergy( int seed, double cpcsratio ) : m_seed(seed), m_cpcsratio(cpcsratio)
{
   srand48( m_seed );
}


int m_seed;
double m_cpcsratio;

private:
TestEnergy(const TestEnergy&);
TestEnergy& operator=(const TestEnergy&);

};

#endif
