//-*-c++-*-
#ifndef TEST_ENERGY_H
#define TEST_ENERGY_H

class TestEnergy
{
public:

TestEnergy( double seed ) : m_seed(seed)
{
}


double m_seed;

private:
TestEnergy(const TestEnergy&);
TestEnergy& operator=(const TestEnergy&);

};

#endif
