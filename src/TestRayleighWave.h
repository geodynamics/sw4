//-*-c++-*-
#ifndef TEST_RAYLEIGH_WAVE_H
#define TEST_RAYLEIGH_WAVE_H

#include <iostream>


class TestRayleighWave
{
public:

TestRayleighWave( double rho, double cs, double cp, int nwl, double xmax );

double m_rho, m_cp, m_cs, m_cr, m_lambda, m_mu, m_omega, m_alpha;

private:
TestRayleighWave(const TestRayleighWave&);
TestRayleighWave& operator=(const TestRayleighWave&);

};

#endif
