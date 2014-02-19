//-*-c++-*-
#ifndef FORCING_TWILIGHT_H
#define FORCING_TWILIGHT_H

#include "Sarray.h"
#include "boundaryConditionTypes.h"
#include <string>

class ForcingTwilight
{
public:

ForcingTwilight( double omega, double c, double phase, double momega, double mphase,
		 double amprho, double ampmu, double amplambda );

//void default_bcs( boundaryConditionType bcs[6] );

double m_omega, m_c, m_phase, m_momega, m_mphase, m_amprho, m_ampmu, m_amplambda;
bool m_use_attenuation;
int m_number_mechanisms;

private:
ForcingTwilight(const ForcingTwilight&);
ForcingTwilight& operator=(const ForcingTwilight&);

};

#endif
