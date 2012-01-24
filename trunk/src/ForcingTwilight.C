#include "ForcingTwilight.h"

//-----------------------------------------------------------------------
ForcingTwilight::ForcingTwilight( double omega, double c, double phase,
				  double momega, double mphase, double amprho,
				  double ampmu, double amplambda )
{
   m_omega  = omega;
   m_c      = c;
   m_phase  = phase;
   m_momega = momega;
   m_mphase = mphase;
   m_amprho = amprho;
   m_ampmu  = ampmu;
   m_amplambda = amplambda;
   m_number_mechanisms = 1;
   m_use_attenuation = false;
}
