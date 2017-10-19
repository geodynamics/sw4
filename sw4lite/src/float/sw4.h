#ifndef SW4_H
#define SW4_H

#define float_sw4 float

enum boundaryConditionType { bStressFree, bDirichlet, bSuperGrid, bPeriodic, bCCInterface, bRefInterface, bAEInterface,
 bProcessor, bNone };

enum timeDep { iRicker, iGaussian, iRamp, iTriangle, iSawtooth, iSmoothWave, iErf, iRickerInt, iVerySmoothBump, iBrune,
 iBruneSmoothed, iGaussianWindow, iLiu, iDiscrete, iDBrune, iDirac, iC6SmoothBump, iDiscrete6moments };

#endif
