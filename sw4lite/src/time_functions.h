//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
#ifndef EW_TIMEFUNCTIONS_H
#define EW_TIMEFUNCTIONS_H
#include "sw4.h"

float_sw4 VerySmoothBump(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_tom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_omom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_tttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_tttom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 VerySmoothBump_ttomom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 RickerWavelet(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerWavelet_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerWavelet_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerWavelet_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerWavelet_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerWavelet_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 RickerInt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerInt_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerInt_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerInt_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerInt_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 RickerInt_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Gaussian(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_tom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_omom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_tttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_tttom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Gaussian_ttomom(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Erf( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Erf_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Erf_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Erf_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Erf_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Erf_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Ramp(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Ramp_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Ramp_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Ramp_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Ramp_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Ramp_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Triangle(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Triangle_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Triangle_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Triangle_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Triangle_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Triangle_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Sawtooth(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Sawtooth_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Sawtooth_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Sawtooth_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Sawtooth_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Sawtooth_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 SmoothWave(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 SmoothWave_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 SmoothWave_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 SmoothWave_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 SmoothWave_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 SmoothWave_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Brune( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Brune_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Brune_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Brune_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Brune_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Brune_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 DBrune( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 DBrune_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 DBrune_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 DBrune_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 DBrune_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 DBrune_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 BruneSmoothed( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 BruneSmoothed_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 BruneSmoothed_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 BruneSmoothed_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 BruneSmoothed_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 BruneSmoothed_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 GaussianWindow( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 GaussianWindow_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 GaussianWindow_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 GaussianWindow_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 GaussianWindow_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 GaussianWindow_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Liu( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Liu_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Liu_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Liu_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Liu_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Liu_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 NullFunc( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Dirac( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_tttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_tom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_omom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_tttom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Dirac_ttomom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 Discrete( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_t( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_tt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_ttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_tttt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_om( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_tom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_omom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_omtt( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_tttom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 Discrete_ttomom( float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

float_sw4 C6SmoothBump(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 C6SmoothBump_t(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 C6SmoothBump_om(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 C6SmoothBump_tt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 C6SmoothBump_ttt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );
float_sw4 C6SmoothBump_omtt(float_sw4 freq, float_sw4 t, float_sw4* par, int npar, int* ipar, int nipar );

#endif
