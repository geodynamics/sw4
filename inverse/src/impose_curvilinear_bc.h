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
#ifndef impose_curvilinear_bc_h
#define impose_curvilinear_bc_h
void
impose_curvilinear_bc(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
		      Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J,
		      double t, boundaryConditionType bcType[] );
void eval_curvilinear_bc_stress(Sarray & a_u, double ** bcForcing, Sarray & a_x, Sarray & a_y, Sarray & a_z,
				Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J);
/* void curvilinear_bc_forcing(double ** bcForcing, int * numberOfBCPoints, Sarray & a_x, Sarray & a_y, Sarray & a_z, */
/* 			    Sarray & a_mu, Sarray & a_lam, Sarray & a_q, Sarray & a_r, Sarray & a_s, Sarray & a_J, */
/* 			    double t, boundaryConditionType bcType[], Forcing * f, */
/* 			    int noof_mechanisms, Sarray* alphaVE, Sarray* muVE, Sarray* lambdaVE ); */

#endif
