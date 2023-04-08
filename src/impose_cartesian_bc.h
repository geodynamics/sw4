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
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General
// Public License"
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
#ifndef impose_cartesian_bc_h
#define impose_cartesian_bc_h
void impose_cartesian_bc(Sarray& u, Sarray& mu, Sarray& lam, double** bcForcing,
                         boundaryConditionType bcType[], sw4_type onesided[],
                         double curlcoeff, double h);

/* void */
/* cartesian_bc_forcing(Sarray & mu, */
/* 		     double ** bcForcing, */
/* 		     sw4_type * numberOfBCPoints, */
/* 		     boundaryConditionType bcType[], */
/* 		     Forcing* forcing, double zmin, double t, double h, sw4_type
 * onesided[], double curlcoeff, */
/* 		     sw4_type numberof_mechanisms, Sarray* AlphaVE, Sarray* MuVE,
 * Sarray* LambdaVE ); */
#endif
