// -*-c++-*-
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
#ifndef SUPERGRID_H
#define SUPERGRID_H

#include "sw4.h"

class SuperGrid 
{

public:
SuperGrid();
void define_taper(bool left, float_sw4 leftStart, bool right, float_sw4 rightEnd, 
		  float_sw4 width );
float_sw4 dampingCoeff(float_sw4 x) const;
float_sw4 stretching( float_sw4 x ) const;
float_sw4 cornerTaper( float_sw4 x ) const;
float_sw4 tw_stretching( float_sw4 x ) const;
float_sw4 get_tw_omega() const {return m_tw_omega;}
void   set_twilight( float_sw4 omega );
void   print_parameters() const;

private:
bool m_left, m_right;
float_sw4 m_x0, m_x1, m_width, m_trans_width, m_const_width;
float_sw4 m_epsL, m_tw_omega;
float_sw4 Psi0(float_sw4 xi) const;
float_sw4 PsiAux(float_sw4 x) const;
float_sw4 PsiDamp(float_sw4 x) const;
float_sw4 linTaper(float_sw4 x) const;

};

#endif
