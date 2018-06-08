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
#ifndef WPP_MATERIALPROPERTY_H
#define WPP_MATERIALPROPERTY_H

#include <string>

class MaterialProperty
{
public:
// -----------------------------------------------------------------
// class MaterialProperty
//
// Defines Vp, Vs, Rho and, optionally, Qp and Qs.
// ------------------------------------------------------------------
MaterialProperty(int id, double vp0, double vp1, double vp2, double vs0, double vs1, double vs2, 
		 double rho0, double rho1, double rho2, double qp, double qs );

void setSqrtCoefficients( double vp1o2, double vs1o2, double rho1o2 );

int m_materialID;
double m_vp0, m_vp1, m_vp2, m_vp1o2, m_vs0, m_vs1, m_vs2, m_vs1o2, m_rho0, m_rho1, m_rho2, m_rho1o2, m_qp, m_qs;

private:
MaterialProperty();
};
#endif
