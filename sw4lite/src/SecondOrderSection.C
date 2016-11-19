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
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "SecondOrderSection.h"

using namespace std;

SecondOrderSection::SecondOrderSection()
{
  for (int q=0; q<3; q++)
  {
    m_n.m_c[q] = 0;
    m_d.m_c[q] = 0;
  }
  
} // end default constructor

SecondOrderSection::SecondOrderSection(Polynomial &nom, Polynomial &denom)
{
  for (int q=0; q<3; q++)
  {
    m_n.m_c[q] = nom.m_c[q];
    m_d.m_c[q] = denom.m_c[q];
  }
  
} // end constructor

float_sw4 SecondOrderSection::numer(unsigned int q)
{
  return m_n.m_c[q];
}

float_sw4 SecondOrderSection::denom(unsigned int q)
{
  return m_d.m_c[q];
}

