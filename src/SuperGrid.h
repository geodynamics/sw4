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

class SuperGrid 
{

public:
   SuperGrid();
   void define_taper(bool left, double leftStart, bool right, double rightEnd, 
                     double width );
   double dampingCoeff(double x) const;
   double stretching( double x ) const;
   double cornerTaper( double x ) const;
   double tw_stretching( double x ) const;
   double get_tw_omega() const {return m_tw_omega;}
   void   set_twilight( double omega );
   void   print_parameters() const;
   void set_eps( double new_eps );

private:
   bool m_left, m_right;
   double m_x0, m_x1, m_width, m_trans_width, m_const_width;
   double m_epsL, m_tw_omega;
   double Psi0(double xi) const;
   double PsiAux(double x) const;
   double PsiDamp(double x) const;
   double linTaper(double x) const;
};

#endif
