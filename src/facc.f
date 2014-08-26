!  SW4 LICENSE
! # ----------------------------------------------------------------------
! # SW4 - Seismic Waves, 4th order
! # ----------------------------------------------------------------------
! # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
! # Produced at the Lawrence Livermore National Laboratory. 
! # 
! # Written by:
! # N. Anders Petersson (petersson1@llnl.gov)
! # Bjorn Sjogreen      (sjogreen2@llnl.gov)
! # 
! # LLNL-CODE-643337 
! # 
! # All rights reserved. 
! # 
! # This file is part of SW4, Version: 1.0
! # 
! # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
! # 
! # This program is free software; you can redistribute it and/or modify
! # it under the terms of the GNU General Public License (as published by
! # the Free Software Foundation) version 2, dated June 1991. 
! # 
! # This program is distributed in the hope that it will be useful, but
! # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
! # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
! # conditions of the GNU General Public License for more details. 
! # 
! # You should have received a copy of the GNU General Public License
! # along with this program; if not, write to the Free Software
! # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
      subroutine facc(x,y,z,t,om,c,ph,omm,phm,amprho,ampmu,amplambda,cre
     #a_par)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision om
      doubleprecision c
      doubleprecision ph
      doubleprecision omm
      doubleprecision phm
      doubleprecision amprho
      doubleprecision ampmu
      doubleprecision amplambda
      doubleprecision crea_par(3)

      doubleprecision acc(3)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t14
      doubleprecision t19
      doubleprecision t22
      doubleprecision t30
      doubleprecision t4
      doubleprecision t5
      doubleprecision t7

        t1 = c*t
        t4 = sin(om*(x-t1))
        t5 = om**2
        t7 = c**2
        t10 = sin(om*y+ph)
        t14 = sin(om*z+ph)
        acc(1) = -t4*t5*t7*t10*t14
        t19 = sin(om*x+ph)
        t22 = sin(om*(y-t1))
        acc(2) = -t19*t22*t5*t7*t14
        t30 = sin(om*(z-t1))
        acc(3) = -t19*t10*t30*t5*t7
        crea_par(1) = acc(1)
        crea_par(2) = acc(2)
        crea_par(3) = acc(3)
        return
        return
      end
