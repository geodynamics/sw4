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
      subroutine ffsurfz(x,y,z,t,omega,c,phase,momega,mphase,amprho,ampm
     #u,amplambda,crea_par)
      doubleprecision x
      doubleprecision y
      doubleprecision z
      doubleprecision t
      doubleprecision omega
      doubleprecision c
      doubleprecision phase
      doubleprecision momega
      doubleprecision mphase
      doubleprecision amprho
      doubleprecision ampmu
      doubleprecision amplambda
      doubleprecision crea_par(3)

      doubleprecision forces(3)
      doubleprecision t10
      doubleprecision t13
      doubleprecision t15
      doubleprecision t16
      doubleprecision t19
      doubleprecision t2
      doubleprecision t20
      doubleprecision t21
      doubleprecision t23
      doubleprecision t24
      doubleprecision t28
      doubleprecision t29
      doubleprecision t3
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t37
      doubleprecision t38
      doubleprecision t43
      doubleprecision t44
      doubleprecision t49
      doubleprecision t54
      doubleprecision t56
      doubleprecision t6
      doubleprecision t60
      doubleprecision t62
      doubleprecision t65
      doubleprecision t9

        t2 = momega*x+mphase
        t3 = cos(t2)
        t6 = sin(momega*y+mphase)
        t9 = momega*z+mphase
        t10 = sin(t9)
        t13 = ampmu*(3+t3*t6*t10)
        t15 = omega*x+phase
        t16 = cos(t15)
        t19 = omega*y+phase
        t20 = sin(t19)
        t21 = c*t
        t23 = omega*(z-t21)
        t24 = sin(t23)
        t28 = omega*(x-t21)
        t29 = sin(t28)
        t32 = omega*z+phase
        t33 = cos(t32)
        t34 = t33*omega
        forces(1) = t13*(t16*omega*t20*t24+t29*t20*t34)
        t37 = sin(t15)
        t38 = cos(t19)
        t43 = omega*(y-t21)
        t44 = sin(t43)
        forces(2) = t13*(t37*t38*omega*t24+t37*t44*t34)
        t49 = cos(t23)
        t54 = sin(t2)
        t56 = cos(t9)
        t60 = cos(t28)
        t62 = sin(t32)
        t65 = cos(t43)
        forces(3) = 2*t13*t37*t20*t49*omega+amplambda*(2+t54*t6*t56)*(t6
     #0*omega*t20*t62+t37*t65*omega*t62+t37*t20*t49*omega)
        crea_par(1) = forces(1)
        crea_par(2) = forces(2)
        crea_par(3) = forces(3)
        return
        return
      end
