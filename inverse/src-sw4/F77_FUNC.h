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
#ifndef F77_FUNC_H
#define F77_FUNC_H
/*--------------------------------------------------------------------------*/

/* F77_FUNC - attempt a uniform naming of FORTRAN functions which 
 *          - gets around loader naming conventions
 *          - F77_FUNC(foo, FOO)(x, y, z)
 */

#ifndef F77_FUNC

#ifdef CMAKE_FC_FUNC
#include "FC.h"
#define F77_FUNC FC_GLOBAL
#else /* CMAKE_FC_FUNC */

/* MACOSX predefines __APPLE__ */
#ifdef __APPLE__
#  define F77_FUNC(x, X) x##_
#endif

# ifdef ANSI_F77
#  define F77_FUNC(x, X)  X
# endif /* ANSI_F77 */

# ifndef __GNUC__

#  ifdef __xlC__
#   define F77_FUNC(x, X)  x
#  endif /* IBM XL compiler */

#  ifdef HPUX
#   define F77_FUNC(x, X)  x
#  endif /* HPUX */

# endif /* __GNUC__ */

# ifndef F77_FUNC
#  define F77_FUNC(x, X)  x ## _
# endif /* F77_FUNC */

#endif /* CMAKE_FC_FUNC */

#endif /* F77_FUNC */

#endif /* F77_FUNC_H */

