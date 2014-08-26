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
#ifndef REQUIRE_H
#define REQUIRE_H

//---------------------------------------------------------------------------
//
// Require.h -- Simplyfied design by contract tools.
//
//---------------------------------------------------------------------------
#include <mpi.h>
#include <iostream>
#include <string>
#include <cmath>

//----------------------------------------------------------------------------
//                            REQUIRE & ASSERT -- Preconditions
//----------------------------------------------------------------------------

#ifdef ASSERT2
#undef ASSERT2
#endif

#ifdef REQUIRE2
#undef REQUIRE2
#endif

#ifdef VERIFY2
#undef VERIFY2
#endif

#ifdef CHECK_INPUT
#undef CHECK_INPUT
#endif

#ifdef CHECK_INPUT2
#undef CHECK_INPUT2
#endif

// This macro is used both for optimized and non-optimized code
#define CHECK_INPUT(x, msg) \
if (!(x)) { \
  int myRank; \
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank); \
  std::cout << "Fatal input error: " << msg << std::endl;	\
  MPI_Abort( MPI_COMM_WORLD, 1 );\
}

#define CHECK_INPUT2(x, msg, fp)			\
if (!(x)) { \
  fp << "Fatal input error: " << msg << std::endl;	\
}

// these macros are also used both for optimized and non-optimized code
#define DBC_ASSERTION(x, msg, kind) \
if (!(x)) { \
  int myRank; \
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank); \
  if (myRank==0){ \
     std::cout << kind << ": " << msg << std::endl;	\
     std::cout << "...at line " << __LINE__ <<		\
	" of file " << __FILE__ << "." << std::endl;	\
  }\
  MPI_Abort( MPI_COMM_WORLD, 1 );\
}
#define REQUIRE2(x, msg) DBC_ASSERTION(x, msg, "Precondition violated")
#define ASSERT2(x, msg) DBC_ASSERTION(x, msg, "Assertion violated")
#define VERIFY2(x, msg) DBC_ASSERTION(x, msg, "Verification failed");

//----------- Define one-argument forms
#ifdef ASSERT
#undef ASSERT
#endif

#ifdef REQUIRE
#undef REQUIRE
#endif

#ifdef VERIFY
#undef VERIFY
#endif

#define ASSERT(x) ASSERT2(x, #x)
#define REQUIRE(x) REQUIRE2(x, #x)
#define VERIFY(x) VERIFY2(x, #x)

#endif // REQUIRE_H

