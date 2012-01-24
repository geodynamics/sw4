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

// This macro is used both for optimized and non-optimized code
#define CHECK_INPUT(x, msg) \
if (!(x)) { \
  int myRank; \
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank); \
  std::cout << "Fatal input error: " << msg << std::endl;	\
  MPI_Abort( MPI_COMM_WORLD, 1 );\
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

