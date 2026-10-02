#ifndef __SW4_CALIPER_ANNOTATIONS__
#define __SW4_CALIPER_ANNOTATIONS__

#if defined(ENABLE_CALIPER)

#include <caliper/cali.h>

#define SW4_MARK_FUNCTION CALI_CXX_MARK_FUNCTION

#define SW4_MARK_BEGIN(tag) CALI_MARK_BEGIN((tag))

#define SW4_MARK_END(tag) CALI_MARK_END((tag))

#else
// Stubs for removing Caliper dependencies

#define SW4_MARK_FUNCTION

#define SW4_MARK_BEGIN(tag)

#define SW4_MARK_END(tag)

#endif

#endif
