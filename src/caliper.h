#ifndef __SW4_CALIPER_ANNOTATIONS__
#define __SW4_CALIPER_ANNOTATIONS__

#define SW4_CPU_WARN std::cout<<"WARNING "<<__func__<<" running on CPU \n"

#if defined(ENABLE_CALIPER)

#include <caliper/cali.h>
#include <caliper/cali-manager.h>

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
