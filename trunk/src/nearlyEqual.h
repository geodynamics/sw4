#ifndef NEARLYEQUAL_H
#define NEARLYEQUAL_H
namespace dbc
{

// Lock for assertions
// bool assertionLock (void);
// void assertionUnLock(void);

// nearlyEqual returns true if x and y approximately equal
template<typename T, typename U>
inline bool nearlyEqual(const T& x, 
                 const U& y, 
                 double relativeTolerance = 1.0e-5,
                 double absoluteTolerance = 1.0e-12)
{
   return (std::fabs(x-y) <=
           (std::fabs(y) * relativeTolerance + absoluteTolerance)
          );
}
} // namespace dbc
#endif
