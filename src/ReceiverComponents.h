#ifndef SW4_RECEIVER_COMPONENTS_H
#define SW4_RECEIVER_COMPONENTS_H
#include "TimeSeries.h"
#include <vector>

namespace sw4 {
struct ReceiverComponent { const char* name; const char* suffix; };
inline std::vector<ReceiverComponent> receiver_components(TimeSeries::receiverMode mode,
                                                          bool cartesian)
{
   switch(mode) {
   case TimeSeries::Displacement:
      return cartesian ? std::vector<ReceiverComponent>{{"X","x"},{"Y","y"},{"Z","z"}}
                       : std::vector<ReceiverComponent>{{"EW","e"},{"NS","n"},{"UP","u"}};
   case TimeSeries::Velocity:
      return cartesian ? std::vector<ReceiverComponent>{{"Vx","xv"},{"Vy","yv"},{"Vz","zv"}}
                       : std::vector<ReceiverComponent>{{"Vew","ev"},{"Vns","nv"},{"Vup","uv"}};
   case TimeSeries::Div: return {{"Div","div"}};
   case TimeSeries::Curl: return {{"Curlx","curlx"},{"Curly","curly"},{"Curlz","curlz"}};
   case TimeSeries::Strains:
      return {{"Uxx","xx"},{"Uyy","yy"},{"Uzz","zz"},{"Uxy","xy"},{"Uxz","xz"},{"Uyz","yz"}};
   case TimeSeries::DisplacementGradient:
      return {{"DUXDX","duxdx"},{"DUXDY","duxdy"},{"DUXDZ","duxdz"},
              {"DUYDX","duydx"},{"DUYDY","duydy"},{"DUYDZ","duydz"},
              {"DUZDX","duzdx"},{"DUZDY","duzdy"},{"DUZDZ","duzdz"}};
   }
   return {};
}
inline bool receiver_vector(TimeSeries::receiverMode mode)
{ return mode==TimeSeries::Displacement || mode==TimeSeries::Velocity || mode==TimeSeries::Curl; }
inline const char* receiver_unit(TimeSeries::receiverMode mode)
{ return mode==TimeSeries::Displacement ? "m" : mode==TimeSeries::Velocity ? "m/s" : "1"; }
}
#endif
