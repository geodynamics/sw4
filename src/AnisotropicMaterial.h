#ifndef SW4_ANISOTROPICMATERIAL_H
#define SW4_ANISOTROPICMATERIAL_H

#include "Sarray.h"

class AnisotropicMaterial {
 public:
  AnisotropicMaterial() { mCoversAllPoints = false; };
  virtual void set_material_properties(std::vector<Sarray>& rho,
                                       std::vector<Sarray>& c) = 0;
  bool coversAllPoints() { return mCoversAllPoints; };

 protected:
  bool mCoversAllPoints;
};

#endif
