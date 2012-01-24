#ifndef MATERIALDATA_H
#define MATERIALDATA_H

#include "Sarray.h"

class MaterialData
{
public:
MaterialData();
virtual void set_material_properties(std::vector<Sarray> & rho, std::vector<Sarray> & cs,
                                     std::vector<Sarray> & cp,
                                     std::vector<Sarray>& xis, std::vector<Sarray>& xip) =0;

// a block can be ignored if it is followed by another block which covers all its grid points

bool coversAllPoints(){return mCoversAllPoints;};

protected:
bool mCoversAllPoints;
};

#endif


