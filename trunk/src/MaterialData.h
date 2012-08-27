//-*-c++-*-
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

virtual int get_material_pt( double x, double y, double z, double& rho, double& cs, double& cp,
			     double& qs, double& qp )=0;

protected:
bool mCoversAllPoints;
};

#endif


