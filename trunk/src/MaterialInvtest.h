//-*-c++-*-
#ifndef MATERIALINVTEST_H
#define MATERIALINVTEST_H

#include "MaterialData.h"

class EW;

class MaterialInvtest : public MaterialData
{
 public:
   MaterialInvtest( EW* a_ew, int nr );
   void set_material_properties(std::vector<Sarray> & rho, std::vector<Sarray> & cs,
				std::vector<Sarray> & cp,
				std::vector<Sarray>& xis, std::vector<Sarray>& xip);
 private:
   int m_nr;
   EW *mEW; 
};

#endif

