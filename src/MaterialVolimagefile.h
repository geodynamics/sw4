//-*-c++-*-
#ifndef MATERIALVOLIMAGEFILE_H
#define MATERIALVOLIMAGEFILE_H

#include "MaterialData.h"
#include <string>

class EW;

class MaterialVolimagefile : public MaterialData
{
 public:
   MaterialVolimagefile( EW* a_ew, bool rhomula, std::string path, std::string rhofile, 
	 std::string mufile, std::string lambdafile, std::string qpfile, std::string qsfile );

   void set_material_properties(std::vector<Sarray> & rho, std::vector<Sarray> & cs,
				std::vector<Sarray> & cp,
				std::vector<Sarray>& xis, std::vector<Sarray>& xip);
 private:
   std::string m_path, m_rho, m_mu, m_lambda, m_qp, m_qs;
   bool m_rhomula;
   EW *mEW; 
};

#endif
