#ifndef MOPT_H
#define MOPT_H

class EW;

class Mopt
{
 private:
   EW* m_ew;
   int m_myrank;
   double m_rhoscale, m_muscale, m_lambdascale, m_misfitscale;
   void badOption(string name, char* option) const;
   bool startswith(const char begin[], char *line);
   void processMaterialParCart( char* buffer );
   void processMrun( char* buffer );
   void processMscalefactors( char* buffer );
 public:
   Mopt( EW* a_ew );
   bool parseInputFileOpt( std::string filename );
   void get_scalefactors( double& rhoscale, double& muscale,
			  double& lambdascale, double& fscale );
   int m_opttest;
   MaterialParameterization *m_mp;   
};

#endif
