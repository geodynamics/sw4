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
   void processLBFGS( char* buffer );
   void processNLCG( char* buffer );
 public:
   Mopt( EW* a_ew );
   bool parseInputFileOpt( std::string filename );
   void get_scalefactors( double& rhoscale, double& muscale,
			  double& lambdascale, double& fscale );
   int m_opttest;
   int m_maxit, m_maxsubit, m_nbfgs_vectors, m_optmethod, m_ihess_guess;
   bool m_dolinesearch, m_fletcher_reeves, m_wolfe, m_mcheck, m_output_ts;
   double m_tolerance;

   MaterialParameterization *m_mp;   
};

#endif
