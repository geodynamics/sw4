#include <cstring>

#include "MaterialParCartesian.h"

#include "EW.h"

#include "Mopt.h"

//-----------------------------------------------------------------------
Mopt::Mopt( EW* a_ew )
{
   m_ew = a_ew;
   m_opttest = 1;

   m_rhoscale    = 1;
   m_muscale     = 1;
   m_lambdascale = 1;
   m_misfitscale = 1;
   MPI_Comm_rank( MPI_COMM_WORLD, &m_myrank );
}

//-----------------------------------------------------------------------
bool Mopt::parseInputFileOpt( std::string filename )
{
   char buffer[256];
   ifstream inputFile;
   MPI_Barrier(MPI_COMM_WORLD);

   inputFile.open(filename.c_str());
   if (!inputFile.is_open())
   {
      if (m_myrank == 0)
	 cerr << endl << "ERROR OPENING INPUT FILE: " << filename << endl << endl;
      return false;
   }
   while (!inputFile.eof())
   {    
      inputFile.getline(buffer, 256);
      if( strlen(buffer) > 0 && !(startswith("#",buffer)||startswith("\n",buffer)||startswith("\r",buffer) ) )
      {
	 if( startswith("mparcart", buffer) )
	    processMaterialParCart( buffer );
         else if( startswith("mrun",buffer) )
	    processMrun(buffer);
         else if( startswith("mscalefactors",buffer) )
	    processMscalefactors(buffer);
      }
   }
   inputFile.close();
   MPI_Barrier(MPI_COMM_WORLD);
// wait until all processes have read the input file
   if( m_ew->getVerbosity() >=3 && m_myrank == 0 )
      cout << "********Done reading the input file*********" << endl;
   return true;
}

//-----------------------------------------------------------------------
void Mopt::badOption(string name, char* option) const
{
   if (m_myrank == 0)
      cout << "\tWarning: ignoring " << name << " line option '" << option << "'" << endl;
}

//-----------------------------------------------------------------------
bool Mopt::startswith(const char begin[], char *line)
{
  int lenb = strlen(begin);

  // We ignore any preceeding whitespace
  while (strncmp(line, " ", 1) == 0 || strncmp(line, "\t", 1) == 0)
    line++;
    
  if (strncmp(begin, line, lenb) == 0)
     return true;
  else 
     return false;
}

//-----------------------------------------------------------------------
void Mopt::processMaterialParCart( char* buffer )
{
   int nr=1;
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mparcart", token) == 0,
	       "ERROR: not an mparcart line: " << token);
   token = strtok(NULL, " \t");

   int nx=3, ny=3, nz=3, init=0;
   char file[256];

   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("type=",token) )
      {

      }
      else if( startswith("nx=",token) )
      {
         token += 3;
         nx = atoi(token);
      }
      else if( startswith("ny=",token) )
      {
         token += 3;
         ny = atoi(token);
      }
      else if( startswith("nz=",token) )
      {
         token += 3;
         nz = atoi(token);
      }
      else if( startswith("init=",token) )
      {
         token += 5;
	 if( strcmp(token,"0")==0 )
            init = 0;
	 else if( strcmp(token,"sineperturbation")==0 )
	    init = 1;
	 else if( strcmp(token,"interpolate")== 0 )
	    init = 3;
	 else
	 {
            strncpy(file,token,256);
	    init = 2;
	 }
      }
      else
      {
	 badOption("mparcart", token);
      }
      token = strtok(NULL, " \t");
   }
   m_mp = new MaterialParCartesian( m_ew, nx, ny, nz, init, file );
}

//-----------------------------------------------------------------------
void Mopt::processMrun( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mrun", token) == 0,
	       "ERROR: not an mrun line: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("task=",token) )
      {
         token += 5;
	 if( strcmp(token,"minvert") == 0 )
	    m_opttest = 1;
	 else if( strcmp(token,"gradtest") == 0 )
	    m_opttest = 2;
	 else if( strcmp(token,"hesstest") == 0 )
	    m_opttest = 3;
	 else if( strcmp(token,"func1d") == 0 )
	    m_opttest = 4;
	 else if( strcmp(token,"func2d") == 0 )
	    m_opttest = 5;
	 else if( strcmp(token,"forward") == 0 )
	    m_opttest = 6;
	 else
	    cout << "ERROR: mrun task=" << token << " not recognized " << endl;
      }
      else
         badOption("mrun",token);
      token = strtok(NULL," \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::processMscalefactors( char* buffer )
{
   char* token = strtok(buffer, " \t");
   CHECK_INPUT(strcmp("mscalefactors", token) == 0,
	       "ERROR: not an mscalefactors line: " << token);
   token = strtok(NULL, " \t");
   while (token != NULL)
   {
      // while there are tokens in the string still
      if (startswith("#", token) || startswith(" ", buffer))
	// Ignore commented lines and lines with just a space.
	 break;
      else if( startswith("rho=",token) )
      {
         token += 4;
	 m_rhoscale = atof(token);
      }
      else if( startswith("mu=",token) )
      {
         token += 3;
	 m_muscale = atof(token);
      }
      else if( startswith("lambda=",token) )
      {
         token += 7;
	 m_lambdascale = atof(token);
      }
      else if( startswith("misfit=",token) )
      {
         token += 7;
	 m_misfitscale = atof(token);
      }
      else
         badOption("mscalefactors",token);
      token = strtok(NULL," \t");
   }
}

//-----------------------------------------------------------------------
void Mopt::get_scalefactors( double& rhoscale, double& muscale,
			     double& lambdascale, double& fscale )
{
   rhoscale = m_rhoscale;
   muscale = m_muscale;
   lambdascale = m_lambdascale;
   fscale = m_misfitscale;
}
