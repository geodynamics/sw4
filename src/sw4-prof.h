#ifndef SW4_PROFILE
#define SW4_PROFILE

#include <string>
#include <vector>

class SW4Prof0
{
 public:
   SW4Prof0();
   virtual void time_stamp( std::string label );
   virtual void time_stamp(std::string label, int nr );
   virtual void flush();
};

class SW4Prof : public SW4Prof0
{
 public:
   SW4Prof(std::string path = ".");
   void time_stamp( std::string label );
   void time_stamp(std::string label, int nr );
   void flush();
 private:
   std::vector<std::string> m_labels;
   std::vector<double> m_times;
   int m_myrank, m_nprocs;
   double m_t0;
   std::string m_fname;
   static int m_instances;
};

extern SW4Prof0* sw4_profile;

#endif
