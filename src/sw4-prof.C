#include <mpi.h>

#include "sw4-prof.h"

#include <fstream>
#include <iostream>
#include <sstream>

int SW4Prof::m_instances = 0;

// Define an instance of a pointer to SW4Prof
SW4Prof0* sw4_profile;

SW4Prof::SW4Prof(std::string path) {
  if (m_instances == 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nprocs);

    std::stringstream s;
    if (path != ".") s << path;
    s << "sw4prof" << m_myrank << ".txt";
    m_fname = s.str();
    m_instances = 1;

    MPI_Barrier(MPI_COMM_WORLD);
    // Set zero-time
    m_t0 = MPI_Wtime();
  } else
    std::cout << "SW4Prof: Error, trying to create a second instance"
              << std::endl;
}

void SW4Prof::time_stamp(std::string label) {
  m_times.push_back(MPI_Wtime() - m_t0);
  m_labels.push_back(label);
}

void SW4Prof::time_stamp(std::string label, int nr) {
  m_times.push_back(MPI_Wtime() - m_t0);
  std::stringstream s;
  s << label << nr;
  m_labels.push_back(s.str());
}

void SW4Prof::flush() {
  std::ofstream oobj(m_fname.c_str(), std::ofstream::trunc);
  oobj << "Task " << m_myrank << " out of " << m_nprocs << std::endl;
  oobj << m_times.size() << std::endl;
  oobj << m_t0 << " Profiling start time" << std::endl;
  for (int i = 0; i < m_times.size(); i++)
    oobj << m_times[i] << " " << m_labels[i] << std::endl;
  oobj.close();
}

//-----------------------------------------------------------------------
SW4Prof0::SW4Prof0(){};
void SW4Prof0::time_stamp(std::string label){};
void SW4Prof0::time_stamp(std::string label, int nr){};
void SW4Prof0::flush(){};
