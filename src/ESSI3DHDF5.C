//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// # Produced at the Lawrence Livermore National Laboratory.
// #
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// #
// # LLNL-CODE-643337
// #
// # All rights reserved.
// #
// # This file is part of SW4, Version: 1.0
// #
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// #
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991.
// #
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details.
// #
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#include "mpi.h"
#include "EW.h"

#include "ESSI3DHDF5.h"

#ifdef USE_HDF5
#include "hdf5.h"
#endif

using std::cout;

ESSI3DHDF5* ESSI3DHDF5::nil=static_cast<ESSI3DHDF5*>(0);

//-----------------------------------------------------------------------
ESSI3DHDF5::ESSI3DHDF5(const std::string& filename, int (&global)[3],
    int (&window)[6], bool ihavearray, int precision) :
  m_start_cycle(-1),
  m_end_cycle(-1),
  m_filename(filename),
  m_ihavearray(ihavearray),
  m_precision(precision)
{
#ifdef USE_HDF5
  bool debug = false;
  for (int d=0; d < 3; d++)
  {
    m_global[d] = global[d];
    m_window[2*d] = window[2*d]; // lo, relative to global
    m_window[2*d+1] = window[2*d+1]; // hi
  }
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if (debug)
  {
    char msg[1000];
    sprintf(msg, "Rank %d global size = %d %d %d\n",
        myRank, m_global[0], m_global[1], m_global[2]);
    cout << msg;
    sprintf(msg, "Rank %d ihavearray = %s\n", myRank,
        (m_ihavearray) ? "yes" : "no");
    cout << msg;
    sprintf(msg, "Rank %d index range = (%d:%d %d:%d %d:%d)\n",
        myRank, m_window[0], m_window[1], m_window[2], m_window[3],
        m_window[4], m_window[5]);
    cout << msg;
    cout.flush();
  }

  // Time step is 1st comp
  m_slice_dims[0] = 1; // write one time step
  m_global_dims[0] = 1; // Extendible dimension for cycle
  m_window_dims[0] = 1; // write one time step
  m_cycle_dims[0] = 1; // this will become cycle-1


  for (int d=1; d < 4; d++)
  {
    m_window_dims[d] = m_window[2*(d-1)+1] - m_window[2*(d-1)] + 1;
    m_global_dims[d] = m_global[d-1];
    m_cycle_dims[d]  = m_global[d-1];
  }

  // 1 i,j column per slice
  m_slice_dims[1] = 1;
  m_slice_dims[2] = 1;
  m_slice_dims[3] = m_global[2];

  // MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//-----------------------------------------------------------------------
ESSI3DHDF5::~ESSI3DHDF5()
{

}

void ESSI3DHDF5::create_file()
{
#ifdef USE_HDF5
  m_mpiprop_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(m_mpiprop_id, H5FD_MPIO_INDEPENDENT);

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  hid_t prop_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(prop_id, comm, info);
  m_file_id = H5Fcreate( const_cast<char*>(m_filename.c_str()),
      H5F_ACC_TRUNC, H5P_DEFAULT, prop_id);
  if (m_file_id < 0)
  {
    cout << "Could not open hdf5 file: " << m_filename << endl;
    MPI_Abort(comm,m_file_id);
  }
  H5Pclose(prop_id);
#endif
}

void ESSI3DHDF5::write_header(double h, double (&lonlat_origin)[2], double az,
  double (&origin)[3], int cycle, double t, double dt)
{
#ifdef USE_HDF5
  /* bool debug=true; */
  bool debug=false;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // Add all the metadata to the file
  if (debug && (myRank == 0))
     cout << "Writing hdf5 metadata for file: " << m_filename << endl;

  // Save the cycle for later
  m_start_cycle = cycle;

  // Global grid lat,lon origin
  hsize_t dim=2;
  hid_t dataspace_id = H5Screate_simple(1,&dim,NULL);
  hid_t dataset_id = H5Dcreate2 (m_file_id,
      "Grid lon-lat origin", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
     m_mpiprop_id, &lonlat_origin);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // Grid azimuth
  dim=1;
  dataspace_id = H5Screate_simple(1,&dim,NULL);
  dataset_id = H5Dcreate2 (m_file_id, "Grid azimuth", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, &az);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // grid spacing
  dim=1;
  dataspace_id = H5Screate_simple(1,&dim,NULL);
  dataset_id = H5Dcreate2 (m_file_id, "ESSI xyz grid spacing", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, &h);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // x,y origin
  dim=3;
  dataspace_id = H5Screate_simple(1,&dim,NULL);
  dataset_id = H5Dcreate2 (m_file_id, "ESSI xyz origin", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, &origin);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // t start
  dim=1;
  dataspace_id = H5Screate_simple(1,&dim,NULL);
  dataset_id = H5Dcreate2 (m_file_id, "time start", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, &t);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // dt
  dim=1;
  dataspace_id = H5Screate_simple(1,&dim,NULL);
  dataset_id = H5Dcreate2 (m_file_id, "timestep", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, &dt);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  if (debug && (myRank == 0))
     cout << "Writing hdf5 metadata done" << endl;
#endif
}

void ESSI3DHDF5::write_topo(double* window_array)
{
#ifdef USE_HDF5
  bool debug=false;
  /* bool debug=true; */
  MPI_Comm comm = MPI_COMM_WORLD;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if (debug && (myRank == 0))
     cout << "Writing hdf5 z coordinate: " << m_filename << endl;

  // Create the data space and chunk for z coordinate
  char msg[1000];
  if (debug)
  {
    sprintf(msg, "Rank %d creating z dataspace size = %d %d %d\n", myRank,
        (int) m_global_dims[0], (int) m_global_dims[1], (int) m_global_dims[2]);
    cout << msg;
  }
  /*
  if (m_ihavearray)
    sprintf(msg, "Rank %d creating z dataspace size = %d %d %d\n", myRank,
        (int) m_global_dims[0], (int) m_global_dims[1], (int) m_global_dims[2]);
  else
    sprintf(msg, "Rank %d exiting z write\n", myRank);
  */
  hsize_t z_dims = 3;
  // Modify dataset creation properties to enable chunking
  hid_t prop_id = H5Pcreate (H5P_DATASET_CREATE);
  herr_t ierr;
  /* ierr = H5Pset_chunk (prop_id, z_dims, m_slice_dims); */

  // Create dataset for z field
  hid_t dataspace_id = H5Screate_simple(z_dims, &m_global_dims[1], NULL);
  hid_t dataset_id = H5Dcreate2 (m_file_id, "z coordinates", H5T_IEEE_F64LE,
      dataspace_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
  ierr = H5Pclose(prop_id);

  if (debug)
  {
    sprintf(msg, "Rank %d creating chunk/slice size = %d %d %d\n", myRank,
        (int) m_slice_dims[1], (int) m_slice_dims[2], (int) m_slice_dims[3]);
    cout << msg;
  }

  /*
  if (debug)
     cout << "Rank " << myRank <<
       " in write_topo_real: starting slab loop" << endl;
  sprintf(msg, "Rank %d looping over hyperslab = (%d:%d %d:%d %d:%d)\n",
      myRank, m_window[0], m_window[1], m_window[2], m_window[3],
      m_window[4], m_window[5]);
  cout << msg;
  */

  hsize_t start[3]={0,0,0};
  start[0] = m_window[0]; // i window offset
  start[1] = m_window[2]; // j window offset
  start[2] = m_window[4]; // k local index lo relative to global
  // We should write data
  if (debug)
  {
    sprintf(msg, "Rank %d selecting z hyperslab = %d %d %d\n",
        myRank, (int) start[0], (int) start[1], (int) start[2]);
    cout << msg;
    cout.flush();
  }
  hid_t window_id = H5Screate_simple(z_dims, &m_window_dims[1], NULL);
  ierr = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL,
      &m_window_dims[1], NULL);
  if (ierr < 0)
  {
    cout << "Error from z H5Sselect_hyperslab" << endl;
    MPI_Abort(comm,ierr);
  }
  if (debug)
  {
    char msg[1000];
    sprintf(msg, "Writing z array Rank %d\n", myRank);
    cout << msg;
    cout.flush();
  }

  ierr = H5Dwrite(dataset_id, H5T_IEEE_F64LE, window_id, dataspace_id,
      m_mpiprop_id, window_array);
  if (ierr < 0)
  {
    cout << "Error from z H5Dwrite " << endl;
    MPI_Abort(comm,ierr);
  }

  ierr = H5Sclose(dataspace_id);
  ierr = H5Dclose(dataset_id);
  /*
  MPI_Barrier(comm);
  if (debug && (myRank == 0))
     cout << "In write_topo_real: done writing data" << endl;
  */

  /*
  MPI_Barrier(comm);
  */
  if (debug && (myRank == 0))
     cout << "Done writing hdf5 z coordinate: " << m_filename << endl;
#endif
}

void ESSI3DHDF5::init_write_vel(int ntimestep)
{
#ifdef USE_HDF5
  /* bool debug=true; */
  bool debug=false;
  MPI_Comm comm = MPI_COMM_WORLD;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  hid_t dtype = H5T_NATIVE_DOUBLE;
  if (m_precision == 4) 
      dtype = H5T_NATIVE_FLOAT;

  if (debug)
  {
    char msg[1000];
    sprintf(msg, "Rank %d creating vel dataspaces size = %d %d %d\n", myRank,
        (int) m_global_dims[0], (int) m_global_dims[1], (int) m_global_dims[2]);
    cout << msg;
  }
  // Create the extendible velocity data space and chunk
  hsize_t num_vel=3;
  hsize_t num_dims=4;
  // m_xvel_dataspace_id = H5Screate_simple(num_dims, m_global_dims, dims);
  hid_t prop_id = H5Pcreate (H5P_DATASET_CREATE);
  // Tang
  /* herr_t ierr = H5Pset_chunk (prop_id, num_dims, m_slice_dims); */
  /* if (myRank == 0) { */
  /*     printf("Chunk sizes:\n"); */
  /*     for (int i = 0; i < num_dims; i++) */ 
  /*         printf("%llu  ", m_slice_dims[i]); */
  /* } */
  /* printf("\n"); */
  if (ntimestep > 0) 
    m_cycle_dims[0]  = ntimestep + 1;
  else
    printf("Error with m_ntimestep=%d!\n", ntimestep);

  /* if (myRank == 0) { */
  /*     printf("m_cycle_dims sizes:\n"); */
  /*     for (int i = 0; i < num_dims; i++) */ 
  /*         printf("%llu  ", m_cycle_dims[i]); */
  /* } */
  /* printf("\n"); */

  for (int c=0; c < num_vel; c++)
  {
    m_vel_dataspace_id[c] =
      H5Screate_simple(num_dims, m_cycle_dims, m_cycle_dims);
      /* H5Screate_simple(num_dims, m_cycle_dims, m_global_dims); */
    // Modify dataset creation properties to enable chunking
    char var[100];
    sprintf(var, "vel_%d ijk layout", c);
    m_vel_dataset_id[c] = H5Dcreate2 (m_file_id, var, dtype,
        m_vel_dataspace_id[c], H5P_DEFAULT, prop_id, H5P_DEFAULT);
  }
  H5Pclose(prop_id);
  if (debug )
     cout << "Rank" << myRank << "Done creating vel" << endl;
#endif
}

void ESSI3DHDF5::write_vel(double* window_array, int comp, int cycle)
{
#ifdef USE_HDF5
  /* bool debug=true; */
  bool enable_timing = false;
  bool debug=false;
  MPI_Comm comm = MPI_COMM_WORLD;
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  herr_t ierr;
  double write_time_start, write_time;
  int write_size = 1;
  m_end_cycle = cycle; // save for header for later when we close the file

  for (int i = 0; i < 4; i++) 
      write_size *= m_window_dims[i];
  
  if (write_size == 0) 
      return;
  
  if (enable_timing) 
    write_time_start = MPI_Wtime();
  // Loop over slices and write a new time step

  m_cycle_dims[0] = cycle - m_start_cycle + 1;
  /* ierr = H5Dset_extent(m_vel_dataset_id[comp], m_cycle_dims); */
  /* if (ierr < -1) */
  /* { */
  /*   cout << "Error from H5Dset_extent " << endl; */
  /*   MPI_Abort(comm,ierr); */
  /* } */
  hsize_t vel_dims=4;
  hid_t window_id = H5Screate_simple(vel_dims, m_window_dims, NULL);
  // MPI_Barrier(comm);
  char msg[1000];
  /*
  if (debug)
     cout << "Rank " << myRank <<
       " in write_vel: starting slab loop" << endl;
  sprintf(msg, "Rank %d writing vel hyperslab = (%d:%d %d:%d %d:%d %d)\n",
      myRank, m_window[0], m_window[1], m_window[2], m_window[3],
      m_window[4], m_window[5], );
  cout << msg;
  */

  m_vel_dataspace_id[comp] = H5Dget_space(m_vel_dataset_id[comp]);
  hsize_t start[4]={0,0,0,0};
  start[0] = cycle - m_start_cycle;
  start[1] = m_window[0]; // local index lo relative to global
  start[2] = m_window[2]; // local index lo relative to global
  start[3] = m_window[4]; // local index lo relative to global
  if (debug)
  {
    sprintf(msg, "Rank %d writing vel dataset index = %d %d %d %d\n",
        myRank,
        (int) start[0], (int) start[1], (int) start[2], (int) start[3]);
    cout << msg;
  }
  ierr = H5Sselect_hyperslab(m_vel_dataspace_id[comp], H5S_SELECT_SET,
            start, NULL, m_window_dims, NULL);
  if (ierr < 0)
  {
    cout << "Error from vel H5Sselect_hyperslab" << endl;
    MPI_Abort(comm,ierr);
  }
  if (debug)
  {
    char msg[1000];
    sprintf(msg, "Writing vel array Rank %d\n", myRank);
    cout << msg;
    cout.flush();
  }
  ierr = H5Dwrite(m_vel_dataset_id[comp], H5T_IEEE_F64LE, window_id,
      m_vel_dataspace_id[comp], m_mpiprop_id, window_array);
  if (ierr < 0)
  {
    cout << "Error from vel H5Dwrite " << endl;
    MPI_Abort(comm,ierr);
  }
  if (debug )
     cout << "Rank" << myRank << "Done writing vel" << endl;
  // H5Sclose(slice_id);
  
  if (enable_timing) {
    write_time = MPI_Wtime() - write_time_start;
    if (cycle % 100  == 0) {
      printf("rank=%d, cycle=%d, comp=%d, write size=%fMB, write time=%f\n", 
              myRank, cycle, comp, write_size/1048576.0, write_time);
      fflush(stdout);
    }
  }

#endif
}

void ESSI3DHDF5::close_file()
{
#ifdef USE_HDF5
  // Updated header with final start,end cycle
  hsize_t dim=2;
  hid_t dataspace_id = H5Screate_simple(1,&dim,NULL);
  hid_t dataset_id = H5Dcreate2 (m_file_id, "cycle start, end", H5T_STD_I32LE,
      dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int cycles[2] = {m_start_cycle, m_end_cycle};
  herr_t ierr = H5Dwrite(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL,
      m_mpiprop_id, cycles);
  ierr = H5Dclose(dataset_id);
  ierr = H5Sclose(dataspace_id);

  // Close each velocity dataspace, then dataset
  for (int comp=0; comp < 3; comp++)
  {
    H5Sclose(m_vel_dataspace_id[comp]);
    H5Dclose(m_vel_dataset_id[comp]);
  }
  H5Pclose(m_mpiprop_id);
  H5Fclose(m_file_id);
#endif
}

