//-*-c++-*-
#ifndef EW_WPPPIO_H
#define EW_WPPPIO_H

class Comminfo
{
public:
   Comminfo();
   ~Comminfo();
   void print( int recv );
   bool m_has_values;
   int m_steps;
   int*  m_ncomm; 
   int** m_comm_id;
   int** m_comm_index[6];
   int  m_maxbuf;
   int* m_ilow;
   int* m_jlow;
   int* m_klow;
   int* m_niblock;
   int* m_njblock;
   int* m_nkblock;
};

class Parallel_IO
{
public:
   Parallel_IO( int iwrite, int pfs, int globalsizes[3], int localsizes[3],
	    int starts[3], int nptsbuf=1000000, int padding=0 );
   void write_array( int* fid, int nc, void* array, off_t pos0, char* type );
   void read_array( int* fid, int nc, double* array, off_t pos0, char* typ );
			      
   void print( );
   void begin_sequential( MPI_Comm comm );
   void end_sequential( MPI_Comm comm );
   int proc_zero();
   int i_write() const {return m_iwrite==1;}
   void writer_barrier( );
   int n_writers() const {return m_nwriters;}
private:
   void init_pio( int iwrite, int pfs, int ihave_array=-1 );
   void init_array( int globalsizes[3], int localsizes[3], 
		    int starts[3], int nptsbuf, int padding=0 );

   int m_iwrite, m_nwriters, m_parallel_file_system;
   int m_csteps;
   int* m_writer_ids;
   int ni, nj, nk, nig, njg, nkg, oi, oj, ok;

   MPI_Comm m_write_comm; 
   MPI_Comm m_data_comm;

   Comminfo m_isend;
   Comminfo m_irecv;
};

#endif
