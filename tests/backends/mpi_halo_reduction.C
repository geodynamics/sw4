#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <cstdio>
#include <vector>

// Isolate the halo/reduction sequence from SW4's scientific operators.
int main(int argc, char** argv)
{
   MPI_Init(&argc, &argv);
   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   std::vector<double> work(1000);
#pragma omp parallel for
   for (int i=0; i<1000; ++i) work[i]=i;
   int cpu=sched_getcpu();
   std::vector<int> cpus(size);
   MPI_Allgather(&cpu, 1, MPI_INT, cpus.data(), 1, MPI_INT, MPI_COMM_WORLD);
   // This on-demand probe is intentionally a single-node, single-thread test.
   for (int i=0; i<size; ++i)
      for (int j=0; j<i; ++j)
         if (cpus[i]==cpus[j]) {
            if (!rank) std::fprintf(stderr, "FAIL: ranks %d and %d share logical CPU %d\n", i, j, cpus[i]);
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
   int dims[2]={0,0}, periods[2]={0,0};
   MPI_Dims_create(size, 2, dims);
   MPI_Comm cart;
   MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
   int left, right, down, up;
   MPI_Cart_shift(cart, 0, 1, &left, &right);
   MPI_Cart_shift(cart, 1, 1, &down, &up);
   std::vector<double> send(300, rank), recv(300);
   for (int iteration=0; iteration<10000; ++iteration) {
#pragma omp parallel for
      for (int i=0; i<1000; ++i) work[i]=rank+iteration+i;
      MPI_Sendrecv(send.data(), 300, MPI_DOUBLE, right, 345, recv.data(), 300,
                   MPI_DOUBLE, left, 345, cart, MPI_STATUS_IGNORE);
      MPI_Sendrecv(send.data(), 300, MPI_DOUBLE, left, 346, recv.data(), 300,
                   MPI_DOUBLE, right, 346, cart, MPI_STATUS_IGNORE);
      MPI_Sendrecv(send.data(), 300, MPI_DOUBLE, up, 347, recv.data(), 300,
                   MPI_DOUBLE, down, 347, cart, MPI_STATUS_IGNORE);
      MPI_Sendrecv(send.data(), 300, MPI_DOUBLE, down, 348, recv.data(), 300,
                   MPI_DOUBLE, up, 348, cart, MPI_STATUS_IGNORE);
      double local=rank+10.*iteration, result=-1;
      MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, cart);
      if (result!=size-1+10.*iteration) {
         std::fprintf(stderr, "FAIL: rank=%d iteration=%d reduction=%.17g\n", rank, iteration, result);
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
   if (!rank) std::puts("PASS: distinct logical CPUs and 10000 halo/reduction iterations");
   MPI_Comm_free(&cart);
   MPI_Finalize();
}
