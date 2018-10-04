#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv){

  //compile with mpicc -o mpi mpi.c
  //mpiexec -n 4 ./mpi

  long long int Ninside = 0; // number of random points inside 1/4 circle
  long long int Ntests = 1000000000;
  long long n;
  
  //int MPI_Barrier( MPI_Comm comm );
  
  //STEP A: Initialize the MPI process
  MPI_Init(&argc, &argv);

  //STEP B and C: Find the rank and size of the MPI process
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //STEP D: sets the random number seed to rank
  srand48(rank);

  //STEP E: Send and Receive
  long long int sum;
  MPI_Reduce(&Ninside,&sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  
  double estpi = 0;

  srand48(12345);
  
  for(n=0;n<Ntests;++n){
    double x = drand48();
    double y = drand48();
    
    if(x*x+y*y<1){
      ++Ninside;
    }
  }

  estpi = 4.*(Ninside/(double)Ntests);

  printf("estPi = %lf\n", estpi);

  //STEP A: Finalize the MPI process
  MPI_Finalize();

  return 0;
}
