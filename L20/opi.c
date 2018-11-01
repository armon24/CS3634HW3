#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv){

  //2A - Completed above in an #include statement..........................
  
  //2B.....................................................................

  double stime = omp_get_wtime();
  
  int threadCount = atoi(argv[1]);
  omp_set_num_threads(threadCount);

  long long int Ninside = 0; // number of random points inside 1/4 circle
  long long int Ntests = 1000000000;

  //long int* array;
  //array = calloc(threadCount, sizeof(long int));

   double estpi = 0;

    

  //2C and 2D .............................................................

  #pragma omp parallel num_threads(threadCount) reduction(+:Ninside)                               //stuff in this scope gets executed by all OpenMP threads
  {
    int rank = omp_get_thread_num();
    long long n;

     struct drand48_data buff;
       
    //srand48_r( rank, array + rank );
    srand48_r(12345*rank, &buff);

    #pragma omp for
    for(n=0;n<Ntests;++n)
    {
      double x;
      drand48_r(&buff, &x);
      
      double y;
      drand48_r(&buff, &y);
      
      if(x*x+y*y<1)
      {
	++Ninside;
      }
    }
    
  }

  estpi = 4.*(Ninside/(double)Ntests);

  double time = omp_get_wtime()-stime;
  printf("estPi = %lf\n", estpi);
  printf("time = %lf\n", time);

  return 0;
}
