/* 

To compile:

   gcc -O3 -o mandelbrot mandelbrot.c -lm

To create an image with 4096 x 4096 pixels

    ./mandelbrot 4096 4096 

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

int writeMandelbrot(const char *fileName, int width, int height, float *img, int minI, int maxI);

// Q2a: add include for CUDA header file here:

#define MXITER 1000

typedef struct {
  
  double r;
  double i;
  
}complex_t;

// return iterations before z leaves mandelbrot set for given c
int testpoint(complex_t c){

  int iter;
  complex_t z;
  double temp;
  
  z = c;
  
  for(iter=0; iter<MXITER; iter++){  
    temp = (z.r*z.r) - (z.i*z.i) + c.r;
    
    z.i = z.r*z.i*2. + c.i;
    z.r = temp;
    
    if((z.r*z.r+z.i*z.i)>4.0){
      return iter;
    }
  }
  return iter; 
}

// perform Mandelbrot iteration on a grid of numbers in the complex plane
// record the  iteration counts in the count array
__global__ void mandelbrot(int Nre, int Nim, complex_t cmin, complex_t dc, float *count){ 

  int xthreadIndex = threadIdx.x;
  int xblockIndex  = blockIdx.x;
  int xthreadCount = blockDim.x;

  int ythreadIndex = threadIdx.y;
  int yblockIndex  = blockIdx.y;
  int ythreadCount = blockDim.y;

  int nX = xthreadIndex + xthreadCount*xblockIndex;
  int nY = ythreadIndex + ythreadCount*yblockIndex;

  complex_t c;

  c.r = cmin.r + dc.r*nX;
  c.i = cmin.i + dc.i*nY;
      
  count[nX+nY*Nre] = (float) testpoint(c);

}

int main(int argc, char **argv){

  // to create a 4096x4096 pixel image 
  // usage: ./mandelbrot 4096 4096 

  int Nre = (argc==3) ? atoi(argv[1]): 4096;
  int Nim = (argc==3) ? atoi(argv[2]): 4096;

  int N = Nre*Nim;
  // Q2b: set the number of threads per block and the number of blocks here:

  //HOST
  float *h_a = (float*) malloc(N*sizeof(float));

  //DEVICE
  float *c_a;

  //CUDAMALLOC
  cudaMalloc(&c_a, N*sizeof(float));

  //CUDAMEMCPY
  cudaMemcpy(c_a, h_a, N*sizeof(float), cudaMemcpyHostToDevice);

  dim3 TPB(16,16,1); // Bx * By threads in thread-block
  dim3 BPG((Nre+15)/16,(Nim+15)/16,1); // Gx * Gy grid of thread-blocks


  // storage for the iteration counts
  float *count;
  count = (float*) malloc(Nre*Nim*sizeof(float));

  // Parameters for a bounding box for "c" that generates an interesting image
  const float centRe = -.759856, centIm= .125547;
  const float diam  = 0.151579;

  complex_t cmin; 
  complex_t cmax;
  complex_t dc;

  cmin.r = centRe - 0.5*diam;
  cmax.r = centRe + 0.5*diam;
  cmin.i = centIm - 0.5*diam;
  cmax.i = centIm + 0.5*diam;

  //set step sizes
  dc.r = (cmax.r-cmin.r)/(Nre-1);
  dc.i = (cmax.i-cmin.i)/(Nim-1);

  clock_t start = clock(); //start time in CPU cycles

  // compute mandelbrot set
  mandelbrot <<< BPG, TPB >>> (Nre, Nim, cmin, dc, count); 

  cudaMemcpy(h_a, c_a, N*sizeof(float), cudaMemcpyDeviceToHost);
  
  // copy from the GPU back to the host here

  clock_t end = clock(); //start time in CPU cycles
  
  // print elapsed time
  printf("elapsed = %f\n", ((double)(end-start))/CLOCKS_PER_SEC);

  // output mandelbrot to ppm format image
  printf("Printing mandelbrot.ppm...");
  writeMandelbrot("mandelbrot.ppm", Nre, Nim, count, 0, 80);
  printf("done.\n");

  free(count);

  exit(0);
  return 0;
}  


/* Output data as PPM file */
void saveppm(const char *filename, unsigned char *img, int width, int height){

  /* FILE pointer */
  FILE *f;
  
  /* Open file for writing */
  f = fopen(filename, "wb");
  
  /* PPM header info, including the size of the image */
  fprintf(f, "P6 %d %d %d\n", width, height, 255);

  /* Write the image data to the file - remember 3 byte per pixel */
  fwrite(img, 3, width*height, f);

  /* Make sure you close the file */
  fclose(f);
}



int writeMandelbrot(const char *fileName, int width, int height, float *img, int minI, int maxI){

  int n, m;

  unsigned char *rgb   = (unsigned char*) calloc(3*width*height, sizeof(unsigned char));
  
  for(n=0;n<height;++n){
    for(m=0;m<width;++m){
      int id = m+n*width;

      int I = (int) (768*sqrt((double)(img[id]-minI)/(maxI-minI)));
      
      // change this to change palette
      if(I<256)      rgb[3*id+2] = 255-I;
      else if(I<512) rgb[3*id+1] = 511-I;
      else if(I<768) rgb[3*id+0] = 767-I;
      else if(I<1024) rgb[3*id+0] = 1023-I;
      else if(I<1536) rgb[3*id+1] = 1535-I;
      else if(I<2048) rgb[3*id+2] = 2047-I;

    }
  }

  saveppm(fileName, rgb, width, height);

  free(rgb);
}


