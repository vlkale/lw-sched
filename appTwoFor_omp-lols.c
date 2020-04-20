/* -- Libraries for C/C++  -- */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <iostream>
using namespace std;

/* -- library for parallelization of code -- */
// #include <pthread.h>

#ifdef MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/* variables to support parallelization of code */

#define DEFAULT_NUMTHREADS 16 // default constant and set value for the number of threads
int numThreads;

/* -- Application specific #defines and variables -- */
int numOuterIters;
int probSize;

#define DEFAULT_NUMOUTERITERS 1
#define DEFAULT_PROBSIZE 500 // Default values based on architecture. For a proper test, the problem size ought to be such that data goes out of cache.                                                               
#define MAX_NUMITERS 100000000
#define MAX_PROBSIZE 16384 // Default values based on architecture. For a proper test, the problem size ought to be such that data goes out of cache.                                                               
/* -- Debugging -- */

#define VERBOSE 1                                                                                                                                  

/* --  Performance Measurement -- */
double totalTime = 0.0;
FILE* myfile;// output file for experimental data                                                                                                                                      

/* --  Hardware Profiling -- */
// #define USE_PAPI                                                                                                                                                                                       
 #ifdef USE_PAPI
 #include <papi.h>  // Leave PAPI out for now to make fully portable. Some platforms don't have the L2 and L3 cache misses available. TODO: need a way to check the counters in config or programmatically.   
 #endif

/* --Library for scheduling strategy and variables and macros associated with the library -- */
#include "vSched.h"

// in the below macros, strat is how we specify the library                                    
#define FORALL_BEGIN(strat, s,e, start, end, tid, numThds )  loop_start_ ## strat (s,e ,&start, &end, tid, numThds);  do {
#define FORALL_END(strat, start, end, tid)  } while( loop_next_ ## strat (&start, &end, tid));

int main (int argc, char** argv );
int i4_min ( int i1, int i2 );
void timestamp ( );

/******************************************************************************/

int main (int argc, char** argv )

/******************************************************************************/
/*
  Purpose

    MAIN is the main program for MANDELBROT_OPENMP.

  Discussion:

    MANDELBROT_OPENMP computes an image of the Mandelbrot set.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Original Modified:

    03 September 2012

  Revision Modified: 
    19 April 2020

  Author:

    John Burkardt

  Revision Author: 

     Vivek Kale 

  Purpose: 
   The revision author uses the original code to demonstrate and test the use low-overhead loop scheduling strategy through library vSched

  Local Parameters:

    Local, int COUNT_MAX, the maximum number of iterations taken
    for a particular pixel.
*/
{
  int m = 500;
  int n = 500;
  
  int b[m][n];
  int c;
  int count[m][n];
  int count_max = 2000;
  int g[m][n];
  int i;
  int j;
  int jhi;
  int jlo;
  int k;
  char const *output_filename = "mandelbrot_openmp.ppm"; // use const to allow for c++ string conversion
  FILE *output_unit;
  int r[m][n];
  double wtime;
  double x_max =   1.25;
  double x_min = - 2.25;
  double x;
  double x1;
  double x2;
  double y_max =   1.75;
  double y_min = - 1.75;
  double y;
  double y1;
  double y2;
  
  #ifdef MPI 
  MPI_Init(&argc, &argv);
  #endif

  //for vSched 
  int numThreads = DEFAULT_NUMTHREADS;
  int threadNum;

  double static_fraction = 0.5;
  double constraint = 0.1;
  int chunk_size = 4;
  
  
  if(argc <= 2) // if user fails to put in minimum args, which are for application domain specific for this test
    {
      //  char userReplyDefault;

      cout << "Usage: testAppTwo_omp_{loopSched} [probSize] [numOuterIters] [chunk_size] <static_fraction> <constraint>" << endl;
      cout << "Usage(cont'd): where {loopSched} is the implementation strategy or library you use, e.g., vSched's low-overhead scheduling, OpenMP rtl's low-overhead scheduling" << endl;

      //      cout << "Use defaults? [y/N]" << endl ;
      //  cin << varName;
      cout << "continuing with default problem size and app parameters." << endl;
      probSize = DEFAULT_PROBSIZE;
      numOuterIters = DEFAULT_NUMOUTERITERS; 
    }

  //   else // required args

   else
    {
      probSize = atoi(argv[1]);
      numOuterIters = atoi(argv[2]);
    }
  if (argc > 3) chunk_size = atoi(argv[3]);  
  if (argc > 4) static_fraction = atof(argv[4]);
  if (argc > 5) constraint = atof(argv[5]);


  // set number of threads to input 
  //omp_set_num_threads(numThreads);

#pragma omp parallel
  {

    #pragma omp master 
    cout << "Number of threads is : " << omp_get_num_threads() << endl;
    
#ifdef USE_VSCHED
    numThreads = omp_get_num_threads();
#endif
  }

#ifdef USE_VSCHED
  vSched_init(numThreads);
#endif

  timestamp ( );
  printf ( "\n" );
  printf ( "MANDELBROT_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  Create an ASCII PPM image of the Mandelbrot set.\n" );
  printf ( "\n" );
  printf ( "  For each point C = X + i*Y\n" );
  printf ( "  with X range [%g,%g]\n", x_min, x_max );
  printf ( "  and  Y range [%g,%g]\n", y_min, y_max );
  printf ( "  carry out %d iterations of the map\n", count_max );
  printf ( "  Z(n+1) = Z(n)^2 + C.\n" );
  printf ( "  If the iterates stay bounded (norm less than 2)\n" );
  printf ( "  then C is taken to be a member of the set.\n" );
  printf ( "\n" );
  printf ( "  An ASCII PPM image of the set is created using\n" );
  printf ( "    M = %d pixels in the X direction and\n", m );
  printf ( "    N = %d pixels in the Y direction.\n", n );
    
  omp_sched_t schedule;
  int chunksize = chunk_size;  
  omp_get_schedule(&schedule, &chunksize);  
  printf ( "OpenMP schedule: OMP_SCHEDULE=%d\t%d\n", (int) schedule, chunksize);
  
   // chunk start and end indices for vSched runtime to set and retrieve 
  int startInd;
  int endInd;
  // upper bound of loop that gets parallelized by OpenMP, used by vSched
  int probSize; 
  probSize = m; // set to loop bound of loop that gets parallelized by OpenMP
 
  wtime = omp_get_wtime ( );
  
/*
  Carry out the iteration for each pixel, determining COUNT.
*/

  for (int iter = 0; iter < numOuterIters; iter++)
    {
# pragma omp parallel						   \
  shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) \
  private ( i, j, k, x, x1, x2, y, y1, y2 )
{
  // # pragma omp for

   //TODO: figure out how to better template the scheduling strategy name in the below lines of code for both library implmentations
   // # pragma omp for schedule(user:statdynstaggered, &lr) collapse(2)  // prototype UDS placeholder
   // The first parameter is the loop scheduling strategy.
  threadNum = omp_get_thread_num();
  numThreads = omp_get_num_threads();
#ifdef USE_VSCHED

  setCDY(static_fraction, constraint, chunk_size); // set parameter of scheduling strategy for vSched
  FORALL_BEGIN(statdynstaggered, 0, probSize, startInd, endInd, threadNum, numThreads)   
#ifdef VERBOSE
   if(VERBOSE==1) printf("Thread [%d] : iter = %d \t startInd = %d \t  endInd = %d \t\n", threadNum,iter, startInd, endInd);
#endif
#else
  #ifdef VERBOSE
  if(VERBOSE==1) printf("Thread [%d] : iter = %d executing a chunk \n", threadNum,iter);
#endif
  #pragma omp for schedule(guided, chunk_size)
#endif
  { // add parens to show comparison with forall macro


  for ( i = 0; i < m; i++ )
  {
    y = ( ( double ) (     i - 1 ) * y_max   
        + ( double ) ( m - i     ) * y_min ) 
        / ( double ) ( m     - 1 );
    
    for ( j = 0; j < n; j++ )
    {
      x = ( ( double ) (     j - 1 ) * x_max   
          + ( double ) ( n - j     ) * x_min ) 
          / ( double ) ( n     - 1 );

      count[i][j] = 0;

      x1 = x;
      y1 = y;

      for ( k = 1; k <= count_max; k++ )
      {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2 * x1 * y1 + y;

        if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
        {
          count[i][j] = k;
          break;
        }
        x1 = x2;
        y1 = y2;
      }

      if ( ( count[i][j] % 2 ) == 1 )
      {
        r[i][j] = 255;
        g[i][j] = 255;
        b[i][j] = 255;
      }
      else
      {
        c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt ( 
          ( ( double ) ( count[i][j] ) / ( double ) ( count_max ) ) ) ) ) );
        r[i][j] = 3 * c / 5;
        g[i][j] = 3 * c / 5;
        b[i][j] = c;
      }
    }
  }
  
     } //End parallelized for block

#ifdef USE_VSCHED
 FORALL_END(statdynstaggered, startInd, endInd, threadNum)
#endif


   
   } // end parallel

    } // end outer iter loop 

wtime = omp_get_wtime ( ) - wtime;
printf ( "\n" );
printf ( "  Time = %g seconds.\n", wtime );
/*
  Write data to an ASCII PPM file.
*/
  output_unit = fopen ( output_filename, "wt" );
  
  fprintf ( output_unit, "P3\n" );
  fprintf ( output_unit, "%d  %d\n", n, m );
  fprintf ( output_unit, "%d\n", 255 );

  for ( i = 0; i < m; i++ )
  {
    for ( jlo = 0; jlo < n; jlo = jlo + 4 )
    {
      jhi = i4_min ( jlo + 4, n );
      for ( j = jlo; j < jhi; j++ )
      {
        fprintf ( output_unit, "  %d  %d  %d", r[i][j], g[i][j], b[i][j] );
      }
      fprintf ( output_unit, "\n" );
    }
  }

  fclose ( output_unit );
  printf ( "\n" );
  printf ( "  Graphics data written to \"%s\".\n", output_filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MANDELBROT_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  #ifdef USE_VSCHED
  vSched_finalize(numThreads);
  #endif 

  #ifdef MPI
  MPI_Finalize(MPI_COMM_WORLD);
  #endif 

  return 0;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
