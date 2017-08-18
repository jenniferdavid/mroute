/* -------------------------------------------------------------------- 
   
     filename: random.cc  	program: mroute  	project: mroute
     written at Theoretical Physics 2, Lund University
  									
     this file contains the random stuff definitions.
  									
   -------------------------------------------------------------------- */

#include "random.hh"

#include <sys/time.h>
#include <math.h>

extern "C" int random(void);                    /*     long random();    */
extern "C" int srandom( unsigned int seed );   /* void srandom(int seed);    */




double Rnd() {			// !! change to something from limits.h
  return ( random() % 1073741824 ) / 1073741824.0 ;
}



double Grand() {
  return sqrt( -2 * log(Rnd())) * cos( M_PI * Rnd() );
}



int Randomize(int seed/*=-1*/) {
  
  struct timeval tp; /* { long tv_sec; long tv_usec; } */

  if (seed==-1) {
    gettimeofday(&tp,0);
    seed = tp.tv_sec;
  }
  srandom(seed);

  return seed;
}



