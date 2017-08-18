/*

  filename: mroute.hh

  */

#ifndef mroute_hh
#define mroute_hh

#define FALSE 0
#define TRUE 1

//#define DEBUG TRUE
#define DEBUG FALSE

class Manswer;
class Mreq; 
class Net;
class Option;

typedef void Alg (Mreq&, Manswer&, const Option&);

// problem solving algorithms
Alg bellman;
Alg NN_alg;
Alg solveexact;
Alg seq_bellman;

void plot(int,Net*,Mreq*,Manswer*);


// includes needed for the STL stuff
// STL headers
#include <list>
typedef std::list<int>	Path;
typedef std::list<Path>	ManyPaths;

#endif
