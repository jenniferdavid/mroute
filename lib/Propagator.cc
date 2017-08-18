/********************************/
/*      Propagator.cc           */
/*      Oct -96 Martin  	*/
/********************************/

#include "Propagator.hh"

#define small 1.e-13
#define onebysmall 1.e13
#define onebysmall3 1.e39

using namespace std;
Propagator::Propagator(int n, PottsNeuron* v)
{
  int warning = 0;
  int i,j;
  n_n = n;
  P = new double*[n_n];
  for (i = 0; i <n_n; i++)
    P[i] = new double[n_n];
  
  // init P[][] as if v == 0
  for (i = 0; i < n_n; i++) 
    {
      for (j = 0; j < n_n; j++) 
	P[i][j] = 0;
      P[i][i] = 1;
    }

  // use the update function to get the "real" 1/(1-v)
  for (i = 0; i < n_n; i++) 
    warning += update(v[i]);
}

Propagator::~Propagator()
{
  if (P)
    {
      for (int i = 0; i < n_n; i++) 
	delete [] P[i];
      delete [] P;
    }
}

int Propagator::update(const PottsNeuron& dv) // subroutine to run, r=req#, i=node#
{
  double R[n_n];
  double L = 0;
  double aux;
  int i,j,k;
  int warning = 0;
  int row = dv.label();
  
  for (j = 0; j < n_n; j++)
    {
      R[j] = 0;
      for (k = 0; k < dv.n(); k++)
	R[j] += dv[k] * P[dv.ne(k)][j];	// - bec: (1-Vnew) - (1-Vold) = -dv
    }
  
  //  cout << " row = " << row << endl;
  double lambda = 1. - R[row];	
  
  // to avoid overflow & core dump if small numbers
  if ( lambda < small ) // Warns also for negative lambda! /BS
    {
      warning++;
      lambda = small; // => wrong P[][] ! /BS
      cerr << "[Propagator]: Avoiding OF?" << endl;
      if( lambda < 0. )
	cerr << "Negative lambda!" << endl;
    }
  if ( lambda > 1/small ) // To avoid underflow !!!??
    {
      cerr << "[Propagator]: Avoiding UF?" << endl;
      warning++;
      lambda = 1/small; // => wrong P[][] ! /BS
    }
  
  for (i = 0; i < n_n; i++) 
    {
      L = P[i][row] / lambda;
      aux = L * R[i];
      if( aux > onebysmall )
	{
	  cerr << "WARNING [Propagator]: Pii large (i = " << i << ")" << endl;
	  if( aux > onebysmall3 )
	    {
	      cerr << "Pii very large -- exiting" << endl;
	      exit(-1);
	    }
	}
      for (j = 0; j < n_n; j++) 
	P[i][j] += L * R[j];
    }
  return warning;
}

// friends

ostream& operator<<(ostream& s,const Propagator& P)
{
  int i,j;	// loop dummy

  s.setf(ios::fixed);
  s.precision(3);

  s << P.rlabel<< " P[][]:\n";
  for (i = 0; i < P.n_n; i++)
    {
      for (j = 0; j < P.n_n; j++)
	s << "\t" << P.P[i][j];
      s << "\n";
    }
  return s;
}
