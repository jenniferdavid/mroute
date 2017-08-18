/********************************/
/*      Propagator.hh           */
/*          Oct -96 Martin	*/
/********************************/

#ifndef Propagator_hh
#define Propagator_hh

#include <stdlib.h>
#include <iostream>
#include "PottsNeuron.hh"


class Propagator
{
public:
  Propagator(int n, PottsNeuron* v);
  ~Propagator();
  int 		update(const PottsNeuron& dv);// subroutine to run, r=req#, i=node#
  inline int 	n() const;		// returns the n in the nxn P matrix
  inline int 	label() const;		// returns the label (index) of this neuron 
  inline double g(int i , int j) const; // returns element i,j
  inline int 	s(int i , int j, double value); // sets element i,j
  friend std::ostream& operator<<(std::ostream& s,const Propagator& P);
private:
  double** P;
  int n_n;
  int rlabel;
};

// inline functions:
inline int Propagator::n() const { return n_n; }
inline int Propagator::label() const { return rlabel; }
inline double Propagator::g(int i , int j) const 
{ 
  if ( i < n_n && i >=0 && j < n_n && j >= 0) // Remove check later !!
    return P[i][j]; 

  std::cerr << "ERROR: P index out of bounds" << std::endl;
  exit(-1);
}
inline int Propagator::s(int i , int j, double value) // sets element i,j
{
  if ( i < n_n && i >=0 && j < n_n && j >= 0) // Remove check later !!
    P[i][j] = value;
  
  return 0;
  std::cerr << "ERROR: P index out of bounds" << std::endl;
  exit(-1);
}
#endif
