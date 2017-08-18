/********************************/
/*      PottsNeuron.hh		*/
/*          Oct -96 Martin	*/
/********************************/

#ifndef PottsNeuron_hh
#define PottsNeuron_hh

#include <iostream>
#include <stdlib.h>

class PottsNeuron
{
public:
  PottsNeuron(int n_components, int* neighbourvec, int* linkvec, int labeltmp = -1);
  PottsNeuron(PottsNeuron& PN, int labeltmp = -1);
  PottsNeuron();
  ~PottsNeuron();

  inline int	n() const;		// returns the number of components, i.e. #local links
  inline int	ne(const int i) const;	// global index of neighbour neuron i
  inline int	link(int i) const;	// link index of the link to ne(i)
  inline int	label() const;	// returns the label (index) of this neuron 

  inline double& operator[] (const int i) const;
  PottsNeuron&	operator= (const PottsNeuron& PN);
  PottsNeuron&	operator+= (const PottsNeuron& );
  PottsNeuron&	operator-= (const PottsNeuron& );   
  
private:
  int this_index;
  int* neighbour_index;
  int* link_index;
  double* value;		// the components
  int n_c;		// the number of components
};
// help functions, that do not need to be "close" to the class PottsNeuron 
std::ostream& operator<<(std::ostream& s,PottsNeuron& PN); 
PottsNeuron operator+(PottsNeuron& L, PottsNeuron& R);
PottsNeuron operator-(PottsNeuron& L, PottsNeuron& R);
// do not want as a friend bec want to use the public funct. as the rest of the program must

// inline functions:
inline int PottsNeuron::label() const { return this_index; }
inline int PottsNeuron::n() const {  return n_c; }
inline int PottsNeuron::ne(const int i) const
{ 
  //  if ( i > n_c || i < 0)	// remove when it seems to work
  //  cerr << "Error:PottsNeuron::ne(i), argument i out of range" << endl;
  //if ( neighbour_index[i] == -1)	// remove when it seems to work
  //  cerr << "Error:PottsNeuron::ne(i) No global neighbour index set" << endl;
  return neighbour_index[i]; 
}

inline int PottsNeuron::link(int i) const
{ 
  //if ( link_index[i] == -1)	// remove when it seems to work
  //  cerr << "Error: No link index set" << endl;
  return link_index[i]; 
}


inline double& PottsNeuron::operator[] (const int i) const
{
  // if ( i >= n_c || i < 0)	// remove when it seems to work
  //  cerr << "Error: index out of bounds, this_index="<<this_index<<" n_c-1="<<n_c-1<<" i="<<i<<endl;
  return  value[i];
}
#endif

