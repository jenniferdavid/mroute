/********************************/
/*      PottsNeuron.cc		*/
/*       Oct -96 Martin		*/
/********************************/

#include "PottsNeuron.hh"
#include "random.hh"
#include <iostream>
using namespace std;

PottsNeuron::PottsNeuron(int n_components, int* neighbourvec, int* linkvec, int labeltmp/* = -1*/)
{
  n_c = n_components;
  this_index = labeltmp;
  neighbour_index = new int[n_c];
  link_index = new int[n_c];
  value = new double[n_c];
  double tmp = 1./n_c;
  for (int i = 0; i < n_c; i++)
    {
      value[i] = tmp + 0.02*(Rnd() - 0.5);	// +-1 % noice
      neighbour_index[i] = neighbourvec[i];	
// save above and below somewhere else and use a pointer instead, a lot of unneccery inf. transported around (grows with n_r). !!
      link_index[i] = linkvec[i];
    }
}

PottsNeuron::PottsNeuron(PottsNeuron& PN, int labeltmp/* = -1*/)
{
  n_c = PN.n_c;
  if (labeltmp != -1)
    this_index = labeltmp;
  else
    this_index = PN.this_index;
  neighbour_index = new int[n_c];
  link_index = new int[n_c];
  value = new double[n_c];
  for (int i = 0; i < n_c; i++)
    {
      value[i] = PN.value[i];
      neighbour_index[i] = PN.neighbour_index[i] ;
      link_index[i] = PN.link_index[i];
    }
}

PottsNeuron::PottsNeuron()
{
  //nedded for the contruction: 
  // PottsNeuron v[10];
  // v[3] = PottsNeuron(....... != tom);
  // all needed to avoid pointers and reason for that is want 
  // be able to add as v[5] = v[3] + v[2] , i.e. not *(v[3]) + ...
  n_c = -1;
  this_index = -1;
  neighbour_index = NULL;
  link_index = NULL;
  value = NULL;
}

PottsNeuron::~PottsNeuron()
{
  delete [] value;
  delete [] neighbour_index;
  delete [] link_index;
}

PottsNeuron& PottsNeuron::operator=(const PottsNeuron& PN)
{
  delete [] value;
  delete [] neighbour_index;
  delete [] link_index;
  n_c = PN.n_c;
  this_index = PN.this_index;
  neighbour_index = new int[n_c];
  link_index = new int[n_c];
  value = new double[n_c];
  for (int i = 0; i < n_c; i++)
    {
      value[i] = PN.value[i];
      neighbour_index[i] = PN.neighbour_index[i] ;
      link_index[i] = PN.link_index[i];
    }
  return *this;
}

PottsNeuron& PottsNeuron::operator+=(const PottsNeuron& PN)
{
//   int flag = 0;
//   if ( n_c != PN.n_c)
//     flag = 1;
//   for (int i = 0; i < n_c; i++)
//     {
//       if (neighbour_index[i] != PN.neighbour_index[i] ||
// 	  link_index[i] != PN.link_index[i])
// 	flag = 1;
//     }
//   if (flag)
//     {
//       cerr<<"PottsNeuron: Error: object do not";
//       cerr<<" have the same neighbours and/or links or is not the same size";
//       cerr << " Nothing done" << endl;
//     }
//   else
    for (int i = 0; i < n_c; i++)
      value[i] += PN.value[i];
  return *this;
}

PottsNeuron& PottsNeuron::operator-=(const PottsNeuron& PN)
{
//   int flag = 0;
//   if ( n_c != PN.n_c)
//     flag = 1;
//   for (int i = 0; i < n_c; i++)
//     {
//       if (neighbour_index[i] != PN.neighbour_index[i] ||
// 	  link_index[i] != PN.link_index[i])
// 	flag = 1;
//     }
//   if (flag)
//     {
//       cerr<<"PottsNeuron: Error: object do not";
//       cerr<<" have the same neighbours and/or links or is not the same size";
//       cerr << " Nothing done" << endl;
//     }
//   else
    for (int i = 0; i < n_c; i++)
      value[i] -= PN.value[i];
  return *this;
}

std::ostream& operator<<(std::ostream& s,PottsNeuron& PN)
{
  int i;	// loop dummy
  
  s.setf(std::ios::fixed);
  s.precision(4);
  
  s <<"PN "<< PN.label() <<"  neighbour link \t value:\n";
  for (i = 0; i < PN.n(); i++)
    s << "\t" << PN.ne(i) << "\t " << PN.link(i) << "\t " << PN[i] << "\n";
//   for (i = 0; i < PN.n(); i++)
//     s << "\t" << PN[i] ;
//   s << "\t";
//   for (i = 0; i < PN.n(); i++)
//     s << "\t"<< PN.ne(i);
//   s << "\t";
//   for (i = 0; i < PN.n(); i++)
//     s << "\t"<< PN.link(i);
  return s;
}

// help functions, that do not need to be "close" to the class PottsNeuron

PottsNeuron operator+(PottsNeuron& L,PottsNeuron& R)
{
  PottsNeuron tmp(L);
  tmp += R;
  return (tmp);
}

PottsNeuron operator-(PottsNeuron& L, PottsNeuron& R)
{
  PottsNeuron tmp(L);
  tmp -= R;
  return (tmp);
}
