// rmreq.hh
// reduced mreq class
// BS, nov -96

#include "rmreq.hh"
#include <iostream>
#include "stdio.h"
#include <cstdlib>

using namespace std;

// *********************************************
RMreq::RMreq(RNet& RN_in, Mreq& MR_in)
  : name(MR_in.get_name() + ":R")
  , RN(&RN_in)
  , OMR(&MR_in)
{
  int i, j, nn, nr;
  int a, b, len;

  nMR = RN->get_n_S();
  MR = new Mreq [nMR];

  nn = RN->get_M()->get_nn();
  nr = OMR->get_nr();

  int *path, a_glob, b_glob, a_loc, b_loc, ss;
  int nrr[nMR], arr[nMR][nr], brr[nMR][nr];

  for( i = 0; i < nMR; i++)
    nrr[i] = 0;

  nleg = new int [nr];
  iMRleg = new int * [nr];
  iRleg = new int * [nr];
  for( i = 0; i < nr; i++)
    {
      a = OMR->get_startnode(i);
      b = OMR->get_endnode(i);
      len = RN->get_NSpath(a, b, path); // skapar path[]; must be explicitly deleted
      if( path[0] != a || path[len-1] != b )
	{
	  cerr << "\7RMreq::RMreq(): erroneous path <1>" << endl;
	  exit(-1);
	}

      nleg[i] = (len - 1) / 2;
      iMRleg[i] = new int [nleg[i]];
      iRleg[i] = new int [nleg[i]];

      for( j = 0; j < nleg[i]; j++)
	{
	  a_glob = path[2 * j];
	  ss = path[2 * j + 1] - nn;
	  b_glob = path[2 * j + 2];
	  if( ss != RN->get_SNN(a_glob,b_glob) )
	    {
	      cerr << "\7RMreq::RMreq(): erroneous path <2>" << endl;
	      exit(-1);
	    }
	  a_loc = RN->get_local(a_glob, ss);
	  b_loc = RN->get_local(b_glob, ss);

	  arr[ss][nrr[ss]] = a_loc;
	  brr[ss][nrr[ss]] = b_loc;

	  iMRleg[i][j] = ss;
	  iRleg[i][j] = nrr[ss];	  

	  nrr[ss]++;
	}
      delete [] path;
    }

  // OK, make subproblems MR
  for( i = 0; i < nMR; i++)
    MR[i].init("", *RN->get_S(i), *OMR->get_op(), nrr[i], (int *) arr[i], (int *) brr[i]);

  // Also mark original mreq as reduced, so entropy etc. can be computed via subreqs?
}


// *********************************************
RMreq::~RMreq()
{
  delete [] MR;
  for(int i = 0; i < nMR; i++)
    {
      delete [] iMRleg[i];
      delete [] iRleg[i];
    }
  delete [] iMRleg;
  delete [] iRleg;
  delete [] nleg;
}

// ****************************************
ostream& operator<< (ostream& s, RMreq& n)
{
  int i;	// loop dummy

  s << n.name << endl;
  s << " acting on " << n.RN->get_name() << endl;
  s << " based on " << n.OMR->get_name() << endl;
  s << " having " << n.nMR << " parts:\n";
  for (i = 0; i < n.nMR; i++)
    s << " MR[" << i << "]: " << n.MR[i];
  return s;
}

