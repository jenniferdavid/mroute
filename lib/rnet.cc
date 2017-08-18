/********************************/
/*	rnet.cc		*/
/*	Extension needed for	*/
/*	reducing problem.	*/
/*	Bo S 10/10-96		*/
/********************************/

#include <iostream>
#include <cstdlib>
#include "net.hh"
#include "itoa.hh"
#include "rnet.hh"

using namespace std;

// *********************************************
RNet::~RNet()
{
if( DEBUG ) cerr << "Entering Rnet::~RNet()" << endl;
//  delete &H;
  if(S) delete [] S;
if( DEBUG ) cerr << "Leaving Rnet::~RNet()" << endl;
}

// *********************************************
RNet::RNet(Net & net) : M(&net)
{
  int i, j, k, l, m, n, i1, i2;
  
  int nL = net.get_nl();
  int nN = net.get_nn();

  int n_newnode, newnode[nN], used[nN];
  int mult;

  int nS;
  int *iNS[nN], nNS[nN];
  int *iSN[nN], nSN[nN];

  int iNL[nL][2], iSL[nL];
  int nLN[nN], iLN[nN][nL];

  int nC = 0, iPN[nN][nN];

  n_S = 0; // data member
  name = net.get_name() + ":R";

  // --------------------------

  if( (i = net.isConnected()) < 2) // 0 disconn, 1 semi-conn, 2 conn.
    {
      if(i == 1)
	cerr << "\7Full net only semi-connected!" << endl;
      if( i == 0 )
	cerr << "\7Full net not connected!" << endl;
      cerr << "Exiting." << endl;
      exit(-1);
    }
  
  if( !net.isNormal() )
    {
      cerr << "\7RNet::RNet(const Net*): Net abnormal! Exiting." << endl;
      exit(-1); //this = 0;
    }

  // setup nLN[], iLN[][], iNL[]
  
  for( i = 0; i < nN; i++)
    nLN[i] = 0;

  for( i = 0; i < nL; i++)
    {
      i1 = net.get_link(i)->get_node(0);
      i2 = net.get_link(i)->get_node(1);
      iNL[i][0] = i1;
      iNL[i][1] = i2;
      iLN[i1][nLN[i1]++] = i;
    }

  // CHECK FOR NODE REDUCIBILTY IN FULL NET !!!!!
  // (skip link reduction)
  //      __.__
  // e.g. \/ \/
  //
  // => find those nodes that if removed disconnects the net
  //   for each node N:
  //     check if net\N connected

// loop over all nodes to check if reducible
  if( nN == 1 ) // net = 1 free node
    {
      cerr << "Net has a single node" << endl;
      exit(-1);
    }

  // > 1 node in net:
  for( i1 = 0; i1 < nN; i1++)
    {
      // flag all nodes as unused
      for( j = 0; j < nN; j++)
	used[j] = 0;

      // have i1 (checked)
      // flag??
      n_newnode = 0;
      used[newnode[n_newnode++] = i1] = -1;
      mult = 0;
      while( n_newnode < nN )
	{
	  mult++;
	  // pick an unused ref. node i2 in subnet i
	  i2 = -1;
	  for( j = 0; j < nN; j++)
	    if( !used[j] )
	      {
		i2 = j;
		used[newnode[n_newnode++] = i2] = mult;
		break;
	      }
	  if( i2 == -1 )
	    {
	      cerr << "Error (RNet::RNet(...)): i2 not found!" << endl;
	      exit(-1);
	    }
	  
	  for( j = n_newnode - 1 /* don't use the old ones! */; j < n_newnode; j++)
	    {
	      k = newnode[j];
	      // get neighbours
	      for( l = 0; l < nLN[k] ; l++)
		{
		  m = iNL[iLN[k][l]][1];
		  if( !used[m] )
		    {
		      newnode[n_newnode++] = m;
		      used[m] = mult;
		    }
		}
	    }
	} // while
      
      if( mult == 1 ) // not reducible
	continue;

      // Still here? Reducible.
      
      used[i1] = 1 - mult;
      for( j = 0; j < nN; j++)
	iPN[j][nC] = used[j] - 1;
      nC++;
    }

  // get patterns nP, iP[nP][nC]
  int nP;
  nP = 1;
  for( i = 0; i < nN; i++)
    {
      for( j = 0; j < nC; j++)
	if( iPN[i][j] < 0 )
	  nP += -1 - iPN[i][j];
    }

  int iP[nP][nC];
  int iP0[nC];

  // find patterns
  nP = 0;
  for( i = 0; i < nN; i++)
    {
      i1 = 0;
      mult = 1;
      for( j = 0; j < nC; j++)
	{
	  iP0[j] = iPN[i][j];
	  if( iP0[j] < 0 )
	    {
	      i1 = j;
	      mult = -iP0[j];
	    }
	}
      // i1 = poss. critical iN; mult = multiplicity
      for( k = 0; k < mult; k++)
	{
	  if( mult > 1 )
	    iP0[i1] = k;
	  // compare to old patterns
	  i2 = 1;
	  for( j = 0; j < nP; j++)
	    {
	      m = 1; // found flag
	      for( l = 0; l < nC; l++)
		if( iP[j][l] != iP0[l] )
		  {
		    m = 0;
		    break;
		  }
	      if( !m )
		continue;
	      // Still here? Found.
	      i2 = 0;
	      break;
	    }
	  if( i2 ) // not found: new pattern
	    {
	      for( j = 0; j < nC; j++)
		iP[nP][j] = iP0[j];
	      nP++;
	    }
	}
    }


  // OK, redo subnet grouping, allow for virtual nodes and links
  for( i = 0; i < nN; i++)
    used[i] = 0;

  // preliminary pattern matching: find nNS[], nSN[]
  for( k = 0; k < nP; k++) // patterns
    {
      nNS[k] = 0;
      // find matches
      for( i = 0; i < nN; i++) // nodes: match pattern?
	{
	  mult = 1;
	  for( j = 0; j < nC; j++)
	    if( iPN[i][j] != iP[k][j] && iPN[i][j] >= 0 ) // no match
	      {
		mult = 0;
		break;
	      }
	  if( mult ) // match
	    {
	      used[i]++;
	      nNS[k]++;
	      //	      iNS[k][nNS[k]++] = i;
	    }
	}
    }
  

  nS = nP;
  
  // define iSN[iN][jS], nSN[iN]
  for( i = 0; i < nN; i++)
    {
      nSN[i] = 0; // = used[i]; increment later
      iSN[i] = new int [used[i]];
    }

  // redefine nS, iNS
  for( i = 0; i < nS; i++)
    iNS[i] = new int[nNS[i]];

  // final pattern matching: do group-node assignment
  for( k = 0; k < nS; k++) // patterns
    {
      nNS[k] = 0;
      // find matches
      for( i = 0; i < nN; i++) // nodes: match pattern?
	{
	  mult = 1;
	  for( j = 0; j < nC; j++)
	    if( iPN[i][j] != iP[k][j] && iPN[i][j] >= 0 ) // no match
	      {
		mult = 0;
		break;
	      }
	  if( mult ) // match
	    {
	      iSN[i][nSN[i]++] = k;
	      iNS[k][nNS[k]++] = i;
	    }
	}
    }

  int nLS[nS], iLS[nS][nL];

  for( i = 0; i < nS; i++)
    nLS[i] = 0;
  
  for( i = 0; i < nL; i++)
    {
      i1 = iNL[i][0];
      i2 = iNL[i][1];
      for( j = 0; j < nSN[i1]; j++)
	{
	  m = -1;
	  n = iSN[i1][j];
	  for( k = 0; k < nSN[i2]; k++)
	    if( n == iSN[i2][k] )
	      {
		m = n;
		break;
	      }
	  if( m >= 0 )
	    break;
	}
      if( m < 0 )
	{
	  cerr << "common S m not found" << endl;
	  exit(-1);
	}
      iSL[i] = m;
      iLS[m][nLS[m]++] = i;
    }

  // OK
  // i) make a super-net, containing the nN nodes + the nS subnets as nodes, and
  //     links connecting the nodes with the subnets where shadows of it appears
  // ii) make each sub-net, with its contained nodes and links 

  n_S = nS;
  // M = &net; is set
  // Net::init(int nN, int nL, int *Nlbl, int (*iNL)[2], double *Lw, double *Lcap);

  // need nNH=nN+nS, nLH, Nlbl[], Llbl[][2], Lw = 0, Lcap = 0) 
  int nLH = 0;
  for( i = 0; i < nN; i++)
    nLH += nSN[i];
  nLH *= 2;
  int ** jNL = new int * [nLH];
  for( i = 0; i < nLH; i++)
    jNL[i] = new int [2];
  int * jN = new int [nN + nS];

  for( i = 0; i < nN; i++)
    jN[i] = i;
  for( i = 0; i < nS; i++)
    jN[i + nN] = i;

  k = 0;
  for( i = 0; i < nN; i++)
    for( j = 0; j < nSN[i]; j++)
      {
	jNL[k][0] = i;
	jNL[k][1] = nN + iSN[i][j];
	k++;
	if( k > nLH )
	  cerr << "\7k > nLH" << endl;
      }
  for( i = 0; i < nS; i++)
    for( j = 0; j < nNS[i]; j++)
      {
	jNL[k][0] = nN + i;
	jNL[k][1] = iNS[i][j];
	k++;
	if( k > nLH )
	  cerr << "\7k > nLH" << endl;
      }
  
  // Shoot!
  H.init(net.get_name() + ":H", nN + nS, nLH, jN, jNL, 0, 0);

  for( i = 0; i < nLH; i++)
    delete [] jNL[i];
  delete [] jNL;
  delete [] jN;
  
  S = new Net[n_S];

  int j1, j2;
  for( i = 0; i < n_S; i++)
    // need nNS[i], nLS[i], iNS[i][], iNL[iLS[i][]][2], Lw[], Lcap[]) 
    {
      jNL = new int* [nLS[i]];
      for( j = 0; j < nLS[i]; j++)
	jNL[j] = new int [2];
      double * Lw = new double [nLS[i]];
      int * Lcap = new int [nLS[i]];
      for( j = 0; j < nLS[i]; j++)
	{
	  l = iLS[i][j];
	  j1 = iNL[l][0];
	  j2 = iNL[l][1];
	  i1 = i2 = -1;
	  for( k = 0; k < nNS[i]; k++)
	    {
	      if( iNS[i][k] == j1 )
		i1 = k;
	      if( iNS[i][k] == j2 )
		i2 = k;
	      if( i1 >= 0 && i2 >= 0 )
		break;
	    }
	  if( i1 == -1 || i2 == -1 )
	    {
	      cerr << "\7i1 or i2 not found" << endl;
	      exit(-1);
	    }
	  jNL[j][0] = i1;
	  jNL[j][1] = i2;
	  Lw[j] = net.get_link(l)->getweight();
	  Lcap[j] = net.get_link(l)->get_capacity();
	}

      // Shoot!
      S[i].init(net.get_name() + ":S" + itoa(i), nNS[i], nLS[i], iNS[i], jNL, Lw, Lcap);

      for( j = 0; j < nLS[i]; j++)
	delete [] jNL[j];
      delete [] jNL;
      delete [] Lw;
      delete [] Lcap;
    }

  // clean up

  for( i = 0; i < nS; i++)
    delete [] iNS[i];

  for( i = 0; i < nN; i++)
    delete [] iSN[i];

  // Also mark original net as reduced, so things can be computed via subnets?
}

// *********************************************
int RNet::get_NSpath(int a, int b, int*& p) const
{
 int i, j, k, l, l1;

 if( a >= M->get_nn() ||
     a < 0 ||
     b >= M->get_nn() ||
     b < 0 )
   {
     cerr << "\7Illegal a, b" << endl;
     return 0;
   }

 int n = H.get_nn();
 int ln[n], nn[n];
 
 // mini-Bellman-Ford
 for( i = 0; i < n; i++)
   {
     ln[i] = n+1;
     nn[i] = i;
   }
 ln[b] = 0;
 int flag = 1;
 while(flag)
   {
     flag = 0;
     for( i = 0; i < n; i++)
       {
	 l = ln[i];
	 for( k = 0; k < H.get_node(i)->get_nl(); k++)
	   {
	     j = H.get_node(i)->node(k);
	     l1 = ln[j] + 1;
	     if( l1 < l )
	       {
		 l = ln[i] = l1;
		 nn[i] = j;
		 flag = 1;
	       }
	   }
       }
   }

 l = ln[a] + 1;
 p = new int[l];
 k = a;
 p[0] = a;
 for( i = 1; i < l; i++)
   {
     k = nn[k];
     p[i] = k;
   }
 return l;
}

// *********************************************
int RNet::get_global(int n, int s) const
{
  if( s >= n_S || s < 0 )
    {
      cerr << "\7Illegal s" << endl;
      return -1;
    }
  if( n >= S[s].get_nn() || n < 0 )
    {
      cerr << "\7Illegal n" << endl;
      return -1;
    }
  return S[s].get_node(n)->getlabel();
}

// *********************************************
int RNet::get_local(int n, int s) const
{
  if( n >= M->get_nn() || n < 0 )
    {
      cerr << "\7Illegal n" << endl;
      return -1;
    }
  if( s >= n_S || s < 0 )
    {
      cerr << "\7Illegal s" << endl;
      return -1;
    }

  int l = -1;
  for( int j = 0; j < S[s].get_nn(); j++)
    if( n == S[s].get_node(j)->getlabel() )
      {
	l = j;
	break;
      }
  return l;
}

// *********************************************
int RNet::get_SNN(int i, int j) const
{
  Node *a, *b;
  int k, l, na, nb, n;

  if( i >= M->get_nn() ||
      i < 0 ||
      j >= M->get_nn() ||
      j < 0 )
    {
      cerr << "\7Illegal i, j" << endl;
      return -1;
    }

  a = H.get_node(i);
  b = H.get_node(j);
  n = -1;
  for( k = 0; k < a->get_nl(); k++)
    {
      na = a->node(k);
      for( l = 0; l < b->get_nl(); l++)
	{
	  nb = b->node(l);
	  if( nb == na )
	    {
	      n = na;
	      break;
	    }
	}
      if( n > -1 )
	break;
    }
  if( n > -1 )
    return H.get_node(n)->getlabel();
  else
    return -1;
}

// *********************************************
ostream& operator<< (ostream& s, const RNet &r)
{
  s << r.name << endl;
  s << " based on: " << r.M->get_name() << endl;
  s << " HyperNet: " << r.H;
  s << " with " << r.n_S << " subnets:" << endl;
  for( int i = 0; i < r.n_S; i++)
    s << " Subnet[" << i << "]: " << r.S[i];    
  return s;
}
