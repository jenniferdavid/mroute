/************************************/
/*      mreq.cc                     */
/************************************/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include "stdio.h"
#include <cstdlib>
#include "mreq.hh"
#include "mroute.hh"
#include "net.hh"
#include "Option.hh"
#include "random.hh"
#include "itoa.hh"

using namespace std;

/***************************************************************************/
Mreq::Mreq() // default (bare) constructor
  : name("<null>"), energybound(-1), solvable(-1)
{
}

/***************************************************************************/
Mreq::Mreq(Option& opt, Net* nettmp)
  : net(nettmp), op(&opt), energybound(-1), solvable(-1)
{
  if (opt.get_randomflag())
    {
      // random Mreq
      n_r = op->get_nr();
      int n_n = net->get_nn();
      startnode = new int[n_r];
      endnode = new int[n_r];
      
      // generate Mreq name
      name = (std::string) "RndMreq" + "_nr" + itoa(n_r);
      if( opt.get_randomflag() == 2 ) // net also random
	name += "_" + nettmp->get_name();
      else // only Mreq random
	{
	  name += "_nn" + itoa(net->get_nn());
	  name += "_seed" + itoa(op->get_seed());
	}
      
      long int redocounter = 0;
      // generate the random Mreq
      for (int i=0;i<n_r;i++)
	{
	  startnode[i]=int(Rnd()*n_n);
	  endnode[i]=int(Rnd()*n_n);
	  if(startnode[i]==endnode[i])
	    {
	      i--;
	      redocounter++;
	      if(redocounter > 1e+06)
		{
		  cerr << "WARNING startnode==endnode for the " ;
		  cerr << redocounter << ":th time, Exiting" << endl;
		  exit(-1);
		}
	    }
	}

      // save it
      write_data(opt.get_requestfile()); // name also written
    }
  else
    {
      // Mreq from file
      //      name += ":RFile_" + opt.get_requestfile();
      // Bullshit:     if( opt.get_randomflag() == 1 )
      //	name += "_" + itoa(opt.get_seed());
      read_data(opt.get_requestfile()); // reads also name
    }

  //set_entropy();
}

/***************************************************************************/
Mreq::~Mreq()
{
  delete[] startnode;
  delete[] endnode;
}

/***************************************************************************/
void Mreq::init(const std::string& tmp_name, Net& tmp_net, Option& opt, int nreq, int *a, int *b)
  // fills up Mreq from default (bare) constructor
  // for a sub_Mreq in a subnet of an RNet
{
  int i;

  name = tmp_name;
  net = &tmp_net;
  op = &opt;
  n_r = nreq;
  startnode = new int [n_r];
  endnode = new int [n_r];
  for( i = 0; i < n_r; i++)
    {
      startnode[i] = a[i];
      endnode[i] = b[i];
    }
//<<<<<<< mreq.cc
  set_entropy();
//=======
//  set_entropy(); // why do for the whole net?
//>>>>>>> 1.25
}

/***************************************************************************/
double _np(int a, int b, int nn, int *used, int *nnei, int **nei, int l, int *freq, int *path) {

  if( DEBUG )
    if( a < 0 || b < 0 || a >= nn || b >= nn )
      cerr << "\7a|b out of bounds" << endl;

  path[l] = a;
  
  if( a == b )
    {
      freq[l]++;

      if( DEBUG )
	{
	  cout << "found l=" << l << endl;
	  cout << "*";
	  for( int i = 1; i < l; i++)
	    cout << " ";
	  cout << 1. << "[" << l << "]";
	  cout << "\tf:";
	  for( int i = 0; i < nn; i++)
	    cout << " " << freq[i];
	  cout << "\tp:";
	  for( int i = 0; i <= l; i++)
	    cout << " " << path[i];
	  cout << endl;
	}
      
      return 1.;
    }

  int aa;
  double x = 0;
  used[a] = 1;
  for( int i = 0; i < nnei[a]; i++)
    if( !used[aa = nei[a][i]] )
      x += _np(aa, b, nn, used, nnei, nei, l + 1, freq, path);
  used[a] = 0;

  if( DEBUG )
    {
      for( int i = 0; i < l; i++)
	cout << " ";
      cout << x << endl;
    }

  return x;
}

/***************************************************************************/
double _np_s(int a, int b, int *used, int *nnei, int **nei) {
  // BS: for large nn, need to log tmp result x
  // ML: Seems to (at least in some cases) to go on forever, a bug?
  //cout << "Entering _np_s " << endl;
  
  /**/
  static int deepcounter = 0;
  static int callcounter = 0;
  static int MLflag = 0;
  int k;

  if( deepcounter == 0 )
    callcounter = 0;
  /**/

  if( a == b )
    {
      /**/
      for( k = 0; k < deepcounter; k++)
	cout << ".";
      cout << ". " << deepcounter << endl;
      /**/

      return 1.;
    }

  if (MLflag)
    return 10000000;

  deepcounter++;
  callcounter++;
  if ( callcounter == 10000000)
    {
      MLflag = 1;
      cerr << "Warning: the entropy seems to take forever to compute,";
      cerr << "so skipping the calc. and just putting it to 10000000 (ML)" << endl;
    }

  int aa, i;
  used[a] = 1;
  double x = 0;
  for( i = 0; i < nnei[a]; i++)
    if( !used[aa = nei[a][i]] )
      x += _np_s(aa, b, used, nnei, nei); 
  used[a] = 0;

  deepcounter--;

  if( deepcounter == 0 )
    cout << " (" << callcounter << " calls)" << endl;
  return x;
}

/***************************************************************************/
double _get_np(int a, int b, int nn, int *nnei, int **nei) {


  int used[nn];
  for( int i = 0; i < nn; i++)
    used[i] = 0;

  if( nn < 3 /* 20 */ )
    {
      int freq[nn];
      int path[nn];
      for( int i = 0; i < nn; i++)
	freq[i] = 0;

      /**/
      cout << a << " -> " << b;
      cout << ", choosing _np():" << endl;
      /**/
      
      return _np(a, b, nn, used, nnei, nei, 0, freq, path);
    }
  else
//<<<<<<< mreq.cc
    {
      /**/
      cout << a << " -> " << b;
      cout << ", choosing _np_s():" << endl;
      /**/
      
      return _np_s(a, b, used, nnei, nei);
    }
//=======
    return 123456789; // Takes to long time (ML)
  //return _np_s(a, b, used, nnei, nei);
//>>>>>>> 1.25
}

/***************************************************************************/
void Mreq::set_entropy() {
  int i, j;

  // prepare info
  int nn = net->get_nn();
  int nnei[nn];
  int *nei[nn];

  for( i = 0; i < nn; i++)
    {
      nnei[i] = net->get_node(i)->get_nl();
      nei[i] = new int [nnei[i]];
      for( j = 0; j < nnei[i]; j++)
	nei[i][j] = net->get_node(i)->node(j);
    }

  //compute

  /**/
  double partentropy, epe;
  int nl = net->get_nl();
  cout << "========= New subnet(" << nn << "," << nl << ";" << n_r << ") =========" << endl;
  /**/

  entropy = 0;
  for( i = 0; i < n_r; i++)
    {
      entropy += partentropy = log(epe = _get_np(startnode[i], endnode[i], nn, nnei, nei));

      /**/
      cout << "part entropy[" << i + 1 << "] = log(" << epe << ") = " << partentropy << endl;
      double q = nl / (double) nn;
      cout << "\t(est: " << (nn-1.) * (log(q) - 1 + 1./q) - .5*log(nn-1.) + .5*log(2 * M_PI) << ")" << endl;
      /**/
    }

  /**/
  cout << "Total entropy = " << entropy << endl << endl;
  /**/

  
  for( i = 0; i < nn; i++)
    delete [] nei[i];
} // of alt. calc_entropy

/***************************************************************************/
void Mreq::read_data(const std::string filename) {

  ifstream input(filename.c_str());
  if (!input) {
    cerr << "Mreq::read_data:1       cannot open input file."
	 << filename << endl;
    exit(-1);
  }

  // skipping until first '@' character
  char ch=' ';
  while (ch!='@')
    input.get(ch);

  input >> name; // NEW

  input >> n_r;
  startnode = new int[n_r];
  endnode = new int[n_r];

  int n_n = net->get_nn();
  // there is no error checking below, should probably be implemented later
  int i;
  for (i = 0; i < n_r; i++)
    {
      int n_end_nodes; // to be compatible with mcast
      input >> n_end_nodes;
      if( n_end_nodes != 1)
	{
	  cerr << "mreq::ERROR can only run mspp, to be compatible with mcast\n";
	  cerr << "the req file must now have the same structure, so if your \n";
	  cerr << "running a old req file add a 1 in the beginning of each line \n";
	  exit(100);
	}
      input >> startnode[i];
      input >> endnode[i];
      if( startnode[i] < 0 ||
	  startnode[i] >= n_n ||
	  endnode[i] < 0 ||
	  endnode[i] >= n_n )
	cerr << "\7Warning: node # out of bounds for Mreq read from file " << filename << endl;
    }
}

/***************************************************************************/
void Mreq::write_data(const std::string filename) const
{
  ofstream output(filename.c_str());
  if (!output) {
    cerr << "Mreq::write_data:1       cannot open output file."
	 << filename << endl;
    exit(-1);
  }
  
  // writing information to output file
  output << "############ Data for the Mreq ###########\n";
  output << "# all characters until the first 'at'\n";
  output << "# character is neglected when this file is\n";
  output << "# read.\n";
  output << "#\n";
  output << "# Mreq name: " << name << endl;
  output << "#\n";
  output << "# format of this file is:\n";
  output << "#\n";
  output << "# mreq_name\n";
  output << "# total_number_of_transmissions\n";
  output << "# 1 from_node \tto_node\n";
  output << "# 1 from_node \tto_node\n";
  output << "# ...........\n";
  output << "# where the 1 is to be compatible with mcast where it is no of end nodes\n";
  output << "# @\n";
  
  output << name << endl; // NEW

  output << n_r << '\n';
  int i;
  for (i = 0; i < n_r; i++) 
    output << "1  " << startnode[i] << '\t' << endnode[i] << '\n'; 
}

/***************************************************************************/
ostream& operator<< (ostream& s,Mreq& n) 
{
  int i;	// loop dummy
  s << n.name << endl;
  s << " acting on " << n.net->get_name() << endl;
  //  s << " with Option " << *n.op << endl;
  s << " Entropy = " << n.entropy << endl;
  s << " having " << n.n_r << " requests:\n";
  for (i = 0; i < n.n_r; i++)
    s << "\t" << n.startnode[i] << '\t' << n.endnode[i] << endl;
  return s;
}
