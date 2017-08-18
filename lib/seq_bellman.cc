/*

  $Id: seq_bellman.cc,v 1.5 1997/08/27 09:18:01 jari Exp $

  */

#include "manswer.hh"
#include "mreq.hh"
#include "mroute.hh"
#include "net.hh"
#include "Option.hh"
#include "random.hh"
#include "Option.hh"
#include <iostream>

using namespace std;

// the Option& op parameter is not used in this function but needed in the
// fuzzy Bellman-Ford algorithm
void seq_bellman(Mreq& req, Manswer& answer, const Option& op)
{
  const Net& net = *req.get_N();
  answer.sname("MA:" + req.get_name() + ":BF");
  
  int n_n = net.get_nn();	// number of nodes
  int n_l = net.get_nl();	// number of links
  double d[n_n];       		// cost of SP i->B
  double inf = 0;      		// Used to init d (inf > sum(weights) )
  double dtest = 0;		// this route cost
  int next[n_n];       		// keeps track of path formed with min D
  int capacity[n_l];		// keeps track of available capacity of links
  int used_request[req.get_nr()]; // keeps track of which request are done

  for (int i=0; i<n_l; i++) {
    inf+=net.get_link(i)->getweight();
    capacity[i]=net.get_link(i)->get_capacity();
  }
  inf *= 2;

  for (int i=0; i<req.get_nr(); i++)
    used_request[i]=0;
  
  Path* path=new Path[answer.gnpath()];	// the result is going to be stored here
  for (int n=0; n<req.get_nr(); n++) {

    int u_request=-1;
    if (op.get_random_sBF()) {
      for (int i=0; i<Rnd()*(req.get_nr()-n); i++) 
	while (used_request[++u_request]);
      used_request[u_request]=1;
    }
    else
      u_request=n;
    int endnode = req.get_endnode(u_request);
    int startnode = req.get_startnode(u_request);
    dtest = 0;
    if (endnode>n_n ||  startnode>n_n)
      cerr << "bellman():1 Neither end or start node could be found "<<n<<endl;

    // initialize d to 'infinity', to make sure that it is too large, and
    // initialize next to something stupid
    for (int i=0; i<n_n; i++) {
      next[i]=-1;
      d[i] = inf;
    }
    d[endnode]=0;// this avoids outgoing connections from end node
    
    // loop until no updates
    int j; int flag=TRUE;
    while (flag) {
      flag=FALSE;
      for (int i=0; i<n_n; i++)		// from i via j to end node (1 leg more)
	for (j=0; j<net.get_node(i)->get_nl(); j++) {
	  dtest=net.get_node(i)->get_link(j)->getweight()
	    + d[net.get_node(i)->node(j)];
	  if ((dtest + 1.e-10 < d[i]) &&
	      (capacity[net.get_node(i)->get_link(j)->getindex()])) {
	    d[i] = dtest;
	    next[i] = net.get_node(i)->node(j);  // set to global index!
	    flag = TRUE;
	  }
	}
    } // end while

    // checking if good solution (i.e starting from startnode and ending at
    // endnode)
    int i=startnode;
    while (i!=endnode && i!=-1)
      i=next[i];
    if (i==endnode) {
      // store the result
      path[n].push_front(startnode);
      i = startnode;
      while (next[i]>-1) {
	for (int k=0; k<net.get_node(i)->get_nl(); k++)
	  if (next[i]==net.get_node(i)->node(k)) {
	    capacity[net.get_node(i)->get_link(k)->getindex()]--;
	    break;
	  }
	path[n].push_front(next[i]);	
	i = next[i];
      }
    }
  } // end loop over all req

  answer.spath(path);
}

