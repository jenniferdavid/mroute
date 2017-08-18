/************************************************/
/*            bellman.cc                        */
/************************************************/

#include "manswer.hh"
#include "mreq.hh"
#include "mroute.hh"
#include "net.hh"
using namespace std;

// the Option& op parameter is not used in this function but needed in the
// fuzzy Bellman-Ford algorithm
void bellman(Mreq& req, Manswer& answer, const Option& op)
{
  const Net& net = *req.get_N();
  answer.sname("MA:" + req.get_name() + ":BF");
  
  int n_n = net.get_nn();	// number of nodes
  int n_l = net.get_nl();	// number of links
  double d[n_n];       		// cost of SP i->B
  double inf = 0;      		// Used to init d (inf > sum(weights) )
  double dtest = 0;		// this route cost
  int next[n_n];       		// keeps track of path formed with min D

  int i;
  for (i = 0; i < n_l; i++) 
    inf+=net.get_link(i)->getweight();
  inf *= 2;
  
  Path* path=new Path[answer.gnpath()];	// the result is going to be stored here
  for (int n = 0; n < req.get_nr(); n++) { 
    int endnode = req.get_endnode(n);
    int startnode = req.get_startnode(n);
    dtest = 0;
    if ( endnode > n_n ||  startnode > n_n )
      cerr << "bellman():1 Either end or start node could be found"<<n<<endl;
    
    // initialize d to 'infinity', to make sure that it is too large
    for (i = 0; i < n_n; i++)
      d[i] = inf;
    d[endnode]=0;// this avoids outgoing connections from end node
    for (i=0;i<n_n;i++)
      next[i]=-1;	      // initialize next to something stupid
    
    // loop until no updates
    int j; int flag=TRUE;
    while (flag) {
      flag=FALSE;
      for (i=0;i<n_n;i++)		// from i via j to end node (1 leg more)
	for (j = 0; j < net.get_node(i)->get_nl(); j++) {
	  dtest=net.get_node(i)->get_link(j)->getweight()
	    + d[net.get_node(i)->node(j)];
	  if (dtest + 1.e-10 < d[i]) {
	    d[i] = dtest;
	    next[i] = net.get_node(i)->node(j);  // set to global index!
	    flag = TRUE;
	  }
	}
    } // end while
      
    // store the result
    path[n].push_front(startnode);
    i = startnode;
    while (next[i]>-1) {
      path[n].push_front(next[i]);	
      i = next[i];
    }
  } // end loop over all req
  answer.spath(path);
}

