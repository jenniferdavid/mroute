/****************************************/
/*      net.cc				*/
/*      net class def	                */
/*      Martin & Jari 2/7-96	        */
/****************************************/

#define DUMP 0 // temporary - should be in mroute.hh

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include "itoa.hh"
#include "mroute.hh"
#include "net.hh"
#include "Option.hh"
#include "random.hh"
#include "stdio.h"

using namespace std;

/***************************************************************************/
Net::Net(Option& op) : nn(op.get_nn()), nl(2*op.get_nl()), normal(-1)
{
  if (op.get_randomflag() == 2) {
    name = (std::string) "RndNet_nn" + itoa(nn) + "_nl" + itoa(nl) + "_seed" + itoa(op.get_seed());
    random_data(op);
    write_data(op.get_networkfile()); // writes also name
  } 
  else
    {
      // name = "NFile_" + op.get_networkfile();
      read_data(op.get_networkfile()); // reads also name
    }

  if( !amNormal() )
    cerr << "\7Warning, net abnormal" << endl;

  if (DEBUG)
    check_net();
}

/***************************************************************************/
Net::Net()
{
  nn = 0;
  nl = 0;
  name = "<null>";
  n = NULL;
  l = NULL;
  normal = -1;
}

/***************************************************************************/
Net::~Net()
{
  if( DEBUG ) cerr << "Entering Net::~Net() " << name << endl;
  if( l ) delete[] l;
  if( n ) delete[] n;
  if( DEBUG ) cerr << "Leaving Net::~Net() " << name << endl;
}

/***************************************************************************/
int Net::checklinkloads(Path* path,const int npath) const
{
  double load[get_nl()];
  for (int j=0;j<get_nl();j++)		// setting 'zero' level
    load[j]= -get_link(j)->get_capacity();
  
  // checking loads
  for (int i=0;i<npath;i++) {
    Path::iterator prev=path[i].begin();
    Path::iterator jj=prev;
    for (jj++;jj!=path[i].end();jj++) {
      int index;
      // the order to convert into links might seem fishy but remember that
      // the paths are inverted, startnode is in the end of path
      for (int k=0;k<get_node(*jj)->get_nl();k++)
	if (get_node(*jj)->node(k)== *prev) {
	  index=get_node(*jj)->get_link(k)->getindex();
	  break;
	}
      if ((++load[index])>0)
	return 0;
      prev=jj;
    }
  }

  return 1;
} // of checklinkloads

/***************************************************************************/
void Net::l2n(Path& linkpath,Path& nodepath) const
{
  Path::iterator i=linkpath.begin();
  // the order to convert into nodes might seem fishy but remember that
  // the paths are inverted, startnode is in the end of path
  nodepath.push_front(get_link(*i)->get_node(0));
  nodepath.push_front(get_link(*i)->get_node(1));
  for (i++;i!=linkpath.end();i++)
    nodepath.push_back(get_link(*i)->get_node(0));
} // of Net::l2n

/***************************************************************************/
void Net::n2l(Path& nodepath,Path& linkpath) const
{
  Path::iterator i=nodepath.begin();
  Path::iterator j=i;
  for (j++;j!=nodepath.end();j++) {
    // the order to convert into links might seem fishy but remember that
    // the paths are inverted, startnode is in the end of path
    for (int k=0;k<get_node(*j)->get_nl();k++)
      if (get_node(*j)->node(k)==*i) {
	linkpath.push_back(get_node(*j)->get_link(k)->getindex());
	break;
      }
    i=j;
  }
} // of Net::n2l

/***************************************************************************/
double Net::pathcost(Path& path) const
{
  double cost=0;
  for (Path::iterator k=path.begin();k!=path.end();k++)
    cost+=get_link(*k)->getweight();

  return cost;
} // of Net::pathcost(Path&)

/***************************************************************************/
void Net::random_data(Option& op)
{
  if (nl<2*(nn-1)) {
    cerr << "Net::Net():1     too few links" << endl;
    cerr << "    links: " << nl << " and nodes: " << nn << endl;
    exit(-1);
  }   
  if (nl<0) {
    cerr << "Net::Net():1     negative # of links: " << nl << endl;
    exit(-1);
  }
  if (nn<0) {
    cerr << "Net::Net():1     negative # of nodes: " << nn << endl;
    exit(-1);
  }
  
  if( nl ) l=new Link[nl];
  if( nn )
    n=new Node[nn];
  else
    {
      name = "<null>";
      return;
    }
  
  int maxload;
  double weight;
    
  // starting by connecting all the nodes --> tree network
  int link=0;
  for (int k=0;k<nn;k++) 
    n[k].init(k,k,0);	// tell all the node their indices; no connections yet
  for (int k=1;k<nn;k++) {
    int node1=int(Rnd()*k);// Take this random connected node (0 ... k-1), and
    // connect to this random non-connected node (k ... nn-1)
    int node2=int(Rnd()*(nn-k));
    // NOTE node2=number of nonconnected nodes to skip (beginning from node[0]
    // before accepting a node, node1 is sometimes one of the nonconnected
    // nodes and should not be counted as nonconnected, fixing this below
    
    // looking for the connected node
    int match = -1; int counter = -1;
    while ( match!=node1 )
      // k==1 needed for the first node, since all nonconnected then
      match += ( n[++counter].get_nl() || k == 1 );
    node1 = counter;
    
    // looking for the nonconnected node
    match=-1; counter=-1;
    while (match!=node2) 
      match+=((!n[++counter].get_nl()) && (node1!=counter));
    node2=counter;
    
    n[node1].addlink(node2,l[link]);
    weight=Rnd();
    if (op.get_fixcapacity())
      maxload=op.get_maxcapacity();
    else
      maxload=(int) (Rnd()*op.get_maxcapacity())+1;
    l[link].init(n[node1],n[node2],link,weight,maxload);
    link++;
    
    n[node2].addlink(node1,l[link]);
    weight=Rnd();
    if (op.get_fixcapacity())
      maxload=op.get_maxcapacity();
    else
      maxload=(int) (Rnd()*op.get_maxcapacity())+1;
    l[link].init(n[node2],n[node1],link,weight,maxload);
    link++;
  }
  
  // inserting the remaining links
  int loop=0;
  int flag=FALSE;
  for (int k=link;k<nl;k+=2) {
    int node1;
    int node2=node1=int(Rnd()*nn);
    while (node1==node2)
      node2=int(Rnd()*nn);
    int j;
    for (j=0;j<n[node1].get_nl();j++)
      if (n[node1].node(j)==node2)
	flag=TRUE;
    loop+=flag;    // loop is a counter of # tries before link could be placed
    if (loop>10000) {
      cerr << "Net::Net():2     caught in an infinite loop, ";
      cerr << "cannot place all links" << endl;
      exit(-1);
    }
    if (!flag) {
      n[node1].addlink(node2,l[link]);
      weight=Rnd();
      if (op.get_fixcapacity())
	maxload=op.get_maxcapacity();
      else
	maxload=(int) (Rnd()*op.get_maxcapacity())+1;
      l[link].init(n[node1],n[node2],link,weight,maxload);
      link++;
      
      n[node2].addlink(node1,l[link]);
      weight=Rnd();
      if (op.get_fixcapacity())
	maxload=op.get_maxcapacity();
      else
	maxload=(int) (Rnd()*op.get_maxcapacity())+1;
      l[link].init(n[node2],n[node1],link,weight,maxload);
      link++;
      
      loop=0;		// reset loop counter
    }
    else {
      k-=2;		// redo these links
      flag=FALSE;
    }
  }
} // of random_data()

/***************************************************************************/
void Net::init(std::string id, int nodes, int links, int *Nlbl, int **iNL, double *Lw, int *Lcap)
	// initiates a given net
{
  int i, j, m;

  name = id;
  nn = nodes;
  nl = links;
  n = new Node[nn];
  l = new Link[nl];

  for( i = 0; i < nn; i++)
    {
      m = 0;
      for( j = 0; j < nl; j++)
	if( iNL[j][0] == i )
	  m++;
      n[i].init(i, Nlbl[i], m);
    }

  if( Lw )
    {
      if( Lcap )
	for( i = 0; i < nl; i++)
	  l[i].init(n[iNL[i][0]], n[iNL[i][1]], i, Lw[i], Lcap[i]);
      else
	for( i = 0; i < nl; i++)
	  l[i].init(n[iNL[i][0]], n[iNL[i][1]], i, Lw[i], 0);
    }
  else
    {
      if( Lcap )
	for( i = 0; i < nl; i++)
	  l[i].init(n[iNL[i][0]], n[iNL[i][1]], i, 0., Lcap[i]);
      else
	for( i = 0; i < nl; i++)
	  l[i].init(n[iNL[i][0]], n[iNL[i][1]], i, 0., 0);
    }

  for( i = 0; i < nn; i++)
    {
      m = 0;
      for( j = 0; j < nl; j++)
	if( iNL[j][0] == i )
	  {
	    n[i].set_link(m,&l[j],iNL[j][1]);
	    m++;
	  }
      if( m != n[i].get_nl() )
	cerr << "Wrong node_nl" << endl;
    }

  if( !amNormal() )
   cerr << "Net " << name << " not normal" << endl;

  if (DUMP)
    write_data( "subnet." + name); // writes also name

}

/***************************************************************************/
void Net::read_data(const std::string filename) {

  ifstream input(filename.c_str());
  if (!input) {
    cerr << "Net::read_data:1       cannot open input file: "
	 << filename << endl;
    exit(-1);
  }

  normal = -1;

  // skipping until first '@' character
  char ch=' ';
  while (ch!='@')
    input.get(ch);

  input >> name; // NEW

  input >> nn;
  input >> nl;
  //  nl*=2;
  n = new Node[nn];
  l = new Link[nl];

  int maxload;
  double weight;
  int link=0;
  int nlinks;
  int node;
  for (int k=0;k<nn;k++) {
    input >> nlinks;
    n[k].init(k,k, nlinks);
    for (int i=0;i<nlinks;i++) {
      input >> node;
      input >> weight;
      input >> maxload;
      n[k].set_link(i,&(l[link]),node);
      l[link].init(n[k],n[node],link,weight,maxload);
      link++;
    }
  }
}

/***************************************************************************/
int Net::isNormal() const // checks if the links come in single pairs,
{		          // with no self-couplings and no multi-links,
			  // and that the net is connected
  if( normal != -1 )
    return normal; // must be sure that net has not been changed:

  cerr << "Not checked: Have to say abnormal" << endl;
  return 0;
}

/***************************************************************************/
int Net::amNormal() // checks if the links come in single pairs,
{		    // no self-couplings, and that it is connected

  if( normal != -1 )
    return normal; // must be sure that net has not been changed:
		   // every non-const member function must reset normal to -1
  if( nn == 0 && nl == 0 )
    {
      normal = 1;
      return normal;
    }

  int i, j, k;
  int conn[nn][nn];
  int err0 = 0;
  int err1 = 0;
  int err2 = 0;
  int err3 = 0;
  int err4 = 0;

  for( i = 0; i < nn; i++)
    for( j = 0; j < nn; j++)
      conn[i][j] = 0;
  
  for( i = 0; i < nn; i++)
    for( k = 0; k < n[i].get_nl(); k++)
      {
	j = n[i].node(k);
	conn[i][j]++;
      }
  
  //  delete [][] conn;
  for( i = 0; i < nn; i++)
    {
      if( conn[i][i] )
	err0++;
      if( conn[i][i] > 1 )
	++err3;
      for( j = 0; j < i; j++)
	{
	  if( conn[i][j] != conn[j][i] )
	    ++err1;
	  if( conn[i][j] && !conn[j][i] || conn[j][i] && !conn[i][j] )
	    ++err2;
	  if( conn[i][j] > 1 )
	    ++err3;
	  if( conn[j][i] > 1 )
	    ++err3;
	}
    }

  if( err0 )
    cerr << "Net with " << err0 << " links connecting node with itself\7" << endl;

  if( err2 )
    {
      cerr << "Net strongly asymmetric: " << err2 << " cases of link i->j but not j->i\7" << endl;
      cerr << "\tIn addition " << err1 - err2 << " cases of different number of links i->j/j->i\7" << endl;
    }
  else if( err1 )
    cerr << "Net weakly asymmetric: " << err1 << " cases of different number of links i->j/j->i\7" << endl;
  if( err3 )
    cerr << "Net with multi-links: " << err3 << " cases of multiple links i->j\7" << endl;

  i = isConnected();

  if( i != 2 )
    err4 = 1;

  if( err0 || err1 || err2 || err3 || err4 )
    normal = 0;
  else
    normal = 1;

  return normal;
}

/***************************************************************************/
int Net::isConnected(const int *const flag/* = 0*/) const
// flag[n_l]: flag[il] = 0 if link il turned off
// returns 0 (disconn), 1 (semi-conn), 2 (conn)
{
  int i, j, i1, i2;
  int nn1, nn2, nn3, nf[nn], used1[nn], used2[nn];
  Node *v;
  Link *l;
  int allflag = 0, found;

  if( nn == 0 )
    return 1;

  if(!flag) allflag = 1; //full net
  else if( !isNormal() )
    cerr << "\7Warning: Possibly abnormal net" << endl;

  // i) Check if all nodes can be reached from node 0
  nn1 = 0;
  for( i = 0; i < nn; i++)
    used1[i] = 0;
  used1[nf[nn1++] = 0]++;
  for( i = 0; i < nn1; i++)
    {
      i1 = nf[i];
      v = &n[i1];
      for( j = 0; j < v->get_nl(); j++)
	{
	  l = v->get_link(j);
	  if( allflag || flag[l->getindex()] )
	    { 
	      i2 = l->get_node(1); // other end
	      if( !used1[i2] )
		used1[nf[nn1++] = i2]++;
	    }
	}
    }
  
  // i) Check if all nodes can reach node 0
  nn2 = 0;
  for( i = 0; i < nn; i++)
    used2[i] = 0;
  used2[nf[nn2++] = 0]++;
  found = 1;
  while(found)
    {
      found = 0;
      for( i1 = 0; i1 < nn; i1++)
	if( !used2[i1] )
	  {
	    v = &n[i1];
	    for( j = 0; j < v->get_nl(); j++)
	      {
		l = v->get_link(j);
		if( allflag || flag[l->getindex()] )
		  {
		    i2 = l->get_node(1); // other end
		    if( used2[i2] )
		      {
			used2[nf[nn2++] = i1]++;
			found++;
			break;
		      }
		  }
	      }
	  }
    }

  if( nn2 == nn && nn1 == nn )
    return 2;	// fully connected
  else
    {
      nn3 = 0;
      for( i = 0; i < nn; i++)
	if( used1[i] || used2[i] )
	  nn3++;
      if( nn3 == nn )
	return 1;	// semi-connected
      else
	return 0;	// disconnected
    }
}

/***************************************************************************/
void Net::write_data(const std::string filename) const
{
  ofstream output(filename.c_str());
  if (!output) {
    cerr << "Net::write_data:1       cannot open output file: " 
	 << filename << endl;
    exit(-1);
  }
  
  output << "########## Data for the Net ############\n";
  output << "# all characters until the first 'at'\n";
  output << "# character is neglected when this file\n";
  output << "#  is read.\n";
  output << "#\n";
  output << "# Net name: " << name << endl;
  output << "#\n";
  output << "# format of this file is:\n";
  output << "#\n";
  output << "# net_name\n";
  output << "# total_number_of_nodes \ttotal_number_of_links\n";
  output << "# number_of_links_connected_to_this_node\n";
  output << "# connects_to_this_node \tweight max-load\n";
  output << "#\n";
  output << "# @\n";
  
  output << name << endl; // NEW

  output << nn << '\t' << nl << '\n';
  int k;
  for (k=0;k<nn;k++) {
    output << n[k].get_nl() << '\n'; 
    for (int i=0;i<n[k].get_nl();i++) {
      output << '\t' << n[k].node(i);
      output << '\t' << n[k].get_link(i)->getweight();
      output << '\t' << n[k].get_link(i)->get_capacity();
      output << '\n';
    }
  }
}

/***************************************************************************/
void Net::check_net() const
  {
    cout << "----- CHECKING THE NET ----\n";
    cout <<"current status\n";      
    cout << *this;
    cout << " \nThe nodes:\n";
    int i;
    for (i = 0; i < nn; i++)
      cout << n[i];
    cout << endl;
    if( isNormal() )
      cout << "Normal net" << endl;
    else
      cout << "Abnormal net" << endl;
    cout << "----- DONE WITH NET CHECKING ----\n";
  }  

/***************************************************************************/
ostream& operator<< (ostream& s, const Net& n)
{
  s << n.name << endl;

  s << "  (" << n.nn << " nodes, " << n.nl << " links)" << endl;

  s << "Nodes:";
  for (int i=0;i<n.nn;i++)
    s << " " << i << " [" << n.n[i].getlabel() << "]";
  s << endl;

  if (DEBUG)
    {
      for (int i=0;i<n.nn;i++)
	s << n.n[i];
    }

  s << "Links:\n";
  for (int i=0;i<n.nl;i++)
    s << "  " << n.l[i] << endl;

  return s;
}
