/****************************************/
/*      node.cc		       	*/
/*      node class def		*/
/*      Martin & Jari 2/7-96	        */
/****************************************/

#include "node.hh"
#include "mroute.hh"
#include "link.hh"
#include <iostream>

using namespace std;
/***************************************************************************/
Node::Node()
/***************************************************************************/
{ 
  index = -1;
  label = -1;
  n_l=0;
  l=NULL;
  node_index=NULL;
}

/***************************************************************************/
Node::~Node()
/***************************************************************************/
{
if( DEBUG ) cerr << "Entering Node::~Node() " << index << ":" << label << endl;
  if( l ) delete [] l;
  if ( node_index ) delete [] node_index;
if( DEBUG ) cerr << "Leaving Node::~Node() " << index << ":" << label << endl;
}

/***************************************************************************/
void Node::addlink( int		node,
		    Link	&link )
/***************************************************************************/
{
  // this function would gain in efficiency if a linked list were introduced

  int* tmp_nodeindex = new int[n_l+1];
  Link** tmp_l = new Link* [n_l+1];
  
  int i;
  for (i = 0; i < n_l; i++)
    {
      tmp_nodeindex[i] = node_index[i];
      tmp_l[i] = l[i];
    }
  tmp_nodeindex[n_l] = node;
  tmp_l[n_l] = &link;
  n_l++;
  delete[] l;
  delete[] node_index;
  l = tmp_l;
  node_index = tmp_nodeindex;
}

/***************************************************************************/
void Node::init( int	tmp_index,
		 int	tmp_label,
		 int	links )
/***************************************************************************/
{
  //
  // cerr << endl;
  // cerr << "Node::init("
  // << tmp_index << ", "
  // << tmp_label << ", "
  // << links << ")" << endl;
  //


  n_l = links;
  l = new Link* [links];
  node_index = new int[links];
  index = tmp_index;
  label = tmp_label;
}

/***************************************************************************/
void Node::set_link( int	i,
		     Link	*link,
		     int	nodeno )
/***************************************************************************/
{
  l[i]=link;
  node_index[i]=nodeno;
}

/***************************************************************************/
ostream& operator<< (ostream& s, const Node& n)
/***************************************************************************/
{
  int i;
  
  cout << "Node " << n.index << " [" << n.label << "] (" << n.n_l << " links)\n";
  if (n.l)
    for (i=0;i<n.n_l;i++)
      cout << "\t(-> " << n.node_index[i] << "): " << *n.l[i] << endl;
  else
    cout << "No links defined\n";
  return s;
}

