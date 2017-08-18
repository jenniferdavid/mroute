/****************************************/
/*      link.cc		       	*/
/*      link class def		*/
/*      Martin & Jari 2/7-96	        */
/****************************************/

#include "mroute.hh"
#include "link.hh"
#include <iostream>
using namespace std;

/***************************************************************************/
Link::Link()
/***************************************************************************/
{ 
  index=-1;
  w = -1;
  capacity = -1;
  node[0] = 0;	// has to be set by using "init (not setnodes which is never defined)" 
  node[1] = 0;
}

/***************************************************************************/
Link::~Link()
/***************************************************************************/
{
if( DEBUG ) cerr << "Entering and leaving Link::~Link() " << index << endl;
}




/***************************************************************************/
void Link::init( const Node & node1,
		 const Node & node2,
		 int linkno,
		 double tmpw,
		 int tmpmaxload )
/***************************************************************************/
{
  //
  // cerr << endl;
  // cerr << "Link::init("
  // << node1.getindex() << ", "
  // << node2.getindex() << ", "
  // << linkno << ", "
  // << tmpw << ", "
  // << tmpmaxload << ")" << endl;
  //

  index = linkno;
  capacity = tmpmaxload;
  node[0] = (Node *) &node1;
  node[1] = (Node *) &node2;
  w = tmpw;
}




/***************************************************************************/
ostream& operator<< (ostream& s, const Link& l)
/***************************************************************************/
{
  //  s.width(2);
  s.setf(ios::fixed);
  s.precision(3);

  s << "Link " << l.index;
  s << "\t(" << l.node[0]->getindex() << " -> " << l.node[1]->getindex() << ")"; 
  s << "\tWeight " << l.w;
  s << "\tCapacity " << l.capacity;
  return s;
}


