/****************************************/
/*      node.hh                        */
/*      node class def.		*/
/*      Martin & Jari 2/6-96	        */
/****************************************/

#ifndef node_hh
#define node_hh

#include <iostream>

class Link;

class Node
{
public:
  Node();
  ~Node();

  // add a link at the end of linklist: very inefficient!
  void 		addlink(int, Link&);
  inline Link*	get_link(int) const;	// returns pointer to a link
  inline int 	getindex() const;	// return index of this node
  inline int 	getlabel() const;	// return label of this node
  inline int	get_nl() const;
  void		init(int, int, int);		// set index, label and # of links
  inline int	node(int) const;	// return connected node index
  void		set_link(int, Link*, int); // place a link in a specified place
  						 // in the list; better than addlink()
  friend std::ostream& operator<<(std::ostream& s, const Node& n);

private:
  int		label;		// to keep track of original for node in sub-nets, etc.
  int		index;		// global index of this node
  int		n_l;		// # of nodes (links) connected to this node
  Link		** l;		// list of links (connections to other nodes)
  int		* node_index;	// neighbour nodes, returns global index
};

inline int	Node::getindex()	const	{ return index; }
inline int	Node::getlabel()	const	{ return label; }
inline int	Node::node(int i)	const	{ return node_index[i]; }
inline Link*	Node::get_link(int i)	const	{ return l[i]; }
inline int	Node::get_nl()		const	{ return n_l; }

#endif




