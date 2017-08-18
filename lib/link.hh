/****************************************/
/*      link.hh                        */
/*      generell link class def       */
/*      Martin & Jari  2/6-96	        */
/****************************************/

#ifndef link_hh
#define link_hh

#include <iostream>
#include "node.hh"

class Link
{
public:
		Link();
		~Link();

  inline int	get_capacity()	const;
  inline int	get_node(int x)	const;
  inline int	getindex()	const;
  inline double	getweight()	const;

  void		init(const Node &, const Node &, int, double, int);
  //  void		setnodes(const Node & node1, const Node & node2);

  friend std::ostream& operator<<(std::ostream&, const Link&);

private:
  int		index;  	// remove this? No needed in fuzzy for 
				// putting the load at the correct link
  Node		* node[2];	// the nodes that this links connects
  double	w;		// cost coefficient
  int		capacity;
};

inline int	Link::get_node(int x)	const	{ return node[x]->getindex(); }
inline int	Link::get_capacity()	const	{ return capacity; }
inline int	Link::getindex()	const	{ return index; }
inline double	Link::getweight()	const	{ return w; }

#endif
