/****************************************/
/*      net.hh                          */
/*      net class def                   */
/****************************************/

#ifndef net_hh
#define net_hh

#include <iostream>
#include <string>

#include "link.hh"
#include "mroute.hh"
#include "node.hh"

class Option;

class Net
{
public:
  Net(Option&);
  Net();
  ~Net();

  void init(std::string id, int nodes, int links, int *Nlbl, int **iNL, double *, int *);	// initiates a given net
  int			checklinkloads(Path*,const int) const;
  inline Link*		get_link(int)	const;
  inline int		get_nl()	const;
  inline Node*		get_node(int)	const;
  inline int		get_nn()	const;
  inline const std::string &	get_name()	const;
  void			l2n(Path& linkpath, Path& nodepath) const;
  void			n2l(Path& nodepath, Path& linkpath) const;
  double		pathcost(Path& path) const;
  int			isNormal()	const;
  int			isConnected(const int *const = 0) const;
  // int * = n_l vector of 0/1, flagging turned on links; NULL => all on
  friend std::ostream& operator<<(std::ostream&, const Net&);

private:
  int		nn;		// # of nodes
  int		nl;		// # of links
  Node		* n;		// the nodes
  Link		* l;		// the links
  int		normal;
  std::string	name;
  int			amNormal(); // not const, since may change data member normal
  void			check_net() const;
  void			random_data(Option&);
  void			read_data(const std::string); // read node-net, links and weights
  // write a data file of used data (readable by above)
  void			write_data(const std::string) const;
};


inline Link* Net::get_link(int i) const
{ return &l[i]; } 


inline int Net::get_nl() const
{ return nl; }


inline Node* Net::get_node(int i) const
{ return &n[i]; }


inline int Net::get_nn() const
{ return nn; }

inline const std::string &	Net::get_name()	const
{ return (const std::string&) name; }

#endif
