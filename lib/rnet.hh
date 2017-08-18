// rnet.hh
// for reducing problem into subproblems
// bs, Oct 1996

#ifndef rnet_hh

#define rnet_hh

#include "net.hh"
#include <string>
// #include "net.hh"
// #include "mreq.hh"

class RNet // the tree of subnets (obtained by reducing a net)
{
public:
  RNet(Net &); // construct from an existing net
  ~RNet();

  inline int		get_n_S() const;
  inline Net *		get_M() const;
  inline Net *		get_H();
  inline Net *		get_S(int) const;
  inline std::string		get_name() const;
  int			get_NSpath(int, int, int*&) const;
  int			get_SNN(int, int) const;
  int			get_local(int n, int s) const;
  int			get_global(int n, int s) const;

  friend std::ostream& operator<<(std::ostream&, const RNet &s);

private:
  std::string	name;
  int		n_S;
  Net		*M;	// the mother net
  Net		H;	// the HyperNet, with a node for each subnet and each orig. node
  Net		*S;	// the subnets
};

inline int		RNet::get_n_S() const {return n_S;}
inline Net *		RNet::get_M() const {return M;}
inline Net *		RNet::get_H() {return &H;}
inline Net *		RNet::get_S(int i) const {return &S[i];}
inline std::string		RNet::get_name() const {return name;}

#endif
