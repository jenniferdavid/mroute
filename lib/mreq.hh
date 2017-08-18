/************************************/
/*      mreq.hh                      */
/************************************/

#ifndef mreq_hh
#define mreq_hh

#include <iostream>
#include <string>

class Net;
class Option;
//class String;

class Mreq
{
public:
  Mreq(Option&, Net* nettmp);
  Mreq();
  ~Mreq();
  void init(const std::string&, Net&, Option&, int nreq, int *a, int *b);
  inline int		get_endnode(const int) const;	// 
  inline double		get_eb() const;			// Energy bound
  inline double		get_entropy() const;
  inline int		isSolvable() const;
  inline void		set_Solvable(int);
  inline Net		*get_N() const;			// returns the net
  inline int 		get_nr() const;			// number of requests
  inline int		get_startnode(const int) const;
  inline const std::string&	get_name() const;
  inline Option*	get_op() const;
  inline void		set_eb(double);
  friend std::ostream& operator<<(std::ostream& s,Mreq& n);


private:
  std::string	name;		// derived from net and request source (file/seed)
  Net		*net;		// The net the requests works on
  Option	*op;		// The options used
  int	       	n_r;
  int		*startnode;	// one for each req.
  int		*endnode;	// one for each req.
  double	energybound;
  double	entropy;
  int		solvable;
  void read_data(const std::string);
  void set_entropy();
  void write_data(const std::string) const;
};


inline double Mreq::get_eb() const
{
  return energybound;
}


inline const std::string& Mreq::get_name() const
{
  return name;
}

inline int Mreq::get_endnode(int n) const
{
return endnode[n];
}


inline double Mreq::get_entropy() const
{
  return entropy;
}


inline int Mreq::isSolvable() const
{
  return solvable;
}


inline void Mreq::set_Solvable(int i)
{
  solvable=i;
}


inline Net* Mreq::get_N() const
{
return net;
}


inline int Mreq::get_nr() const
{
return n_r;
}


inline int Mreq::get_startnode(int n) const
{
return startnode[n];
}

inline Option* Mreq::get_op() const
{
  return op;
}

inline void Mreq::set_eb(const double d)
{
  energybound=d;
}

#endif
