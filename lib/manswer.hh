/************************/
/*	manswer.hh	*/
/************************/

#ifndef manswer_hh
#define manswer_hh

#include <iostream>
#include <string>

#include "mroute.hh"
#include "mreq.hh"

class Option;
class RMreq;

class Manswer
{
public: 
  Manswer(const Option&,RMreq&,const Alg);
  Manswer(const Option&,Mreq&,const Alg);
  ~Manswer();

  void		analyse();
  inline int	gnpath() const;
  inline Path*	gppath() const;
  inline int	isLegal() const;
  inline int	isSolved() const;
  inline void	setSolved(const int);
  inline void	sname(const std::string);
  inline std::string	gname() const;
  inline void	spath(const Path * const);
  inline void	spathenergy(int, double);
  inline double	genergy() const;
  friend std::ostream& operator<<(std::ostream&, const Manswer&);

private:
  std::string	name;
  double*	pathenergy;
  double	totalenergy;
  int		legal;
  int		solved;
  Mreq*		request;
  Net*		net;
  Path*		path;
};


inline int Manswer::gnpath() const
{
  return request->get_nr();
}


inline Path* Manswer::gppath() const
{
  return path;
}


inline int Manswer::isLegal() const
{
  return legal;
}


inline int Manswer::isSolved() const
{
  return solved;
}


inline void Manswer::setSolved(const int i)
{
  solved=i;
}


inline void Manswer::sname(std::string s)
{
  name=s;
}

inline std::string Manswer::gname() const
{
  return name;
}

inline void Manswer::spath(const Path * const p)
{
  path = (Path *) p;
}


inline double Manswer::genergy() const
{
  return totalenergy;
}

inline void Manswer::spathenergy(int i, double e)
{
  pathenergy[i]=e;
}


#endif
