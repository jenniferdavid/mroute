/* -------------------------------------------------------------------- 
   
     filename: Option.hh  	program: mroute  	project: mroute
     written at Theoretical Physics 2, Lund University
  									
     this file contains the Option class definition.
  									
   -------------------------------------------------------------------- */

#ifndef Option_hh
#define Option_hh

#include <iostream>
#include <string>

class Option
{
public:
  Option(int,char*[],char* filename="mroute.hlp");

  inline int		get_fixcapacity() const;
  inline int		get_maxcapacity() const;
  inline std::string		get_networkfile() const;
  inline int		get_nl() const;
  inline int		get_nn() const;
  inline std::string		get_NNfile(int) const;
  inline int		get_NNflag() const;
  inline int		get_nr() const;
  inline int		get_random_sBF() const;
  inline int		get_randomflag() const;
  inline std::string		get_requestfile() const;
  inline int		get_seed() const;
  inline int		get_sx() const;
  inline void		set_NNflag(const int);

  friend std::ostream&	operator<<(std::ostream&, const Option&);
  std::string networkfile;
  std::string NNfile[16];
  std::string requestfile;
  
  
private:
  int fixcapacity;
  int maxcapacity;
  int nl;
  int nn;
  int NNflag;
  int nr;
  int random_sequential_BF;
  int randomflag;
  int seed;
  int sx;
};


inline int Option::get_fixcapacity() const
{
  return fixcapacity;
}


inline int Option::get_maxcapacity() const
{
  return maxcapacity;
}


inline std::string Option::get_networkfile() const
{
  return networkfile;
}


inline int Option::get_nl() const
{
  return nl;
}


inline int Option::get_nn() const
{
  return nn;
}


inline std::string Option::get_NNfile(const int i) const
{
  return NNfile[i];
}


inline int Option::get_NNflag() const
{
  return NNflag;
}


inline int Option::get_nr() const
{
  return nr;
}


inline int Option::get_random_sBF() const
{
  return random_sequential_BF;
}


inline int Option::get_randomflag() const
{
  return randomflag;
}


inline std::string Option::get_requestfile() const
{
  return requestfile;
}


inline int Option::get_seed() const
{
  return seed;
}


inline int Option::get_sx() const
{
  return sx;
}


inline void Option::set_NNflag(const int i)
{
  NNflag=i;
}


#endif
