// rmreq.hh
// reduced mreq class
// BS, nov -96

#ifndef rmreq_hh
#define rmreq_hh

#include "net.hh"
#include "mreq.hh"
#include "rnet.hh"
#include <string>

class RMreq // the set of mreqs induced by a mreq and a RNet
{
public:
  RMreq(RNet&, Mreq&);
  ~RMreq();
  inline const	std::string&	get_name() const;
  inline	RNet*	get_RN() const;
  inline	Mreq*	get_OMR() const;
  inline	int	get_nleg(int) const;
  inline	int	get_iMRleg(int,int) const;
  inline	int	get_iRleg(int,int) const;
  inline	int	get_nMR() const;
  inline	Mreq*	get_MR(int) const;
  // add set/get_entropy() ?
  friend std::ostream& operator<<(std::ostream&, RMreq&);

private:
  std::string name;
  RNet *RN;	// = arg 1: the reduced net
  Mreq *OMR;	// = arg 2: the original problem
  int	nMR;	// # of parts (Mreq's)
  Mreq *MR;	// [] the corresponding subproblems; one for each subnet
  int *nleg;	// nleg[i] = number of legs obtained upon reducing a single request i of OMR
  int **iMRleg; // iMR[i][j] = index k of the subnet MR[k] associated with leg j of request i of OMR
  int **iRleg;	// iR[i][j] = the sub_request of MR[iMR_leg[i][j]] associated with leg j of request i of OMR
};


inline const std::string& RMreq::get_name() const
{
  return name;
}

inline RNet * RMreq::get_RN() const
{
  return RN;
}

inline Mreq * RMreq::get_OMR() const
{
  return OMR;
}

inline int RMreq::get_nleg(int i) const
{
  return nleg[i];
}

inline int RMreq::get_iMRleg(int i, int j) const
{
  return iMRleg[i][j];
}

inline int RMreq::get_iRleg(int i, int j) const
{
  return iRleg[i][j];
}

inline int RMreq::get_nMR() const
{
  return nMR;
}

inline Mreq * RMreq::get_MR(int i) const
{
  return &MR[i];
}





#endif
