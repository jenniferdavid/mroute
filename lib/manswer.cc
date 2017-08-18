/*

  filename: manswer.cc

  */


#include "manswer.hh"
#include "mreq.hh"
#include "mroute.hh"
#include "net.hh"
#include "Option.hh"
#include "rmreq.hh"

using namespace std;

// declaration of output operators, not done in STL
ostream& operator<< (ostream&,const Path&);

// ***********************************************************************
Manswer::Manswer(const Option& op,RMreq& rr,const Alg alg)
  : name(rr.get_OMR()->get_name()+":R"), totalenergy(0), solved(1), path(NULL)
{
  // shouldn't these things be added in the algorithms?
  if (alg==bellman)
    name+="BF";
  //else if (alg==fuzzybellman)
  // name+="FZ";
  else if (alg==solveexact)
    name+="SX";
  else if (alg==NN_alg)
    name+="NN";
  else
    name+="??";

  request=rr.get_OMR();
  net=request->get_N();
  pathenergy=new double[gnpath()];

  // creating the reduced answers
  Manswer* ra[rr.get_nMR()];
  for (int i=0;i<rr.get_nMR();i++)
    ra[i]=NULL;

  int unsolvableflag = 0;
  int refuseflag = 0;
  for (int i=0;i<rr.get_nMR();i++)
    {
      ra[i]=new Manswer(op,*rr.get_MR(i),alg);
      setSolved(isSolved()*ra[i]->isSolved());
      if (alg==solveexact && !rr.get_MR(i)->isSolvable())
	{
	  rr.get_OMR()->set_Solvable(FALSE);
	  unsolvableflag = 1;
	  break;
	}
      if (alg==solveexact && !ra[i]->isLegal())
	{
	  refuseflag = 1;
	  break;
	}
    }

  // merging all answers into *this
  if (!unsolvableflag && !refuseflag)
    //rr.get_OMR()->isSolvable() || (alg!=solveexact)) //see above (break;)
    {
      path=new Path[gnpath()];
      for (int h=0;h<gnpath();h++)
	{
	  int pathisthere=TRUE;
	  for (int i=rr.get_nleg(h)-1;i>-1;i--)
	    {
	      int aaa=rr.get_iMRleg(h,i);
              int bbb=rr.get_iRleg(h,i);
              Path p;
              if(ra[aaa])
                p = ra[aaa]->gppath()[bbb];
	      if (!p.size())
		pathisthere=FALSE;
	      for (Path::iterator j=p.begin();j!=p.end();j++)
		*j=rr.get_RN()->get_global(*j,aaa);
	      if (path[h].size() && p.size())
		p.pop_front();
	      path[h].splice(path[h].end(),p);
	    }
	  if (!pathisthere)
	    path[h].erase(path[h].begin(),path[h].end());
	}
    }
  
  // unnecessary?
  analyse();

  // there are things to delete !!
  for (int i=0;i<rr.get_nMR();i++)
    if( ra[i] )
      delete ra[i];

} // of Manswer(Option&,RMreq&,Alg)

// ***********************************************************************
Manswer::Manswer(const Option& op,Mreq& req,const Alg alg)
  : name("<?????>"), totalenergy(0), solved(0), request(&req), path(NULL)
{
  // name will be set in the Alg's
  net=req.get_N();
  pathenergy=new double[gnpath()];
  for (int i=0;i<gnpath();i++)
    pathenergy[i]=0;
  
  if( !req.get_nr() ) // what if ... ? /BS
    {
      //cerr << "\7MA: Trivial problem! -- nr = 0" << endl;
      legal = TRUE;
      path = new Path[0];
      pathenergy = new double[0];
      setSolved(TRUE);
      return;
    }

  //  if( net->get_nn() < 3 ) // what if ... ? /BS
    //cerr << "\7MA: Trivial problem? -- nn = " << net->get_nn() << endl;
  
  if( !net->get_nl() ) // what if ... ? /BS
    cerr << "\7MA: Weird problem -- nl = 0" << endl;

  alg(*request,*this,op);
  
  analyse();

  // if legal solution obtained, try to lower the energybound for the problem
  if (legal && (req.get_eb()<0 || req.get_eb()>totalenergy))
    req.set_eb(totalenergy);
} // of Manswer(Option&,Mreq&,Alg)

// ***********************************************************************
Manswer::~Manswer()
{
  if( path )
    delete [] path;
  delete [] pathenergy;
}

// ***********************************************************************
void Manswer::analyse()
{
  int allpathsexist=TRUE;

  if (path) {
    totalenergy=0;
    for (int i=0;i<gnpath();i++) {
      if (!path[i].size())
	allpathsexist=FALSE;
      else {
	Path linkpath;
	net->n2l(path[i],linkpath);
	pathenergy[i]=net->pathcost(linkpath);
	totalenergy+=pathenergy[i];
      }
    }
  }
  else
    allpathsexist=FALSE;

  if ((allpathsexist) && (net->checklinkloads(path,gnpath()))) {
    legal=TRUE;
    request->set_Solvable(TRUE);
  }
  else
    legal=FALSE;
}

// ***********************************************************************
ostream& operator<< (ostream& s,const Manswer& a) 
{
  s << a.name << endl;
  s << " based on: " << a.request->get_name() << endl;

  if (!a.path)
    s << " - No solution created." << endl;

  else {
    s << " (" << a.gnpath() << " paths)" << endl;
    for (int i = 0; i < a.gnpath(); i++)
      s << " Answer " << i << "\t cost: " << a.pathenergy[i]
	<< "\tPath (nodes):  " << "\t" << a.path[i] << endl;
    s << "  Total cost: " << a.totalenergy;
    if(a.legal )
      s << " (legal)";
    else
      s << " (illegal)";
    s << endl;
  }

  return s;
}



ostream& operator<<(ostream& s,const Path& p)
{

  for (Path::const_iterator i=p.begin();i!=p.end();i++) 
    s << *i << ' ';

  return s;
}
