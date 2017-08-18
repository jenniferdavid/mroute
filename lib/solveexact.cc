/*
  
  $Id: solveexact.cc,v 1.27 1997/08/07 15:39:28 jari Exp $

  project: mroute

  this file contains the following functions:
	extractpaths
	initload
	solveexact

  */

#define ELC 50000000	// cut value for energy loops
#define LLC 3000000	// cut value for load loops

#include "net.hh"
#include "manswer.hh"
#include "mreq.hh"
#include "mroute.hh"
#include "Option.hh"
#include <math.h>
#include <iostream>

using namespace std;



void extractpaths(Net&,ManyPaths&,Path&,int*,const int,const int);
int initload(const Net&,ManyPaths::iterator*,double**,int,int);
// the ManyPaths::iterator* was a Pix*

void extractpaths(Net& net,ManyPaths& paths,Path& pathsofar,int* nodesflags,
		  const int here,const int end)
{
  // looking for loops
  if (nodesflags[here])
    return;

  pathsofar.push_front(here);
  nodesflags[here]++;
  // checking if end of path reach
  if (here==end) {
    // return the last (end) node to the path
    Path linkpath;
    // converting the nodes list to a links list
    net.n2l(pathsofar,linkpath);
    paths.push_back(linkpath);
    pathsofar.pop_front();
    nodesflags[here]--;
    return;
  }
  
  Node* thisnode=net.get_node(here);
  for (int i=0;i<thisnode->get_nl();i++)
    extractpaths(net,paths,pathsofar,nodesflags,thisnode->node(i),end);
  pathsofar.pop_front();
  nodesflags[here]--;
  return;
} // of extractpaths




int initload(const Net& net,ManyPaths::iterator* reqp,
	     double** load,int nreq,int level)
// the ManyPaths::iterator* was a Pix*
{
  if (!level)
    for (int j=0;j<net.get_nl();j++)		// setting 'zero' level
      load[0][j]= -net.get_link(j)->get_capacity();
    
  // setting load configuration, and checking if load violations occur
  for (int i=level;i<nreq;i++) {
    for (int j=0;j<net.get_nl();j++)
      load[i+1][j]=load[i][j];
    Path path= *(reqp[i]);
    for (Path::iterator k=path.begin();k!=path.end();k++) {
      load[i+1][net.get_link(*k)->getindex()]=
	load[i][net.get_link(*k)->getindex()]+1;
      if (load[i+1][net.get_link(*k)->getindex()]>0)
	return i;	// load violation
    }
  }
  
  return -1;		// no load violations occured
} // of initload



void solveexact(Mreq& request, Manswer& answer, const Option& op)
{
  Net& net = *request.get_N();
  answer.sname("MA:" + request.get_name() + ":SX");

  int nreq=request.get_nr();
  ManyPaths reqpaths[nreq];	// list of links used in a path!!

  // indexmap is used to map the request order with the sorting order below 
  int indexmap[nreq];
  for (int i=0;i<nreq;i++)
    indexmap[i]=0;

  // extracting and sorting all possible paths for all requests
  Path pathsofar;		// list of nodes used in present path
  ManyPaths reqpath;
  int nodesflags[net.get_nn()]; // flags to mark nodes used in present path
  for (int i=0;i<net.get_nn();i++) // used to speed up program
    nodesflags[i]=0;
  for (int i=0;i<nreq;i++) {
    // extract possible paths for this request
    extractpaths(net,reqpath,pathsofar,nodesflags,request.get_startnode(i),
		 request.get_endnode(i));

    // sorting all paths in cost order
    double pathcosts[reqpath.size()];
    double thiscost;
    ManyPaths::iterator pathpix[reqpath.size()];
    int thisindex=0;
    for (ManyPaths::iterator j=reqpath.begin();j!=reqpath.end();j++) {
      thiscost=net.pathcost(*j);
      for (int k=thisindex-1;TRUE;k--)
	if (k==-1) {
	  pathcosts[0]=thiscost;
	  pathpix[0]=j;
	  break;
	}
	else
	  if (thiscost<pathcosts[k]) {
	    pathcosts[k+1]=pathcosts[k];
	    pathpix[k+1]=pathpix[k];
	  }
	  else {
	    pathcosts[k+1]=thiscost;
	    pathpix[k+1]=j;
	    break;
	  }
      thisindex++;
    }

    // transferring the sorted list
    for (unsigned int j=0;j<reqpath.size();j++)
      reqpaths[i].push_back(*pathpix[j]);

    // sorting in descending order of number of paths for a request
    // the next line takes care of the situation when the if statment in 
    // the loop below is FALSE for all loops
    indexmap[i]=i;
    for (int j=0;j<i;j++) {
      if (reqpaths[i].size()>reqpaths[indexmap[j]].size()) {
	for (int k=nreq-1;k>j-1;k--) 
	  indexmap[k]=indexmap[k-1];
	indexmap[j]=i;
	break;
      }
    }

    // clear the temporary list
    reqpath.erase(reqpath.begin(),reqpath.end());
  } // finished extarcting paths and sorting

  // looking for global minimum
  // initalizing path pointers
  ManyPaths::iterator reqp[nreq];
  ManyPaths::iterator minreqp[nreq];
  for (int i=0;i<nreq;i++)
    minreqp[i]=reqp[i]=reqpaths[indexmap[i]].begin();

  // initilizing variables to keep track of costs
  // cost[0] is used but is always 0, introduced to make code symmetric for
  // the costs, i.e. what normally is read cost[0] should be cost[1]
  // all other, except the load matrix, arrays are 'normal'
  double cost[nreq+1];
  cost[0]=0;
  for (int i=0;i<nreq;i++)
    cost[i+1]=net.pathcost(*reqp[i]) + cost[i];
  double mincost=cost[nreq];

  // initilizing loads on the net
  // load[0][*] is used but is always at 'zero' level, introduced to make code
  // symmetric for the loads, i.e. what normally is read load[0][*] should be
  // load[1][*]
  // all other, except the cost vector, arrays are 'normal'
  double* load[nreq+1];
  for (int i=0;i<nreq+1;i++)
    load[i]=new double[net.get_nl()];
  // level keeps track of the first level where loads go out of bounds, or
  // which level to be updated
  int level=initload(net,reqp,load,nreq,0);
  int goodsolution=FALSE;
  int flag=TRUE;
  if (level==-1) {
    level=nreq-1;
    goodsolution=TRUE;
    flag=FALSE; //no need to loop over configurations, this is the best solution
  }
  else { // the illegal initial solution can give small cost, or a bound exists
    mincost=request.get_eb();
    if (mincost>0)
      mincost+=0.001;	// bound exists
    else
      mincost=1e10;	// no known bound, use something large
  }

  // using different cutoffs for loops due to energy and load bounds
  int energyloop=0;
  int loadloop=0;
  // start check configurations
  while (TRUE) {
    while ((reqp[level]==reqpaths[indexmap[level]].end()) || (flag)) {
      // next configuration
      reqp[level]++;
      if (reqp[level]!=reqpaths[indexmap[level]].end()) {

	// check energy
	flag=FALSE;
	for (int j=level;j<nreq;j++) {
	  cost[j+1]=net.pathcost(*(reqp[j])) + cost[j];
	  if (cost[j+1]>mincost) {
	    if (j) {
	      reqp[j]=reqpaths[indexmap[j]].begin();
	      level=j-1;
	    }
	    flag=TRUE;
	    energyloop++;
	    break;
	  }
	}
	
	// check loads if energy was ok
	if (!flag) {
	  level=initload(net,reqp,load,nreq,level);
	  flag=level+1;		// flag becomes FALSE if loads ok!!
	  if (flag)
	    loadloop++;
	}
      }
      else
	if (level) {
	  reqp[level]=reqpaths[indexmap[level]].begin();
	  level--;
	}
	else
	  break; // out of while ((!reqp[level]) || (flag))
      
      if ((energyloop>ELC) || (loadloop>LLC)) {
	cerr << "# loop counter reached cutoff, giving up ..." << endl;
	break;			// reached cutoff
      }
    }
    
    // check if all necessary configurations have been checked
    if ((cost[1]>mincost) || (reqp[0]==reqpaths[indexmap[0]].end()) || 
	(energyloop>ELC) || (loadloop>LLC))
      break;		// only escape from while (TRUE) loop

    // minimum reached
    mincost=cost[nreq];
    goodsolution=TRUE;
    for (int i=0;i<nreq;i++)
      minreqp[i]=reqp[i];

    // reseting counters and flags
    level=nreq-1;
    flag=TRUE;
    energyloop=loadloop=0;
  }

  for (int i=0;i<nreq+1;i++)
    delete[] load[i];

  if ((energyloop<ELC) && (loadloop<LLC))
    answer.setSolved(TRUE);

  if (!goodsolution) {
    request.set_Solvable(FALSE);
    return;
  }
  request.set_Solvable(TRUE);

  // storing the result to a Manswer
  Path* path=new Path[answer.gnpath()];
  for (int i=0;i<nreq;i++)
    net.l2n(*minreqp[i],path[indexmap[i]]);

  answer.spath(path);
} // of solveexact
