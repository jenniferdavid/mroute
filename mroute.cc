/*

  $Id: mroute.cc,v 1.7 1997/08/27 09:14:11 jari Exp $

  */

#define BELLMAN 1
#define SEQ_BELLMAN 1
#define SEQ_BELLMAN_TRIES 1
#define XOUT 1

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <manswer.hh>
#include <mreq.hh>
#include <mroute.hh>
#include <net.hh>
#include <rmreq.hh>
#include <rnet.hh>
#include <Option.hh>

using namespace std;

int main(int argc,char *argv[])
{

  Option op(argc,argv);

  Net net(op);
  Mreq req(op,&net);

  double entropy;

  RNet rnet(net);
  RMreq rreq(rnet,req);

//  Manswer* nanswer;
//  Manswer* sanswer;
  
  entropy=req.get_entropy();


  for (int i=0; i<argc; i++)
    cout << argv[i] << ' ';

  cout << "\tSubnets: " << rnet.get_n_S() << endl;

  for (int i=0; i<rreq.get_nMR(); i++)
    cout << "\tnn: " << rreq.get_MR(i)->get_N()->get_nn()
	 << "\tnl: " << rreq.get_MR(i)->get_N()->get_nl()/2
	 << "\tnr: " << rreq.get_MR(i)->get_nr()
	 << "\t--> " << rreq.get_MR(i)->get_N()->get_nn()*
      rreq.get_MR(i)->get_N()->get_nl()/2*rreq.get_MR(i)->get_nr()
	 << '/' << req.get_N()->get_nn()*req.get_N()->get_nl()/2*req.get_nr()
	 << '=' << (double)rreq.get_MR(i)->get_N()->get_nn()*
      rreq.get_MR(i)->get_N()->get_nl()*rreq.get_MR(i)->get_nr() /
      req.get_N()->get_nn()/req.get_N()->get_nl()/req.get_nr()
	 << endl;

 // exit(0);

  Manswer* banswer=NULL;
  Manswer* nanswer=NULL;
  Manswer* sanswer=NULL;
  Manswer* sbanswer=NULL;
  

  if (op.get_sx()) {
    sanswer = new Manswer(op,rreq,solveexact);
    if (XOUT)
      cout << "SolveExact done." << endl;
  }

  if (op.get_NNflag()>-1) {
    nanswer = new Manswer(op,rreq,NN_alg);
    if (XOUT)
      cout << "NeurAlg done." << endl;
  }

  if (BELLMAN) {
    banswer = new Manswer(op,rreq,bellman);
    if (XOUT)
      cout << "Bellman-Ford done." << endl;
  }

  if (SEQ_BELLMAN) {
    for (int i=0; i<SEQ_BELLMAN_TRIES; i++) {
      sbanswer = new Manswer(op,rreq,seq_bellman);
      if (!sbanswer->isLegal())
	delete sbanswer;
      else
	break;
      cout << "Failed attempt " << i+1 << " of " << SEQ_BELLMAN_TRIES
	   << " tries." << endl;
    }
    if (XOUT)
      cout << "Sequential Bellman-Ford done." << endl;
  }

  // write used seed
  cout << "\n Writing used seed \n";
  cout.width(7);
  cout << op.get_seed();

  // writing network information
  cout << "\n Writing n/w info \n";
  cout.width(7);
  cout << net.get_nn();
  cout.width(7);
  cout << net.get_nl()/2;
  cout.width(7);
  cout << req.get_nr();
  cout.width(7);
  cout.precision(3);
  cout << entropy;

  // the exact result
  if (op.get_sx()) {
    cout.width(7);
    cout.precision(4);
    cout << "\n Writing exact result \n";
    cout << ' ' << sanswer->genergy();
    if (sanswer->isSolved()) {
      cout.width(7);
      cout << sanswer->isLegal();
    }
    else {
      cout.width(7);
      cout << -1;
    }
  }

  // the NeurAlg result
  if (op.get_NNflag()>-1) {
    cout.precision(4);
    cout.width(7);
    cout << "\n Writing NeurAlg result \n";
    cout << ' ' << nanswer->genergy();
    cout.width(7);
    cout << nanswer->isLegal();
  }

  if (BELLMAN) {
    cout.precision(4);
    cout.width(7);
    cout << "\n Writing BF result \n";
    cout << ' ' << banswer->genergy();
    cout.width(7);
    cout << banswer->isLegal();
  }
    
  if (SEQ_BELLMAN) {
    cout.precision(4);
    cout.width(7);
    cout << "\n Writing SEQ BF result \n";
    cout << ' ' << sbanswer->genergy();
    cout.width(7);
    cout << sbanswer->isLegal();
  }
    
  cout << '\n';

    plot(2,&net,&req,sbanswer);
  
  if (XOUT) {
    cout << "# The problem ";
    if (!req.isSolvable())
      cout << "was not ";
    else
      if (req.isSolvable()==1)
	cout << "was ";
      else
	cout << "might be ";
    cout << "solvable." << endl;
  }

  return 0;
}
