/************************************************/
/*     neuralgclass.hh                          */
/* Multi (point to point)			*/
/************************************************/

#ifndef NeuronAlgClass_hh
#define NeuronAlgClass_hh

#include <stdio.h> // needed for popen
#include <iostream>
#include <string>
#include "PottsNeuron.hh"
#include "Propagator.hh"
#include "mroute.hh"
#include "manswer.hh"
#include "mreq.hh"
#include "net.hh"
#include "Option.hh"


class NN
{
public:
  NN(Mreq& req, Manswer& answer, const Option& op);
  // use check_P online in "run" or to check the final status from the main prog.
  double check_P(int req_no);
  ~NN();

private:
  void 	 run(std::string parfile);	// iterate and produce the answer
  //createanswer stores the result in answer
  void 	 createanswer(const std::string&, const Net& net, Manswer& answer);
  int 	 update(int r, PottsNeuron& v, double kT, double Load_coeff, double Pji_coeff);
  double mean_d(PottsNeuron & vi, double * di);
  double loop_cost(int r, int l, int ni); // The strenght of loop formation supression
  double load_cost(int);// The strenght of load supression
  double get_sat();	// calculates the saturation
  double get_psat();	// calculates the saturation
  void   reinit_vPdL();	// reinits v as high T limit and then P,d and L concistently
  double E_w(int r);	// sum of used links for req r (only used as additional info)
  double E_loop(int r); // Sum Pii (count loops)

  int n_n;	// global number of nodes  
  int n_l;	// global number of links  
  int n_r;	// number of requests      
  PottsNeuron** vall; // all the neurons , one vector for each req. 
  Propagator** P; // Propagator, one matrix for each req.([n_r][n_n][n_n])
  int *start, *endn;	// start/end node in a req.
  double** d;	// local average cost of link i->next ([n_r][n_n])
  double* Load;	// the load of all transmissions
  double* maxLoad;// The capacity of a link
  double* w;	// The cost to use a link
  int dummy;	// The index for the dummy-neuron
  double **Elocall;	//saves Eloc only needed for DSOFT = 0, i.e. D hard !!! 
  int DN, DSOFT, PSOFT, EA;	// NN_flags temorary solution !!
  std::string set_NNflags(const Option& op);		//  - " -
};
#endif
