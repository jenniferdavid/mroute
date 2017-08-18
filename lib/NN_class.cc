/************************************************/
/*    NN_class.cc                               */
/* Multi (point to point)			*/
/************************************************/

#include <fstream>
#include <math.h>
#include "NN_class.hh"
#include <iomanip> // needed for setw(int)
#include <string>

using namespace std;
/*****************************************/
NN::NN(Mreq& req, Manswer& answer, const Option& op)
{
  std::string parfile = set_NNflags(op);
  const Net& net = *req.get_N();
  
  int n_n_real = net.get_nn();	// the dummy neuron not included
  int n_l_real = net.get_nl();  // acces to "real" (not to dummy) n_l is nedded in load, ..
  n_n = n_n_real + DN; 		// Same is nn_real if no dummy neuron (DN = 0)
  n_l = n_l_real + DN*2*n_n_real;// all but dummy connects to dummy, dummy connects back.
  n_r = req.get_nr();		// The number of requests
  dummy = n_n - 1;     		// index for the dummy neuron
  // So n_n and n_l includes dummy neuron (if DN) if no dummy is used _real is the same number.
  
  /* start and endpoint of the requests */
  start = new int[n_r];
  endn = new int[n_r];
  for (int i = 0; i < n_r; i++) 
    {
      start[i] = req.get_startnode(i);
      endn[i] = req.get_endnode(i);
    }

  // init and allocate a Potts-neuron for every node in the net
  // when v is created all values is put to  1/#components
  int tmpn[n_n];// could connected to all neurons, the used lenght is determined from tmp_nc
  int tmpl[n_n];// index of link used to reach the neuron neighbour (tmpn) with the same index
  int tmplabel;	// Assign the label of the neuron to the index of the corresponding node
  int tmpn_c;	// used for # outgoing conn. from a neuron, i.e. #of comp in the Potts spin
  int addlinkno = n_l_real; // Index of the extra links added due to the dummy neuron.

  vall = new PottsNeuron*[n_r];	// one set for each req.
  for (int r = 0; r < n_r; r++) 
    vall[r] = new PottsNeuron[n_n];
  
  for (int i = 0; i < n_n_real; i++)
    {
      tmpn_c = net.get_node(i)->get_nl(); 
      tmplabel = net.get_node(i)->getindex();
      if (tmplabel != i)
	std::cerr << "#NN:Error: index not correct (in nodes?)" << std::endl;
      for (int j = 0; j < tmpn_c; j++)
	{
	  tmpn[j] = net.get_node(i)->node(j);
	  tmpl[j] = net.get_node(i)->get_link(j)->getindex();
	}
      if (DN)	// all connected to the dummy neuron
	{
	  int j = tmpn_c;
	  tmpn_c++;	// one more component
	  // tmplabel not changed
	  tmpn[j] = dummy;
	  tmpl[j] = addlinkno++;
	 }
      vall[0][i] = PottsNeuron(tmpn_c,tmpn,tmpl,tmplabel);
    }


  if (DN) // The dummy neuron
    {
      for (int i = 0; i < n_n_real; i++) // connected to all normal neurons
	{ 
	  tmpn[i] = i;
	  tmpl[i] = addlinkno++;
	}
      vall[0][dummy] = PottsNeuron(n_n_real,tmpn,tmpl,dummy);
      if ( addlinkno != n_l )
	std::cerr<<"#NN:error addlinkno="<<addlinkno<<" but n_l="<<n_l<<std::endl;
    }
  

  // Make a copy of vall[0] for each req.
  for (int r = 1; r < n_r; r++)	// n.b. r=0 already done (i.e. vall[0])
    for (int i = 0; i < n_n; i++)
      vall[r][i] = PottsNeuron(vall[0][i]);
  
  // no outgoing paths from the endnode, the endnode is the sink, so have to modify each v[r].
  for (int r = 0; r < n_r; r++) 
    {
      for (int j = 0; j < vall[r][endn[r]].n(); j++) 
	vall[r][endn[r]][j] = 0;
      if (DN) // The dummy node should only connect to the end node
	{
	  for (int i = 0; i < vall[r][dummy].n(); i++) 
	    vall[r][dummy][i] = 0;
	  // n.b. global and local index is same for the dummy neuron.
	  vall[r][dummy][endn[r]] = 1; 
	}
    }
  // init v done
  
  
  // init and allocate P
  P = new Propagator*[n_r];
  for (int r = 0; r < n_r; r++) 
    P[r] = new Propagator(n_n,vall[r]);  
  
  // init and allocate w and maxLoad
  maxLoad = new double[n_l];
  w = new double[n_l];		// the cost to use the link
  double wsum = 0;
  //ordinary links
  for (int i = 0; i < n_l_real; i++) 
    {
      maxLoad[i] = net.get_link(i)->get_capacity();  
      w[i] = net.get_link(i)->getweight();
      wsum += w[i];
    }
  if (DN)  // links to and from dummy node
    for (int i = n_l_real; i < n_l ; i++) 
      {
	maxLoad[i] = n_r+1;	// all should be able to connect to the dummy node
	// move magic numbers to file .. !!! /ML
	w[i] = wsum;  // at more cost than going via all the rest
      }
  
  // init link load
  Load = new double[n_l];
  for (int i = 0; i < n_l; i++) 
    Load[i] = 0;
  for (int r = 0; r < n_r; r++)	// all the transmissions
    for (int i = 0; i < n_n; i++)
      for (int j = 0; j < vall[r][i].n(); j++)
	Load[vall[r][i].link(j)] += P[r]->g(start[r],i)*vall[r][i][j]/P[r]->g(i,i);
  // Mult with P to filter out the the path starting from the start node 
  // from the tree produced by the BF-like alg.. Divide with P[i][i] to 
  // normalize, don't count loop contributions.
  
  
  // Allocate and init d[r][i] (cost from i to endnode (of req. r)) 
  // In ordinary BF, d[i] is init as infty if i != end, else 0 
  // Di = sum_k Pik sum_j Vkj Dkj gives the "real" hight temperature limit.
  d = new double*[n_r];
  for (int r = 0; r < n_r; r++) 
    {
      d[r] = new double[n_n];  
      for (int i = 0; i < n_n; i++) 
	d[r][i] = 0;	
      for (int k = 0; k < n_n; k++)
	{
	  double tmp = 0;
	  for (int l = 0; l < vall[r][k].n(); l++)
	    tmp += vall[r][k][l] * w[vall[r][k].link(l)];
	  for (int i = 0; i < n_n; i++)
	    d[r][i] += tmp * P[r]->g(i,k);
	}
      if(DN)
	d[r][dummy] = 0;	// the only cost should be the link weight
    }
  
  // remove or improve when alg variant is fixed !!!!
  Elocall = new double*[n_r];
  for (int r = 0; r < n_r; r++) 
    {
      Elocall[r] = new double[n_n];
      for (int i = 0; i < n_n; i++) 
	Elocall[r][i] = d[r][i];
    }

  
  // GO!
  run(parfile);
  if(DEBUG)
    {
      std::cout << "P[][]:\n" << *P[0] << std::endl;
      std::cout << "v[][]:\n";
      for (int i = 0; i < n_n; i++) 
	std::cout << vall[0][i];
      std::cout << "i w[]:\n";
      for (int i = 0; i < n_l; i++) 
	std::cout << i << "\t" << w[i] << std::endl;
    }
  createanswer("Ma:" + req.get_name() + ":NN", net,answer);
}

/*****************************************/

void NN::reinit_vPdL()
{// used if the start kT is raised due to a too large step in the sat.
  
  /* reinit v */
  for (int r = 0; r < n_r; r++) 
    for (int i = 0; i < n_n; i++)
      for (int l = 0; l < vall[r][i].n(); l++) 
	vall[r][i][l] = 1./vall[r][i].n();

  // no outgoing paths from the endnode, the endnode is the sink, so have to modify each v[r].
  for (int r = 0; r < n_r; r++) 
    {
      for (int j = 0; j < vall[r][endn[r]].n(); j++) 
	vall[r][endn[r]][j] = 0;
      if (DN) // The dummy node should only connect to the end node
	{
	  for (int i = 0; i < vall[r][dummy].n(); i++) 
	    vall[r][dummy][i] = 0;
	  // n.b. global and local index is same for the dummy neuron.
	  vall[r][dummy][endn[r]] = 1; 
	}
    }

  /* reinit P, by deleting it and recreating it */
  delete [] P;
  P = new Propagator*[n_r];
  for (int r = 0; r < n_r; r++) 
    P[r] = new Propagator(n_n,vall[r]);  
  
  /* reinit d */
  for (int r = 0; r < n_r; r++) 
    {
      for (int i = 0; i < n_n; i++) 
	d[r][i] = 0;	
      for (int k = 0; k < n_n; k++)
	{
	  double tmp = 0;
	  for (int l = 0; l < vall[r][k].n(); l++)
	    tmp += vall[r][k][l] * w[vall[r][k].link(l)];
	  for (int i = 0; i < n_n; i++)
	    d[r][i] += tmp * P[r]->g(i,k);
	}
      if(DN)
	d[r][dummy] = 0;	// the only cost should be the link weight
    }

  // reinit link load
  for (int i = 0; i < n_l; i++) 
    Load[i] = 0;
  for (int r = 0; r < n_r; r++)	// all the transmissions
    for (int i = 0; i < n_n; i++)
      for (int j = 0; j < vall[r][i].n(); j++)
	Load[vall[r][i].link(j)] += P[r]->g(start[r],i)*vall[r][i][j]/P[r]->g(i,i);
  

}

/*****************************************/
void NN::run(std::string parfile)
{
  int r; //  (request) 
  int i; //  ( global node#)
  int j; //  ( global node#)
  int FLAG = TRUE;
  double sat = get_sat();
  double E;
  double oldsat = sat;
  int warningp1 = 0;
   

  /** read parameters **/
  // default values:
  double kT_start = 10.0;
  double kT_stop  = 0.001;
  double kT_swfac = 0.95;
  double sat_stop = .9999;
  double Pji_coeff = .001;
  double Load_coeff = 10.;
  int PLOT = 0;
  int MONITOR = 0;
  ifstream  input(parfile.c_str());
  if (!input) 
    {
      cerr << "#NN: cannot open input file '" << parfile << "'\n";
      cerr << "#    Default values will be used \n";
    }
  else
    {
      char ch=' ';
      while (ch!='@') input.get(ch);
      input >> kT_start;
      ch=' '; while (ch!='@') input.get(ch);
      input >> kT_stop;
      ch=' '; while (ch!='@') input.get(ch);
      input >> kT_swfac;
      ch=' '; while (ch!='@') input.get(ch);
      input >> sat_stop;
      ch=' '; while (ch!='@') input.get(ch);
      input >> Pji_coeff;
      ch=' '; while (ch!='@') input.get(ch);
      input >> Load_coeff;
      ch=' '; while (ch!='@') input.get(ch);
      input >> PLOT;
      ch=' '; while (ch!='@') input.get(ch);
      input >> MONITOR;
      if (DEBUG)
	{
	  cerr << "#kT_start =" << kT_start << endl; 
	  cerr << "#kT_stop  =" << kT_stop  << endl; 
	  cerr << "#kT_swfac =" << kT_swfac << endl;  
	  cerr << "#sat_stop =" << sat_stop << endl;
	  cerr << "#Pji_coeff =" << Pji_coeff << endl;
	  cerr << "#Load_coeff =" << Load_coeff << endl;
	  cerr << "#PLOT =" << PLOT << endl;
	  cerr << "#MONITOR =" << MONITOR << endl;
	}
      input.close();
    }
  double kT_fac = exp( log(kT_swfac) / (n_r * n_n) ); // lower T after every neuron update
  double kT = kT_start * kT_fac;  // * to give the ini conf a own datapoint in sat plot


  /** init plot **/
  FILE* command;
  char* plotfile="convf.plot";
  ofstream* output;
  if (PLOT)
    {
      output = new ofstream(plotfile);
      output->setf(ios::fixed);
      output->precision(4);
      // write ini sat & E to the file
      *output << 1./kT << " " << sat << " " << 0 << " " << 0 << endl;	
    }


  /** start iterating **/
  int nsw = 0; //number of sweeps
  double loadfac;
  while(FLAG)
    {
      if (MONITOR || PLOT )// compute E = \sum_r d[A][B]
	{
	  E = 0;
	  for( r = 0; r < n_r; r++)
	    E += d[r][start[r]];
	}
      if( MONITOR )
	{
	  if( nsw % 30 == 0 )
	    {
	      cout << "#sw \tsat \t\tkT \t\tE  \t\tEw \t\tloopcount";
	      if(DN)
		cout << " \tsumdummy \tmaxdummy";
	      cout  << endl;
	    }

	  sat =  get_sat();
	  double Ew = 0;
	  for (r = 0; r < n_r; r++) 
	    Ew += E_w(r);	      
	  double loopcount = 0;
	  for (r = 0; r < n_r; r++)
	    loopcount += E_loop(r);

	  cout.setf(ios::fixed,ios::floatfield);cout.precision(8);
	  cout << setw(5) << nsw;
	  cout << " \t" << setw(5) << sat;
	  cout << " \t" << setw(5) << kT;
	  cout << " \t" << setw(5) << E;
	  cout << " \t" << setw(5) << Ew;
	  cout << " \t" << setw(5) << loopcount;
	  if (DN) // Print some stuff about the dummy neuron
	    {
	      double sumdummy = 0;
	      double maxdummy = 0;
	      for (r = 0; r < n_r; r++)
		for (i = 0; i < n_n; i++) 
		  if ( i != dummy && i != endn[r])
		    {
		      double tmp = vall[r][i][vall[r][i].n() - 1];
		      sumdummy += tmp;
		      if (maxdummy < tmp)
			maxdummy = tmp;
		    }
	      cout << " \t" << setw(5) << sumdummy ;
	      cout << " \t" << setw(5) << maxdummy ;
	    }
	  cout << endl;
	}

      for (r = 0; r < n_r; r++)	// all the transmissions
	{ 
	  // subtract the old loads for r
	  for (i = 0; i < n_n; i++)
	    {
	      loadfac = P[r]->g(start[r],i) / P[r]->g(i,i);
	      for (j = 0; j < vall[r][i].n(); j++) 
		Load[vall[r][i].link(j)] -= loadfac * vall[r][i][j];
	    }

	  // update neurons for r
	  for (i = 0; i < n_n; i++) 
	    // The endnode (should always be 0, the sink) should not be updated.
	    // do not update the dummy, it's clamped to end
	    if ( i != endn[r] && i != dummy)
	      {
		// update v, d and P
		warningp1 += update(r,vall[r][i],kT,Load_coeff, Pji_coeff);
		kT *= kT_fac;
	      }  
	  
	  // add the new loads for r
	  for (i = 0; i < n_n; i++)
	    {
	      loadfac = P[r]->g(start[r],i) / P[r]->g(i,i);
	      for (j = 0; j < vall[r][i].n(); j++) 
		Load[vall[r][i].link(j)] += loadfac * vall[r][i][j];
	    }

	  // if (PSOFT)	// remove after test !!!
// 	    {
// 	      double check_res = check_P(r);
// 	      if (check_res  > 0.1 * n_n)
// 		{
// 		  cerr << "#NN:run:Check_P: r="<< r << " kT=" << kT; 
// 		  cerr.setf(ios::fixed);
// 		  cerr.precision(5);
// 		  cerr << "\t\t" << "sum_ij( sqr(P_ij_diff) )/n_n=";
// 		  cerr << check_res / (1.* n_n) << endl;
// 		}
// 	    }
	} // end r loop (all requests)
      nsw++;

      
      oldsat = sat;
      sat = get_sat();
      
      
      if (kT < kT_stop || sat > sat_stop)
	FLAG = FALSE;

      if(PLOT)
	*output << 1./kT << " " << sat << " " << E << " " << nsw << endl;

      if (nsw == 1 && (sat-oldsat)/oldsat > 0.1)
	{
	  nsw = 0;
	  reinit_vPdL();
	  cerr << "#NN:Warning temperature raised from " << kT;
	  kT *= 2;
	  cerr << " to "<< kT << "\n";;
	}
    } // end FLAG
  

  if(PLOT) // move up (before end FLAG) for online-plotting
    {
      output->close();
      command = popen("gnuplotbeta","w");
      fprintf(command,"set term x11\n");
      fprintf(command,"set autoscale x\n");
      fprintf(command,"set yr[0:1.1]\n");
      fprintf(command,"set label 'sat' at 1,1.05\n");
      fprintf(command,"plot '%s' u 1:2 w linesp\n",plotfile);
      fflush(command);
      cerr << " DONE sat,  p for print ENTER to quit " << endl;
      if ( getchar() == 'p' )
	{
	  fprintf(command,"set term postscript\n");
	  fprintf(command,"set output 'tmp.ps'\n");
	  fprintf(command,"replot\n");
	  fprintf(command,"!laser tmp.ps\n");
	  fprintf(command,"set term x11\n");
	}
      fprintf(command,"set auto y\n");
      fprintf(command,"set nolabel\n");
      fprintf(command,"set label 'E' at 1,1.05\n");
      fprintf(command,"plot '%s' u 1:3 w linesp\n",plotfile);
      fflush(command);
      cerr << " DONE E,  p for print ENTER to quit " << endl;
      if ( getchar() == 'p' )
	{
	  fprintf(command,"set term postscript\n");
	  fprintf(command,"set output 'tmp.ps'\n");
	  fprintf(command,"replot\n");
	  fprintf(command,"!laser tmp.ps\n");
	  fprintf(command,"set term x11\n");
	}
      pclose(command);
    }

  if (warningp1 > 0)
    {
      cerr << "#NN::update_p 2: Warning " << warningp1;
      cerr << "# numerical problems with 1/(1+lamda)" << endl;
    }
  
} // end run

/*****************************************/
int NN::update(int r, PottsNeuron& v, double kT, double Load_coeff, double Pji_coeff)
{// subroutine to run, r=req#, i=node#, kT=temperature
  double Etest[n_n];
  double min_Etest = HUGE_VAL;
  double arg;
  double sumv=0;
  int l;	// the # of local links (and nodes)
  int ni = v.label(); // neuron index
  //double lk;	// loopkiller 

  PottsNeuron vold(v);  
  double Eloc[v.n()];	// d_ij

  /** calculating updating variables **/
  for (l = 0; l < v.n(); l++)  // n.b. local #links == #nodes
    {
      
      Eloc[l] = (
		 + w[v.link(l)]				// local cost
		 + Load_coeff * load_cost(v.link(l))	// local load-penalty 
		 + Pji_coeff  * loop_cost(r,l,ni)	// local loop penalty
		 );
    }
  
  if (EA)
    {
      double Pref = 0;
      for (l = 0; l < v.n(); l++) 
	{
	  int j = v.ne(l);
	  Pref = P[r]->g(start[r],ni) / ( P[r]->g(ni,ni) - P[r]->g(j,ni) );
	  if ( Pref > 1 ) Pref = 1;
	  Etest[l] = Pref * ( Eloc[l] - d[r][ni] + d[r][j]  );
	  if (Etest[l] < min_Etest)	// use to avoid sumv=0
	    min_Etest = Etest[l];      
	}
    }
  else // this one used in the paper
    {
      //cout << "ni="<<ni<<"Etest:";
      for (l = 0; l < v.n(); l++) 
	{
	  //cout << Eloc[l] <<"+"<< d[r][v.ne(l)]<<" \t ";
	  Etest[l] = ( Eloc[l] + d[r][v.ne(l)]);  
	  if (Etest[l] < min_Etest)	// use to avoid sumv=0
	    min_Etest = Etest[l];      
	}
      //cout << endl;
    }


  /** update v **/
  for (l = 0; l < v.n(); l++) 
    {// Should be enough. Stupid exp() function, dumps around -708.4
      if ((arg = (min_Etest-Etest[l])/kT) < -50. ) 
	v[l] = 0;
      else
	v[l] = exp(arg); // Don't compute twice /BS ((min_Etest - Etest[l])/kT);
      sumv += v[l];
    }
  for (l = 0; l < v.n(); l++) 
    v[l] /= sumv;      
  
  /** update P **/
  int warning = 0;
  if(PSOFT)	// this one used in the paper
    {// P_im = delta_im + sum_j( v_ij * Pjm )
      double pvalue;
      for (int m = 0; m < n_n; m++) 
	{
	  if ( m != ni )
	    pvalue = 0;
	  else
	    pvalue = 1;
	  for (l = 0; l < vall[r][ni].n(); l++)
	    pvalue += vall[r][ni][l] * P[r]->g(vall[r][ni].ne(l) , m );
	  // if ( pvalue > 1.5 ) // !!
	  // 	    {
	  // 	      cerr << "# Warning pvalue = " << pvalue;
	  // 	      cerr << "\t kT=" << kT << "\t ni=" << ni << endl;
	  // 	    }
	  P[r]->s(ni,m,pvalue);
	}
    }
  else 
    warning = P[r]->update(v - vold); 
  

  /** update d_i **/
  if (DSOFT) // this one used in the paper
    {
      d[r][ni] = 0;
      for (l = 0; l < v.n(); l++) 
	{
	  //      d[r][ni] += v[l] * (w[v.link(l)] + d[r][v.ne(l)]);
	  d[r][ni] += v[l] * (Eloc[l] + d[r][v.ne(l)]);
	}
    }
  else
    {
      // update Elocall for this neuron
      Elocall[r][ni] = 0;
      for (l = 0; l < v.n(); l++) 
	Elocall[r][ni] += v[l] * Eloc[l];
      // update all d based on the local energy of all the other
      for (int k = 0; k < n_n; k++) 
       	{
       	  d[r][k] = 0;
       	  for (int j = 0; j < n_n; j++) 
	    d[r][k] += P[r]->g(k,j) * Elocall[r][j];
       	}
    }
  
  
  return warning;
}

// /*****************************************/
// double NN::mean_d(PottsNeuron & vi, double * Eloc)
// { 	// d-hat = sum_j ( v_ij * d_ij )
//   double aux = 0;
//   return aux;
// }

/*****************************************/
double NN::loop_cost(int r, int l, int ni)
{ // I will speed up by handle constants better, inline,  ... later !!! /ML
  static const double small = 1e-15;
  static const double onemsmall = 1 - small;
  static const double lk0 = 1/small - 1; 
  double lk = P[r]->g(vall[r][ni].ne(l),ni)/P[r]->g(ni,ni);	// the "zeroed" Pji
  if ( lk < onemsmall )
    lk = lk/(1-lk); // => the resulting Pji for choice j
  else
    lk = lk0;
  return lk;
}

/*****************************************/
void NN::createanswer(const string& slbl, const Net& net , Manswer& a) // store the result in the answer
{
  a.sname(slbl);

  // the object that receives the answer also has to delete it!      
  int r,j;
  double hardload[n_l];
  for (j = 0; j < n_l; j++)
    hardload[j] = 0;

  Path* path=new Path[n_r];      // this is to be deleted from Manswer
  for (r = 0; r < n_r; r++)
    { 
      int counter = 0;
      int flag=TRUE;
      int winner = start[r]; 
      int old;
      int jariflag=0;
      while (flag) 
	{
	  path[r].push_front(winner);	
	  old = winner;
	  double bsf = 0;
	  for (j = 0; j < vall[r][old].n(); j++)
	    if (vall[r][old][j] > bsf)
	      {
		winner = vall[r][old].ne(j); // need global index
		bsf = vall[r][old][j];
	      }
	  counter++;
	  if (winner == endn[r])
	    flag = FALSE;
	  else if ( winner == dummy )
	    {
	      flag = FALSE;
	      cerr << "#NN: Error: path via dummy formed" << endl;
	      jariflag=TRUE;
	    }
	  else if ( counter >= n_n)
	    {
	      flag = FALSE;
	      cerr << "#\7NN: Error: path not formed" << endl;
	      jariflag=TRUE;
	    }
	} // end while flag

      path[r].push_front(winner);	
      if (jariflag)
	path[r].erase(path[r].begin(),path[r].end());

      a.spath(path);
      
      // extract the produced load
      Path::iterator left=path[r].begin();
      Path::iterator right=left;
      for (right++;right!=path[r].end();right++)
	{ 
	  int iright=*right;
	  for (int l = 0; l < net.get_node(iright)->get_nl(); l++)
	    if ( *left == net.get_node(iright)->node(l) )
	      hardload[net.get_node(iright)->get_link(l)->getindex()] ++;
	  left++;
	}
    } // end r loop 
  
  // Compare the load ANN sees with that it really produced
  double load_diff2 = 0;
  double aux;
  for (j = 0; j < n_l; j++) 
    {
      aux = hardload[j] - Load[j];
      load_diff2 += aux * aux;
    }
  load_diff2 = sqrt(load_diff2);
  if ( load_diff2 > 0.1 )
    {
      cerr << "#NN WARNING The propagator was not on the shell,";
      cerr << "# loads corrected, sqrt(sum( (real-load - netload)^2))=";
      cerr << load_diff2 << endl;
    }

  //return answer;
}

/*****************************************/
double NN::check_P(int req_no)
{
  //  could be used in "run" online or in a main prog. to check the final status
  //  p=1/(1-v) => p-pv-1=0 , used p-vp-1=0 bec. then the sum vp is over local index
  int i,j,k,l;
  double check;
  double sumcheck = 0;
  
  for (i = 0; i < n_n; i++) 
    for (j = 0; j < n_n; j++)
      { 
 	check = P[req_no]->g(i,j);
	for (l = 0; l < vall[req_no][i].n(); l++)
	  {
	    k = vall[req_no][i].ne(l);
	    check -= vall[req_no][i][l] * P[req_no]->g(k,j);
	  }
 	if ( i == j )
 	  check -= 1;
 	check *= check;
 	sumcheck += check;
      }
  return sumcheck;
}    

/*****************************************/
double NN::load_cost(int link)
{
  double F_load_0, F_load_1;
  F_load_0 = Load[link] - maxLoad[link];
  F_load_1 = F_load_0 + 1;
  
  if (F_load_0 < 0)
    F_load_0 = 0;
  
  if (F_load_1 < 0)
    F_load_1 = 0;
  
  return F_load_1 - F_load_0;
}

/*****************************************/
NN::~NN()
{
  int r;
  for (r = 0; r < n_r; r++) 
    {
      delete [] d[r];
      delete [] vall[r];
      delete [] Elocall[r]; // temporary variable solution !!
    }
  delete [] start;
  delete [] endn;
  delete [] d;
  delete [] Elocall;	// !!
  delete [] vall;
  delete [] P;
  
  delete [] Load;
  delete [] maxLoad;

  // Do not delete answer, a pointer to it is passed to the caller
}

/*****************************************/
double NN::get_sat()
{
  int r,i,l;
  double sat = 0;
  for (r = 0; r < n_r; r++) 
    for (i = 0; i < n_n; i++)
      if( i != endn[r] )
	for (l = 0; l < vall[r][i].n(); l++)
	  sat += vall[r][i][l] * vall[r][i][l];
  sat /= n_r * (n_n - 1);
  return sat;
}

/*****************************************/
double NN::get_psat()
{
  int r,i,l;
  double sat = 0;
  double dsat, ddsat;
  double norm = 0;
  double p;

  for (r = 0; r < n_r; r++) 
    for (i = 0; i < n_n; i++)
      if( i != endn[r] )
	{
	  norm += p = P[r]->g(start[r],i);
	  dsat = 0;
	  for (l = 0; l < vall[r][i].n(); l++)
	    {
	      ddsat = vall[r][i][l];
	      dsat += ddsat * ddsat;
	    }
	  sat += p * dsat;
	}
  return sat / norm;
}

/*****************************************/
string NN::set_NNflags(const Option& op)
{
  
  /** default choise of update type **/
  /*
  String parameterfile = op.get_NNfile(op.get_NNflag());
  DN    = 0; // DN=1 if dummy neuron should be used, also used as +1 in vector lenghts
  DSOFT = 1; // else hard update
  PSOFT = 0; // else hard
  EA    = 0; // energy weighted with P_Ai else not
  
  dec-number	DN	DSOFT	PSOFT	EA
  0		0	0	0       0        
  1		0	0	0       1
  2		0	0	1       0
  3		0	0	1       1
  4		0	1	0       0
  5		0	1	0       1
  6		0	1	1       0
  7		0	1	1       1
  8		1	0	0       0
  9		1	0	0       1
  10		1	0	1       0
  11		1	0	1       1
  12		1	1	0       0
  13		1	1	0       1
  14		1	1	1       0 	<= in paper
  15		1	1	1	1
  */
  int nnflag=op.get_NNflag();
  if (nnflag>15) {
    cerr << "# Option::NNflag too large: " << nnflag;
    cerr << "# \n\tForced to be zero." << endl;
    nnflag=0;
  }

  DN = nnflag & 8 ? 1 : 0;
  DSOFT = nnflag & 4 ? 1 : 0;
  PSOFT = nnflag & 2 ? 1 : 0;
  EA = nnflag & 1 ? 1 : 0;
  
  if (DEBUG)
    {
      cout << "# DN = " 		<< DN           << "\t";
      cout << "  DSOFT = " 		<< DSOFT        << "\t";
      cout << "  PSOFT = " 		<< PSOFT        << "\t";
      cout << "  EA = " 		<< EA           << "\t";
      cout << "  parameterfile = " 	<< op.get_NNfile(nnflag) << "\n";
    }
  return op.get_NNfile(nnflag);
}

/*****************************************/
double NN::E_w(int r)
{
  double E = 0;
  int i,j,l;
  for (i = 0; i < n_n; i++) 
    for (l = 0; l < vall[r][i].n(); l++) 
      if ( i != endn[r] && i != dummy)
      {
	j = vall[r][i].ne(l);
	E += P[r]->g(start[r],i) * vall[r][i][l] * w[vall[r][i].link(l)];
      }
  return E;
}

double NN::E_loop(int r)
{
  // OBS not a real energy only a count of "soft" number of loops /ML
  double E = -n_n;
  int i;
  for (i = 0; i < n_n; i++) 
    E += P[r]->g(i,i);
  return E;
}

/********* Paste stuff ********/

/***** Could fastend thing up to clamp links going to startnode to 0 ****/
// no incoming paths to the endnode
//   for (r = 0; r < n_r; r++) 
//     for (i = 0; i < n_n; i++) 
//       for (j = 0; j < vall[r][i].n(); j++) 
// 	if ( vall[r][i].ne(j) == start[r] )
// 	  vall[r][i][j] = 0; 	
