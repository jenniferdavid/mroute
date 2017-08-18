/* -------------------------------------------------------------------- 
   
     filename: Option.cc  	program: mroute  	project: mroute
     written at Theoretical Physics 2, Lund University
  									
     this file contains the following functions:			
  	parnls		analysis of the given parameters to the 	
  			program, if the parameters are accepted then the	
  			variables are changed accordingly		
  									
   -------------------------------------------------------------------- */

#include "Option.hh"

#include <fstream>
#include <stdlib.h>
#include <cstring>
#include "random.hh"

#define FALSE 0
#define TRUE 1

using namespace std;

Option::Option(int argc,char** argv,char* filename/*="mroute.hlp"*/)
  : fixcapacity(FALSE), maxcapacity(6), networkfile("gennet.data"), nl(6),
    nn(4), NNflag(-1), nr(3), random_sequential_BF(FALSE), randomflag(2),
    requestfile("genreq.data"), seed(-1), sx(FALSE)
{
  for (int i=0;i<16;i++) 
  NNfile[i]="fuzzybellman.parameters";
  //NNfile[6]="na-6-parameters";
  //NNfile[14]="na-14-parameters";
  NNfile[6]="NN.parameters";
  NNfile[14]="NN.parameters";


  char ch=' ';		/* 1st check if help option given	*/
  for (int i=1;i<argc;i++)
    if ((!strcmp(argv[i],"-help")) || (!strcmp(argv[i],"-h"))) {
      ifstream	input(filename);
      cout << '\n';
      if (!input.is_open())
        cout << "Option::Option:    Cannot open helpfile: " << filename << '\n';
      else
        while ((ch=input.get())!=EOF)
          cout << ch;
      cout << '\n';
      exit(1);			/* always exit after printing help	*/
    }
  for (int i=1;i<argc;i++) {	/* check the rest of the options	*/
    int ok=FALSE;
    if (!strcmp(argv[i],"-cf")) {
      fixcapacity=!fixcapacity;
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-cm")) {
      if ((++i)<argc)
 	maxcapacity=atoi(argv[i]);
      if (maxcapacity<1) {
 	cerr << "Option::Option:    Erroneous value for maxcapacity: "
	     << maxcapacity << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-gn")) {
      if ((++i)<argc)
	networkfile=argv[i];
      if (!networkfile.length()) {
	cerr << "Option::Option:    Null length of network file name\n";
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-gr")) {
      if ((++i)<argc)
	requestfile=argv[i];
      if (!requestfile.length()) {
	cerr << "Option::Option:    Null length of requests file name\n";
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-na")) {
      if ((++i)<argc)
 	NNflag=atoi(argv[i]);
      if ((NNflag<-1) || (NNflag>16)) {
 	cerr << "Option::Option:    Erroneous value for NNflag: " << NNflag
	     << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-nl"))	{
      if ((++i)<argc)
 	nl=atoi(argv[i]);
      if (nl<1) {
 	cerr << "Option::Option:    Erroneous value for nl: " << nl << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-nn"))	{
      if ((++i)<argc)
 	nn=atoi(argv[i]);
      if (nn<1) {
 	cerr << "Option::Option:    Erroneous value for nn: " << nn << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-nr"))	{
      if ((++i)<argc)
 	nr=atoi(argv[i]);
      if (nr<1) {
 	cerr << "Option::Option:    Erroneous value for nr: " << nr << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-pf")) {
      if ((++i)<argc)
	NNfile[0]=argv[i];
      if (!NNfile[0].length()) {
	cerr << "Option::Option:    Null length of FZ parameter file name\n";
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-pn")) {
      if ((++i)<argc)
	NNfile[1]=argv[i];
      if (!NNfile[1].length()) {
	cerr << "Option::Option:    Null length of NN parameter file name\n";
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-ra"))	{
      if ((++i)<argc)
	randomflag=atoi(argv[i]);
      if ((randomflag<0) || (randomflag>2)) {
	cerr << "Option::Option:    Illegal value for randomflag: "
	     << randomflag << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-rsbf"))	{
      random_sequential_BF=!random_sequential_BF;
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-se"))	{
      if ((++i)<argc)
	seed=atoi(argv[i]);
      if (seed<0) {
	cerr << "Option::Option:    Illegal value for seed: " << seed << '\n';
	cerr << "Option::Option:    Exiting...\n";
	exit(-1);
      }
      ok=TRUE;
    }
    if (!strcmp(argv[i],"-sx")) {
      sx=!sx;
      ok=TRUE;
    }
    if (!ok)
      cout << "Option::Option:    Invalid option: " << argv[i] << '\n';
  }

  // randomize if necessary
  if (randomflag)
    if (seed==-1)
      seed=Randomize();
    else
      Randomize(seed);

}




ostream& operator<< (ostream& s, const Option& o)
{
  s << "\nDepending on the set of used options, the program will of course\n"
    << "behave differently, but be careful, some options disable, or expects\n"
    << "sensible values on, others options.\n"
    << "See file 'mroute.hlp' for more information, or type 'mroute -h' for\n"
    << "a printout of the file to stdout.\n";

  s << "\nThe option configuration is:\n\n";

  s << "Defined file names:\n";
  s << "\tnetwork file:\t\t\t" << o.networkfile << '\n';
  s << "\trequest file:\t\t\t" << o.requestfile << '\n';
  s << '\n';

  s << "Number of nodes:\t" << o.nn << '\n';
  s << "Number of links:\t" << o.nl << '\n';
  s << "Number of requests:\t" << o.nr << '\n';

  s << "\nUsing ";
  if (o.fixcapacity)
    s << "fixed capacity links. The capacity is ";
  else
    s << "random capacity links. The maximum capacity is ";
  s << o.maxcapacity << '\n';

  s << "\nRandom flag: " << o.randomflag << ", NNflag: " << o.NNflag 
    << ", and seed: " << o.seed << '\n';

  switch (o.get_randomflag()) {
  case 0:
    s << "\tBoth network and request data are read from file\n";
    break;
  case 1:
    s << "\tNetwork data is read from file, request data is randomized\n";
    break;
  case 2:
    s << "\tBoth network and request data are randomized\n";
  }

  s << "\nSolveExact flag is set to ";
  if (o.sx)
    s << "TRUE\n";
  else
    s << "FALSE\n";

  s << "\nThe following NA parameter files are used:\n";
  for (int i=0;i<16;i++)
    s << "\tNA[" << i << "]:\t" << o.NNfile[i] << '\n';

  return s;
}
