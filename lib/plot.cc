/*
  filename: plot.cc		project: mroute

  this file contains the following functions:
	embed(int,double**,coord*)
	plot(int,Net*,Mreq*,Manswer*)
	
	*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "random.hh"
#include "net.hh"
#include "manswer.hh"
#include "mreq.hh"

using namespace std;

struct coord {
  double x;    	
  double y;
};

double embed(int,double**,struct coord*);

void plot(int plot,Net* net,Mreq* req,Manswer* answer) {
  
  char*		transfile="net.answer.plot";
  char*		reqfile="net.req.plot";
  char*		plotfile="net.plot";
  struct coord	point[net->get_nn()];
  double**	distance;
  FILE*		command;
  int		i, j;
  ofstream	network(plotfile);

  // extract all the distances between nodes, 0 -> infinity, except for the
  // diagonal elements in which 0 represents distance from i to i
  distance=new double* [net->get_nn()];
  for (i=0;i<net->get_nn();i++)
    distance[i]=new double[net->get_nn()];
  for (i=0;i<net->get_nn();i++) {
    for (j=0;j<net->get_nn();j++)
      distance[i][j]=0;
    for (j=0;j<net->get_node(i)->get_nl();j++)
      distance[i][net->get_node(i)->node(j)]=
	net->get_node(i)->get_link(j)->getweight();
  }

  // embedding the network into two dimensions
  for (i=0;i<net->get_nn();i++)
    point[i].x=point[i].y=0;
  embed(net->get_nn(),distance,point);

  // plotting the network
  for (i=0;i<net->get_nn();i++) {		// prepare file for gnuplot
    for (j=0;j<net->get_node(i)->get_nl();j++) {
      network << point[i].x << '\t' << point[i].y << '\n';
      network << point[net->get_node(i)->node(j)].x << '\t';
      network << point[net->get_node(i)->node(j)].y << '\n';
    }
    network << '\n';
  }
  network << endl;	// flushing buffer, otherwise cannot read below 

  if (plot) {		// plotting the transmissions and marking request nodes
    ofstream	transmission(transfile);
    ofstream	request(reqfile);

    // extracting request nodes
    for (i=0;i<req->get_nr();i++) {
      request << point[req->get_startnode(i)].x << '\t'
	      << point[req->get_startnode(i)].y << '\n';
      request << point[req->get_endnode(i)].x << '\t'
	      << point[req->get_endnode(i)].y << '\n';
    }
    request << endl;

    if (plot==2) {	// extracting the transmissions
      Path* path=answer->gppath();
      for (i=0;i<answer->gnpath();i++) {
	for (Path::iterator k=path[i].begin();k!=path[i].end();k++) 
	  transmission << point[*k].x << '\t'
		       << point[*k].y << '\n';
	transmission << '\n';
      }
      transmission << endl;
    }
  }

  // C stuff follows, need a pipe to a process. Can this be done in C++?
  command = popen("gnuplot","w");
  fprintf(command,"set term x11\n");
  //fprintf(command,"set xr[-.2:%lf]\nset yr[%lf:%lf]\n",
  //0.2+1.0*(n_n-1), -0.2-0.5*n_n, 0.5+0.5*n_n);
  fflush(command);
  if (plot)
    fprintf(command,"plot '%s' u 1:2 w linesp,'%s' u 1:2 w p,'%s' u 1:2 w l\n",
	    plotfile,reqfile,transfile);
  else
    fprintf(command,"plot '%s' u 1:2 w linesp\n",plotfile);
  fflush(command);
  cout << " DONE,  press ENTER to quit " << endl;
  getchar();
  pclose(command);
}


#define nsw0_bs 20
#define nsw_bs 100
#define eps_bs 1.

double embed(int n, double **D, struct coord *x)
{ 
  int i, j, k, l;
  double T[n][n];
  int nne[n];
  int* ne[n];
  double * (d2[n]);

  for( i = 0; i < n; i++)
    nne[i] = 0;

  for( i = 0; i < n; i++)
    {
      T[i][i] = 0;
      for( j = 0; j < i; j++)
	{
	  T[i][j] = 0;
	  if( (D[i])[j] > 0. )
	    {
	      T[i][j] = 1;
	      nne[i]++;
	      nne[j]++;
	    }
	  T[j][i] = T[i][j];
	}
    }

  for( i = 0; i < n; i++)
    {
      ne[i] = new int[nne[i]];
      d2[i] = new double[nne[i]];
    }

  for( i = 0; i < n; i++)
    {
      k = 0;
      for( j = 0; j < n; j++)
	if( T[i][j] )
	  {
	    ne[i][k] = j;
	    d2[i][k] = (D[i])[j] * (D[i])[j];
	    k++;
	  }
    }

  //Randomize(); Randomize is done in Options.cc

  double xaux, yaux;

  xaux = yaux = 0;
  for( i = 0; i < n; i++)
    {
      xaux += x[i].x = Rnd();
      yaux += x[i].y = Rnd();
    }
  xaux /= n;
  yaux /= n;
  for( i = 0; i < n; i++)
    {
      x[i].x -= xaux;
      x[i].y -= yaux;
    }
  
  double coeff, saux2, aux, Fx, Fy, Gx, Gy, detM, Mxx, Mxy, Myy, detMtM, auxmin;
  int flag;

  for( int isw = 0; isw < nsw0_bs; isw++)
    {
      saux2 = 0;
      auxmin = 1.e10;
      for( i = 0; i < n; i++)
	{
	  flag = 0;
	  Fx = Fy = Mxx = Mxy = Myy = 0;
	  for( l = 0; l < nne[i]; l++)
	    {
	      j = ne[i][l];
	      xaux = x[j].x - x[i].x;
	      yaux = x[j].y - x[i].y;
	      aux = xaux * xaux + yaux * yaux;
	      if( aux < auxmin )
		auxmin = aux;
	      if( aux < 1.e-5 )
		{
		  flag = 1;
		  break;
		}
	      aux -= d2[i][l];
	      saux2 += aux * aux;
	      Fx += aux * xaux;
	      Fy += aux * yaux;
	      Mxx += aux + 2 * xaux * xaux;
	      Mxy += 2 * xaux * yaux;
	      Myy += aux + 2 * yaux * yaux;
	    }
	  if( flag )
	    {
	      cerr << "\7  *** r^2 = 0; kicking (" << auxmin<< ") ***" << endl;
	      x[i].x += (.4 * Rnd() - .2);
	      x[i].y += (.4 * Rnd() - .2);
	      continue;
	    }
	  detMtM = Mxx * Mxx + Myy * Myy + 2 * Mxy * Mxy;
	  if( detMtM == 0. )
	    {
	      cerr << "\7  *** detMtM = 0; left as was;  ***\n";
	      continue;
	    }
	  coeff = eps_bs / sqrt(detMtM);
	  x[i].x += coeff * Fx;
	  x[i].y += coeff * Fy;
	}
      if( saux2 < 1.e-20 )
	goto end;
    }

  for( int isw = 0; isw < nsw_bs; isw++)
    {
      saux2 = 0;
      auxmin = 1.e10;
      for( i = 0; i < n; i++)
	{
	  flag = 0;
	  Fx = Fy = Mxx = Mxy = Myy = 0;
	  for( l = 0; l < nne[i]; l++)
	    {
	      j = ne[i][l];
	      xaux = x[j].x - x[i].x;
	      yaux = x[j].y - x[i].y;
	      aux = xaux * xaux + yaux * yaux;
	      if( aux < auxmin )
		auxmin = aux;
	      if( aux < 1.e-10 )
		{
		  flag = 1;
		  break;
		}
	      aux -= d2[i][l];
	      saux2 += aux * aux;
	      Fx += aux * xaux;
	      Fy += aux * yaux;
	      Mxx += aux + 2 * xaux * xaux;
	      Mxy += 2 * xaux * yaux;
	      Myy += aux + 2 * yaux * yaux;
	    }
	  if( flag )
	    {
	      cerr << "\7  *** r^2 = 0; kicking ***" << endl;
	      x[i].x += Rnd();
	      x[i].y += Rnd();
	      continue;
	    }
	  detM = Mxx * Myy - Mxy * Mxy;
	  if( detM == 0. )
	    {
	      cerr << "\7  *** detM = 0; left as was;  ***\n";
	      continue;
	    }
	  Gx = Myy * Fx - Mxy * Fy;
	  Gy = -Mxy * Fx + Mxx * Fy;
	  coeff = eps_bs / detM;
	  x[i].x += coeff * Gx;
	  x[i].y += coeff * Gy;
	}
      if( saux2 < 1.e-20 )
	goto end;
    }

end:
  
  for( i = 0; i < n; i++)
    {
      delete ne[i];
      delete d2[i];
    }

  saux2 *= 2. / (n * (n - 1.)) ;
  saux2 = sqrt(saux2);

  return saux2;
}
