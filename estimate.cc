#include <math.h>
#include <iostream.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int i, nn, nl, nn1;
  double p, q, np, np1, np0, term0, fac0, fac, fac1;
  double npalt;

  if( argc != 3 )
    { cerr << "Usage: " << argv[0] << " nn nl\n"; exit(-1); }

  nn = atoi(argv[1]);
  nl = atoi(argv[2]);

  if( nn < 3 || nl < nn )
    { cerr << "Usage: " << argv[0] << " nn nl\n"; exit(-1); }
  
  cout << "nn = " << nn << ", nl = " << nl << endl;

  nn1 = nn - 1;
  q = 2 * nl / (double) nn;
  p = q / nn1;

  cout << "q = " << q << ", p = " << p << endl;
  
  np = 1. / nn1;
  np1 = sqrt(2. * M_PI / nn1);
  np0 = 0;
  fac0 = p;
  term0 = 1./nn1;
  fac =  q * exp(1. / q) / nn1;
  fac1 =  q * exp(1. / q - 1);
  for( i = 1; i <= nn1; i++)
    {
      np *= fac * i;
      np1 *= fac1;
      term0 *= fac0 * (nn - i);
      np0 += term0;
    }
  cout << "est <np> = " << np0 << " approx " << np << " approx " << np1 << endl;
  cout << "est S = " << log(np0) << " approx " << log(np) << " approx " << log(np1) << endl;


  npalt = 1;
  for( i = 0; i < nl - nn + 1; i++)
    npalt *= 2;

  cout << "alt est <np> <= " << npalt << endl;
  cout << "alt S <= " << log(npalt) << endl;
}
