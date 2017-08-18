#include <stdio.h>
#include <string>
#include "itoa.hh"
#include <iostream>

using namespace std;
string & itoa(int i)
{
  char s[100];
  sprintf(s,"%d",i);
  
  string *S = new string(s);
  return *S;

std::cerr << "itoa: i = " << i << ", S = " << *S << std::endl;
}
