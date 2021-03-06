/*

  $Id: mroute_templates.cc,v 1.3 1997/08/25 17:04:11 jari Exp $

  */

// stl headers
#include <list>

#include "../mroute.hh"

// instantation of the templates we use
template class list<int>;
template class list<Path>;

// this is a part of the template list class but is needed when compiling with
// optimization -0
//template void construct(int *, int const &);

// instantation of additional operators used by the templates above
bool operator!=(Path::const_iterator const &,Path::const_iterator const &);
bool operator!=(Path::iterator const &,Path::iterator const &);
bool operator!=(ManyPaths::iterator const &,ManyPaths::iterator const &);
