#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <boost/foreach.hpp>
#define REP(i,a) for(int i=0;i<(a);i++)
#define foreach BOOST_FOREACH
#define pb push_back

namespace utils {
  int stoi (const std::string str);
  double stod( const std::string str);
}

#endif
