#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include "lprop.h"
#include "graph.h"
extern "C" {
#include <unistd.h>
}
#ifndef HAVE_GETOPT_H
#include <getopt.h>
#endif

#if defined HAVE_GETOPT_H && defined HAVE_GETOPT_LONG
#include <getopt.h>
#else
extern "C" {
#include "getopt.h"
}
#endif

#define PREC_DEFAULT 8
#define EPS_DEFAULT 1.0e-9
#define HELP " [-m iter -e eps -p prec] -i matrix -l label -r output"
#define ITER_DEFAULT 10000

int main (int argc, char** argv)
{
  std::string input_matrix;
  std::string input_labels;
  std::string result;
  unsigned int max_iter = ITER_DEFAULT;
  unsigned int prec     = PREC_DEFAULT;
  double eps = EPS_DEFAULT;

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc,argv,"e:p:m:i:l:r:h")) != -1) {
    switch(opt) {
      case 'm': max_iter = atoi(optarg); break;
      case 'e': eps = atof(optarg); break;
      case 'p': prec = atoi(optarg); break;
      case 'i': input_matrix = std::string(optarg); break;
      case 'l': input_labels = std::string(optarg); break;
      case 'r': result = std::string(optarg); break;
      case 'h': std::cout << argv[0] << HELP << std::endl;
        return EXIT_FAILURE;
      default:
        std::cout << argv[0] << HELP << std::endl; return EXIT_FAILURE;
    }
  }

  try { // validation
    if (input_matrix.empty()) throw std::runtime_error ("no matrix input");
    if (input_labels.empty()) throw std::runtime_error ("no label input");
    if (result.empty()) throw std::runtime_error ("specify output file");
  } catch (const std::exception& e) {
    std::cerr << "ERROR: main() says: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  try { // execute label propagation
    ssl_lprop::LP lprop(eps);
    lprop.read(input_matrix, input_labels);
    lprop.train(max_iter);
    lprop.write(argv[argc-1], prec);
    //    lprop.show(prec);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
