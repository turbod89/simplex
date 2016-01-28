#ifndef H_TOOLS
#define H_TOOLS

#include <iostream>
#include <cstdlib> 
#include <string.h>

namespace Tools {

  // 
  int gcd(int n, int * v, int * c = NULL);
  
  // block operations
  bool leq(int n, int const * const a, int const * const b);
  bool eq(int n, int const * const a, int const * const b );
  void swap(int n, int * const a, int * const b);
  int mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions = false);
}

#endif