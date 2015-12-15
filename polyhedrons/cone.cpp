#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  simplicialPolyhedron A;
  A.read(cin).cone().print(cout);

  return 0;
}
