#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  simplicialPolyhedron A,B,C;
  A.read(cin);
  B.read(cin);
  C = A*B;
  C.print(cout);

  return 0;
}
