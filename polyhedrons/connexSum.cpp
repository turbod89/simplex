#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  simplicialPolyhedron A,B;
  A.read(cin);
  B.read(cin);
  A << B;
  A.print(cout);
  
  
  return 0;
}
