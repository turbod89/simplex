#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  simplicialPolyhedron A;
  A.read(cin).suspension().print(cout);

  return 0;
}
