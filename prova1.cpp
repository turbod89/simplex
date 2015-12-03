#include <iostream>
#include "classes/simplicialChainComplex.h"
using namespace std;

int main(int argc, char *argv[]) {

  simplicialPolyhedron P;
  P.read(cin);
  simplicialChainComplex S(P);

  S.boundaryOperator(S.dim()).print(cout);
  
  sparseMatrix L,D,U,R,Q;
  S.boundaryOperator(S.dim()).LDU_efficient(L,D,U,R,Q);
  L.print(cout);
  D.print(cout);
  U.print(cout);
  R.print(cout);
  Q.print(cout);

  return 0;
}
