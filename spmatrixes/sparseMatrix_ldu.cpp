#include <iostream>
#include "sparseMatrix.h"
using namespace std;

int main (int argc, char * argv[]) {
  
  sparseMatrix M;
  M.read(cin);
  
  sparseMatrix L, D, U, P, Q;
  
  M.LDU_efficient(L,D,U,P,Q);
  L.print(cout);
  D.print(cout);
  U.print(cout);
  P.print(cout);
  Q.print(cout);
    
  return 0;
}
