#include <iostream>
#include "classes/simplicialChainComplex.h"
using namespace std;

int main(int argc, char *argv[]) {

  simplicialPolyhedron P;
  P.read(cin);
  simplicialChainComplex S(P);
  /*cout << "Euler characteristic: " << S.eulerCharacteristic() << endl;
  cout << "Fundamental Class: " << endl;
  S.fundamentalClass().print_full(cout);
  cout << "and its boundary has " << S.boundary(S.dim(),S.fundamentalClass()).length() << " simplexes." << endl;
  cout << "Its support:" << endl;
  S.support(S.dim()-1,S.boundary(S.dim(),S.fundamentalClass())).print(cout);
  
  sparseMatrix L,D,U,R,Q;
  
  S.boundaryOperator(S.dim()).print(cout);
  
  sparseMatrix H_n = S.boundaryOperator(S.dim()).ker().transpose();
  
  for (int i = 0; i < H_n.size(1); i++) {
    cout << "The " << i+1 << "-th maximal class of homology has support: " << endl;
    S.support(S.dim(),H_n[i]).print(cout);
  } */

  return 0;
}
