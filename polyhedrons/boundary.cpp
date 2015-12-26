#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  simplicialPolyhedron A,B;
  A.read(cin);
  int * signs = (int *) malloc((A.dim()+1) * A.length()*sizeof(int));
  B = A.boundary(signs);

  for (int i = 0; i < B.length(); i++)
    cout << " " << signs[i];
  cout << endl;

  B.orientSimplexes(signs);
  B.print(cout);

  return 0;
}
