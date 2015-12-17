#include <iostream>
#include "sparseMatrix.h"
using namespace std;

int main (int argc, char * argv[]) {
  
  sparseMatrix M,N;
  M.read(cin);
  N.read(cin);
  sparseMatrix K = M*N;
  K.print(cout);
  
  return 0;
}
