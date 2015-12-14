#include <iostream>
#include "../classes/sparseMatrix.h"
using namespace std;

int main (int argc, char * argv[]) {
  
  sparseMatrix M;
  M.read(cin);
  sparseMatrix K = M.ker();
  K.print(cout);
  
  return 0;
}
