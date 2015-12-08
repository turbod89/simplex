#include <iostream>
#include "classes/sparseMatrix.h"
using namespace std;

int main (int argc, char * argv[]) {
  
  sparseMatrix M,N;
  
  M.read(cin);
  N.read(cin);
  cout << "Succefully readed two matrixes of dimensions " << M.size(1) << "x" << M.size(2) << " and " << N.size(1) << "x" << N.size(2) << endl;

  cout << "First matrix has " << M.length() << " elements. Showing it:" << endl;
  M.print_octave(cout);
  cout << "List of its values: " << endl;
  int * a = (int *) malloc(M.size(1)*M.size(2)*sizeof(int));
  M.getFullValues(a);
  for (int i = 0; i < M.size(1)*M.size(2); i++)
    cout << " " << a[i];
  cout << endl;
  
  cout << "Second matrix has " << N.length() << " elements. Showing it:" << endl;
  N.print_octave(cout);
  cout << "List of its values: " << endl;
  int * b = (int *) malloc(N.size(1)*N.size(2)*sizeof(int));
  N.getFullValues(b);
  for (int i = 0; i < N.size(1)*N.size(2); i++)
    cout << " " << b[i];
  cout << endl;
  
  if ( M.size(2) == N.size(1)) {
    cout << "Their product" << endl;
    (M*N).print_octave(cout);
  }
  
  return 0;
}
