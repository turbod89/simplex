#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  bool itself = false;

  for ( int i = 1;i < argc; i++ )
    if ( string(argv[i]) == "--itself" || string(argv[i]) == "-d") {
  	  itself = true;
    } else {
  	  cerr << "Unknow argument \""<< string(argv[i]) <<"\"." << endl;
    }

  simplicialPolyhedron A,B,C;
  A.read(cin);
  if (itself)
  	B = A;
  else
  	B.read(cin);

  C = A*B;
  C.print(cout);

  return 0;
}
