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


  simplicialPolyhedron A,B;
  A.read(cin);

  if (itself)
  	B = A;
  else
  	B.read(cin);

  A << B;
  A.print(cout);
  
  
  return 0;
}
