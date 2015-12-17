#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <vector>
#include "sparseMatrix.h"
using namespace std;


int main(int argc, char *argv[]) {
  
  bool reversed = false;
  
  for ( int i = 1; i < argc; i++ )
    if ( string(argv[i]) == "-r") {
      reversed = true;
    }

  if (! reversed) {
	  
    int r,c;
  
    cin >> r >> c;
    int * M = (int *) malloc(r*c*sizeof(int));
    int cnt = 0;
  
    for (int i = 0; i < r*c; i++)
      cin >> M[i];

    sparseMatrix spM(M,r,c);
    spM.print(cout);
    
  } else {
    
    sparseMatrix spM;
    
    spM.read(cin);
        
    spM.print_full(cout);
    
  }
  
  return 0;
}
