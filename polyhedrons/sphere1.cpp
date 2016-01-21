#include <iostream>
#include <cstdlib>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  if (argc < 1) {
    cerr << "S1 [n]" << endl;
    cerr << endl;
    cerr << "\t" << "Description: returns a S^1 trianguled with n simplexes."<< endl;
    cerr << "\t" << "             Of curse, n >= 3. By default, n = 3."<< endl;
    cerr << endl;
    return 1;
  }
  
  int n = 3;

  if (argc > 1)
    n = max(3,atoi(argv[1]));
  
  for (int i = 2;i < argc; i++ ) {
	  cerr << "Unknow argument \""<< string(argv[i]) <<"\"." << endl;
	}

  // creating S^1_n
  
  int * a = (int *) malloc(2*n*sizeof(int));
  for (int i = 0; i < n; i++) {
    a[2*i] = i;
    a[2*i+1] = (i+1)%n;
  }
  
  simplicialPolyhedron A(1,n,a);
  A.print(cout);

  return 0;
}
