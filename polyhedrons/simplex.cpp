#include <iostream>
#include <cstdlib>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char * argv[]) {

  if (argc < 2) {
    cerr << "simplex d [-n n]" << endl;
    cerr << endl;
    cerr << "\t" << "Description: returns n simplexes of dimension d."<< endl;
    cerr << "\t" << "             By defalut, n = 1."<< endl;
    cerr << endl;
    return 1;
  }
  
  int n = 1, d = atoi(argv[1]);
  
  for ( int i = 2;i < argc; i++ )
    if ( string(argv[i]) == "-n") {
      i++;
      n = atoi(argv[i]);
    } else {
	  cerr << "Unknow argument \""<< string(argv[i]) <<"\"." << endl;
	}

  // creating simplex
  
  int * a = (int *) malloc((d+1)*n*sizeof(int));
  
  for (int i=0; i<(d+1)*n; i++)
    a[i] = i;
  
  simplicialPolyhedron A(d,n,a);
  A.print(cout);

  return 0;
}
