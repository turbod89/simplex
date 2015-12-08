#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <fstream>
#include <vector>
using namespace std;

void print_matrix(int * const M, int rows, int cols) {
  for (int i = 0 ; i < rows; i++) {
    for (int j = 0; j < cols; j++)
      cout << " " << M[i*cols+j];
    cout << endl;
  }
}

void random_sparse_matrix_generator(int * &M, int rows, int cols, int k, int maximun, int minimun) {
  int n = k*ceil(sqrt(cols*rows));
  M = (int *) calloc(cols*rows, sizeof(int));
  for (int i = 0; i < n; i++) {
	int r = rand()%rows, c = rand()%cols, v = rand()%(maximun-minimun+1) + minimun;	
    M[cols*r + c] = v;
  }
}


int main(int argc, char *argv[]) {

  if (argc < 3) {
    cout << "spmatrix r c [-c constant] [-d] [-M maximun] [-m minimun]" << endl;
    cout << endl;
    cout << "\t" << "Description: returns a matrix with O(sqrt(r*c)) non-null elements"<< endl;
    cout << "\t" << "             of r rows and c columns."<< endl;
    cout << endl;
    cout << "\t" << "-c matrix has almost constant*sqrt(r*c). By default constant = 1. "<< endl;
    cout << "\t" << "-d outs matrix dimensions before out matrix. "<< endl;
    return 1;
  }
  
  srand(time(NULL));
  
  int * M;
  int r = atoi(argv[1]), c = atoi(argv[2]), k = 1, maximunValue = 100, minimunValue = -100;
  
  for ( int i = 3;i < argc; i++ )
    if ( string(argv[i]) == "-c") {
      i++;
      k = atoi(argv[i]);
    } else if ( string(argv[i]) == "-M") {
      i++;
      maximunValue = atoi(argv[i]);
    } else if( string(argv[i]) == "-m") {
      i++;
      minimunValue = atoi(argv[i]);
    }
  
  random_sparse_matrix_generator(M,r,c,k,maximunValue, minimunValue);
 
  cout << r << " " << c << endl; 
  
  print_matrix(M,r,c);
  
  return 0;
}
