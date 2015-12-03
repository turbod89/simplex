#include <iostream>
#include <time.h>
#include "classes/simplicialChainComplex.h"
using namespace std;

int main(int argc, char *argv[]) {

  sparseMatrix M,L1,D1,U1,P1,Q1;
  sparseMatrix L2,D2,U2,P2,Q2;
  
  M.read(cin);
  
  cout << "Matrix read succefully" << endl;
  
  clock_t tStart1 = clock();
  M.LDU(L1,D1,U1,P1,Q1);
  clock_t tFinish1 = clock();
  cout << "Full method: " << (double)(tFinish1 - tStart1)/CLOCKS_PER_SEC << "s" << endl;

  clock_t tStart2 = clock();
  M.LDU_efficient(L2,D2,U2,P2,Q2);
  clock_t tFinish2 = clock();
  cout << "Efficient method: " << (double)(tFinish2 - tStart2)/CLOCKS_PER_SEC << "s" << endl;
  
  

  return 0;
}
