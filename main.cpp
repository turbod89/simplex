#include <iostream>
#include <ctime>
#include "simplicialChainComplex.h"
using namespace std;

float timeMark(clock_t & c) {
  float t = 1.0 * (clock() - c) / CLOCKS_PER_SEC;
  c = clock();
  return t;
}


int main(int argc, char *argv[]) {

  clock_t c_start = clock();

  simplicialPolyhedron P,Q;
  P.read(cin);
cerr << "Read polyhedron: " << timeMark(c_start) << "s" << endl;
  simplicialChainComplex S(P);
cerr << "Read and inflate simplicial Chain: " << timeMark(c_start) << "s" << endl;
  cout << "# SIMPLICIAL CHAIN COMPLEX" << endl;
  S.print(cout);
  
  cout << "# EULER CHARACTERISTIC" << endl;
  cout << S.eulerCharacteristic() << endl;

  sparseMatrix  H[S.dim()+1];

  for (int i = S.dim(); i >= 0 ; i--) {
    H[i] = S.getHomology(i);
cerr << "Find homology classes of grade " << i << ": " << timeMark(c_start) << "s" << endl;
    cout << "# HOMOLOGY CLASSES OF GRADE " << i << endl;
    if (H[i].size(1) == 1 ) {
      H[i].print(cout);
    } else if (H[i].size(1) > 1) {
      for (int j = 0; j < H[i].size(1); j++) {
        H[i][j].print(cout);
      }
    }
  }

  sparseMatrix cH[S.dim()+1];

  for (int i = S.dim(); i >= 0 ; i--) {
    cH[i] = S.getCohomology(i);
cerr << "Find cohomology classes of grade " << i << ": " << timeMark(c_start) << "s" << endl;
    cout << "# COHOMOLOGY CLASSES OF GRADE " << i << endl;
    if (cH[i].size(1) == 1 ) {
      cH[i].print(cout);
    } else if (cH[i].size(1) > 1) {
      for (int j = 0; j < cH[i].size(1); j++) {
        cH[i][j].print(cout);
      }
    }
  }

  cout << endl << endl << endl;

  for (int i = S.dim(); i > 0 ; i--) {
    cout << "# PD OF COHOMOLOGY CLASSES OF GRADE " << i << endl;
    sparseMatrix L,D,U,Q,P;
    H[S.dim()-i].transpose().LDU_efficient(L,D,U,P,Q);

    for (int j = 0; j < cH[i].size(1) ; j++) {
      sparseMatrix C;
      if (cH[i].size(1) == 1)
        C = cH[i];
      else
        C = cH[i][j];

      sparseMatrix PD_C = S.flat(i,C);


      PD_C = S.d_ldu[S.dim()-i].P * S.d_ldu[S.dim()-i].L.LComplementary(S.d_ldu[S.dim()-i].P.transpose()*PD_C.transpose());

      sparseMatrix X = L.LXeqY(P.transpose()*PD_C);
      X.print_octave(cout);
    }

  }



  //delete[] H;
  //delete[] cH;

  return 0;
}
