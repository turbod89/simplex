#include <iostream>
#include "simplicialChainComplex.h"
using namespace std;

int main(int argc, char *argv[]) {

  simplicialPolyhedron P,Q;
  P.read(cin);
  simplicialChainComplex S(P);
  
  cout << "# SIMPLICIAL CHAIN COMPLEX" << endl;
  S.print(cout);
  
  cout << "# EULER CHARACTERISTIC" << endl;
  cout << S.eulerCharacteristic() << endl;

  for (int i = S.dim(); i >= 0 ; i--) {
    sparseMatrix H = S.getHomology(i);
    cout << "# HOMOLOGY CLASSES OF GRADE " << i << endl;
    if (H.size(1) == 1 ) {
      H.print(cout);
    }
    else if (H.size(1) > 1) {
      for (int j = 0; j < H.size(1); j++)
        H[j].print(cout);
    }
  }

  for (int i = S.dim(); i >= 0 ; i--) {
    sparseMatrix H = S.getCohomology(i);
    cout << "# COHOMOLOGY CLASSES OF GRADE " << i << endl;
    if (H.size(1) == 1 ) {
      H.print(cout);
      cout << "PD " << endl;
      S.flat(i,H).print(cout);
    }
    else if (H.size(1) > 1) {
      for (int j = 0; j < H.size(1); j++) {
        H[j].print(cout);
        cout << "PD " << endl;
        S.flat(i,H[j]).print(cout);
      }
    }
  }

 /*
  sparseMatrix dd0 = S.boundaryOperator(1).transpose();
  sparseMatrix dd1 = S.boundaryOperator(2).transpose();

  sparseMatrix L0, D0, U0, P0, Q0;
  sparseMatrix L1, D1, U1, P1, Q1;

  dd0.LDU_efficient(L0,D0,U0,P0,Q0);
  dd1.LDU_efficient(L1,D1,U1,P1,Q1);

/* homology

  sparseMatrix Q1P0KerL0t = Q1*P0*(L0.transpose().ker());
//Q1P0KerL0t.print(cerr);
  sparseMatrix ImU1t = U1.transpose();
//ImU1t.print(cerr);
  sparseMatrix X = Q1.transpose()*ImU1t.LComplementary(Q1P0KerL0t);

  //P0KerL0t.print_octave(cout);
  //dd0.transpose().ker().print_octave(cout);
  for (int i = 0; i < X.size(2); i++)
      X.transpose()[i].print(cout);


*/

/* cohomology

  sparseMatrix P0tQ1tKerU1 = P0.transpose()*Q1.transpose()*(U1.ker());
//Q1P0KerL0t.print(cerr);
  sparseMatrix ImL0 = L0;
//ImU1t.print(cerr);
  sparseMatrix X = P0*ImL0.LComplementary(P0tQ1tKerU1);

  //P0KerL0t.print_octave(cout);
  //dd0.transpose().ker().print_octave(cout);
  for (int i = 0; i < X.size(2); i++)
      X.transpose()[i].print(cout);


*/


/*
  sparseMatrix Kerd1 = dd0.transpose().ker();
  sparseMatrix P0tKer = P0.transpose()*L0.transpose().ker();
  cout << "dd0 = ";
  dd0.print_octave(cout);
  cout << "L0 = ";
  L0.print_octave(cout);
  cout << "D0 = ";
  D0.print_octave(cout);
  cout << "U0 = ";
  U0.print_octave(cout);
  cout << "P0 = ";
  P0.print_octave(cout);
  cout << "Q0 = ";
  Q0.print_octave(cout);
  //Kerd1.print_octave(cout);
  //P0tKer.print_octave(cout);
*/

  /*
  cout << "Fundamental Class: " << endl;
  S.fundamentalClass().print_full(cout);
  cout << "and its boundary has " << S.boundary(S.dim(),S.fundamentalClass()).length() << " simplexes." << endl;
  cout << "Its support:" << endl;
  S.support(S.dim()-1,S.boundary(S.dim(),S.fundamentalClass())).print(cout);
  
  sparseMatrix L,D,U,R,Q;
  
  S.boundaryOperator(S.dim()).print(cout);
  
  sparseMatrix H_n = S.boundaryOperator(S.dim()).ker().transpose();
  
  for (int i = 0; i < H_n.size(1); i++) {
    cout << "The " << i+1 << "-th maximal class of homology has support: " << endl;
    S.support(S.dim(),H_n[i]).print(cout);
  } */

  return 0;
}
