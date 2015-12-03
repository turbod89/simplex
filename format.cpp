#include <iostream>
#include <cstdlib>
#include <string>
#include "classes/simplicialChainComplex.h"
using namespace std;

int main(int argc, char* argv[]) {
  
  string name_simplexes = "simplexes";
  string name_boundary_operator = "d";
  string name_fundamental_class = "Gamma";
  
  
  int n;
  cin >> n; //nombre de matrius de simplexs. Es a dir, dim + 1


  // simplexs  
  cout << "# name: " << name_simplexes << endl;
  cout << "# type: cell" << endl;
  cout << "# rows: 1" << endl;
  cout << "# columns: " << n << endl;
  
  for (int i = 0 ; i < n ; i++) {
    int r, c;
    cin >> r >> c;
    cout << "# name: <cell-element>" << endl;
    cout << "# type: matrix" << endl;
    cout << "# rows: " << r << endl;
    cout << "# columns: " << c << endl;
    
    int a;
    for (int i = 0 ; i < r; i++) {
      for (int j = 0; j < c; j++) {
        cin >> a;
        cout << " " << a;
      }
      cout << endl;
    }
    cout << endl << endl << endl;
  }
  
  // operadors bora
  cout << "# name: " << name_boundary_operator << endl;
  cout << "# type: cell" << endl;
  cout << "# rows: 1" << endl;
  cout << "# columns: " << n-1 << endl;
  
  for (int i = 0 ; i < n-1 ; i++) {
    sparseMatrix d;
    d.read(cin);
    
    cout << "# name: <cell-element>" << endl;
    cout << "# type: matrix" << endl;
    cout << "# rows: " << d.size(1) << endl;
    cout << "# columns: " << d.size(2) << endl;
    d.print_full(cout);

    cout << endl << endl << endl;
  }
  
  // classe fonamental
  
  cin >> n; // llargada de la classe fonamental
  cout << "# name: " << name_fundamental_class << endl;
  cout << "# type: matrix" << endl;
  cout << "# rows: " << n << endl;
  cout << "# columns: 1" << endl;
  
  int a;
  for (int i = 0; i < n; i++) {
    cin >> a;
    cout << " " << a << endl;
  }
  
  cout << endl << endl << endl;
    
  
  return 0;
}
