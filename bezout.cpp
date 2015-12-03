#include <iostream>
#include <cstdlib>
using namespace std;

int gcd_extended(int n, int * v, int * c) {

  int indexMin = -1;
  int min = 0;
  int numNonZeros = 0;
  
  for (int i = 0; i < n ; i++)
    if (v[i] != 0)
      numNonZeros++;
  
  if (numNonZeros == 0)
    return 0;
  
  for (int i = 0; i < n ; i++)
    if (v[i] != 0 && (min == 0 || v[i] < min) ) {
      indexMin = i;
      min = v[i];
    }
  
  if (numNonZeros == 1) {
    
    for (int i = 0; i < n ; i++)
      if (i != indexMin)
        c[i] = 0;
      else
        c[i] = 1;
        
    return v[indexMin];
    
  } else if (min < 0) {
    // cas amb negatius
    
    int * signs = (int *) malloc(n*sizeof(int));
    int * v2 = (int *) malloc(n*sizeof(int));

    for (int i = 0; i < n ; i++)
      if (v[i] < 0) {
        v2[i] = -v[i];
        signs[i] = -1;
      } else if (v[i] > 0) {
        v2[i] = v[i];
        signs[i] = 1;
      } else {
        v2[i] = 0;
        signs[i] = 0;
      }
      
    int gcd = gcd_extended(n,v2,c);
    
    for (int i = 0; i < n ; i++)
      if ( v[i] != 0)
        c[i] *= signs[i];
    
    return gcd;

  } else {
    // cas general
    int * r = (int *) malloc(n*sizeof(int));

    for (int i = 0; i < n ; i++)
      if (v[i] == 0) {
        r[i] = 0;
      } else if (i != indexMin) {
        r[i] = v[i]%v[indexMin];
      } else {
        r[i] = v[i];
      }
  
    int gcd = gcd_extended(n,r,c);
 
    for (int i = 0; i < n ; i++)
      if (i != indexMin) {
        c[indexMin] -= c[i]*((v[i]-r[i])/v[indexMin]);
      }
      
    return gcd;  
  }
      
  return 0;

}

void print(ostream& out , int r, int c, int * A) {
  out << r << " " << c << endl;
  for (int i = 0 ; i < r ; i++) {
    for (int j = 0; j < c; j++) 
      out << " " << A[i*c+j];
    out << endl;
  }
}

void read (istream& in, int r, int c, int * A) {
  for (int i=0; i <r; i++)
    for (int j = 0; j < c ; j++)
      cin >> A[i*c+j];
}

void sumColums(int r, int c, int k, int * A, int * coeffs) {

  int * col = (int *) calloc(r,sizeof(int));
  
  for (int j = 0 ; j < c; j++)
    if (coeffs[j] != 0)
      for (int i = 0; i<r; i++)
        col[i] += coeffs[j]*A[i*c+j];
  
  for (int i = 0 ; i < r ; i++)
    A[i*c+k] = col[i];
  
}

void sumRows(int r, int c, int k, int * A, int * coeffs) {

  int * row = (int *) calloc(c,sizeof(int));
  
  for (int i = 0; i < r; i++)
    if (coeffs[i] != 0)
      for (int j = 0 ; j < c ; j++)
        row[j] += coeffs[i]*A[i*c+j];
        
  for (int j = 0; j < c; j++)
    A[k*c+j] = row[j];

}

void swapCols(int r, int c, int * A, int k, int l) {

  if (k == l)
    return;

  for (int i = 0; i < r; i++) {
    int a = A[i*c+k];
    A[i*c+k] = A[i*c+l];
    A[i*c+l] = a;
  }
}

void swapRows(int r, int c, int * A, int k, int l) {

  if (k == l)
    return;

  for (int j = 0; j < c; j++) {
    int a = A[k*c+j];
    A[k*c+j] = A[l*c+j];
    A[l*c+j] = a;
  }
}

void decomposition(int n, int * M, int * L, int * R) {
  
  for (int i = 0 ; i < n-1; i++) {
    // prepare values to range bezout
    int * v = (int *) malloc((n-i)*(n-i)*sizeof(int));
    int * coeffs = (int *) malloc((n-i)*(n-i)*sizeof(int));
    
    for (int j = 0 ; j < n-i; j++)
      for (int k = 0; k < n-i; k++)
        v[j*(n-i) + k] = M[(j+i)*n + (k+i)];
    
    // bezout
    int gcd = gcd_extended(n-i,v,coeffs);
    
    if (gcd == 0)
      return;
    
    // choose pivot
    int r = -1, c = -1;
    for (int j = 0 ; j < n-i; j++)
      for (int k = 0; k < n-i; k++)
        if ( coeffs[j*(n-i)+k] != 0 && (r<0 || c<0 || abs(coeffs[r*(n-i)+c]) > abs(coeffs[j*(n-i)+k]))) {
          r = j;
          c = k;
        }
    
    r += n-i;
    c += n-i;
    
    // appli bezout
    
    //
    swapCols(n,n,M,n-i,c);
    swapCols(n,n,R,n-i,c);
    swapRows(n,n,M,n-i,r);
    swapRows(n,n,L,n-i,r);
    
  }
}


int main(int argc, char *argv[]) {
  int n;
  cin >> n;
  int * v = (int *) malloc(n*sizeof(int));
  int * c = (int *) malloc(n*sizeof(int));
  for (int i = 0; i < n; i++)
    cin >> v[i];
  
  int gcd = gcd_extended(n,v,c);
  cout << gcd << endl;
  for (int i = 0; i < n ; i++)
    cout << " " << c[i];
  cout << endl;
}
