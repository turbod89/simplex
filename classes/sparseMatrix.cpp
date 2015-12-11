#include "sparseMatrix.h"

sparseMatrix::sparseMatrix(){
  this->values = NULL;
  this->cols = NULL;
  this->rows = NULL;
}

sparseMatrix::sparseMatrix(int r, int c, const int * const rows, const int * const cols, const int * const values) {
  this->numRows = r;
  this->numCols = c;
  this->rows = (int *) malloc((r+1)*sizeof(int));
  memcpy(this->rows,rows,(r+1)*sizeof(int));

  this->values = (int *) malloc(this->length()*sizeof(int));
  this->cols = (int *) malloc(this->length()*sizeof(int));
  memcpy(this->values,values,this->length()*sizeof(int));
  memcpy(this->cols,cols,this->length()*sizeof(int));
}

sparseMatrix::sparseMatrix(int r , int c) {
  this->numRows = r;
  this->numCols = c;
  this->values = NULL;
  this->cols = NULL;
  this->rows = (int *) malloc((r+1)*sizeof(int));
  
  for (int i = 0 ; i < this->numRows + 1; i++)
    this->rows[i] = 0;

}

sparseMatrix::sparseMatrix(int * M, int r , int c) {

  int cnt = 0;
  for (int i = 0; i < r*c; i++)
    if (M[i] != 0)
      cnt++;

    this->numCols = c;
    this->numRows = r;
    this->values = (int *) malloc(cnt*sizeof(int));
    this->cols = (int *) malloc(cnt*sizeof(int));
    this->rows = (int *) malloc((r+1)*sizeof(int));

    this->rows[0] = 0;
    for (int i = 0; i < r; i++) {
      this->rows[i+1] = this->rows[i];
      for (int j = 0; j < c; j++) {
        if (M[i*c + j] != 0) {
          this->cols[this->rows[i+1]] = j;
          this->values[this->rows[i+1]] = M[i*c + j];
          this->rows[i+1]++;
        }
      }
    }

}

sparseMatrix::sparseMatrix(const sparseMatrix& M){
  *this = M;
}

sparseMatrix::~sparseMatrix(){
  free(this->values);
  free(this->rows);
  free(this->cols);
}

sparseMatrix& sparseMatrix::operator=(const sparseMatrix& M) {
  this->numCols = M.numCols;
  this->numRows = M.numRows;

  this->values = (int *) malloc(M.length()*sizeof(int));
  this->cols = (int *) malloc(M.length()*sizeof(int));
  this->rows = (int *) malloc((this->numRows+1)*sizeof(int));
  memcpy(this->values,M.values,M.length()*sizeof(int));
  memcpy(this->cols,M.cols,M.length()*sizeof(int));
  memcpy(this->rows,M.rows,(this->numRows+1)*sizeof(int));
  
  return *this;
}

const int * sparseMatrix::getValues() const {
  return this->values;
}

int * sparseMatrix::getValues() {
  return this->values;
}

const int * sparseMatrix::getRows() const {
  return this->rows;
}

int * sparseMatrix::getRows() {
  return this->rows;
}

const int * sparseMatrix::getCols() const {
  return this->cols;
}

int * sparseMatrix::getCols() {
  return this->cols;
}

int sparseMatrix::size(int i) const {
  if (i == 1)
    return this->numRows;
  else if (i == 2)
    return this->numCols;
  
  return 0;
}

inline int sparseMatrix::length() const {
  return this->rows[this->numRows];
}

const sparseMatrix& sparseMatrix::decompose(int ** values, int ** cols, int ** rows) const {
  *values = (int *) malloc(this->length()*sizeof(int));
  *cols = (int *) malloc(this->length()*sizeof(int));
  *rows = (int *) malloc((this->numRows+1)*sizeof(int));
  memcpy(*values,this->values,this->length()*sizeof(int));
  memcpy(*cols,this->cols,this->length()*sizeof(int));
  memcpy(*rows,this->rows,(this->numRows+1)*sizeof(int));
  return *this;
}
sparseMatrix& sparseMatrix::read_index_format(int numRows, int numCols, int length, int * rows, int * cols, int * values, bool sorted) {
  
  if ( !sorted) {
    // sorting
    int ** a = (int **) malloc(2*sizeof(int *));
    a[0] = rows;
    a[1] = values;
    this->mergeSort(length, cols, 2, a);
    a[0] = cols;
    this->mergeSort(length, rows, 2, a);
  }
  
  int * crows = (int *) malloc((numRows + 1)*sizeof(int));
  
  crows[0] = 0;
  for (int i = 0, v = 0; i < numRows && v < length; v++) {
      while ( i < rows[v]) {
        crows[i+2] = crows[i+1];
        i++;
      }
      crows[i+1]++;
  }
  
  *this = sparseMatrix(numRows,numCols,crows,cols,values);
  
  return *this;
}
const sparseMatrix& sparseMatrix::index_format(int * rows, int * cols, int * values) const {

  if (rows != NULL)
    for (int i = 0, r = 0; i < this->length(); i++) {
      while (this->rows[r+1] <= i)
        r++;
      rows[i] = r;
    }

  if (cols != NULL)  
    memcpy(values,this->values,this->length()*sizeof(int));
  if (values != NULL)
    memcpy(cols,this->cols,this->length()*sizeof(int));
    
  return *this;
}


sparseMatrix& sparseMatrix::read(istream& in) {

  in >> this->numRows >> this->numCols;
  
  this->rows = (int *) malloc((this->numRows + 1)*sizeof(int));

  for (int i = 0; i < this->numRows + 1; i++)
    in >> this->rows[i];
    
  this->values = (int *) malloc(this->length()*sizeof(int));
  this->cols = (int *) malloc(this->length()*sizeof(int));


  for (int i = 0; i < this->length(); i++)
    in >> this->cols[i];
  for (int i = 0; i < this->length(); i++)
    in >> this->values[i];

  return *this;
  
}

const sparseMatrix& sparseMatrix::print(ostream& out) const {

  out << this->numRows << " " << this->numCols << endl;

  for (int i = 0; i < this->numRows + 1; i++)
    out << " " << this->rows[i];
  out << endl;

  for (int i = 0; i < this->length(); i++)
    out << " " << this->cols[i];
  out << endl;

  for (int i = 0; i < this->length(); i++)
    out << " " << this->values[i];
  out << endl;

  return *this;
}

const sparseMatrix& sparseMatrix::print_full(ostream& out) const {

  out << this->numRows << " " << this->numCols << endl;

  for (int i = 0, nextValue = 0; i < this->numRows; i++) {
    for ( int j = 0; j < this->numCols ; j++)
      if (nextValue < this->rows[i+1] && this->cols[nextValue] == j) {
        out << " " << this->values[nextValue];
        nextValue++;
      } else {
        out << " 0";
      }
    out << endl;
  }

  return *this;
}

const sparseMatrix& sparseMatrix::print_octave(ostream& out) const {

  out  << "[";

  for (int i = 0, nextValue = 0; i < this->numRows; i++) {
    for ( int j = 0; j < this->numCols ; j++)
      if (nextValue < this->rows[i+1] && this->cols[nextValue] == j) {
        out << " " << this->values[nextValue];
        nextValue++;
      } else {
        out << " 0";
      }
    if (i < this->numRows -1)
      out << ";" << endl;
  }
  
  out << "]" << endl;

  return *this;
}

const sparseMatrix& sparseMatrix::getFullValues(int * a) const {
  
  ////////////////
  //
  // write on the array the values of the full matrix
  //
  ////////////////
  
  for (int i = 0, nextValue = 0; i < this->numRows; i++) {
    for ( int j = 0; j < this->numCols ; j++)
      if (nextValue < this->rows[i+1] && this->cols[nextValue] == j) {
        a[i*this->numCols + j] = this->values[nextValue];
        nextValue++;
      } else {
        a[i*this->numCols + j] = 0;
      }
  }
  
  return (*this);
}

void sparseMatrix::swap(int &a, int &b) const {
  int c = a;
  a = b;
  b = c;
}

int sparseMatrix::binary_search_position(int length, const int * const v, int needle) const {
  
  ////////////////////
  //
  // binary search, return the position where element is placed if it
  // exists or the position where it SHOULD be placed
  //
  ////////////////////
  
  if (length <= 0)
    return 0;
  
  int i = length/2;
  if ( v[i] == needle)
    return i;
  else if ( v[i] > needle)
    return binary_search_position(i,v, needle);
  else if ( v[i] < needle)
    return i + 1 + binary_search_position(length - i - 1,&v[i+1], needle); 
}

int sparseMatrix::binary_search(int length, const int * const v, int needle) const {
  
  ////////////////////
  //
  // binary search, return the position where element is placed if it
  // exists or -1 otherwhise
  //
  ////////////////////
  
  
  int i = this->binary_search_position(length,v, needle);
  if (i == length)
    return -1;
  if (v[i] == needle)
    return i;
  return -1;
}

void sparseMatrix::mergeSort(int length, int * v, int auxLength, int * const * const aux) const {
  
  ////////////////////
  //
  // merge sort, sort v.
  //
  // Moreover, it also sorts auxLength number of arrays in same order
  // it did it with v.
  //
  ////////////////////
  
  if (length <= 1)
    return;

  int l1 = (int) (length/2),
      l2 = length - l1;
  int * a = (int *) malloc(l1*sizeof(int)),
      * b = (int *) malloc(l2*sizeof(int));
  memcpy(a,v,l1*sizeof(int));
  memcpy(b,&v[l1],l2*sizeof(int));

  int ** aux1 = NULL, **aux2 = NULL;

  if (auxLength > 0 && aux != NULL) {
    aux1 = (int **) malloc(auxLength*sizeof(int *));
    aux2 = (int **) malloc(auxLength*sizeof(int *));
    for ( int i = 0; i < auxLength; i++) {
      aux1[i] = (int *) malloc(l1 * sizeof(int));
      aux2[i] = (int *) malloc(l2 * sizeof(int));
      memcpy(aux1[i],aux[i],l1 * sizeof(int));
      memcpy(aux2[i],&(aux[i])[l1],l2 * sizeof(int));
    }
  }

  mergeSort(l1, a, auxLength, aux1);
  mergeSort(l2, b, auxLength, aux2);

  int i = 0, j = 0;
  while( i < l1 || j < l2)
    if (    (i < l1 && j < l2 && b[j] < a[i])
         || (j < l2 && i == l1) ) {
      v[i+j] = b[j];
      j++;
    } else {
      v[i+j] = a[i];
      i++;
    }
    
  if (auxLength > 0 && aux != NULL) {
    i = 0;
    j = 0;
    
    while( i < l1 || j < l2)
      if (    (i < l1 && j < l2 && b[j] < a[i])
         || (j < l2 && i == l1) ) {
        for (int k = 0; k <auxLength; k++)  
          (aux[k])[i+j] = (aux2[k])[j];
        j++;
      } else {
        for (int k = 0; k <auxLength; k++)  
          (aux[k])[i+j] = (aux1[k])[i];
        i++;
      }
  }
}

int sparseMatrix::gcd(int n, int * v, int * c) const {

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
    
    if ( c != NULL)
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
      
    int g = this->gcd(n,v2,c);
    
    if ( c != NULL)
      for (int i = 0; i < n ; i++)
        if ( v[i] != 0)
          c[i] *= signs[i];
    
    return g;

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
  
    int g = this->gcd(n,r,c);
 
    if ( c != NULL)
      for (int i = 0; i < n ; i++)
        if (i != indexMin) {
          c[indexMin] -= c[i]*((v[i]-r[i])/v[indexMin]);
        }
      
    return g;  
  }
      
  return 0;

}

int sparseMatrix::LDU_full(int n, int m, int * const M, int * const L, int * const D, int * const U, int * const rowPerm, int * const colPerm) const {
  
  // Nota: La matriu L esta transposta, per conveniencia del calcul
  
  // inicialitzem
  int acum = 1;
  int RANK_MAX = min(n,m);
  bool M_is_zero = true;
  int k = 0;
  
  for (int i = 0; i < n*RANK_MAX; i++)
    L[i] = 0;
    
  for (int i = 0; i < m*RANK_MAX; i++)
    U[i] = 0;
  
  for (int i = 0; i < RANK_MAX; i++)
    D[i] = 0;
  
  for (int i = 0; i < m; i++)
    colPerm[i] = i;
  
  for (int i = 0; i < n; i++)
    rowPerm[i] = i;
  
  // comprovem que no es la matriu 0
  int c = k;
  int r = k;
  
    for (int i = k; i < n && M_is_zero; i++)
      for (int j = k; j < m && M_is_zero; j++)
        if ( M[i*m + j] != 0) {
          M_is_zero = false;
          r = i;
          c = j;
        }
  
  while (k < RANK_MAX && !M_is_zero) {

    // permutem files i columnes
    int aux = rowPerm[k];
    rowPerm[k] = rowPerm[r];
    rowPerm[r] = aux;
    
    aux = colPerm[k];
    colPerm[k] = colPerm[c];
    colPerm[c] = aux;
    
    for (int j = 0; j < RANK_MAX ; j++) {
      aux = L[j*n + k];  // recordem: L esta trasposta
      L[j*n + k] = L[j*n + r];
      L[j*n + r] = aux;
    }

    for (int i= 0; i < RANK_MAX ; i++) {
      aux = U[i*m + k];
      U[i*m + k] = U[i*m + c];
      U[i*m + c] = aux;
    }

    for (int j = 0 ; j < m ; j++) {
      aux = M[k*m+j];
      M[k*m+j] = M[r*m+j];
      M[r*m+j] = aux;
    }
    
    for (int i = 0 ; i < n ; i++) {
      aux = M[i*m+k];
      M[i*m+k] = M[i*m+c];
      M[i*m+c] = aux;
    }
    
    // algoritme
    for (int j = k; j < m; j++)
      U[k*m + j] = M[k*m + j];
      
    for (int i = k; i < n; i++)
      L[k*n + i] = M[i*m + k]; // L esta transposta

    int g = this->gcd(n - k, &L[k*n + k]);
    int d = M[k*m + k]/g;

    for (int i = k; i < n; i++)
      L[k*n + i] = M[i*m + k]/g; // L esta transposta
    
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) {
        M[i*m + j] = d*M[i*m + j] - L[k*n + i]*U[k*m+j]; // L esta trasposta
      }
    
    acum *= d;
    D[k] = acum;


    // seguent iteracio
    k++;
    
    // comprovem si la matriu es 0
    // i donem les coordenades del
    // primer element no nul
   
   
    M_is_zero = true;
    for (int i = k; i < n && M_is_zero; i++)
      for (int j = k; j < m && M_is_zero; j++)
        if ( M[i*m + j] != 0) {
          M_is_zero = false;
          r = i;
          c = j;
        }
        
  }
  
  
  return k;
}

int sparseMatrix::minAbs(int n, int const * const v) const {
    int j = -1;
    
    for (int i = 0; i < n && (j < 0 || v[j] != 0) ; i++)
        if (j < 0 || abs(v[i]) < abs(v[j]))
            j = i;
    
    return j;
}

int sparseMatrix::maxAbs(int n, int const * const v) const {
    int j = -1;
    
    for (int i = 0; i < n ; i++)
        if (j < 0 || abs(v[i]) > abs(v[j]))
            j = i;
    
    return j;
}


int sparseMatrix::minAbsNotZero(int n, int const * const v) const {
    int j = -1;
    
    for (int i = 0; i < n ; i++)
        if (v[i] != 0 && (j < 0 || abs(v[i]) < abs(v[j])))
            j = i;
    
    return j;
}


void sparseMatrix::multiply(int a, int n, int * const v) const {
    for (int i = 0 ; i < n ; i++)
    v[i] *= a;
}

int sparseMatrix::multiply(int n, int const * const v, int const * const w) const {
    int s = 0;
    for (int i = 0; i < n ; i++)
        s += v[i]*w[i];
    return s;
}

int sparseMatrix::multiplyRows(int n1, int const * const c1, int const * const v1, int n2, int const * const c2, int const * const v2) const {

  int s = 0;
  int k = 0, l = 0;
  while (k < n1 && l < n2)
    if (c1[k] == c2[l]) {
      s += v1[k]*v2[l];
      k++;
      l++;
    } else if (c1[k] > c2[l]) {
      l++;
    } else if (c1[k] < c2[l]) {
      k++;
    }
    return s;
}

int sparseMatrix::sumRows(int n1, int const * const c1, int const * const v1, int n2, int const * const c2, int const * const v2, int * const c3, int * const v3) const {
  
  int numValues = 0;
  for (int i = 0, j = 0 ; i < n1 || j < n2;)
    if ( i >= n1 || (j < n2 && c1[i] > c2[j]) ) {
      v3[numValues] = v2[j];
      c3[numValues] = c2[j];
      j++;
      numValues++;
    } else if ( j >= n2  || (i < n1 && c1[i] < c2[j]) ) {
      v3[numValues] = v1[i];
      c3[numValues] = c1[i];
      i++;
      numValues++;
    } else if ( c1[i] == c2[j] && v1[i] + v2[j] != 0) {
      v3[numValues] = v1[i] + v2[j];
      c3[numValues] = c1[i];
      i++;
      j++;
      numValues++;
    } else if ( c1[i] == c2[j] && v1[i] + v2[j] == 0) {
      i++;
      j++;
    }

  return numValues;
}


inline int sparseMatrix::numValuesInRow(int row) const {

  ////////////////
  //
  // return the number of non null elements in a row
  //
  ////////////////
  
  return this->rows[row+1] - this->rows[row];

}

sparseMatrix& sparseMatrix::swapRows(int swap1, int swap2) {
  
  if (swap1 > swap2)
    return this->swapRows(swap2,swap1);
  if (swap1 == swap2)
    return *this;
  if (swap2 >= this->numRows || swap1 < 0)
    return *this;
    
  int numValuesRow1 = this->numValuesInRow(swap1),
      numValuesRow2 = this->numValuesInRow(swap2);
      
  int row1Starts = this->rows[swap1],
      row2Starts = this->rows[swap2];
      
  if (numValuesRow1 > numValuesRow2) {
  
    int * valuesRow1 = (int *) malloc(numValuesRow1*sizeof(int)),
        * colsRow1 = (int *) malloc(numValuesRow1*sizeof(int));
    memcpy(valuesRow1,&this->values[row1Starts],numValuesRow1*sizeof(int));
    memcpy(colsRow1,&this->cols[row1Starts],numValuesRow1*sizeof(int));

    // cpy row2 on row1
    for (int i = row1Starts; i < row1Starts + numValuesRow2 ; i++){
      this->values[i] = this->values[i - row1Starts + row2Starts];
      this->cols[i] = this->cols[i - row1Starts + row2Starts];    
    }
    
    // move inner rows
    for (int i = row1Starts + numValuesRow2; i < row2Starts - numValuesRow1 + numValuesRow2; i++) {
      this->values[i] = this->values[i + numValuesRow1 - numValuesRow2];
      this->cols[i] = this->cols[i + numValuesRow1 - numValuesRow2];
    }

    // copy row1 on row2
    for (int i = row2Starts - numValuesRow1 + numValuesRow2; i < row2Starts+ numValuesRow2; i++) {
      this->values[i] = valuesRow1[i - row2Starts - numValuesRow2 + numValuesRow1];
      this->cols[i] = colsRow1[i - row2Starts - numValuesRow2 + numValuesRow1];
    }
    
  } else if (numValuesRow1 < numValuesRow2) {
  
    int * valuesRow2 = (int *) malloc(numValuesRow2*sizeof(int)),
        * colsRow2 = (int *) malloc(numValuesRow2*sizeof(int));
    memcpy(valuesRow2,&this->values[row2Starts],numValuesRow2*sizeof(int));
    memcpy(colsRow2,&this->cols[row2Starts],numValuesRow2*sizeof(int));

    // cpy row1 on row2
    for (int i = row2Starts + numValuesRow2 - numValuesRow1; i < row2Starts + numValuesRow2 ; i++){
      this->values[i] = this->values[i - row2Starts - numValuesRow2 + numValuesRow1 + row1Starts];
      this->cols[i] = this->cols[i - row2Starts - numValuesRow2 + numValuesRow1 + row1Starts];    
    }
    
    // move inner rows
    for (int i = row2Starts + numValuesRow2 - numValuesRow1 - 1; i >= row1Starts + numValuesRow2 ; i--) {
      this->values[i] = this->values[i + numValuesRow1 - numValuesRow2];
      this->cols[i] = this->cols[i + numValuesRow1 - numValuesRow2];
    }

    // copy row2 on row1
    for (int i = row1Starts; i < row1Starts+ numValuesRow2; i++) {
      this->values[i] = valuesRow2[i - row1Starts];
      this->cols[i] = colsRow2[i - row1Starts];
    }
  } else {
    int * valuesRow1 = (int *) malloc(numValuesRow1*sizeof(int)),
        * colsRow1 = (int *) malloc(numValuesRow1*sizeof(int));
    memcpy(valuesRow1,&this->values[row1Starts],numValuesRow1*sizeof(int));
    memcpy(colsRow1,&this->cols[row1Starts],numValuesRow1*sizeof(int));

    // cpy row2 on row1
    for (int i = row1Starts; i < row1Starts + numValuesRow2 ; i++){
      this->values[i] = this->values[i - row1Starts + row2Starts];
      this->cols[i] = this->cols[i - row1Starts + row2Starts];    
    }

    // copy row1 on row2
    for (int i = row2Starts ; i < row2Starts+ numValuesRow1; i++) {
      this->values[i] = valuesRow1[i - row2Starts ];
      this->cols[i] = colsRow1[i - row2Starts ];
    }
    
  }
  
  // update rows

  for (int i = swap1 + 1 ; i <= swap2 ; i++)
    this->rows[i] += numValuesRow2 - numValuesRow1; 

  return *this;

}

sparseMatrix& sparseMatrix::swapCols(int col1, int col2) {

  if (col1 > col2)
    return this->swapCols(col2,col1);
  if (col1 == col2)
    return *this;
  if (col2 >= this->numCols || col1 < 0)
    return *this;

  for (int i = 0 ; i < this->numRows; i++) {
    
    int nr = this->numValuesInRow(i);      
    int colPos1 = this->binary_search_position(nr, &this->cols[this->rows[i]], col1);
      
    if ( colPos1 >= nr || this->cols[this->rows[i] + colPos1] != col1) { // col1 does not exist
        int colPos2 = this->binary_search(nr, &this->cols[this->rows[i]], col2);
        if (colPos2 >= 0) { // col 2 exists
          int colValue2 = this->values[this->rows[i] + colPos2];
          for (int j = colPos2; j > colPos1; j--) {
            this->values[this->rows[i] + j] = this->values[this->rows[i] + j - 1];
            this->cols[this->rows[i] + j] = this->cols[this->rows[i] + j - 1];
          }
          this->cols[this->rows[i] + colPos1] = col1;
          this->values[this->rows[i] + colPos1] = colValue2;

        }
        
      } else { // col1 exists
          
        int colPos2 = this->binary_search_position(nr, &this->cols[this->rows[i]], col2);
        
        if ( colPos2 >= nr || this->cols[this->rows[i] + colPos2] != col2) { // col2 does not exist
          int colValue1 = this->values[this->rows[i] + colPos1];
          for ( int j = colPos1 ; j < colPos2 - 1 ; j++) {
            this->cols[this->rows[i] + j] = this->cols[this->rows[i] + j + 1];
            this->values[this->rows[i] + j] = this->values[this->rows[i] + j + 1];
          }
          this->cols[this->rows[i] + colPos2 -1 ] = col2;
          this->values[this->rows[i] + colPos2 - 1] = colValue1;
          
        } else { // col2 exists
          int a = this->values[this->rows[i] + colPos1];
          this->values[this->rows[i] + colPos1] = this->values[this->rows[i] + colPos2];
          this->values[this->rows[i] + colPos2] = a;
        }
      }
    }
      
  return *this;
}

sparseMatrix& sparseMatrix::eye(int r, int c) {
  this->numRows = r;
  this->numCols = c;
  int numValues = min(r,c);
  this->values = (int *) malloc(numValues*sizeof(int));
  this->cols = (int *) malloc(numValues*sizeof(int));
  this->rows = (int *) malloc((this->numRows+1)*sizeof(int));
  
  for (int i = 0; i < numValues; i++) {
    this->values[i] = 1;
    this->cols[i] = i;
  }
  
  for (int i = 0; i < this->numRows + 1; i++)
    if (i < numValues)
      this->rows[i] = i;
    else
      this->rows[i] = numValues;
  
  return *this;
}

sparseMatrix& sparseMatrix::eye(int n) {
  return eye(n,n);
}

sparseMatrix sparseMatrix::transpose() const {

  // to index format
  int * rows = (int *) malloc(this->length()*sizeof(int));
  this->index_format(rows);
  
  // swaping & to rcompresed format
  sparseMatrix M;
  M.read_index_format(this->numCols,this->numRows,this->length(),this->cols,rows,this->values);
  
  return M;
  

}

sparseMatrix& sparseMatrix::multiplyByTransposed(sparseMatrix const & M1,sparseMatrix const & M2) {
  
  // TODO : Per a que necesito fer una copia?
  
M1.print_octave(cerr);
M2.print_octave(cerr);
  
  int numRows1 = M1.numRows, numCols1 = M1.numCols, *rows1, *cols1, *values1;
  int numRows2 = M2.numRows, numCols2 = M2.numCols, *rows2, *cols2, *values2;
  
  M1.decompose(&values1, &cols1, &rows1);
  M2.decompose(&values2, &cols2, &rows2);

  int numNonNullRows1 = 0, numNonNullRows2 = 0;
  for (int i = 0; i<numRows1; i++)
    if (M1.numValuesInRow(i) > 0)
      numNonNullRows1++;
  for (int i = 0; i<numRows2; i++)
    if (M2.numValuesInRow(i) > 0)
      numNonNullRows2++;
  
  this->numRows = M1.size(1);
  this->numCols = M2.size(1); // since we multiply by its transposed
  this->values = (int *) malloc(this->length()*sizeof(int));
  this->cols = (int *) malloc(this->length()*sizeof(int));
  this->rows = (int *) malloc((this->numRows+1)*sizeof(int));

  
  this->rows[0] = 0;
  for (int i = 0; i < numRows1; i++) {
    this->rows[i+1] = this->rows[i];
    for (int j = 0; j < numRows2; j++) {
      int s = this->multiplyRows(M1.numValuesInRow(i),&cols1[rows1[i]],&values1[rows1[i]],M2.numValuesInRow(j),&cols2[rows2[j]],&values2[rows2[j]]);
      if (s != 0) {
        this->values[this->rows[i+1]] = s;
        this->cols[this->rows[i+1]] = j;
        this->rows[i+1]++;
      }
    }
  }
   
  // TODO: resize cols and values arrays
  
  return *this;
 
}

sparseMatrix sparseMatrix::operator*(const sparseMatrix& M2) const {
  sparseMatrix M1;
  M1.multiplyByTransposed(*this,M2.transpose());
  return M1;
}

sparseMatrix sparseMatrix::operator*(int a) const {
  
  if (a == 0)
    return sparseMatrix(this->numRows, this->numCols);
    
  sparseMatrix M(*this);
  for (int i = 0; i < M.length() ; i++)
    M.values[i] *= a;
    
  return M;
   
}

sparseMatrix sparseMatrix::operator+(const sparseMatrix& M) const {

  if ( M.numCols != this->numCols || this->numRows != M.numRows)
    cerr << "Error: In operation A + B dimensions does not match: A is " << this->numRows << "x" << this->numCols << " and B is " << M.numRows << "x" << M.numCols << endl;
  
  if (this->numCols < M.numCols)
    return M + (*this);
  
  int * values = (int *) malloc((this->length()+M.length())*sizeof(int));
  int * cols = (int *) malloc((this->length()+M.length())*sizeof(int));
  int * rows = (int *) malloc((max(this->numRows,M.numRows) + 1)*sizeof(int));
  
  int numValues = 0;
  
  for (int i = 0; i < this->numRows; i++) {
    rows[i] = numValues;
    numValues += this->sumRows(this->numValuesInRow(i), &this->cols[this->rows[i]], &this->values[this->rows[i]],
                               M.numValuesInRow(i), &M.cols[M.rows[i]], &M.values[M.rows[i]],
                               &cols[numValues], &rows[numValues]);
  }
  
  rows[this->numRows] = numValues;
  
  return sparseMatrix(this->numRows,this->numCols,rows,cols,values);
}

sparseMatrix sparseMatrix::operator[](int row) const {

  if (this->numRows > 1 && this->numValuesInRow(row) > 0) {
    // returning row
    int * rows = (int *) malloc(2*sizeof(int));
    rows[0] = 0;
    rows[1] = this->numValuesInRow(row);
    return sparseMatrix(1,this->numCols, rows, &(this->cols[this->rows[row]]), &(this->values[this->rows[row]]));
  } else if ( this->numRows > 1 && this->numValuesInRow(row) <= 0 ) {
    // row of zeros
    int * rows = (int * ) calloc(2,sizeof(int));
    return sparseMatrix(1,this->numCols,rows, NULL,NULL);
  } else if ( this->numRows == 1 && this->length() > 0 ) {
    int pos = this->binary_search(this->length(), this->cols, row);
    if (pos < 0) {
      // return 1x1 matrix with 0 value
      int * rows = (int *) calloc(2,sizeof(int));
      return sparseMatrix(1,1,rows,NULL,NULL);  
    } else {
      // return 1x1 matrix with the value of the column
      int cols = 0;
      int values = this->values[pos];
      int * rows = (int *) malloc(2*sizeof(int));
      rows[0] = 0;
      rows[1] = 1;
      return sparseMatrix(1,1,rows,&cols,&values);      
    }
      
  }
  
  // empty matrix
  int a = 0;
  return sparseMatrix(0,0,&a,NULL,NULL);
}

int sparseMatrix::operator()(int row, int col) const {

  if (this->numValuesInRow(row) > 0) {
    int pos = this->binary_search(this->numValuesInRow(row), &(this->cols[this->rows[row]]), col);
    if (pos < 0)
      return 0;
    else
      return this->values[this->rows[row] + pos];
  }  
   
  return 0;    
}

const sparseMatrix& sparseMatrix::LDU_efficient(sparseMatrix& L, sparseMatrix& D, sparseMatrix& U, sparseMatrix& rowPerm, sparseMatrix& colPerm) const {
  
  ///////////////////////////////
  //
  //  If call this matrix M, with dimensions RxC,
  //  above matrices will be matrices so that
  //
  //  P*M*Q = L*D^{-1}*U
  //
  // for P_{RxR} and Q_{CxC} permutation, L_{RxK} lower triangular,
  // D_{KxK} diagonal, and U_{KxC} upper triangular entire matrices,
  // where K is the rank of M.
  //
  // Moreover, if d_i are the diagonal elements of D,
  // d_i divides d_j for all j >= i.
  //
  ////////////////////////////////////
  
  sparseMatrix M(*this);
  
  int RANK_MAX = min(M.numCols, M.numRows);
  
  L = sparseMatrix(RANK_MAX,M.numRows); // L is transposed for convenience
  U = sparseMatrix(RANK_MAX,M.numCols);
  int * seq = (int *) malloc((max(this->numCols,this->numRows) + 1)*sizeof(int));
  int * ones = (int *) malloc((max(this->numCols,this->numRows) + 1)*sizeof(int));
  
  for (int i = 0 ; i < max(this->numCols,this->numRows) + 1; i++) {
    seq[i] = i;
    ones[i] = 1;
  }
  
  rowPerm = sparseMatrix(this->numRows, this->numRows,  seq,seq,ones);
  colPerm = sparseMatrix(this->numCols, this->numCols,  seq,seq,ones);
  
  int * diagonal = (int *) malloc(RANK_MAX*sizeof(int));
  
  int acum = 1;
  int k = 0;
  
  while (k < RANK_MAX && M.length() > 0) {
    // take the non-zero element
    int r = 0, c = 0;
    while (r < M.numRows && M.rows[r] < 0)
      r++; 
    c = M.cols[M.rows[r]];
    // put it on (k,k)
    M.swapCols(k,c).swapRows(k,r);
    L.swapCols(k,r);
    U.swapCols(k,c);
    rowPerm.swapCols(k,r);
    colPerm.swapRows(k,c);
  
    sparseMatrix M_trans = M;
    M_trans.transpose();
  
    //algorimth
    int M_numValuesInRow = M.numValuesInRow(k);
    U.values = (int *) realloc(U.values,(U.length() + M_numValuesInRow)*sizeof(int));
    U.cols = (int *) realloc(U.cols,(U.length() + M_numValuesInRow)*sizeof(int));
    U.rows[k] = U.length();
    for (int i = 0; i < M_numValuesInRow; i++) {
      U.cols[U.rows[k] + i] = M.cols[M.rows[k] + i];
      U.values[U.rows[k] + i] = M.values[M.rows[k] + i];
    }

    int M_numValuesInCol = M_trans.numValuesInRow(k);
    int g = this->gcd(M_numValuesInCol,&M_trans.values[M_trans.rows[k]]);
if (g == 0) {
  cerr << "ERROR, g = 0" << endl;
  cerr << M_numValuesInCol << endl;
  for (int l=0; l < M_numValuesInCol ; l++)
    cerr << " " << M_trans.values[M_trans.rows[k]+l];
  cerr << endl;
}
    int d = M.values[M.rows[k]] / g;

    L.values = (int *) realloc(L.values,(L.length() + M_numValuesInCol)*sizeof(int));
    L.cols = (int *) realloc(L.cols,(L.length() + M_numValuesInCol)*sizeof(int));
    L.rows[k] = L.length();
    for (int i = 0; i < M_numValuesInCol; i++) {
      L.cols[L.rows[k] + i] = M_trans.cols[M.rows[k] + i];
      L.values[L.rows[k] + i] = M_trans.values[M.rows[k] + i]/g;
    }
  
    acum *= d;
    diagonal[k] = acum;
    
    sparseMatrix A = L[k]*(-1);
    
    A.transpose();
    
    M = M*d + A*U[k];
    
    // iteration

    k++;
  
  }
  
  L.numRows = k;
  U.numRows = k;
  D = sparseMatrix(k,k,seq,seq,diagonal);

  L.transpose();
  return *this;
}

sparseMatrix sparseMatrix::ker() const {
  sparseMatrix L,D,U,P,Q;
  this->LDU_efficient(L,D,U,P,Q);
  int rank = D.size(1);
  int dim = this->numCols - rank;
  
  int * v = (int *) malloc(dim*(rank+1)*sizeof(int));
  int * c = (int *) malloc(dim*(rank+1)*sizeof(int));
  int * r = (int *) malloc(dim*sizeof(int));

  int nonNullValuesInRow = U.numValuesInRow(rank-1) - 1;
  
  int cnt_nonNullValues = 0;
  int cnt_rows = 0;

  for (int i = rank, nv = 0 ; i < this->numCols; i++) {
    int * coeffs = (int *) malloc((rank + 1)*sizeof(int));
    int * col_coeffs = (int *) malloc((rank + 1)*sizeof(int));
    for (int j = 0; j < rank; j++)
      col_coeffs[j] = j;
    col_coeffs[rank] = i;
    
    int * actualRow = &U.rows[rank-1];
    int * cols = &U.cols[*actualRow];
    int * values = &U.values[*actualRow];
    
    if ( nv >= nonNullValuesInRow || i < cols[nv + 1]) {
      coeffs[rank-1] = 0;
    } else if ( i == cols[nv+1]) {
      coeffs[rank - 1] = -values[nv+1];
      nv++;
    } else {
      // no pot passar mai
cerr << "Violation assertion: In sparseMatrix::ker()" << endl;
    }
    
    coeffs[rank] = values[0];
    
    int g = this->gcd(2,&coeffs[rank-1]);
    coeffs[rank-1] /= g;
    coeffs[rank] /= g;
    
    for (int j = rank - 2; j >= 0 ; j--) {
      actualRow = &U.rows[j];
      cols = &U.cols[*actualRow];
      values = &U.values[*actualRow];
      int numValues = U.numValuesInRow(j);
      
      int * coeffs2 = (int *) malloc(2*sizeof(int));
      if (numValues > 1)
        coeffs2[0] = -1*this->multiplyRows(rank-j,&col_coeffs[j+1],&coeffs[j+1],numValues-1,&cols[1],&values[1]);
      else
        coeffs2[0] = 0;
      
      coeffs2[1] = values[0];
      
      g = this->gcd(2,coeffs2);
      coeffs[j] = coeffs2[0]/g;
      this->multiply(coeffs2[1]/g,rank-j,coeffs);
      
    }
    
    // copy to global
    
    r[cnt_rows] = cnt_nonNullValues;
    cnt_rows++;
    for (int j = 0; j < rank + 1 ; j++)
      if ( coeffs[j] != 0) {
        c[cnt_nonNullValues] = col_coeffs[j];
        v[cnt_nonNullValues] = coeffs[j];
        cnt_nonNullValues++;
      }
  }
  
  
  sparseMatrix V(cnt_rows,U.numCols,r,c,v),W;
  W.multiplyByTransposed(V,Q).transpose();
  return W;
}


/////////////////////////////////////////////////////////////////////////
//
//  Outside class
//

sparseMatrix operator*(int a, const sparseMatrix& M) {
  return M*a;
}
