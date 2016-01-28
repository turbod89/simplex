#include "sparseMatrix.h"

sparseMatrix::sparseMatrix() {
  this->values = NULL;
  this->cols = NULL;
  this->rows = NULL;
}

sparseMatrix::sparseMatrix(int r, int c, const int * const rows, const int * const cols, const int * const values) {
  this->numRows = r;
  this->numCols = c;

  this->rows = (int *) malloc((r+1)*sizeof(int));
  memcpy(this->rows,rows,(r+1)*sizeof(int));

  this->cols = (int *) malloc(this->length()*sizeof(int));
  memcpy(this->cols,cols,this->length()*sizeof(int));
  this->values = (int *) malloc(this->length()*sizeof(int));
  memcpy(this->values,values,this->length()*sizeof(int));
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
  this->values = NULL;
  this->cols = NULL;
  this->rows = NULL;
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
  if (this->rows && !(this->rows == NULL))
    free(this->rows);
    //this->rows = (int *) realloc(this->rows,(this->numRows+1)*sizeof(int));
  //else
    this->rows = (int *) malloc((this->numRows+1)*sizeof(int));
  memcpy(this->rows,M.rows,(this->numRows+1)*sizeof(int));

  if (this->values && this->values != NULL)
    free(this->values);
    //this->values = (int *) realloc(this->values,M.length()*sizeof(int));
  //else  
    this->values = (int *) malloc(M.length()*sizeof(int));

  if (this->cols && this->cols != NULL)
    free(this->cols);
    //this->cols = (int *) realloc(this->cols,M.length()*sizeof(int));
  //else
    this->cols = (int *) malloc(M.length()*sizeof(int));
  memcpy(this->values,M.values,M.length()*sizeof(int));
  memcpy(this->cols,M.cols,M.length()*sizeof(int));
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

const sparseMatrix& sparseMatrix::decompose(int ** values, int ** cols, int ** rows) const {
  *values = (int *) malloc(this->length()*sizeof(int));
  *cols = (int *) malloc(this->length()*sizeof(int));
  *rows = (int *) malloc((this->numRows+1)*sizeof(int));
  memcpy(*values,this->values,this->length()*sizeof(int));
  memcpy(*cols,this->cols,this->length()*sizeof(int));
  memcpy(*rows,this->rows,(this->numRows+1)*sizeof(int));
  return *this;
}
sparseMatrix& sparseMatrix::read_index_format(int numRows, int numCols, int length, int const * rows, int const * cols, int const * values, bool sorted) {
  
  int * rows2 = (int *) malloc(length*sizeof(int));
  int * cols2 = (int *) malloc(length*sizeof(int));
  int * values2 = (int *) malloc(length*sizeof(int));

  memcpy(rows2,rows,length*sizeof(int));
  memcpy(cols2,cols,length*sizeof(int));
  memcpy(values2,values,length*sizeof(int));

  if ( !sorted) {
    // sorting
    
    int ** a = (int **) malloc(2*sizeof(int *));
    a[0] = rows2;
    a[1] = values2;
    this->mergeSort(length, cols2, 2, a);
    a[0] = cols2;
    this->mergeSort(length, rows2, 2, a);
    free(a);
  }

  int * crows = (int *) malloc((numRows + 1)*sizeof(int));
  
  crows[0] = 0;
  for (int i = 0, v = 0; i < numRows; i++) {
    crows[i+1] = crows[i];
    while (v < length && i == rows2[v]) {
      v++;
      crows[i+1]++;
    }
  }

  *this = sparseMatrix(numRows,numCols,crows,cols2,values2);

  //free
  free(crows);
  free(rows2);
  free(cols2);
  free(values2);

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

sparseMatrix& sparseMatrix::removeZeros() {
  
  for (int v =0, r = 0, numValues = this->length(), numZeros = 0; v < numValues || r < this->numRows;) {
    while (r < this->numRows && (v >= numValues || this->rows[r+1] - numZeros <= v)) {
      r++;
      this->rows[r] -= numZeros;
    }

    if (v < numValues) {    
      this->values[v] = this->values[v+numZeros];
      this->cols[v] = this->cols[v+numZeros];
      if (this->values[v] == 0) {
        numZeros++;
        numValues--;
      } else
        v++;
    }
  }
  
  this->cols = (int *) realloc(this->cols,this->length()*sizeof(int));
  this->values = (int *) realloc(this->values,this->length()*sizeof(int));

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

  //free

  if (auxLength > 0 && aux != NULL) {
    for ( int i = 0; i < auxLength; i++) {
      free(aux1[i]);
      free(aux2[i]);
    }
    free(aux1);
    free(aux2);
  }



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

    int g = Tools::gcd(n - k, &L[k*n + k]);
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
if (nr < 0) {
cerr << "num values in row: " << nr<<endl;
for (int l = 0; l <nr; l++)
  cerr << " " << this->cols[this->rows[i]+l];
cerr << endl;
cerr << " " << col2 << endl;
cerr << "This row is the " << i << "-th row of the matrix" << endl;
this->print(cerr);
}
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

sparseMatrix& sparseMatrix::deleteRows(int n, int * const v) {

  if (n <= 0)
    return *this;

  this->mergeSort(n, v);

  int numValues = this->rows[this->numRows];
  int erasedValues = 0;
  int erasedRows = 0;
  for (int i = 0 ; i < this->numRows - erasedRows ;) {

    if (erasedRows < n && i == v[erasedRows]-erasedRows) {
      int numValuesInRow = this->rows[i+1+erasedRows] - this->rows[i+erasedRows];
      erasedValues += numValuesInRow;
      erasedRows++;
    } else {
      int numValuesInRow = this->rows[i+1+erasedRows] - this->rows[i+erasedRows];
      this->rows[i] = this->rows[i+erasedRows] - erasedValues;
      for (int j = this->rows[i]; erasedValues > 0 && j < this->rows[i] + numValuesInRow ; j++) {
        this->cols[j] = this->cols[j+erasedValues];
        this->values[j] = this->values[j+erasedValues];
      }
      i++;
    }
  }

  this->numRows -= erasedRows;
  this->rows = (int *) realloc(this->rows,(this->numRows + 1)*sizeof(int));
  this->cols = (int *) realloc(this->cols,(numValues - erasedValues)*sizeof(int));
  this->values = (int *) realloc(this->values,(numValues - erasedValues)*sizeof(int));
  return *this;
}

sparseMatrix& sparseMatrix::deleteRow(int i) {
  return this->deleteRows(1,&i);
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
  free(rows);
  return M;
}

sparseMatrix& sparseMatrix::multiplyByTransposed(sparseMatrix const & M1,sparseMatrix const & M2) {
  
  // conta cota superior del nombre de valors
  int numNonNullRows1 = 0, numNonNullRows2 = 0;
  for (int i = 0; i<M1.size(1); i++)
    if (M1.numValuesInRow(i) > 0)
      numNonNullRows1++;
  for (int i = 0; i<M2.size(1); i++)
    if (M2.numValuesInRow(i) > 0)
      numNonNullRows2++;
  
  this->numRows = M1.size(1);
  this->numCols = M2.size(1); // since we multiply by its transposed
  int numEstimatedValues = numNonNullRows1*numNonNullRows2;
  this->values = (int *) malloc(numEstimatedValues*sizeof(int));
  this->cols = (int *) malloc(numEstimatedValues*sizeof(int));
  this->rows = (int *) malloc((this->numRows+1)*sizeof(int));

  
  this->rows[0] = 0;
  for (int i = 0; i < M1.size(1); i++) {
    this->rows[i+1] = this->rows[i];
    for (int j = 0; j < M2.size(1); j++) {
      int s = this->multiplyRows(M1.numValuesInRow(i),&M1.cols[M1.rows[i]],&M1.values[M1.rows[i]],M2.numValuesInRow(j),&M2.cols[M2.rows[j]],&M2.values[M2.rows[j]]);
      if (s != 0) {
        this->values[this->rows[i+1]] = s;
        this->cols[this->rows[i+1]] = j;
        this->rows[i+1]++;
      }
    }
  }

  this->cols = (int *) realloc(this->cols,this->length()*sizeof(int));
  this->values = (int *) realloc(this->values,this->length()*sizeof(int));
  
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
  
  rows[0] = 0;
  for (int i = 0; i < this->numRows; i++) {
    rows[i+1] = rows[i] + this->sumRows(this->numValuesInRow(i), &this->cols[this->rows[i]], &this->values[this->rows[i]],
                               M.numValuesInRow(i), &M.cols[M.rows[i]], &M.values[M.rows[i]],
                               &cols[rows[i]], &values[rows[i]]);
  }

  sparseMatrix N(max(this->numRows,M.numRows),this->numCols,rows,cols,values);

  free(values);
  free(rows);
  free(cols);

  return N;

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

void sparseMatrix::LDU_pivotChoosing(const sparseMatrix * M, const sparseMatrix * Mt, int & r, int & c) const {

  r = 0;  
  c = 0;

  int min_c = 0, min_c_index = 0;
  
  // Check column with minimun number of values
  while (c < Mt->numRows && min_c != 1 ) {
    int nc = Mt->numValuesInRow(c);
      if (min_c == 0 || (nc < min_c && nc > 0) ) {
        min_c = nc;
        min_c_index = c;
      }
      c++;
    }
  
  // If this number is > 2 we check for rows
  if ( min_c > 2 ) {
    
    // check
    int min_r = 0, min_r_index = 0;
    while (r < M->numRows && min_r != 1) {
      int nr = M->numValuesInRow(r);
      if (min_r == 0 || (nr < min_r && nr >0)) {
        min_r = nr;
        min_r_index = r;
      }
      r++;
    }
    
    // is this way better?
    if (min_r > 1) {
        c = min_c_index;
        r = Mt->cols[Mt->rows[c]];
      } else {
        r = min_r_index;
        c = M->cols[M->rows[r]];
      }
      
    } else {
      c = min_c_index;
      r = Mt->cols[Mt->rows[c]];
    }

    return;
}

void sparseMatrix::LDU_permutations(sparseMatrix * M, sparseMatrix * Mt, sparseMatrix * L, sparseMatrix * U, sparseMatrix * rowPerm, sparseMatrix * colPerm, int k, int r, int c) const {

  M->swapCols(k,c).swapRows(k,r);
  Mt->swapCols(k,r).swapRows(k,c);
  L->swapCols(k,r);
  U->swapCols(k,c);
  rowPerm->swapCols(k,r);
  colPerm->swapRows(k,c);

  return;

}

void sparseMatrix::LDU_calculation_dM_LU(sparseMatrix * M, sparseMatrix * Mt, const sparseMatrix * L, const sparseMatrix * U, int k, int d) const {

    int * lu_rows = (int *) malloc(max((M->numRows + 1),M->numCols +1) *sizeof(int));
    int * lu_cols = (int *) malloc(L->numValuesInRow(k)*U->numValuesInRow(k)*sizeof(int));
    int * lu_values = (int *) malloc(L->numValuesInRow(k)*U->numValuesInRow(k)*sizeof(int));
    
    lu_rows[0] = 0;
    for (int i = 0, v = 0, r = -1; i < L->numValuesInRow(k); i++) {
      while ( r < L->cols[L->rows[k]+i]) {
        r++;
        lu_rows[r+1] = lu_rows[r];
      }
      for (int j = 0; j < U->numValuesInRow(k); j++, v++) {
        lu_values[v] = -L->values[L->rows[k]+i]*U->values[U->rows[k]+j];
        lu_cols[v] = U->cols[U->rows[k] + j];
        lu_rows[r+1]++;
      }
    }
    for (int r = L->cols[L->rows[k] + L->numValuesInRow(k) -1]+1; r < M->numRows; r++)
      lu_rows[r+1] = lu_rows[r];

    //*M = (*M)*d + *mLU;
/*
    for (int i = 0; i < M->rows[M->numRows] ; i++)
      M->values[i] *= d;

    *M = M->operator+(*mLU);
*/


    int * M_rows = (int *) malloc((M->numRows - k + 1)*sizeof(int));
    int * M_cols = (int *) malloc((M->rows[M->numRows] - M->rows[k])*sizeof(int));
    int * M_values = (int *) malloc((M->rows[M->numRows] - M->rows[k])*sizeof(int));

    memcpy(M_rows,&(M->rows[k]),(M->numRows - k + 1)*sizeof(int));
    memcpy(M_cols,&(M->cols[M->rows[k]]),(M->rows[M->numRows] - M->rows[k])*sizeof(int));
    memcpy(M_values,&(M->values[M->rows[k]]),(M->rows[M->numRows] - M->rows[k])*sizeof(int));

    for (int i = 0; i < (M->rows[M->numRows] - M->rows[k]) ; i++)
      M_values[i] *= d;

    for (int i = k; i < M->numRows; i++) {
      int a = this->sumRows(M_rows[i+1 - k] - M_rows[i-k], &M_cols[M_rows[i-k] - M_rows[0]], &M_values[M_rows[i-k] - M_rows[0]],
                      lu_rows[i+1]-lu_rows[i], &lu_cols[lu_rows[i]], &lu_values[lu_rows[i]],
                      &M->cols[M->rows[i]], &M->values[M->rows[i]]);
      M->rows[i+1] = M->rows[i] + a;
    }

    free(M_rows);
    free(M_cols);
    free(M_values);

    /////////////////////////
    // transpose

    int * lut_rows = lu_rows;
    int * lut_cols = lu_cols;
    int * lut_values = lu_values;
    
    lut_rows[0] = 0;
    for (int i = 0, v = 0, r = -1; i < U->numValuesInRow(k); i++) {
      while ( r < U->cols[U->rows[k]+i]) {
        r++;
        lut_rows[r+1] = lut_rows[r];
      }
      for (int j = 0; j < L->numValuesInRow(k); j++, v++) {
        lut_values[v] = -U->values[U->rows[k]+i]*L->values[L->rows[k]+j];
        lut_cols[v] = L->cols[L->rows[k] + j];
        lut_rows[r+1]++;
      }
    }
    for (int r = U->cols[U->rows[k] + U->numValuesInRow(k) -1]+1; r < Mt->numRows; r++)
      lut_rows[r+1] = lut_rows[r];

    int * Mt_rows = (int *) malloc((Mt->numRows - k + 1)*sizeof(int));
    int * Mt_cols = (int *) malloc((Mt->rows[Mt->numRows] - Mt->rows[k])*sizeof(int));
    int * Mt_values = (int *) malloc((Mt->rows[Mt->numRows] - Mt->rows[k])*sizeof(int));

    memcpy(Mt_rows,&(Mt->rows[k]),(Mt->numRows - k + 1)*sizeof(int));
    memcpy(Mt_cols,&(Mt->cols[Mt->rows[k]]),(Mt->rows[Mt->numRows] - Mt->rows[k])*sizeof(int));
    memcpy(Mt_values,&(Mt->values[Mt->rows[k]]),(Mt->rows[Mt->numRows] - Mt->rows[k])*sizeof(int));

    for (int i = 0; i < (Mt->rows[Mt->numRows] - Mt->rows[k]) ; i++)
      Mt_values[i] *= d;

    for (int i = k; i < Mt->numRows; i++) {
      int a = this->sumRows(Mt_rows[i+1 - k] - Mt_rows[i-k], &Mt_cols[Mt_rows[i-k] - Mt_rows[0]], &Mt_values[Mt_rows[i-k] - Mt_rows[0]],
                      lut_rows[i+1]-lut_rows[i], &lut_cols[lut_rows[i]], &lut_values[lut_rows[i]],
                      &Mt->cols[Mt->rows[i]], &Mt->values[Mt->rows[i]]);
      Mt->rows[i+1] = Mt->rows[i] + a;
    }

    free(Mt_rows);
    free(Mt_cols);
    free(Mt_values);



    free(lu_rows);
    free(lu_cols);
    free(lu_values);

    return;

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
  
  // set variables
  sparseMatrix * M = new sparseMatrix(*this);
  sparseMatrix * M_trans = new sparseMatrix(M->transpose());

  M->cols = (int *) realloc(M->cols,(M->numRows * M->numCols)*sizeof(int));
  M->values = (int *) realloc(M->values,(M->numRows * M->numCols)*sizeof(int));
  M_trans->cols = (int *) realloc(M_trans->cols,(M_trans->numRows * M_trans->numCols)*sizeof(int));
  M_trans->values = (int *) realloc(M_trans->values,(M_trans->numRows * M_trans->numCols)*sizeof(int));

  int RANK_MAX = min(M->numCols, M->numRows);
  
  L = sparseMatrix(RANK_MAX,M->numRows); // L is transposed for convenience
  L.numRows = 0;
  L.rows[0] = 0;
  L.cols = (int *) realloc(L.cols,(RANK_MAX*M->numRows - RANK_MAX*(RANK_MAX-1)/2)*sizeof(int));
  L.values = (int *) realloc(L.values,(RANK_MAX*M->numRows - RANK_MAX*(RANK_MAX-1)/2)*sizeof(int));
  U = sparseMatrix(RANK_MAX,M->numCols);
  U.numRows = 0;
  U.rows[0] = 0;
  U.cols = (int *) realloc(U.cols,(RANK_MAX*M->numCols - RANK_MAX*(RANK_MAX-1)/2)*sizeof(int));
  U.values = (int *) realloc(U.values,(RANK_MAX*M->numCols - RANK_MAX*(RANK_MAX-1)/2)*sizeof(int));
  
  int * seq = (int *) malloc((max(this->numCols,this->numRows) + 1)*sizeof(int));
  int * ones = (int *) malloc((max(this->numCols,this->numRows) + 1)*sizeof(int));
  
  for (int i = 0 ; i < max(this->numCols,this->numRows) + 1; i++) {
    seq[i] = i;
    ones[i] = 1;
  }
  
  rowPerm = sparseMatrix(this->numRows, this->numRows,seq,seq,ones);
  colPerm = sparseMatrix(this->numCols, this->numCols,seq,seq,ones);
  
  int * diagonal = (int *) malloc(RANK_MAX*sizeof(int));
  
  // start algorimth
  
  int acum = 1;
  int k = 0;

  while (k < RANK_MAX && M->length() > 0) {
        
    // pivot choosing
    int r = 0, c = 0;
    this->LDU_pivotChoosing(M,M_trans,r,c);
    
    // put it on (k,k)
    this->LDU_permutations(M,M_trans,&L,&U,&rowPerm,&colPerm,k,r,c);

    // U matrix
    int M_numValuesInRow = M->numValuesInRow(k);
    U.numRows++;
    U.rows[k+1] = U.rows[k] + M_numValuesInRow;
    memcpy(&U.cols[U.rows[k]], &M->cols[M->rows[k]],M_numValuesInRow*sizeof(int));
    memcpy(&U.values[U.rows[k]], &M->values[M->rows[k]],M_numValuesInRow*sizeof(int));

    // diagonal
    int M_numValuesInCol = (*M_trans).numValuesInRow(k);
    int g = Tools::gcd(M_numValuesInCol,&(*M_trans).values[(*M_trans).rows[k]]);
    int d = M->values[M->rows[k]] / g;

    if (d < 0) {
      d *= -1;
      g *= -1;
    }

    // L matrix
    L.numRows++;
    L.rows[k+1] = L.rows[k] + M_numValuesInCol;
    memcpy(&L.cols[L.rows[k]],&(*M_trans).cols[(*M_trans).rows[k]],M_numValuesInCol*sizeof(int));
    for (int i = 0; i < M_numValuesInCol; i++) {
      L.values[L.rows[k] + i] = (*M_trans).values[(*M_trans).rows[k] + i]/g;
    }
    
    // Diagonal

    acum *= d;
    diagonal[k] = acum;

    // M := dM - LU, Mt := M^t
    this->LDU_calculation_dM_LU(M,M_trans,&L,&U,k,d);
    
    // iteration

    k++;
  
  }
  
  L.numRows = k;
  U.numRows = k;
  D = sparseMatrix(k,k,seq,seq,diagonal);


  delete M;
  delete M_trans;

  L.cols = (int *) realloc(L.cols,L.rows[L.numRows]*sizeof(int));
  L.values = (int *) realloc(L.values,L.rows[L.numRows]*sizeof(int));
  U.cols = (int *) realloc(U.cols,U.rows[U.numRows]*sizeof(int));
  U.values = (int *) realloc(U.values,U.rows[U.numRows]*sizeof(int));

  L = L.transpose();
  return *this;
}

sparseMatrix sparseMatrix::ker() const {
  sparseMatrix L,D,U,P,Q;
  this->LDU_efficient(L,D,U,P,Q);
  int rank = D.size(1);
  int dim = this->numCols - rank;
  int * v = (int *) malloc(dim*(rank+1)*sizeof(int));
  int * c = (int *) malloc(dim*(rank+1)*sizeof(int));
  int * r = (int *) malloc((dim+1)*sizeof(int));

  int nonNullValuesInRow = U.numValuesInRow(rank-1) - 1;
  
  int cnt_nonNullValues = 0;
  int cnt_rows = 0;
  r[0] = 0;
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
cerr << "Violation of assertion: In sparseMatrix::ker()" << endl;
    }
    
    coeffs[rank] = values[0];
    
    int g = Tools::gcd(2,&coeffs[rank-1]);
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
      
      g = Tools::gcd(2,coeffs2);
      coeffs[j] = coeffs2[0]/g;
      this->multiply(coeffs2[1]/g,rank-j,&coeffs[j+1]);

      free(coeffs2);            
    }
    
    // copy to global
    
    cnt_rows++;
    for (int j = 0; j < rank + 1 ; j++)
      if ( coeffs[j] != 0) {
        c[cnt_nonNullValues] = col_coeffs[j];
        v[cnt_nonNullValues] = coeffs[j];
        cnt_nonNullValues++;
      }

    r[cnt_rows] = cnt_nonNullValues;

    //free
    free(coeffs);
    free(col_coeffs);
  }

  sparseMatrix V(cnt_rows,U.numCols,r,c,v),W;

  //free
  free(r);
  free(c);
  free(v);

  return W.multiplyByTransposed(V,Q).transpose();
}

sparseMatrix sparseMatrix::LXeqY(const sparseMatrix &Y) const {
  // THIS IS MATRIX L
  // we will assum this matrix is lower triangular
  // with non null values on the diagonal

  int numRows = Y.numCols;
  int numCols = this->numCols;
  int * r = (int *) malloc((numRows+1)*sizeof(int));
  int * c = (int *) malloc((numRows*numCols)*sizeof(int));
  int * v = (int *) malloc((numRows*numCols)*sizeof(int));
//cerr << this->size(1) << " x " << this->size(2) << " " << Y.size(1) << " x " << Y.size(2)  << endl;
  r[0] = 0;
  for (int row = 0; row < numRows; row++) {
    int * local_c = &c[r[row]];
    int * local_v = &v[r[row]];

    bool isThisRowPossible = true;
    int nvr = 0; // number of non null values on this row
    for (int col = 0; isThisRowPossible && col < numCols; col++) {
      //
      // equation B + x_k l_kk = Y_col,row
      //
      // get independent term B
      int B = 0;
      for (int i = 0, j = 0; i < nvr && j < this->numValuesInRow(col);)
        if ( local_c[i] == this->cols[this->rows[col] + j]) {
          B += local_v[i]*this->values[this->rows[col] + j];
          i++;
          j++;
        } else if (local_c[i] < this->cols[this->rows[col] + j]) {
          i++;
        } else if (local_c[i] > this->cols[this->rows[col] + j]) {
          j++;
        } else {
          // no pot ocorrer
        }

      // GET A
      int A = Y(col,row) - B;
      // GET l_kk
      int l_kk = this->values[this->rows[col] + this->numValuesInRow(col)-1];
      // is this row possible?
      isThisRowPossible = l_kk == 1 || l_kk == -1 || A%l_kk == 0;
      // X_row,col
      int X = A/l_kk;
      // allocate if corresponds
      if (X != 0) {
        local_c[nvr] = col;
        local_v[nvr] = X;
//cerr << "(" << row << ", " << col << ", " << X << ")" << endl;
        nvr++;
      }

    }

    // now check excedent conditions

    for (int i = this->numCols ; isThisRowPossible && i < Y.numRows; i++) {
      int a = this->multiplyRows(nvr, local_c, local_v, this->numValuesInRow(i), &this->cols[this->rows[i]],&this->values[this->rows[i]] );
      isThisRowPossible = (a == Y(i,row));
    }

    // check if this row is a solution
    if (isThisRowPossible) {
      r[row+1] = r[row] + nvr;
    } else {
      r[row+1] = r[row];
    }

  }

  c = (int *) realloc(c,r[numRows]*sizeof(int));
  v = (int *) realloc(v,r[numRows]*sizeof(int));

  sparseMatrix W(numRows, numCols, r,c,v);

  //free
  free(c);
  free(v);

  return W.transpose();
}

sparseMatrix sparseMatrix::LComplementary(const sparseMatrix &Y) const {
  // THIS IS MATRIX L
  // we will assum this matrix is lower triangular
  // with non null values on the diagonal

  sparseMatrix L(this->numCols,this->numCols,this->rows,this->cols,this->values);
  sparseMatrix Y2(this->numCols,Y.size(2),Y.getRows(),Y.getCols(),Y.getValues());
  sparseMatrix X = L.LXeqY(Y2);

//cerr << this->size(1) << " x " << this->size(2) << " " << X.size(1) << " x " << X.size(2)  << endl;


  sparseMatrix A = (*this)*X +Y*(-1);
  // Aqui es on triga
  sparseMatrix l,d,u,p,q;
  A.LDU_efficient(l,d,u,p,q);

  return p*l;

}


/////////////////////////////////////////////////////////////////////////
//
//  Outside class
//

sparseMatrix operator*(int a, const sparseMatrix& M) {
  return M*a;
}
