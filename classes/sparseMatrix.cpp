#include "sparseMatrix.h"

sparseMatrix::sparseMatrix(){
  this->values = NULL;
  this->cols = NULL;
  this->rows = NULL;
}

sparseMatrix::sparseMatrix(int r, int c, int v,const int * const values,const int * const cols,const int * const rows) {
  this->numRows = r;
  this->numCols = c;
  this->numValues = v;
  this->values = (int *) malloc(v*sizeof(int));
  this->cols = (int *) malloc(v*sizeof(int));
  this->rows = (int *) malloc(r*sizeof(int));
  memcpy(this->values,values,v*sizeof(int));
  memcpy(this->cols,cols,v*sizeof(int));
  memcpy(this->rows,rows,r*sizeof(int));
}

sparseMatrix::sparseMatrix(int r , int c) {
  this->numRows = r;
  this->numCols = c;
  this->numValues = 0;
  this->values = NULL;
  this->cols = NULL;
  this->rows = (int *) malloc(r*sizeof(int));
  
  for (int i = 0 ; i < this->numRows; i++)
    this->rows[i] = -1;

}

sparseMatrix::sparseMatrix(int * M, int r , int c) {

  int cnt = 0;
  for (int i = 0; i < r*c; i++)
    if (M[i] != 0)
      cnt++;

    this->numCols = c;
    this->numRows = r;
    this->numValues = cnt;
    this->values = (int *) malloc(cnt*sizeof(int));
    this->cols = (int *) malloc(cnt*sizeof(int));
    this->rows = (int *) malloc(r*sizeof(int));

    cnt = 0;
    for (int i = 0; i < r; i++) {
      this->rows[i] = -1;
      for (int j = 0; j < c; j++) {
        if (M[i*c + j] != 0) {
          if (this->rows[i] == -1)
            this->rows[i] = cnt;
          this->cols[cnt] = j;
          this->values[cnt] = M[i*c + j];
          cnt++;
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
  this->numValues = M.numValues;

  this->values = (int *) malloc(this->numValues*sizeof(int));
  this->cols = (int *) malloc(this->numValues*sizeof(int));
  this->rows = (int *) malloc(this->numRows*sizeof(int));
  memcpy(this->values,M.values,this->numValues*sizeof(int));
  memcpy(this->cols,M.cols,this->numValues*sizeof(int));
  memcpy(this->rows,M.rows,this->numRows*sizeof(int));
  
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

int sparseMatrix::length() const {
  return this->numValues;
}

const sparseMatrix& sparseMatrix::decompose(int &r, int &c, int &v, int ** values, int ** cols, int ** rows) const {
  r = this->numRows;
  c = this->numCols;
  v = this->numValues;
  *values = (int *) malloc(v*sizeof(int));
  *cols = (int *) malloc(v*sizeof(int));
  *rows = (int *) malloc(r*sizeof(int));
  memcpy(*values,this->values,v*sizeof(int));
  memcpy(*cols,this->cols,v*sizeof(int));
  memcpy(*rows,this->rows,r*sizeof(int));
  return *this;
}

const sparseMatrix& sparseMatrix::format2(int &r, int &c, int &v, int ** rows, int ** cols, int ** values) const {
  r = this->numRows;
  c = this->numCols;
  v = this->numValues;
  *values = (int *) malloc(v*sizeof(int));
  *cols = (int *) malloc(v*sizeof(int));
  *rows = (int *) malloc(v*sizeof(int));

  if (rows != NULL) {
    // from format1 to format2
    // O(min(numValues,numRows))
    int lastValue = this->numValues;
    for (int i = this->numRows - 1 ; i >= 0 && lastValue > 0; i--)
      if (this->rows[i] >= 0) {
        for (int j = this->rows[i] ; j < lastValue; j++)
          *rows[j] = i;
        lastValue = this->rows[i];
      }
  }

  if (cols != NULL)  
    memcpy(*values,this->values,v*sizeof(int));
  if (values != NULL)
    memcpy(*cols,this->cols,v*sizeof(int));
    
  return *this;
}


sparseMatrix& sparseMatrix::read(istream& in) {

  in >> this->numRows >> this->numCols >> this->numValues;
  
  this->values = (int *) malloc(this->numValues*sizeof(int));
  this->cols = (int *) malloc(this->numValues*sizeof(int));
  this->rows = (int *) malloc(this->numRows*sizeof(int));

  for (int i = 0; i < this->numValues; i++)
    in >> this->values[i];

  for (int i = 0; i < this->numValues; i++)
    in >> this->cols[i];

  for (int i = 0; i < this->numRows; i++)
    in >> this->rows[i];

  return *this;
  
}

const sparseMatrix& sparseMatrix::print(ostream& out) const {

  out << this->numRows << " " << this->numCols << " " << this->numValues <<endl;

  for (int i = 0; i < this->numValues; i++)
    out << " " << this->values[i];
  out << endl;

  for (int i = 0; i < this->numValues; i++)
    out << " " << this->cols[i];
  out << endl;

  for (int i = 0; i < this->numRows; i++)
    out << " " << this->rows[i];
  out << endl;

  return *this;
}

const sparseMatrix& sparseMatrix::print_full(ostream& out) const {

  out << this->numRows << " " << this->numCols << endl;

  for (int i = 0; i < this->numRows; i++) {
    if ( this->rows[i] < 0 ) {
      for (int j = 0 ; j < this->numCols; j++)
        out << " 0";
    } else {
      int nextValue = 0;
      int nr = this->numValuesInRow(i);
      for ( int j = 0; j < this->numCols ; j++)
        if (nextValue < nr && this->cols[this->rows[i] + nextValue] == j) {
          out << " " << this->values[this->rows[i] + nextValue];
          nextValue++;
        } else {
          out << " 0";
        }
    }
    out << endl;
  }

  return *this;
}

const sparseMatrix& sparseMatrix::get(int * a) const {
  
  for (int i = 0, v = 0; i < this->numRows; i++) {
      
    if ( this->rows[i] < 0 ) {
        for (int j = 0 ; j < this->numCols; j++)
        a[v] = 0;
          
      } else {
        int nextValue = 0;
        int nr = this->numValuesInRow(i);
        for ( int j = 0; j < this->numCols ; j++)
            if (nextValue < nr && this->cols[this->rows[i] + nextValue] == j) {
                a[v] = this->values[this->rows[i] + nextValue];
            v++;
                nextValue++;
            } else {
              a[v] = 0;
          v++;
            }
      }

  }
  
  return (*this);
}

int sparseMatrix::rowStartPosition(int row) const {
  
  if (this->rows[row] >= 0)
    return this->rows[row];
 
  int i = row;
  
  while (i >= 0 && this->rows[i]<0)
    i--;
  
  if ( i < 0)
    return 0;

  return this->rows[i] + this->numValuesInRow(i);
 
}

int sparseMatrix::binary_search_position(const int * const v, int length, int needle) const {
  
  if (length <= 0)
    return 0;
  
  int i = length/2;
  if ( v[i] == needle)
    return i;
  else if ( v[i] > needle)
    return binary_search_position(v,i, needle);
  else if ( v[i] < needle)
    return i + 1 + binary_search_position(&v[i+1], length - i - 1, needle); 
}

int sparseMatrix::binary_search(const int * const v, int length, int needle) const {
  int i = this->binary_search_position(v,length, needle);
  if (i == length)
    return -1;
  if (v[i] == needle)
    return i;
  return -1;
}

void sparseMatrix::mergeSort(int * v, int length, int * const * const aux, int auxLength) const {
  
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

  mergeSort(a, l1, aux1, auxLength);
  mergeSort(b, l2, aux2, auxLength);

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


int sparseMatrix::numValuesInRow(int row) const {
  if (this->rows[row] < 0 )
    return 0;
  int i = row+1;
  while(i < this->numRows && this->rows[i] < 0)
    i++;
  if (i == this->numRows)
    return this->numValues - this->rows[row];
  
  return this->rows[i] - this->rows[row];

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
      
  int row1Starts = this->rowStartPosition(swap1),
      row2Starts = this->rowStartPosition(swap2);
      
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

  for (int i = swap1 + 1 ; i < swap2 ; i++)
    if ( this->rows[i] >= 0)
      this->rows[i] += numValuesRow2 - numValuesRow1; 
  
  if (numValuesRow1 == 0 && numValuesRow2 > 0) {
    this->rows[swap2] = -1;
    this->rows[swap1] = row1Starts;
  } else if (numValuesRow1 > 0 && numValuesRow2 == 0) {
    this->rows[swap1] = -1;
    this->rows[swap2] = row2Starts + numValuesRow2 - numValuesRow1;
  } else if (numValuesRow1 > 0 && numValuesRow2 > 0) {
    this->rows[swap2] += numValuesRow2 - numValuesRow1;
  }

  return *this;

}

sparseMatrix& sparseMatrix::swapCols(int col1, int col2) {

  if (col1 > col2)
    return this->swapCols(col2,col1);
  if (col1 == col2)
    return *this;
  if (col2 >= this->numCols || col1 < 0)
    return *this;

  for (int i = 0 ; i < this->numRows; i++) 
    if (this->rows[i] >= 0) {
        
      int nr = this->numValuesInRow(i);      
      int colPos1 = this->binary_search_position(&this->cols[this->rows[i]], nr, col1);
      
      if ( colPos1 >= nr || this->cols[this->rows[i] + colPos1] != col1) { // col1 does not exist
        int colPos2 = this->binary_search(&this->cols[this->rows[i]], nr, col2);
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
          
        int colPos2 = this->binary_search_position(&this->cols[this->rows[i]], nr, col2);
        
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

sparseMatrix& sparseMatrix::deleteCols(int n, int * cols) {
  
  int numElementsDeleted = 0;
  
  this->mergeSort(cols,n);
  for (int i = 0; i < this->numRows; i++)
    if (this->rows[i] >= 0) {
      int numValuesInThisRow = this->numValuesInRow(i);
      this->rows[i] -= numElementsDeleted;
      for (int j = 0, k = 0; j < numValuesInThisRow;)
        if (k < n && this->cols[this->rows[i]+j] == cols[k]) {
          // erase element
          for ( int l = this->rows[i]+j+1; l < this->numValues-numElementsDeleted; l++) {
            this->values[l-1] = this->values[l];
            this->cols[l-1] = this->cols[l];
          }
          // update counters
          numValuesInThisRow--;
          numElementsDeleted++;
          // advance 
          k++;
        } else if (k >= n || this->cols[this->rows[i]+j] < cols[k]) {
          this->cols[this->rows[i]+j] -= k;
          j++;
        } else if (k < n && this->cols[this->rows[i]+j] > cols[k]) {
          k++;
    }
        
      if (numValuesInThisRow == 0)
        this->rows[i] = -1;
    }

  this->numValues -= numElementsDeleted;
  this->numCols -= n;
  return *this;
}

sparseMatrix& sparseMatrix::deleteCol(int col) {
  return this->deleteCols(1,&col);
}

int sparseMatrix::emptyRowsToBottom(int * rowPerm) {
    
    // copy of rowPerm
    int * v = (int *) malloc(this->numRows * sizeof(int));
    for (int i = 0 ; i < this->numRows ; i++)
        v[i] = rowPerm[i];
    
    // counting empty rows
    int numEmptyRows = 0;
    for (int i = 0; i < this->numRows; i++)
        if (this->rows[i] < 0 )
            numEmptyRows++;
    if (numEmptyRows == 0)
        return 0;
    
    // index of empty rows
    int * emptyRows = (int *) malloc(numEmptyRows*sizeof(int));

    for (int i = 0, cnt = 0; i + cnt < this->numRows; )
        if (this->rows[i+cnt] < 0 ) {
            emptyRows[cnt] = i+cnt;
            cnt++;
        } else {
            this->rows[i] = this->rows[i+cnt];
            rowPerm[i] = v[i+cnt];
             i++;
        }
    
    // setting empty rows to bottom
    for (int i = 0; i < numEmptyRows; i++) {
        rowPerm[this->numRows - numEmptyRows + i] = v[emptyRows[i]];
        this->rows[this->numRows - numEmptyRows + i] = -1;
    }

    return numEmptyRows;
}

sparseMatrix& sparseMatrix::eye(int r, int c) {
  this->numRows = r;
  this->numCols = c;
  this->numValues = min(r,c);
  this->values = (int *) malloc(this->numValues*sizeof(int));
  this->cols = (int *) malloc(this->numValues*sizeof(int));
  this->rows = (int *) malloc(this->numRows*sizeof(int));
  
  for (int i = 0; i < this->numValues; i++) {
    this->values[i] = 1;
    this->cols[i] = i;
  }
  
  for (int i = 0; i <this->numRows; i++)
    if (i < this->numValues)
      this->rows[i] = i;
    else
      this->rows[i] = -1;
  
  return *this;
}

sparseMatrix& sparseMatrix::eye(int n) {
  return eye(n,n);
}

sparseMatrix& sparseMatrix::transpose() {
  // from format1 to format2
  // O(min(numValues,numRows))
  int * rows = this->rows;
  this->rows = (int *) malloc(this->numValues*sizeof(int));
  int lastValue = this->numValues;
  for (int i = this->numRows - 1 ; i >= 0 && lastValue > 0; i--)
    if (rows[i] >= 0) {
      for (int j = rows[i] ; j < lastValue; j++)
        this->rows[j] = i;
      lastValue = rows[i];
    }

  // swap cols and rows O(1)
  int * a = this->rows;
  this->rows = this->cols;
  this->cols = a;
  int b = this->numRows;
  this->numRows = this->numCols;
  this->numCols = b;

  // sort by rows O(numValues^2)  
  int **aux = (int **) malloc(2*sizeof(int *));
  aux[0] = this->values;
  aux[1] = this->cols;
  
  this->mergeSort(this->rows, this->numValues,aux,2);

  // from format2 to format1
  // O(max(numValues,numRows))

  rows = this->rows;
  this->rows = (int *) malloc(this->numRows*sizeof(int));
  int row = -1;
  for ( int i = 0; i < this->numValues; i++)
    if (row != rows[i]) {
      for ( int j = row + 1; j < rows[i]; j++)
        this->rows[j] = -1;
      this->rows[rows[i]] = i;
      row = rows[i];
    }
  for (int i = row + 1; i < this->numRows ; i++)
    this->rows[i] = -1;

  return *this;

}

sparseMatrix& sparseMatrix::multiplyByTransposed(sparseMatrix const & M1,sparseMatrix const & M2) {
  
  int numRows1, numCols1, numValues1, *rows1, *cols1, *values1;
  int numRows2, numCols2, numValues2, *rows2, *cols2, *values2;
  
  M1.decompose(numRows1, numCols1, numValues1, &values1, &cols1, &rows1);
  M2.decompose(numRows2, numCols2, numValues2, &values2, &cols2, &rows2);

  int numNonNullRows1 = 0, numNonNullRows2 = 0;
  for (int i = 0; i<numRows1; i++)
    if (rows1[i] >= 0)
      numNonNullRows1++;
  for (int i = 0; i<numRows2; i++)
    if (rows2[i] >= 0)
      numNonNullRows2++;
  
  this->numRows = M1.size(1);
  this->numCols = M2.size(1); // since we multiply by its transposed
  this->numValues = numNonNullRows1*numNonNullRows2; // at most
  this->values = (int *) malloc(this->numValues*sizeof(int));
  this->cols = (int *) malloc(this->numValues*sizeof(int));
  this->rows = (int *) malloc(this->numRows*sizeof(int));

  
  int nonNullValues = 0;
  for (int i = 0; i < numRows1; i++) {
    this->rows[i] = -1;
    if (rows1[i] >= 0) {
      int nonNullValuesInRow = 0;
      for (int j = 0; j < numRows2; j++)
        if ( rows2[j] >= 0) {
          int s = this->multiplyRows(M1.numValuesInRow(i),&cols1[rows1[i]],&values1[rows1[i]],M2.numValuesInRow(j),&cols2[rows2[j]],&values2[rows2[j]]);
          if (s != 0) {
            this->values[nonNullValues + nonNullValuesInRow] = s;
            this->cols[nonNullValues + nonNullValuesInRow] = j;
            nonNullValuesInRow++;
          }
        }
      if (nonNullValuesInRow > 0) {
        this->rows[i] = nonNullValues;
        nonNullValues += nonNullValuesInRow;
      }
    }
  }
  
  // update number of non null values
  this->numValues = nonNullValues;
  
  return *this;
 
}

sparseMatrix sparseMatrix::operator*(sparseMatrix M2) const {
  sparseMatrix M1;
  M2.transpose();
  M1.multiplyByTransposed(*this,M2);
  return M1;
}

sparseMatrix sparseMatrix::operator*(int a) const {
  
  if (a == 0)
    return sparseMatrix(this->numRows, this->numCols);
    
  sparseMatrix M(*this);
  for (int i = 0; i < M.numValues ; i++)
    M.values[i] *= a;
    
  return M;
   
}

sparseMatrix sparseMatrix::operator+(const sparseMatrix& M) const {

  // TODO

  if ( M.numCols != this->numCols || this->numRows != M.numRows)
    cerr << "Error: In operation A + B dimensions does not match: A is " << this->numRows << "x" << this->numCols << " and B is " << M.numRows << "x" << M.numCols << endl;
  
  if (this->numCols < M.numCols)
    return M + (*this);
  
  int * values = (int *) malloc((this->numValues+M.numValues)*sizeof(int));
  int * cols = (int *) malloc((this->numValues+M.numValues)*sizeof(int));
  int * rows = (int *) malloc(max(this->numRows,M.numRows)*sizeof(int));
  
  int numValues = 0;
  
  for (int i = 0; i < this->numRows; i++)
    if (this->numValuesInRow(i) == 0 && M.numValuesInRow(i) == 0)
      rows[i] = -1;
    else if ( this->numValuesInRow(i) == 0) {
      rows[i] = numValues;
      int n = M.numValuesInRow(i);
      for (int j = M.rows[i]; j < M.rows[i] + n; j++, numValues++) {
        cols[numValues] = M.cols[j];
        values[numValues] = M.values[j];
      }
    } else if ( M.numValuesInRow(i) == 0) {
      rows[i] = numValues;
      int n = this->numValuesInRow(i);
      for (int j = this->rows[i]; j < this->rows[i] + n; j++, numValues++) {
        cols[numValues] = this->cols[j];
        values[numValues] = this->values[j];
      }
    } else {
      rows[i] = numValues;
      int n1 = this->numValuesInRow(i);
      int n2 = M.numValuesInRow(i);
      for (int j1 = this->rows[i], j2 = M.rows[i]; j1 < this->rows[i] + n1 || j2 < M.rows[i] + n2;)
        if ( j1 == this->rows[i] + n1) {
          cols[numValues] = M.cols[j2];
          values[numValues] = M.values[j2];
          j2++;
          numValues++;
        } else if ( j2 == M.rows[i] + n2) {
          cols[numValues] = this->cols[j1];
          values[numValues] = this->values[j1];
          j1++;
          numValues++;
        } else if ( this->cols[j1] > M.cols[j2]) {
          cols[numValues] = M.cols[j2];
          values[numValues] = M.values[j2];
          j2++;
          numValues++;
        } else if ( this->cols[j1] < M.cols[j2]) {
          cols[numValues] = this->cols[j1];
          values[numValues] = this->values[j1];
          j1++;
          numValues++;
        } else if (this->cols[j1] == M.cols[j2] && this->values[j1] + M.values[j2] != 0) {
          cols[numValues] = this->cols[j1];
          values[numValues] = this->values[j1] + M.values[j2];
          j1++;
          j2++;
          numValues++;
        } else {
          j1++;
          j2++;
        }
        
        if (numValues == rows[i]) // all row cancel out
          rows[i] = -1;
    }
  
  return sparseMatrix(this->numRows,this->numCols,numValues,values,cols,rows);
}

sparseMatrix sparseMatrix::operator[](int row) const {

  if (this->numRows > 1 && this->numValuesInRow(row) > 0) {
    // returning row
    int rows = 0;
    sparseMatrix M(1,this->numCols, this->numValuesInRow(row), &(this->values[this->rows[row]]), &(this->cols[this->rows[row]]), &rows);
    return M;
  } else if ( this->numRows > 1 && this->numValuesInRow(row) <= 0 ) {
    // row of zeros
    int rows = -1;
    sparseMatrix M(1,this->numCols,0,0,0, &rows);
    return M;  
  } else if ( this->numRows == 1 && this->numValues > 0 ) {
    int pos = this->binary_search(this->cols, this->numValues, row);
    if (pos < 0) {
      // return 1x1 matrix with 0 value
      int rows = -1;
      sparseMatrix M(1,1,0,0,0, &rows);
      return M;  
    } else {
      // return 1x1 matrix with the value of the column
      int rows = 0;
      int cols = 0;
      int values = this->values[pos];
      sparseMatrix M(1,1,1,&values,&cols, &rows);
      return M;      
    }
      
  }
  
  // empty matrix
  sparseMatrix M(0,0,0,0,0,0);
  return M;  
}

int sparseMatrix::operator()(int row, int col) const {

  if (this->numValuesInRow(row) > 0) {
    int pos = this->binary_search(&(this->cols[this->rows[row]]), this->numValuesInRow(row), col);
    if (pos < 0)
      return 0;
    else
      return this->cols[this->rows[row] + pos];
  }  
   
  return 0;    
}

sparseMatrix& sparseMatrix::vcat(const sparseMatrix& M) {
    
    if (M.size(2) != this->numCols)
        return *this;
    
    this->values = (int*) realloc(this->values, (this->numValues + M.length())*sizeof(int));
    this->cols = (int*) realloc(this->cols, (this->numValues + M.length())*sizeof(int));
    this->rows = (int*) realloc(this->rows, (this->numRows + M.size(1))*sizeof(int));
    
    int numCols,numRows,numValues;
    int  * newValues = &(this->values[this->numValues]), * newCols = &(this->cols[this->numValues]),* newRows = &(this->rows[this->numRows]);
    M.decompose(numRows,numCols,numValues,&newValues,&newCols, &newRows);
    
    for (int i = 0; i < numRows; i++)
        if (this->rows[this->numRows + i] >= 0)
          this->rows[this->numRows + i] += this->numValues;
          
    
    this->numValues += numValues;
    this->numRows += numRows;
    
    return *this;
}


sparseMatrix& sparseMatrix::dcat(const sparseMatrix& M) {
    
    this->values = (int*) realloc(this->values, (this->numValues + M.length())*sizeof(int));
    this->cols = (int*) realloc(this->cols, (this->numValues + M.length())*sizeof(int));
    this->rows = (int*) realloc(this->rows, (this->numRows + M.size(1))*sizeof(int));
    
    int numCols,numRows,numValues;
    int  * newValues = &(this->values[this->numValues]), * newCols = &(this->cols[this->numValues]),* newRows = &(this->rows[this->numRows]);
    M.decompose(numRows,numCols,numValues,&newValues,&newCols, &newRows);
    
    for (int i = 0; i < numRows; i++)
        if (this->rows[this->numRows + i] >= 0)
          this->rows[this->numRows + i] += this->numValues;

    for (int i = 0 ; i < numValues; i++)
        this->cols[this->numValues + i] += this->numCols;    
    
    this->numValues += numValues;
    this->numCols += numCols;
    this->numRows += numRows;
    
    return *this;
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
  int * seq = (int *) malloc(max(this->numCols,this->numRows)*sizeof(int));
  int * ones = (int *) malloc(max(this->numCols,this->numRows)*sizeof(int));
  
  for (int i = 0 ; i <max(this->numCols,this->numRows); i++) {
    seq[i] = i;
    ones[i] = 1;
  }
  
  rowPerm = sparseMatrix(this->numRows, this->numRows, this->numRows, ones,seq,seq);
  colPerm = sparseMatrix(this->numCols, this->numCols, this->numCols, ones,seq,seq);
  
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
    U.values = (int *) realloc(U.values,(U.numValues + M_numValuesInRow)*sizeof(int));
    U.cols = (int *) realloc(U.cols,(U.numValues + M_numValuesInRow)*sizeof(int));
    U.rows[k] = U.numValues;
    U.numValues += M_numValuesInRow;
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

    L.values = (int *) realloc(L.values,(L.numValues + M_numValuesInCol)*sizeof(int));
    L.cols = (int *) realloc(L.cols,(L.numValues + M_numValuesInCol)*sizeof(int));
    L.rows[k] = L.numValues;
    L.numValues += M_numValuesInCol;
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
  D = sparseMatrix(k,k,k,diagonal,seq,seq);

  L.transpose();
  return *this;
}

const sparseMatrix& sparseMatrix::LDU(sparseMatrix& L, sparseMatrix& D, sparseMatrix& U, sparseMatrix& rowPerm, sparseMatrix& colPerm) const {
  
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
  
  // ho fem de manera poc efficient
  // passant a matrius enteres i despres tornant
  
  int * M = (int *) malloc(this->numRows*this->numCols*sizeof(int));
  int max_rank = min(this->numRows, this->numCols);
  int * l = (int *) malloc(this->numRows*max_rank*sizeof(int));
  int * u = (int *) malloc(max_rank*this->numCols*sizeof(int));
  int * d = (int *) malloc(max_rank*sizeof(int));
  int * rowperm = (int *) malloc(this->numRows*sizeof(int));
  int * colperm = (int *) malloc(this->numCols*sizeof(int));
  int * seq = (int *) malloc(max(this->numCols,this->numRows)*sizeof(int));
  int * ones = (int *) malloc(max(this->numCols,this->numRows)*sizeof(int));
  
  this->get(M);
  for (int i = 0 ; i <max(this->numCols,this->numRows); i++) {
    seq[i] = i;
    ones[i] = 1;
  }
  
  int k = this->LDU_full(this->numRows,this->numCols,M,l,d,u,rowperm,colperm);
  U = sparseMatrix(u,k,this->numCols);
  L = sparseMatrix(l,k,this->numRows);    
  L.transpose();
  D = sparseMatrix(k,k,k,d,seq,seq);
  rowPerm = sparseMatrix(this->numRows, this->numRows, this->numRows, ones,rowperm,seq);
  colPerm = sparseMatrix(this->numCols, this->numCols, this->numCols, ones,colperm,seq);
  
  
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
  
  
  sparseMatrix V(cnt_rows,U.numCols,cnt_nonNullValues,v,c,r),W;
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
