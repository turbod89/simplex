#include "simplicialChainComplex.h"

bool simplicialChainComplex::leq(int n, int const * const a, int const * const b ) const {
	for (int i = 0 ; i < n ; i++)
		if (a[i] != b[i])
		   return a[i] < b[i];
	return true;
}

bool simplicialChainComplex::eq(int n, int const * const a, int const * const b ) const {
	for (int i = 0 ; i < n ; i++)
		if (a[i] != b[i])
		   return false;
	return true;
}

void simplicialChainComplex::swap(int n, int * const a, int * const b) const {
	int c;
	for ( int i = 0 ; i < n ; i++) {
		c = a[i];
		a[i] = b[i];
		b[i] = c;
	}
}

int simplicialChainComplex::mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions) const{
  
  if (m < 2)
    return m ;
  
  int m1 = m/2, m2 = m - m1;
  int * a = (int *) malloc(n*m1*sizeof(int));
  int * b = (int *) malloc(n*m2*sizeof(int));
  memcpy(a,A,n*m1*sizeof(int));
  memcpy(b,&A[n*m1],n*m2*sizeof(int));

  m1 = this->mergeSortBlocks(n, m1, a, deleteRepetitions);
  m2 = this->mergeSortBlocks(n, m2, b, deleteRepetitions);
  
  // merge
  
  int repetitions = 0;
  
  for (int i = 0 , j = 0; i < m1 || j < m2;)
    if (deleteRepetitions && i < m1 && j < m2 && eq(n,&a[n*i],&b[n*j]) ) {
      memcpy(&A[n*(i+j-repetitions)],&a[n*i],n*sizeof(int));
      i++;
      j++;
      repetitions++;
    } else if ( i < m1 && (j >= m2 || leq(n,&a[n*i],&b[n*j]) )) {
      memcpy(&A[n*(i+j-repetitions)],&a[n*i],n*sizeof(int));
      i++;
    } else {
      memcpy(&A[n*(i+j-repetitions)],&b[n*j],n*sizeof(int));
      j++;
    }
  
  repetitions += m - m1 - m2;
  
  if (repetitions > 0)
    A = (int *) realloc(A,n*(m-repetitions)*sizeof(int));
  
  return m - repetitions;
}

int simplicialChainComplex::bubbleSort(int n, int * const v) const {

  if (n < 1)
    return 0;

  int sign = 1;
  
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n ; j++)
      if ( v[i] > v[j] ) {
        int a = v[i];
        v[i] = v[j];
        v[j] = a;
        sign *= -1;
      }
  
  return sign;
}

void simplicialChainComplex::vcat(int n1, int n2, int m, int * const C, int const * const A, int const * const B) const {
  for (int j = m; j > 0; j--) {
    memcpy(&C[(n1+n2)*(j-1)],&A[n1*(j-1)], n1*sizeof(int));
    memcpy(&C[(n1+n2)*(j-1) + n1],&B[n2*(j-1)], n2*sizeof(int));
  }
}


void simplicialChainComplex::vsplit(int n, int m, int * A, int k, int const * const v, int * const B) const {

  int delay = 0;
  for (int j = 0 ; j < m ; j++) {
    int c = 0;
    for (int i = 0; i < n; i++) {
      while (c < k && v[c] < i)
        c++;
        
      if (c < k && i == v[c]) {
        A[j*n + i - delay] = A[j*n + i];
      } else {
        if (B != NULL)
          B[delay] = A[j*n + i];
        delay++;
      }
    }
  }
}





simplicialChainComplex& simplicialChainComplex::inflate(const simplicialPolyhedron& P) {
  this->n = P.dim()+1;
  this->P = (simplicialPolyhedron *) malloc(this->n*sizeof(simplicialPolyhedron));
  this->d = (sparseMatrix *) malloc((this->n-1)*sizeof(sparseMatrix));
  this->orientation = (int *) malloc(P.length()*sizeof(int));
  for (int i = 0; i < P.length() ; i++)
    this->orientation[i] = 1;
  
  // skeleton

  this->P[this->n-1] = P;
  this->P[this->n-1].simplifySimplexes(this->orientation);
  
  for (int i = this->n-1; i > 0 ; i--)
    this->P[i - 1] = this->P[i].skeleton();
    
  // boundary operators

  for (int i = 0; i < this->n - 1; i++) {
    int numValues = (i+2)*this->P[i+1].length();
    int * values = (int *) malloc(numValues*sizeof(int));
    int * cols = (int *) malloc(numValues*sizeof(int));
    int * rows = (int *) malloc((this->P[i+1].length()+1)*sizeof(int));
    
    for (int j = 0; j < this->P[i+1].length(); j++) {
      
      int * signs = (int *) malloc((i+2)*sizeof(int));
      for (int k = 0; k < (i+2) ; k++)
        signs[k] = 1;
        
      simplicialPolyhedron S = this->P[i+1][j].boundary(signs);
      memcpy(&values[(i+2)*j],signs,(i+2)*sizeof(int));
      
      this->P[i].binarySearch(S,&cols[(i+2)*j]);
      for (int k = 0 ; k < S.length() ; k++)
        if (cols[(i+2)*j + k] < 0 ) {
          cerr << "Error!! facet not found at skeleton" << endl;
          cerr << "facet:" << endl;
          S[k].print(cerr);
          cerr << "of simplex:" << endl;
          this->P[i+1][j].print(cerr);
          cerr << "searched in:" << endl;
          this->P[i].print(cerr);
        }
    }

    for (int j = 0; j <= this->P[i+1].length(); j++)
      rows[j] = (i+2)*j;

    this->d[i] = sparseMatrix (this->P[i+1].length(),this->P[i].length(),rows,cols,values);
    this->d[i] = this->d[i].transpose();
  }

  return *this;
}

simplicialPolyhedron simplicialChainComplex::deflate(int * sign) const {
  if (sign != NULL) {
    for (int i = 0; i < this->P[this->n-1].length(); i++)
      sign[i] = this->orientation[i];
  }
  
  return this->P[n-1];
}

simplicialChainComplex::simplicialChainComplex() {
  this->n = 0;
  this->P = NULL;
  this->d = NULL;
  this->orientation = NULL;
}

simplicialChainComplex::~simplicialChainComplex() {
  free(this->P);
  free(this->d);
  free(this->orientation);
}

simplicialChainComplex::simplicialChainComplex(const simplicialPolyhedron& P) {
  this->inflate(P);
}

simplicialChainComplex::simplicialChainComplex(const simplicialChainComplex& CC) {
  *this = CC;
}

simplicialChainComplex& simplicialChainComplex::operator=(const simplicialChainComplex& CC) {
  this->n = CC.n;
  this->P = (simplicialPolyhedron *) malloc(this->n*sizeof(simplicialPolyhedron));
  this->d = (sparseMatrix *) malloc((this->n-1)*sizeof(sparseMatrix));
  this->orientation = (int *) malloc(CC.P[CC.n - 1].length()*sizeof(int));
  
  for (int i = 0; i < this->n ; i++)
    this->P[i] = simplicialPolyhedron(CC.P[i]);
    
  for (int i = 0; i < this->n -1 ; i++)
    this->d[i] = sparseMatrix(CC.d[i]);
  
  for (int i = 0; i < CC.P[CC.n -1].length(); i++)
    this->orientation[i] = CC.orientation[i];
  
  return *this;
}

simplicialChainComplex& simplicialChainComplex::read(istream& in) {

  in >> this->n;
  this->P = (simplicialPolyhedron *) malloc(this->n*sizeof(simplicialPolyhedron));
  this->d = (sparseMatrix *) malloc((this->n-1)*sizeof(sparseMatrix));
  
  for (int i = 0; i < this->n ; i++)
    this->P[i].read(in);

  for (int i = 0; i < this->n - 1  ; i++)
    this->d[i].read(in);
  
  int n;
  in >> n;
  for (int i = 0; i < n; i++)
    in >> this->orientation[i];

  return *this;
} 

const simplicialChainComplex& simplicialChainComplex::print(ostream & out) const {
  
  out << this->n << endl;
  
  for (int i = 0; i < this->n ; i++)
    this->P[i].print(out);
  
  for (int i = 0; i < this->n - 1; i++)
    this->d[i].print(out);
  
  out << this->P[this->n-1].length() << endl;
  for (int i = 0; i < this->P[this->n-1].length(); i++)
    out << " " << this->orientation[i];
  out << endl;

  
  return *this;
}




int simplicialChainComplex::dim() const {
  return this->n - 1;
}

int simplicialChainComplex::length(int i) const {
  return this->P[i].length();
}

int simplicialChainComplex::eulerCharacteristic() const {
  int s = 0, sign = 1;
  for (int i = 0; i < this->n; i++) {
    s += sign*(this->P[i].length());
    sign *= -1;
  }
    
  return s;
}

sparseMatrix simplicialChainComplex::fundamentalClass() const {
	int * values = (int *) malloc(this->P[this->n -1].length() * sizeof(int));
	int * rows = (int *) malloc(this->P[this->n -1].length() * sizeof(int));
	int * cols = (int *) malloc(2*sizeof(int));
	
  cols[0] = 0;
	for (int i = 0; i < this->P[this->n -1].length(); i++) {
		values[i] = this->orientation[i];
		rows[i] = i;
	}
  cols[1] = this->P[this->n -1].length();
	
	return sparseMatrix(1,this->P[n-1].length(),cols,rows,values);
	
}

simplicialPolyhedron& simplicialChainComplex::operator[](int i) {

  if (i < 0 || i >= n)
    cerr << "Value " << i << " is out of range 0 to max dim = " << this->n-1 << endl;

  return this->P[i];
}
  
const simplicialPolyhedron& simplicialChainComplex::operator[](int i) const {
  if (i < 0 || i >= n)
    cerr << "Value " << i << " is out of range 0 to max dim = " << this->n-1 << endl;
  return this->P[i];
}

sparseMatrix& simplicialChainComplex::boundaryOperator(int i) {
  if (i < 1 || i >= n)
    cerr << "Value " << i << " is out of range 1 to max dim = " << this->n-1 << endl;

  return this->d[i-1];
}

const sparseMatrix& simplicialChainComplex::boundaryOperator(int i) const {
  if (i < 1 || i >= n)
    cerr << "Value " << i << " is out of range 1 to max dim = " << this->n-1 << endl;

  return this->d[i-1];
}

sparseMatrix simplicialChainComplex::adjacencyMatrix(int i, int j) const {

  if (i > j)
    return this->adjacencyMatrix(j,i).transpose();

  if (i < 0 || i >= this->n)
    cerr << "Value (first) =" << i << " is out of range 0 to max dim = " << this->n-1 << endl;

  if (j < 0 || j >= this->n)
    cerr << "Value (second) =" << j << " is out of range 0 to max dim = " << this->n-1 << endl;
    
  if (i == j) {
    sparseMatrix M;
    M.eye(this->P[i].length());
  }
  
  sparseMatrix M(this->d[j-1]);
  for (int l = 0 ; l < M.length() ; l++)
    M.getValues()[l] = 1;
    
  if (i == j-1)
    return M; // we save one matrix multiplication
  
  M = this->adjacencyMatrix(i,j-1)*M;
 
  for (int l = 0 ; l < M.length() ; l++)
    M.getValues()[l] = 1;
  
  return M;

}

sparseMatrix simplicialChainComplex::boundary(int i,const sparseMatrix& M) const {
  if (i < 1 || i >= n)
    return sparseMatrix(0,0);

  sparseMatrix A;
  return A.multiplyByTransposed(M,this->d[i-1]);
}

simplicialPolyhedron simplicialChainComplex::support(int i, const sparseMatrix& M) const {
  if (i < 1 || i >= this->n)
    cerr << "Value " << i << " is out of range 0 to max dim = " << this->n-1 << endl;
    
  if (M.size(2) != this->P[i].length() || M.size(1) != 1)
    cerr << "Chain does not correspond to a " << i << "-chain" << endl;
    
  int * A = (int *) malloc((i+1)*M.length()*sizeof(int));
  
  int * values , *rows, *cols;
  int numCols = M.size(2), numRows = M.size(1), numValues = M.length();
  M.decompose(&values,&cols,&rows);
  
  for (int j = 0; j < numValues; j++)
    memcpy(&A[(i+1)*j],&(this->P[i].values()[(i+1)*cols[j]]),(i+1)*sizeof(int));
  
  return simplicialPolyhedron(i,numValues,A);
  
}

/*sparseMatrix simplicialChainComplex::cup(int k, const sparseMatrix& M, int l, const sparseMatrix& N) const {
  
  if ( k < 0 || l < 0 || k+l > this->n -1 )
    return M(0,0);
  
  sparseMatrix A = this->;

}
*/

sparseMatrix simplicialChainComplex::flat(int k, const sparseMatrix & M) const {

  int * A = (int *) malloc((k+1)*M.length()*sizeof(int));
  
  // copy k-simplexes to A
  for ( int i = 0; i < M.length(); i++)
    memcpy(&A[(k+1)*i],&this->P[k].values()[(k+1)*M.getCols()[i]],(k+1)*sizeof(int));

  simplicialPolyhedron Q(k,M.length(),A);
  free(A);
  
  // search them on maximal simplex
  
  int * start = (int *) malloc(Q.length()*sizeof(int));
  int * end = (int *) malloc(Q.length()*sizeof(int));
  
  this->P[this->dim()].subSearch(Q,start,end);
  
  // get complementary

  int l = 0;
  for (int i = 0; i < Q.length(); i++)
    l += end[i] - start[i];
  A = (int *) malloc((this->dim()-k+1)*l*sizeof(int));
  int * signs = (int *) malloc(l*sizeof(int));
  
  for (int i = 0, cnt = 0 ; i < Q.length(); i++)
    for (int j = start[i]; j < end[i]; j++, cnt++) {
      signs[cnt] = M.getValues()[i] * this->orientation[j];
      memcpy(&A[(this->dim()-k+1)*cnt], &this->P[this->dim()].values()[(this->dim()+1)*j + k], (this->dim()-k+1)*sizeof(int));
    }
  
  simplicialPolyhedron R(this->dim() - k, l, A);
  free(A);
  free(start);
  free(end);
  R.simplifySimplexes(signs);
  
  // search them on (dim - k) - simplexes
  
  start = (int *) malloc(R.length()*sizeof(int));
  end = (int *) malloc(R.length()*sizeof(int));
  
  this->P[this->dim() - k].subSearch(R,start,end);
  
  // put'em on a matrix
  
  l = 0;
  for (int i = 0; i < Q.length(); i++)
    l += end[i] - start[i];
  
  int * rows = (int *) malloc(2*sizeof(int));
  rows[0] = 0;
  rows[1] = l;
  int * cols = (int *) malloc(l*sizeof(int));
  int * values = (int *) malloc(l*sizeof(int));
  
  for (int i = 0, cnt = 0 ; i < R.length(); i++)
    for (int j = start[i]; j < end[i]; j++, cnt++) {
      cols[cnt] = j;
      values[cnt] = signs[i];
    }
    
  sparseMatrix N(1,this->P[this->dim() - k].length(), rows, cols, values);
  N.removeZeros();
  
  return N;

}

sparseMatrix simplicialChainComplex::getHomology(int i) const {

  if (i == this->dim()) {

    return this->d[i-1].ker().transpose();

  } else if (i == 0) {

    sparseMatrix L1, D1, U1, P1, Q1;

    this->d[0].transpose().LDU_efficient(L1,D1,U1,P1,Q1);

    sparseMatrix Q1P0KerL0t = Q1;
    sparseMatrix ImU1t = U1.transpose();
    sparseMatrix X = Q1.transpose()*ImU1t.LComplementary(Q1P0KerL0t);

    return X.transpose();

  } else if ( i < this->dim() && i > 0) {
  
    sparseMatrix L0, D0, U0, P0, Q0;
    sparseMatrix L1, D1, U1, P1, Q1;

    this->d[i-1].transpose().LDU_efficient(L0,D0,U0,P0,Q0);
    this->d[i].transpose().LDU_efficient(L1,D1,U1,P1,Q1);

    sparseMatrix Q1P0KerL0t = Q1*P0*(L0.transpose().ker());
    sparseMatrix ImU1t = U1.transpose();
    sparseMatrix X = Q1.transpose()*ImU1t.LComplementary(Q1P0KerL0t);

    return X.transpose();

  }

  return sparseMatrix();

}