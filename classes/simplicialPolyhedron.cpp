#include "simplicialPolyhedron.h"

bool simplicialPolyhedron::leq(int n, int const * const a, int const * const b ) const {
	for (int i = 0 ; i < n ; i++)
		if (a[i] != b[i])
		   return a[i] < b[i];
	return true;
}

bool simplicialPolyhedron::eq(int n, int const * const a, int const * const b ) const {
	for (int i = 0 ; i < n ; i++)
		if (a[i] != b[i])
		   return false;
	return true;
}

void simplicialPolyhedron::swap(int n, int * const a, int * const b) const {
	int c;
	for ( int i = 0 ; i < n ; i++) {
		c = a[i];
		a[i] = b[i];
		b[i] = c;
	}
}

int simplicialPolyhedron::mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions) const{
  
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

int simplicialPolyhedron::bubbleSort(int n, int * const v) const {

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

void simplicialPolyhedron::vcat(int n1, int n2, int m, int * const C, int const * const A, int const * const B) const {
  for (int j = m; j > 0; j--) {
    memcpy(&C[(n1+n2)*(j-1)],&A[n1*(j-1)], n1*sizeof(int));
    memcpy(&C[(n1+n2)*(j-1) + n1],&B[n2*(j-1)], n2*sizeof(int));
  }
}


void simplicialPolyhedron::vsplit(int n, int m, int * A, int k, int const * const v, int * const B) const {

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






simplicialPolyhedron::simplicialPolyhedron() {
  this->n = 0;
  this->m = 0;
  this->A = NULL;
}

simplicialPolyhedron::simplicialPolyhedron(const simplicialPolyhedron& P) {
  this->n = P.dim() + 1;
  this->m = P.length();
  this->A = (int *) malloc(this->n*this->m*sizeof(int));
  memcpy(this->A,P.values(),this->n*this->m*sizeof(int));
}

simplicialPolyhedron::simplicialPolyhedron(int n, int m, int const * const A) {
  this->n = n + 1;
  this->m = m;
  this->A = (int *) malloc((n+1)*m*sizeof(int));
  memcpy(this->A,A,(n+1)*m*sizeof(int));
}

simplicialPolyhedron::~simplicialPolyhedron() {
  free(this->A);
}


simplicialPolyhedron& simplicialPolyhedron::operator=(const simplicialPolyhedron& P) {
  this->n = P.dim() + 1;
  this->m = P.length();
  this->A = (int *) realloc(this->A,this->n*this->m*sizeof(int));
  memcpy(this->A,P.values(),this->n*this->m*sizeof(int));
  return *this;
}







int simplicialPolyhedron::dim() const {
  return this->n - 1;
}

int simplicialPolyhedron::length() const {
  return this->m;
}

const int * simplicialPolyhedron::values() const {
  return this->A;
};

simplicialPolyhedron& simplicialPolyhedron::read(istream & in) {
  int n, m;
  in >> n >> m;
  n--;
  int * A = (int *) malloc((n+1)*m*sizeof(int));
  
  for (int i = 0; i < n+1; i++)
    for (int j = 0; j < m; j++)
      in >> A[ (n+1)*j + i];

  this->n = n + 1;
  this->m = m;
  this->A = (int *) malloc((n+1)*m*sizeof(int));
  memcpy(this->A,A,(n+1)*m*sizeof(int));

  return *this;

}

const simplicialPolyhedron& simplicialPolyhedron::print(ostream & out) const {
  
  out << this->n << " " << this->m << endl;
  
  for (int i = 0 ; i < this->n; i++) {
    for (int j = 0; j < this->m; j++)
      out << " " << this->A[j*this->n + i];
    out << endl;
  }
  
  return *this;
}



simplicialPolyhedron& simplicialPolyhedron::sortSimplexes(int * signs) {
  if (signs != NULL) 
    for (int i = 0; i < m ; i++)
      signs[i] *= this->bubbleSort(this->n,&(this->A[n*i]));
  else
    for (int i = 0; i < m ; i++)
      this->bubbleSort(this->n,&(this->A[n*i]));
      
  return *this;
}


simplicialPolyhedron& simplicialPolyhedron::simplifySimplexes(int * coeffs) {

  if ( this->n <= 0)
    return *this;
    
  if (coeffs == NULL) {
    coeffs = (int *) malloc(this->m*sizeof(int));
    for (int i = 0; i < this->m ; i++)
      coeffs[i] = 1;
  }

  int * p = (int *) malloc(this->m*sizeof(int));
  for (int i = 0; i < this->m; i++)
    p[i] = i;

  int * Ap = (int *) malloc((this->n+1)*this->m*sizeof(int));
  
  // sort simplexes  
  this->sortSimplexes(coeffs);
  // sort among them and obtain
  // the corresponding permutation
  this->vcat(this->n,1,this->m,Ap,this->A,p);
  this->mergeSortBlocks(this->n+1,this->m,Ap);

  int * v = (int *) malloc(this->n*sizeof(int));
  for (int i = 0; i < this->n; i++)
    v[i] = i;
  
  this->vsplit(this->n+1,this->m,Ap,this->n,v,p);
    
  // count number of different simplexes
  int m2 = 0; // number of DIFERENT simplexes
  for (int i = 0; i < this->m; i++)
    if ( i == 0 || !this->eq(this->n,&Ap[this->n*(i-1)],&Ap[this->n*i]) )
      m2++;
  
  // group by

  int * c = (int *) malloc(this->m*sizeof(int));
  memcpy(c,coeffs,this->m*sizeof(int));
  
  this->A = (int *) realloc(this->A,this->n*m2*sizeof(int));
  if (m2 < this->m)
    coeffs = (int *) realloc(coeffs,m2*sizeof(int));
  int cnt = -1;
  for (int i = 0; i < this->m; i++)
    if ( i == 0 || !this->eq(n,&Ap[this->n*(i-1)],&Ap[this->n*i]) ) {
      cnt++;
      memcpy(&A[n*cnt],&Ap[this->n*i],this->n*sizeof(int));
      memcpy(&coeffs[1*cnt],&c[p[i]],1*sizeof(int));
    } else {
      coeffs[1*cnt] += c[p[i]];
    } 
  // delete 0's

  cnt = 0;
  
  for (int i = 0; i < m2; i++)
    if ( coeffs[i] != 0 ) {
      if ( cnt != i) {
        memcpy(&A[this->n*cnt],&A[this->n*i],this->n*sizeof(int));
        memcpy(&coeffs[1*cnt],&coeffs[1*i],1*sizeof(int));
      }
      cnt++;
    }
        
  int m3 = cnt;
  this->A = (int *) realloc(this->A,this->n*m3*sizeof(int));
  if (m3 < m2)
    coeffs = (int *) realloc(coeffs,1*m3*sizeof(int));

  this->m = m3;

  return *this;
  
}


simplicialPolyhedron simplicialPolyhedron::boundary(int * coeffs2, int * coeffs) const{

  if (this->n == 1) {
    return simplicialPolyhedron();
  }

  if (coeffs == NULL) {
    coeffs = (int *) malloc(this->m*sizeof(int));
    for (int i = 0; i < this->m ; i++)
      coeffs[i] = 1;
  }
  
  if (coeffs2 == NULL) {
    coeffs2 = (int *) malloc(this->n*this->m*sizeof(int));
    for (int i = 0; i < this->m*this->n ; i++)
      coeffs2[i] = 1;
  }
    
  simplicialPolyhedron res((this->n)-1,this->m,this->A);
  res.simplifySimplexes(coeffs);
  int n = res.dim()+1, m = res.length();
  const int * A = res.values();  
  
  int * B = (int *) malloc((n-1)*n*m*sizeof(int));
  int * v = (int *) malloc(n*sizeof(int));
  
  for (int i = 0; i < n - 1; i++)
    v[i] = i + 1;
  
  int sign = 1;
  
  for (int i = 0; i < n ; i++) {
    int * A2 = (int *) malloc(n*m*sizeof(int));
    memcpy(A2,A,n*m*sizeof(int));
    this->vsplit(n,m,A2,n-1,v,NULL);
    memcpy(&B[(n-1)*m*i],A2,(n-1)*m*sizeof(int));
    for (int j = 0; j < m; j++)
      coeffs2[m*i + j] = coeffs[j]*sign;
    if (i < n-1 )
      v[i]--;
    sign *= -1;
    free(A2);
  }

  res = simplicialPolyhedron((n-1)-1,n*m,B);
  res.simplifySimplexes(coeffs2);
  free(B);
  free(v);
  return res;

}

simplicialPolyhedron simplicialPolyhedron::skeleton() const{

  if (this->n == 1) {
    return simplicialPolyhedron();
  }
  
  simplicialPolyhedron res((this->n)-1,this->m,A);
  
  res.simplifySimplexes();
  int n = res.dim()+1, m = res.length();
  const int * A = res.values();
  
  
  int * B = (int *) malloc((n-1)*n*m*sizeof(int));
  
  int * v = (int *) malloc(n*sizeof(int));
  for (int i = 0; i < n - 1; i++)
    v[i] = i + 1;
  
  for (int i = 0; i < n ; i++) {
    int * A2 = (int *) malloc(n*m*sizeof(int));
    memcpy(A2,A,n*m*sizeof(int));
    this->vsplit(n,m,A2,n-1,v,NULL);
    memcpy(&B[(n-1)*m*i],A2,(n-1)*m*sizeof(int));
    if (i < n-1 )
      v[i]--;
  }

  res = simplicialPolyhedron((n-1)-1,n*m,B);
  res.sortSimplexes(); // ordenem *sense coefficients* per evitar cancelacions
  res.simplifySimplexes();  // eliminem duplicats
  
  return res;

}



simplicialPolyhedron simplicialPolyhedron::operator[](int k) const {
  return simplicialPolyhedron(this->n-1,1,&(this->A[k*n]));
}

const simplicialPolyhedron& simplicialPolyhedron::binarySearch(const simplicialPolyhedron& P, int * const v, int * a, int * b) const {

  if (this->n == 0 || this->dim() != P.dim() )
    return *this;
  
  if ( a == NULL || b == NULL) {
    a = (int *) malloc(P.length()*sizeof(int));
    for (int i = 0; i < P.length(); i++)
      a[i] = 0;
    b = (int *) malloc(P.length()*sizeof(int));
    for (int i = 0; i < P.length(); i++)
      b[i] = this->m;
    return this->binarySearch(P,v,a,b);
  }
  
  bool finish = true;
  for (int i = 0; i < P.length() && finish ; i++)
    finish = finish && a[i] >= b[i];
  if (finish)
    return *this;
  
  const int * A = P.values();
  int * c = (int *) malloc(P.length()*sizeof(int));
  for (int i = 0; i < P.length(); i++)
    c[i] = (a[i]+b[i])/2;

  for (int i = 0; i < P.length(); i++) {
	  const int * v1 = &(this->A[this->n*c[i]]);
	  const int * v2 = &(A[this->n*i]);
    if (a[i] >= b[i]) {
      // vam trobar la solucio en una iteracio
      // anterior. No fem res
    } else if ( this->eq(this->n,v1,v2) ) {
  	  // solucio trobada
  	  v[i] = c[i];
  	  a[i] = c[i];
  	  b[i] = c[i];
  	} else if (a[i] == c[i]) {
  	  // solucio no existent
  	  v[i] = -1;
  	  a[i] = c[i];
  	  b[i] = c[i];	  
  	} else if ( this->leq(this->n,v1,v2) ) {
  	  //ens restringim a la dreta
  	  a[i] = c[i];
  	} else {
  	  //ens restringim a l'esquerra
  	  b[i] = c[i];
  	}
  
  }

  return this->binarySearch(P,v,a,b);
}

simplicialPolyhedron& simplicialPolyhedron::remove(int n, int * I) {
  
  this->bubbleSort(n,I);
  
  int numRemoved = 0;
  
  for (int i = 0, k = 0; i < this->m && k < n; )
    if (I[k] == i) {
      for (int j = this->n*(i+1); j < this->n*this->m; j++)
        this->A[j-this->n] = this->A[j];
      numRemoved++;
      k++;
      i++;
    } else if ( I[k] > i) {
      i++;
    } else if (I[k] < i) {
      k++;
    }
    
    this->A = (int *) realloc(this->A,this->n*(this->m - numRemoved)*sizeof(int));
    this->m -= numRemoved;
  
  return *this;
}

simplicialPolyhedron& simplicialPolyhedron::remove(int i) {
  return this->remove(1,&i);
}
