/*
 *  Class: simplicialChainComplex
 *  Author: Daniel Torres
 *  Filename: simplicialChainComplex.h
 * 
 * 
 * 
 *  Notes: 
 * 
 *
 */

#ifndef H_SIMPLICIALCHAINCOMPLEX
#define H_SIMPLICIALCHAINCOMPLEX

#include <iostream>
#include <cstdlib> 
#include <fstream>
#include <string.h>
#include "simplicialPolyhedron.h"
#include "sparseMatrix.h"
using namespace std;

class simplicialChainComplex {

  private:
  
  simplicialPolyhedron * P;
  sparseMatrix * d;
  int * orientation;
  int n;

  private:
  
  bool leq(int n, int const * const a, int const * const b) const;
  bool eq(int n, int const * const a, int const * const b ) const;
  void swap(int n, int * const a, int * const b) const;
  int mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions = false) const;
  int bubbleSort(int n, int * const v) const;
  void vcat(int n1, int n2, int m, int * const C, int const * const A, int const * const B) const;
  void vsplit(int n, int m, int * A, int k, int const * const v, int * const B = NULL) const;
    
  public:

  simplicialChainComplex& inflate(const simplicialPolyhedron& P);
  simplicialPolyhedron deflate(int * sign = NULL) const;
  
  simplicialChainComplex();
  ~simplicialChainComplex();
  simplicialChainComplex(const simplicialPolyhedron& P);
  simplicialChainComplex(const simplicialChainComplex& CC);
  simplicialChainComplex& operator=(const simplicialChainComplex& CC);
  
  int length(int i) const;
  int dim() const;
  simplicialChainComplex& read(istream& in); 
  const simplicialChainComplex& print(ostream & out) const;

  int eulerCharacteristic() const;
  sparseMatrix fundamentalClass() const;
  
  simplicialPolyhedron& operator[](int i);
  const simplicialPolyhedron& operator[](int i) const;
  
  sparseMatrix& boundaryOperator(int i);
  const sparseMatrix& boundaryOperator(int i) const;
  
  sparseMatrix boundary(int i,const sparseMatrix& M) const;
  sparseMatrix adjacencyMatrix(int i, int j) const;
  
  simplicialPolyhedron support(int i, const sparseMatrix& M) const;
  
  //sparseMatrix cup(int k, const sparseMatrix& M, int l, const sparseMatrix& N) const;
  
  sparseMatrix flat(int i, const sparseMatrix & M) const;
  
};

#endif
	
