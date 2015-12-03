/*
 *  Class: simplicialPolyhedron
 *  Author: Daniel Torres
 *  Filename: simplicialPolyhedron.h
 * 
 * 
 * 
 *  Notes: 
 * 
 *
 */

#ifndef H_SIMPLICIALPOLYHEDRON
#define H_SIMPLICIALPOLYHEDRON

#include <iostream>
#include <cstdlib> 
#include <fstream>
#include <string.h>
using namespace std;

class simplicialPolyhedron {

  private:
  
  int * A;
  int n;
  int m;

  private:
  
  bool leq(int n, int const * const a, int const * const b) const;
  bool eq(int n, int const * const a, int const * const b ) const;
  void swap(int n, int * const a, int * const b) const;
  int mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions = false) const;
  int bubbleSort(int n, int * const v) const;
  void vcat(int n1, int n2, int m, int * const C, int const * const A, int const * const B) const;
  void vsplit(int n, int m, int * A, int k, int const * const v, int * const B = NULL) const;
  
  public:
  
  simplicialPolyhedron();
  simplicialPolyhedron(const simplicialPolyhedron& A);
  simplicialPolyhedron(int n, int m, int const * const A);
  ~simplicialPolyhedron();
  simplicialPolyhedron& operator=(const simplicialPolyhedron& A);
  
  int length() const;
  int dim() const;
  const int * values() const;
  simplicialPolyhedron& read(istream & in);
  const simplicialPolyhedron& print(ostream & out) const;
  
  simplicialPolyhedron& sortSimplexes(int * signs = NULL);
  simplicialPolyhedron& simplifySimplexes(int * coeffs = NULL);
  simplicialPolyhedron boundary(int * coeffs2 = NULL, int * coeffs = NULL) const;
  simplicialPolyhedron skeleton() const;
  
  simplicialPolyhedron operator[](int k) const;
  const simplicialPolyhedron& binarySearch(const simplicialPolyhedron& P, int * const v, int * a = NULL, int * b = NULL) const;
  simplicialPolyhedron& remove(int n, int * simplexes);
  simplicialPolyhedron& remove(int i);

  
};

#endif
	
