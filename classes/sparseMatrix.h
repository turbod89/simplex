/*
/*
 *  Class: sparseMatrix
 *  Author: Daniel Torres
 *  Filename: sparseMatrix.h
 * 
 * 
 * 
 *  Notes: 
 * 
 *
 */

#ifndef H_SPARSEMATRIX
#define H_SPARSEMATRIX

#include <iostream>
#include <cstdlib>
#include <string.h>
using namespace std;

class sparseMatrix {

  private:

    int numCols;
    int numRows;
    int * values;
    int * cols;
    int * rows;
    
    void swap(int &a, int &b) const;
    int binary_search_position(int length, const int * const v, int needle) const;
    int binary_search(int length, const int * const v, int needle) const;
    void mergeSort(int length, int * v,int auxLength = 0, int * const * const aux = NULL) const;
    int gcd(int n, int * v, int * c = NULL) const;
    int LDU_full(int n, int m, int * const M, int * const L, int * const D, int * const U, int * const rowPerm, int * const colPerm) const;
    int minAbs(int n, int const * const v) const;
    int maxAbs(int n, int const * const v) const;
    int minAbsNotZero(int n, int const * const v) const;
    void multiply(int a, int n, int * const v) const;
    int multiply(int n, int const * const v, int const * const w) const;
    int multiplyRows(int n1, int const * const c1, int const * const v1, int n2, int const * const c2, int const * const v2) const;
    int sumRows(int n1, int const * const c1, int const * const v1, int n2, int const * const c2, int const * const v2, int * const c3, int * const v3) const;

    
  public:
  
    sparseMatrix();
    sparseMatrix(int r, int c,const int * const rows,const int * const cols,const int * const values);
    sparseMatrix(int r, int c);
    sparseMatrix(int * M, int r, int c);
    sparseMatrix(const sparseMatrix& M);
    ~sparseMatrix();
    sparseMatrix& operator=(const sparseMatrix& M);
    const int * getValues() const;
    int * getValues();
    const int * getRows() const;
    int * getRows();
    const int * getCols() const;
    int * getCols();
    int size(int i) const;
    inline int length() const;
    const sparseMatrix& decompose(int ** values, int ** cols, int ** rows) const;
    sparseMatrix& read_index_format(int numRows, int numCols, int length, int const * rows, int const * cols, int const * values, bool sorted = false);
    const sparseMatrix& index_format(int * rows = NULL, int * cols = NULL, int * values = NULL) const;
    sparseMatrix& removeZeros();
    sparseMatrix& read(istream& in);
    const sparseMatrix& print(ostream& out) const;
    const sparseMatrix& print_full(ostream& out) const;
    const sparseMatrix& print_octave(ostream& out) const;
    const sparseMatrix& getFullValues(int * a) const;
    inline int numValuesInRow(int row) const;
    sparseMatrix& swapRows(int row1, int row2);
    sparseMatrix& swapCols(int col1, int col2);
    sparseMatrix& deleteRows(int n, int * const v);
    sparseMatrix& deleteRow(int i);
    sparseMatrix& eye(int r, int c);
    sparseMatrix& eye(int n);
    sparseMatrix transpose() const;
    sparseMatrix& multiplyByTransposed(sparseMatrix const & M1,sparseMatrix const & M2);
    sparseMatrix operator*(const sparseMatrix & M) const;
    sparseMatrix operator*(int a) const;
    sparseMatrix operator+(const sparseMatrix& M) const;
    sparseMatrix operator[](int row) const;
    int operator()(int row, int col) const;
    const sparseMatrix& LDU_efficient(sparseMatrix& L, sparseMatrix& D, sparseMatrix& U, sparseMatrix& rowPerm, sparseMatrix& colPerm) const;
    sparseMatrix ker() const;
    sparseMatrix LXeqY(const sparseMatrix &Y) const;
    sparseMatrix LComplementary(const sparseMatrix &Y) const;
};

inline int sparseMatrix::length() const {
  return this->rows[this->numRows];
}

inline int sparseMatrix::numValuesInRow(int row) const {

  ////////////////
  //
  // return the number of non null elements in a row
  //
  ////////////////
  
  return this->rows[row+1] - this->rows[row];

}

#endif
    
