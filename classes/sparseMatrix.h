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
    int numValues;
    int * values;
    int * cols;
    int * rows;
    
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

    
  public:
  
    sparseMatrix();
    sparseMatrix(int r, int c, int v,const int * const values,const int * const cols,const int * const rows);
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
    int length() const;
    const sparseMatrix& decompose(int ** values, int ** cols, int ** rows) const;
    const sparseMatrix& format2(int ** rows = NULL, int ** cols = NULL, int ** values = NULL) const;
    sparseMatrix& read(istream& in);
    const sparseMatrix& print(ostream& out) const;
    const sparseMatrix& print_full(ostream& out) const;
    const sparseMatrix& getFullValues(int * a) const;
    inline int numValuesInRow(int row) const;
    sparseMatrix& swapRows(int row1, int row2);
    sparseMatrix& swapCols(int col1, int col2);
    sparseMatrix& deleteCols(int n, int * cols);
    sparseMatrix& deleteCol(int col);
    sparseMatrix& eye(int r, int c);
    sparseMatrix& eye(int n);
    int emptyRowsToBottom(sparseMatrix& rowPerm);
    sparseMatrix& transpose();
    sparseMatrix& multiplyByTransposed(sparseMatrix const & M1,sparseMatrix const & M2);
    sparseMatrix operator*(sparseMatrix M) const;
    sparseMatrix operator*(int a) const;
    sparseMatrix operator+(const sparseMatrix& M) const;
    sparseMatrix operator[](int row) const;
    int operator()(int row, int col) const;
    sparseMatrix& vcat(const sparseMatrix& M);
    sparseMatrix& dcat(const sparseMatrix& M);
    const sparseMatrix& LDU_efficient(sparseMatrix& L, sparseMatrix& D, sparseMatrix& U, sparseMatrix& rowPerm, sparseMatrix& colPerm) const;
    const sparseMatrix& LDU(sparseMatrix& L, sparseMatrix& D, sparseMatrix& U, sparseMatrix& rowPerm, sparseMatrix& colPerm) const;
    sparseMatrix ker() const;
};

#endif
    
