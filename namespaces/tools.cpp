#include "tools.h"

namespace Tools {
  
  int gcd(int n, int * v, int * c) {

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
        
      int g = Tools::gcd(n,v2,c);
      
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
    
      int g = Tools::gcd(n,r,c);
   
      if ( c != NULL)
        for (int i = 0; i < n ; i++)
          if (i != indexMin) {
            c[indexMin] -= c[i]*((v[i]-r[i])/v[indexMin]);
          }
        
      return g;  
    }
        
    return 0;
  }

  ////////////////////////////////////////
  //
  // Block operations
  //
  ////////////////////////////////////////

  bool leq(int n, int const * const a, int const * const b ) {
    for (int i = 0 ; i < n ; i++)
      if (a[i] != b[i])
         return a[i] < b[i];
    return true;
  }

  bool eq(int n, int const * const a, int const * const b ) {
    for (int i = 0 ; i < n ; i++)
      if (a[i] != b[i])
         return false;
    return true;
  }

  void swap(int n, int * const a, int * const b) {
    int c;
    for ( int i = 0 ; i < n ; i++) {
      c = a[i];
      a[i] = b[i];
      b[i] = c;
    }
  }




  int mergeSortBlocks(int n, int m, int * A, bool deleteRepetitions) {
  
    if (m < 2)
      return m ;
    
    int m1 = m/2, m2 = m - m1;
    int * a = (int *) malloc(n*m1*sizeof(int));
    int * b = (int *) malloc(n*m2*sizeof(int));
    memcpy(a,A,n*m1*sizeof(int));
    memcpy(b,&A[n*m1],n*m2*sizeof(int));

    m1 = Tools::mergeSortBlocks(n, m1, a, deleteRepetitions);
    m2 = Tools::mergeSortBlocks(n, m2, b, deleteRepetitions);
    
    // merge
    
    int repetitions = 0;
    
    for (int i = 0 , j = 0; i < m1 || j < m2;)
      if (deleteRepetitions && i < m1 && j < m2 && Tools::eq(n,&a[n*i],&b[n*j]) ) {
        memcpy(&A[n*(i+j-repetitions)],&a[n*i],n*sizeof(int));
        i++;
        j++;
        repetitions++;
      } else if ( i < m1 && (j >= m2 || Tools::leq(n,&a[n*i],&b[n*j]) )) {
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



}