#include <iostream>
#include <cstdlib>
#include "simplicialPolyhedron.h"
#include <string>
#include <vector>
#include <math.h>
using namespace std;

class Manifolds {

  public:

  int n;
  int NUM_MAX;
  int * dim;
  int * num;
  //simplicialPolyhedron * P;
  vector<string> name;
  
  Manifolds() {
    n = 0;
    NUM_MAX = 0;
    this->dim = NULL;
    num = NULL;
    //P = nullptr;
  }
  
  Manifolds(int NUM_MAX) {
    
    this->n = 0;
    this->NUM_MAX = NUM_MAX;
    this->dim = (int *) malloc(NUM_MAX*sizeof(int));
    this->num = (int *) malloc(NUM_MAX*sizeof(int));
    this->name = vector<string>(NUM_MAX);
    //this->P = (int *) malloc(NUM_MAX*sizeof(simplicialPolyhedron));
  
  }
  
  void swapManifolds(int i, int j) {
    swap(this->dim[i], this->dim[j]);
    swap(this->num[i], this->num[j]);
    string a = this->name[i];
    this->name[i] = this->name[j];
    this->name[j] = a;
    
  }
  
  bool le(int i, int j) const {
    if (this->dim[i] != this->dim[j])
      return this->dim[i] < this->dim[j];
    return this->num[i] < this->num[j];
  }
  
  int find (int k, int start = 0, int end = -1) {
  
    if (end < 0)
      end = this->n;

    if (end <= start)
      return -1;
      
    int c = (end + start)/2;
    
    if (this->le(k,c))
      return this->find(k,start,c);
    else if (this->le(c,k))
      return this->find(k,c+1,end);
    else
      return c;
      
    return -1;
  
  }
  
  int insert(int num, int dim, string name) {

    if (this->n >= this->NUM_MAX)
      return -1;

    this->num[this->n] = num;
    this->dim[this->n] = dim;
    this->name[this->n] = name;
    
    int k = this->find(this->n);
    if (k >= 0)
      return k;
    
    int i = this->n;
    while (i > 0 &&  this->le(i,i-1) ) {
      this->swapManifolds(i,i-1);
      i--;
    }

    this->n++;
    return true;
      
  }
  
  const Manifolds& print(ostream & out) const {
  
    for (int i = 0; i < this->n ; i++)
      out << this->num[i] << " " << this->dim[i] << " " << this->name[i] << endl;
  
  }
  
};

long long int choose(int n, int m) {
  
  int k = max(n,m);
  n = min(n,m);  
  m = k;
  n = min(n,m-n);
  long long int res = 1;

  for (int i = 1; i <= n; i++)
    res = res*(m-n+i)/i;
  
  return res;
  
}

int crossAll(Manifolds & list, int j, int MAX_NUM, int MAX_DIM) {
  
  int l = list.n;
  int b = -1;

  for (int i = 0; i < l; i++) {
    if (list.num[i]*list.num[j]*choose(list.dim[j],list.dim[j]+list.dim[i]) <= MAX_NUM && list.dim[i]+list.dim[j] <= MAX_DIM)
      b =  list.insert(list.num[i]*list.num[j]*choose(list.dim[j],list.dim[j]+list.dim[i]),list.dim[i]+list.dim[j], list.name[i]+" x "+list.name[j]);
    if (list.dim[j] == list.dim[i])
      if (list.num[i]+list.num[j]-2 <= MAX_NUM && list.dim[i] <= MAX_DIM)
        b =  list.insert(list.num[i]+list.num[j]-2,list.dim[j], list.name[i]+" # "+list.name[j]);
  }

  return b;  
}

void F (int a, int b) {

  if ( a < b + 2 || b < 1) {
    cout << "IMPOSSIBLE" << endl;
  } else if ( (a-2)%b == 0) {
    if ((a-2)/b > 1)
      cout << "S^" << b << "_" << b+2 << " # ... " << (a-2)/b << " times ... # " << " S^" << b << "_" << b+2;
    else
      cout << "S^" << b << "_" << b+2;
  

}

int main (int argc, char *argv[]) {
  
  int num;
  
  cin >> num;
  
  Manifolds list(num);
  list.insert(3,1,"S1_3");
  
  int last = 0;
  while (list.n < list.NUM_MAX) {
    if (last < 0 || last >= list.n)
      last = 0;
    last = crossAll(list,last,1024,10);
  }
  
  
  list.print(cout);
  
  return 0;
}
