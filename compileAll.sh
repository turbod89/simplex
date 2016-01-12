#!/bin/sh

echo "Compiling headers";
g++ classes/*.h

#echo "Compiling simplicial polyhedrons tools";
#g++ classes/*.cpp polyhedrons/times.cpp -I classes -o times
#g++ classes/*.cpp polyhedrons/boundary.cpp -I classes -o boundary 
#g++ classes/*.cpp polyhedrons/cone.cpp -I classes -o cone
#g++ classes/*.cpp polyhedrons/suspension.cpp -I classes -o suspension
#g++ classes/*.cpp polyhedrons/simplex.cpp -I classes -o simplex
#g++ classes/*.cpp polyhedrons/connexSum.cpp -I classes -o sum
#g++ classes/*.cpp polyhedrons/sphere1.cpp -I classes -o s1

#echo "Compiling sparse Matrixes tools";
#g++ classes/*.cpp spmatrixes/sparseMatrix_generator.cpp -I classes -o spmatrix
#g++ classes/*.cpp spmatrixes/sparseMatrix_convert.cpp -I classes -o spconvert
#g++ classes/*.cpp spmatrixes/sparseMatrix_ker.cpp -I classes -o spker
#g++ classes/*.cpp spmatrixes/sparseMatrix_tester.cpp -I classes -o sptest
#g++ classes/*.cpp spmatrixes/sparseMatrix_ldu.cpp -I classes -o spldu
#g++ classes/*.cpp spmatrixes/sparseMatrix_multiply.cpp -I classes -o spmultiply

echo "Compiling main";
g++ classes/*.cpp main.cpp -I classes -o main
g++ classes/*.cpp cellDecomposition.cpp -I classes -o cellDecomposition
#g++ classes/*.cpp manifolds_generator.cpp -I classes -o manifols_generator

