#!/bin/sh

usage() {
    cat <<- _EOF_


      compileAll [-a | --all | -m | --matrices | -p | --polyhedrons ]

                  Compile main function and optionally programs to work with matrices
                  or polyhedrons.


_EOF_

}

polyhedrons=0;
matrices=0;

while [ "$1" != "" ]; do
    case $1 in
        -a | --all )            polyhedrons=1
                                matrices=1
                                ;;
        -m | --matrices )       matrices=1
                                ;;
        -p | --polyhedrons )    polyhedrons=1
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

echo "Compiling headers";
g++ classes/*.h namespaces/*.h -I namespaces -I classes

if [ "$polyhedrons" = "1" ]; then
  echo "Compiling simplicial polyhedrons tools";
  g++ classes/*.cpp polyhedrons/times.cpp -I classes -o times
  g++ classes/*.cpp polyhedrons/boundary.cpp -I classes -o boundary 
  g++ classes/*.cpp polyhedrons/cone.cpp -I classes -o cone
  g++ classes/*.cpp polyhedrons/suspension.cpp -I classes -o suspension
  g++ classes/*.cpp polyhedrons/simplex.cpp -I classes -o simplex
  g++ classes/*.cpp polyhedrons/connexSum.cpp -I classes -o sum
  g++ classes/*.cpp polyhedrons/sphere1.cpp -I classes -o s1
fi

if [ "$matrices" = "1" ]; then
  echo "Compiling sparse Matrixes tools";
  g++ classes/*.cpp spmatrixes/sparseMatrix_generator.cpp -I classes -o spmatrix
  g++ classes/*.cpp spmatrixes/sparseMatrix_convert.cpp -I classes -o spconvert
  g++ classes/*.cpp spmatrixes/sparseMatrix_ker.cpp -I classes -o spker
  g++ classes/*.cpp spmatrixes/sparseMatrix_tester.cpp -I classes -o sptest
  g++ classes/*.cpp spmatrixes/sparseMatrix_ldu.cpp -I classes -o spldu
  g++ classes/*.cpp spmatrixes/sparseMatrix_multiply.cpp -I classes -o spmultiply
fi

echo "Compiling main";
g++ -g namespaces/*.cpp classes/*.cpp main.cpp -I classes -I namespaces -o main
#g++ classes/*.cpp cellDecomposition.cpp -I classes -o cellDecomposition

