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
g++ classes/*.h namespaces/*.h -I classes -I namespaces -std=c++11

if [ "$polyhedrons" = "1" ]; then
  echo "Compiling simplicial polyhedrons tools";
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/times.cpp -I classes -I namespaces -o times
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/boundary.cpp -I classes -I namespaces -o boundary 
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/cone.cpp -I classes -I namespaces -o cone
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/suspension.cpp -I classes -I namespaces -o suspension
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/simplex.cpp -I classes -I namespaces -o simplex
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/connexSum.cpp -I classes -I namespaces -o sum
  g++ namespaces/*.cpp classes/*.cpp polyhedrons/sphere1.cpp -I classes -I namespaces -o s1
fi

if [ "$matrices" = "1" ]; then
  echo "Compiling sparse Matrixes tools";
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_generator.cpp -I classes -I namespaces -o spmatrix
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_convert.cpp -I classes -I namespaces -o spconvert
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_ker.cpp -I classes -I namespaces -o spker
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_tester.cpp -I classes -I namespaces -o sptest
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_ldu.cpp -I classes -I namespaces -o spldu
  g++ namespaces/*.cpp classes/*.cpp spmatrixes/sparseMatrix_multiply.cpp -I classes -I namespaces -o spmultiply
fi

echo "Compiling main";
g++ -g namespaces/*.cpp classes/*.cpp main.cpp -I classes -I namespaces -std=c++11 -o main

