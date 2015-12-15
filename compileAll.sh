#!/bin/sh

g++ classes/*.h

g++ classes/*.cpp polyhedrons/times.cpp -I classes -o times
g++ classes/*.cpp polyhedrons/boundary.cpp -I classes -o boundary 
g++ classes/*.cpp polyhedrons/cone.cpp -I classes -o cone
g++ classes/*.cpp polyhedrons/suspension.cpp -I classes -o suspension

g++ classes/*.cpp main.cpp -I classes -o main

