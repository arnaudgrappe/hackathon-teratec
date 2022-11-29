#!/bin/bash
g++ -O3 stencil.cxx -o stencil4 -fopenmp
mpicxx -O3 stencil2.cxx -o stencil5 -fopenmp
