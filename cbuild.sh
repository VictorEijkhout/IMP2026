#!/bin/bash

module load gcc/15 catch2/3.16
rm -rf build
cmake -S . -B build -DCMAKE_CXX_COMPILER=g++
cmake --build build
ctest --test-dir build
