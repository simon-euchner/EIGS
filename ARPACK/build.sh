#!/bin/sh

### Run this script to compile and link ARPACK into a C library

CMPL="gfortran -fPIC -Wall -pedantic -fcheck=all -c"
LINK="gfortran -shared ./OBJ.D/* -o"

echo "Compiling ARPACK"
cd ./SRC.D/
$CMPL ./*
mv *.o ../OBJ.D/
cd ../UTIL.D/
$CMPL ./*
mv *.o ../OBJ.D/
cd ../ICB.D/
cd ./ICB.D/
$CMPL ./icbz.f
$CMPL ./icbd.f
mv ./icbz.o ../OBJ.D/
mv ./icbd.o ../OBJ.D/
cd ../
echo "Linking shared library"
$LINK ./LIB.D/libarpack.so
echo "Building process complete"
