#!/bin/bash 

BUILD_DIR=../Debug
FILE=dp1d_vok

mkdir -p $BUILD_DIR
cd ../fortran
gfortran -fno-range-check $FILE.f -o $BUILD_DIR/$FILE

cd $BUILD_DIR
./$FILE
