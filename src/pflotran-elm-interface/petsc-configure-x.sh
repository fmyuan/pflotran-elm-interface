#!/bin/sh -v

PATH=/usr/local/bin:/usr/bin:/bin
PATH=/usr/local/CMake.App/Contents/bin:$PATH
PATH=/usr/local/autotools/bin:$PATH
PATH=/usr/local/gcc-clang-darwin/gcc-x-clang/bin:$PATH
PATH=/usr/local/gcc-clang-darwin/openmpi-4.x-gcc/bin:$PATH

PACKAGE_ROOT=/usr/local/gcc-clang-darwin
MPI_CURRENT=openmpi-4.x-gcc
HDF5_CURRENT=hdf5-1.10-openmpi
NC4_CURRENT=netcdf-4.x-hdf5-openmpi

MPI_DIR=$PACKAGE_ROOT/$MPI_CURRENT 
HDF5_DIR=$PACKAGE_ROOT/$HDF5_CURRENT
NC4_DIR=$PACKAGE_ROOT/$NC4_CURRENT

echo $MPI_DIR
echo $HDF5_DIR


BLASLAPACK_LIB_DIR=$BLASLAPACK_DIR
CMAKE_DIR=/usr/local/CMake.app/Contents

# zlib library and include might not be in same directory in many systems
ZLIB_DIR_LIB=/usr/local/gcc-clang-darwin/zlib-1.2.12/lib
ZLIB_DIR_INC=/usr/local/gcc-clang-darwin/zlib-1.2.12/include

PETSC_SOURCE_DIR=./
INSTALL_DIR=$PACKAGE_ROOT/petsc-x-noopt
PETSC_ARCH=arch-darwin-noopt

cd ./

./configure \
           --prefix=$INSTALL_DIR \
           --with-clean=1 --with-c2html=0 --with-x=0 \
           --with-ssl=0 --with-debugging=0 --with-valgrind=0 \
           --with--cxx-dialect=C++11 \
           --with-shared-libraries=1 --with-debugging=0 --with-precision=double \
           --with-index-size=32 --with-memalign=16 --with-64-bit-indices=0 \
           --with-mpi-dir=$MPI_DIR --known-mpi-shared-libraries=0 --with-mpi=1 \
           --with-blas-lapack-dir=$BLASLAPACK_LIB_DIR \
           --with-zlib-lib=$ZLIB_DIR_LIB/libz.dylib --with-zlib-include=$ZLIB_DIR_INC \
           --with-cmake-dir=$CMAKE_DIR \
           --with-hdf5-dir=$HDF5_DIR --download-hdf5-fortran-bindings=yes \
           --download-sowing=yes \
           --download-metis=yes \
           --download-parmetis=yes \
           --download-mumps=yes \
           --download-scalapack=yes \
           --download-superlu=yes \
           --download-supperlu-dist=yes \
           --download-hypre=yes \
           LIBS=" -L$ZLIB_DIR/lib -lz -lm" \
           PETSC_DIR=$PETSC_SOURCE_DIR \
           PETSC_ARCH=$PETSC_ARCH \
           COPTFLAGS=" -fPIC -O0 " \
           FCOPTFLAGS="-fPIC -O0 -fallow-argument-mismatch" \
           CXXOPTFLAGS=" -fPIC -O0 " \
           FOPTFLAGS=" -O0 "


#           --with-opencl=1 --with-opencl-include=/System/Library/Frameworks/OpenCL.framework/Headers \
#           --with-opencl-lib="-framework opencl" \
#           --with-viennacl=1 --with-viennacl-dir=/usr/local/viennacl-1.7

