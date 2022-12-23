cat ../makefile | sed -e 's/MYFLAGS = -I./MYFLAGS = -I. -fcheck=all -g /' > ./makefile
make -f ./makefile  SRC_DIR='../' PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH rtest -j12
#rm makefile
