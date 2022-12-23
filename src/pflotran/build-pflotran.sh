cat ../makefile | sed -e 's/MYFLAGS = -I./MYFLAGS = -I. -fcheck=all -g /' > ./makefile
make -f ./makefile  SRC_DIR='../' PETSC_DIR=/home/defukuy/software/petsc/ PETSC_ARCH=gnu-c-debug rtest -j12
rm makefile
