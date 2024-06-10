MYFLAGS =

CFLAGS = 
FFLAGS =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all : vecgetarray_test petsc_zero_pivot_test_f petsc_test_f ctx_test \
      petsc_test_c petsc_test_kspsolve_c petsc_test_kspsolve_f \
      petsc_zero_pivot_test_c petsc_test_mat_f petsc_test_mat

vecgetarray_test : vecgetarrayf90_test.o
	${FLINKER}   -o vecgetarrayf90_test vecgetarrayf90_test.o ${PETSC_LIB} ${LIBS} 
petsc_zero_pivot_test_f : petsc_zero_pivot_test_f.o
	${FLINKER}   -o petsc_zero_pivot_test_f petsc_zero_pivot_test_f.o ${PETSC_LIB} ${LIBS} 

petsc_test_f : petsc_test_f.o
	${FLINKER}   -o petsc_test_f petsc_test_f.o ${PETSC_LIB} ${LIBS} 

ctx_test : ctx_test.o
	${FLINKER}   -o ctx_test ctx_test.o ${PETSC_LIB} ${LIBS} 

petsc_test_c: petsc_test_c.o
	${CLINKER}   -o petsc_test_c petsc_test_c.o ${PETSC_LIB} ${LIBS} 

petsc_test_kspsolve_c: petsc_test_kspsolve_c.o
	${CLINKER}   -o petsc_test_kspsolve_c petsc_test_kspsolve_c.o ${PETSC_LIB} ${LIBS} 

petsc_test_kspsolve_f: petsc_test_kspsolve_f.o
	${CLINKER}   -o petsc_test_kspsolve_f petsc_test_kspsolve_f.o ${PETSC_LIB} ${LIBS} 

petsc_zero_pivot_test_c: petsc_zero_pivot_test_c.o
	${CLINKER}   -o petsc_zero_pivot_test_c petsc_zero_pivot_test_c.o ${PETSC_LIB} ${LIBS} 

petsc_test_mat: petsc_test_mat.o
	${CLINKER}   -o petsc_test_mat petsc_test_mat.o ${PETSC_LIB} ${LIBS} 

petsc_test_mat_f: petsc_test_mat_f.o
	${CLINKER}   -o petsc_test_mat_f petsc_test_mat_f.o ${PETSC_LIB} ${LIBS} 

%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<
