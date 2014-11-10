# Makefile for PFLOTRAN coupled with CLM

pflotran_src = ./
common_src   = ./
pflotran_test_dir = ../../regression_tests

MYFLAGS = -I.

CFLAGS   =
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = libepic.a

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif


# Begin Source Block

epic_obj = \
	${common_src}EPIC.o \
	${common_src}EPICConstants.o \
	${common_src}EPICData.o \
  ${common_src}EPICMesh.o

# End Source Block

lib_epic_obj = \
	${epic_obj}

pflotran_interface_test_obj = ${pflotran_src}pflotran_interface_main.o

libpflotran.a : $(lib_pflotran_obj)
	$(AR) $(ARFLAGS) $@ $(lib_pflotran_obj)

libepic.a : $(lib_epic_obj)
	$(AR) $(ARFLAGS) $@ $(lib_epic_obj)

libpflotranchem.a : $(util_obj) $(shared_mode_aux_obj) $(chem_obj)
	$(AR) $(ARFLAGS) $@ $(util_obj) $(shared_mode_aux_obj) $(chem_obj)

pflotran_interface : libpflotran.a $(pflotran_interface_test_obj)
	${FLINKER} -o $@ libpflotran.a $(pflotran_interface_test_obj) ${PETSC_LIB} ${LIBS} 

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# NOTE(bja, 2013-07) dependencies for clm interface. Keep these separate so that
# the long list below is easier to update from the main makefile.
clm_pflotran_interface_data.o :
EPIC.o : EPICConstants.o EPICData.o EPICMesh.o
EPICConstants.o :
EPICData.o : EPICConstants.o
EPICMesh.o : EPICConstants.o

mapping.o : option.o pflotran_constants.o input_aux.o

