# Makefile for PFLOTRAN coupled with CLM

pflotran_src = ./
common_src   = ./
pflotran_test_dir = ../../regression_tests

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

TEST_MANAGER = ../../regression_tests/regression_tests.py
PFLOTRAN = pflotran_interface

TEST_OPTIONS =

# make PERFORMANCE=true check
ifdef PERFORMANCE
	TEST_OPTIONS += --check-performance
endif

ifdef BACKTRACE
	TEST_OPTIONS += --backtrace
endif

ifdef TIMEOUT
	TEST_OPTIONS += --timeout $(TIMEOUT)
endif

ifneq ($(strip $(MPIEXEC)),)
	TEST_OPTIONS += --mpiexec $(MPIEXEC)
endif

STANDARD_CFG = \
	regression_tests/clm_pflotran.cfg

MYFLAGS = -I.
MYFLAGS += ${FC_DEFINE_FLAG}PROCESS_MODEL
MYFLAGS += ${FC_DEFINE_FLAG}CLM_PFLOTRAN
MYFLAGS += ${FC_DEFINE_FLAG}GEOMECH
#MYFLAGS += ${FC_DEFINE_FLAG}LEGACY_SATURATION_FUNCTION

ifdef codecov
  MYFLAGS += -fprofile-arcs -ftest-coverage
  LIBS += -lgcov
endif

# if compiling pflotran_rxn or libpflotranchem.a, need to toggle 
# PFLOTRAN_RXN_FLAG 
ifneq (,$(filter libpflotranchem.a pflotran_rxn,$(MAKECMDGOALS)))
  PFLOTRAN_RXN_FLAG := 1
endif

ifdef uo2
  MYFLAGS += -DCOMPUTE_INTERNAL_MASS_FLUX
endif

ifdef glenn
#  MYFLAGS += -DXINGYUAN_BC
  MYFLAGS += -DYE_FLUX
endif

ifdef solid_solution
  MYFLAGS += -DSOLID_SOLUTION
endif

ifdef radioactive_decay
  MYFLAGS += -DRADIOACTIVE_DECAY
endif

ifdef mualem_spline
  MYFLAGS += -DMUALEM_SPLINE
endif

ifdef petsc_api_3_3
  MYFLAGS += -DHAVE_PETSC_API_3_3
endif

ifdef have_petscpcmg_h
	MYFLAGS += -DHAVE_PETSCPCMG_H
endif

ifdef dbl
  MYFLAGS += -DDOUBLE_LAYER
endif

ifdef pcl
  MYFLAGS += -DPCL -DMUALEM_SPLINE
endif

ifdef condnr
  MYFLAGS += -DCONDNR
endif

ifdef perm
  MYFLAGS += -DPERM
endif

ifdef tab
  MYFLAGS += -DTABLE
endif

ifdef temp
  MYFLAGS += -DTEMP_DEPENDENT_LOGK
endif

ifdef remove_sat
  MYFLAGS += -DREMOVE_SATURATION
endif

ifdef compressibility
  MYFLAGS += -DUSE_COMPRESSIBILITY
endif

ifdef no_vapor_diffusion
  MYFLAGS += -DNO_VAPOR_DIFFUSION
endif

ifdef debug_gen
  MYFLAGS += -DDEBUG_GENERAL
endif 

ifdef coll
  MYFLAGS += -DCOLL
endif

ifdef dont_use_wateos
  MYFLAGS += -DDONT_USE_WATEOS
endif

ifdef have_hdf5
  MYFLAGS += -I$(HDF5_INCLUDE) -I$(HDF5_LIB) -DPETSC_HAVE_HDF5
endif

ifdef gpu
  MYFLAGS += -DCHUNK
endif

ifdef mfd
  MYFLAGS += -DDASVYAT
endif

ifdef dasvyat_debug
  MYFLAGS += -DDASVYAT_DEBUG
endif

ifdef dasvyat_test_cut
  MYFLAGS += -DDASVYAT_TEST_CUT
endif

ifdef ugrid_debug
  MYFLAGS += -DUGRID_DEBUG
endif

ifdef amanzi_bgd
  MYFLAGS += -DAMANZI_BGD
endif

ifdef store_flowrates
  MYFLAGS += -DSTORE_FLOWRATES
endif

ifdef mfd_ugrid
  MYFLAGS += -DMFD_UGRID
endif

ifdef clm_offline
  MYFLAGS += -DCLM_OFFLINE
endif

ifdef clm_pflotran
  MYFLAGS += -DCLM_PFLOTRAN -g
endif

ifdef map_debug
  MYFLAGS += -DMAP_DEBUG
endif

UPDATE_PROVENANCE=0
ifdef provenance
  UPDATE_PROVENANCE=1
endif

CFLAGS   =
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = pflotran pflotran_rxn libpflotranchem.a make.log pflotran_provenance.F90 \
  libpflotran.a pflotran_interface

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif

# Set this accordingly on your platform
# SCORPIO_DIR=${HOME}/soft/scorpio
ifdef scorpio
  LIBS += -L${SCORPIO_DIR}/lib -lscorpio
  MYFLAGS += -DSCORPIO
  MYFLAGS += -DSCORPIO_WRITE 
  MYFLAGS += -I${SCORPIO_DIR}/include
endif

ifdef tau
FC = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
FLINKER = /ccs/proj/geo002/tau/tau-2.19.1/craycnl/bin/tau_f90.sh -optCompInst -tau_makefile=$(TAU_MAKEFILE)
endif

# Begin object file block

include ../pflotran/pflotran_object_files.txt

clm_interface_obj = \
	${common_src}pflotran_model.o \
	${common_src}clm_pflotran_interface_data.o \
	${common_src}mapping.o

# End object file block

lib_pflotran_obj = \
	${clm_interface_obj} \
	$(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(geomech_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) $(srcsink_obj) \
	$(io_obj) $(misc_obj) $(solver_obj) $(well_obj) $(geophysics_obj) \
	$(deriv_obj)

pflotran_interface_test_obj = ${pflotran_src}pflotran_interface_main.o

pflotran_rxn_obj = ${pflotran_src}pflotran_rxn.o

remove_linked_files :
	./remove_linked_files.sh

link_files :
	./link_files.sh

libpflotran.a : $(lib_pflotran_obj)
	$(AR) $(ARFLAGS) $@ $(lib_pflotran_obj)

libpflotranchem.a : $(util_obj) $(shared_mode_aux_obj) $(chem_obj)
	$(AR) $(ARFLAGS) $@ $(util_obj) $(shared_mode_aux_obj) $(chem_obj)

pflotran_interface : libpflotran.a $(pflotran_interface_test_obj)
	${FLINKER} -o $@ $(pflotran_interface_test_obj) libpflotran.a ${PETSC_LIB} ${LIBS} 

pflotran_rxn : remove_linked_files link_files libpflotranchem.a $(pflotran_rxn_obj)
	${FLINKER} -o $@ $(pflotran_rxn_obj) -L. -lpflotranchem ${PETSC_LIB} ${LIBS} 

pflotran_provenance.F90 : FORCE
ifeq ($(UPDATE_PROVENANCE),1)
	@-rm pflotran_provenance.o pflotran_provenance_module.mod pflotran_provenance.F90
	$(PYTHON) ./pflotran-provenance.py \
		--pflotran-fflags "${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}" \
		--petsc-config "${CONFIGURE_OPTIONS}"
else
	@if test ! -f pflotran_provenance.F90; then cp pflotran_no_provenance.F90 pflotran_provenance.F90; fi
endif

# developer level regression testing
test : standard

standard :
	$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(STANDARD_CFG)

clean-tests :
	-find . -type f -name '*.testlog'    -print0 | xargs -0 rm
	-find . -type f -name '*.out'        -print0 | xargs -0 rm
	-find . -type f -name '*.tec'        -print0 | xargs -0 rm
	-find . -type f -name '*.regression' -print0 | xargs -0 rm
	-find . -type f -name '*.stdout'     -print0 | xargs -0 rm
	-find . -type f -name '*.old'        -print0 | xargs -0 rm

# user lever tests to verify pflotran is built correctly
check : FORCE
	@if [ -d $(pflotran_test_dir) ]; then \
		$(MAKE) --directory=$(pflotran_test_dir) $(MAKECMDGOALS); \
	fi

# null rule to force things to happen
FORCE :

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# NOTE(bja, 2013-07) dependencies for clm interface. Keep these separate so that
# the long list below is easier to update from the main makefile.
clm_pflotran_interface_data.o :
mapping.o : input_aux.o option.o pflotran_constants.o string.o
pflotran_model.o : characteristic_curves.o characteristic_curves_common.o \
                   clm_pflotran_interface_data.o connection.o coupler.o \
                   discretization.o factory_geomechanics.o \
                   factory_pflotran.o communicator_aux.o \
                   factory_subsurface.o factory.o \
                   field.o global_aux.o grid.o \
                   grid_unstructured.o grid_unstructured_aux.o \
                   grid_unstructured_cell.o input_aux.o logging.o mapping.o \
                   material.o simulation_multirealization.o option.o patch.o \
                   pflotran_constants.o realization_base.o \
                   realization_subsurface.o richards.o \
                   richards_aux.o saturation_function.o simulation_base.o \
                   simulation_subsurface.o \
                   string.o \
                   th.o th_aux.o timestepper_base.o units.o utility.o \
                   variables.o waypoint.o
richards.o : clm_pflotran_interface_data.o

###############################################################################
# Dependencies
# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
###############################################################################
include ../pflotran/pflotran_dependencies.txt
