# Makefile for running standard regression tests
#
# include the make variables from PETSc so we can use PYTHON and do
# some conditional testing, e.g. test if unstructured mesh is available.
#
# If PYTHON defined from petsc is not correct, override it on the
# command line with: make PYTHON=python3.3 test

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


TEST_MANAGER = regression_tests.py
PFLOTRAN = ../src/pflotran/pflotran

TEST_OPTIONS =

# make PERFORMANCE=true check
ifdef PERFORMANCE
	TEST_OPTIONS += --check-performance
endif

ifdef DRYRUN
	TEST_OPTIONS += --dry-run
endif

ifdef UPDATE
	TEST_OPTIONS += --update
endif

ifdef CHECK
	TEST_OPTIONS += --check-only
endif

ifdef NEW
	TEST_OPTIONS += --new-tests
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

#
# standard tests that are run to verify pflotran is built correctly
#
STANDARD_CFG = general/general.cfg
	#default/anisothermal/anisothermal.cfg \
	#default/column/column.cfg \
	#default/eos/eos.cfg \
	#default/materials/soil_compressibility/soil_compressibility.cfg \
	#default/restart/restart.cfg \
	#default/condition/condition.cfg \
	#default/condition/dataset/dataset.cfg \
	#general/general.cfg


#STANDARD_CFG += \
#	surfaceflow/inlinesurface/inlinesurface.cfg

#
# domain specific problems
#

FLOW_CFG = 
#	default/condition/condition.cfg \
#	default/eos/eos.cfg \
#	default/infiltrometer/infiltrometer.cfg \
#	default/steady/steady.cfg \
#	general/general.cfg \
#	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg


test : standard

standard :
ifeq ($(strip $(PARMETIS_LIB)),)
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif
	$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(STANDARD_CFG)


standard_parallel :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard_parallel \
		--config-files $(STANDARD_CFG)

check : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(CHECK_CFG)

flow :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(FLOW_CFG) --suite flow



#
# Try to cleanup the regression test output.
#
# NOTE: do NOT add a "*.h5" glob. The .h5 extension is used for both
# input and output files, and a simple glob can not distinguish
# between them. Manually add rm commands for each file that needs to be
# removed.
#
clean-tests :
	-find . -type f -name '*.testlog' -print0 | xargs -0 rm
	-find . -type f -name '*.out' -print0 | xargs -0 rm
	-find . -type f -name '*.tec' -print0 | xargs -0 rm
	-find . -type f -name '*.regression' -print0 | xargs -0 rm
	-find . -type f -name '*.stdout' -print0 | xargs -0 rm
	-find . -type f -name '*.old' -print0 | xargs -0 rm
	-find . -type f -name '*~' -print0 | xargs -0 rm
	-find . -type f -name 'tmp-restart-*' -print0 | xargs -0 rm
	-find . -type f -name '*.chk' -print0 | xargs -0 rm
	-find . -type f -name '*.pyc' -print0 | xargs -0 rm
	-find . -type f -name 'sigma*.txt' -print0 | xargs -0 rm 
	-find . -type f -name '*.dpd' -print0 | xargs -0 rm
	-find . -type f -name 'e4d*.log' -print0 | xargs -0 rm
	-find . -type d -name '__pycache__' -print0 | xargs -0 rm -r
	-rm -f general/ictest1D.h5
	-rm -f general/ictest1D-mas.dat
	-rm -f general/ictest2D.h5
	-rm -f general/ictest2D-mas.dat
	-rm -f default/543/543_flow_and_tracer_dbase.h5
	-rm -f ascem/batch/calcite-kinetics-int.dat
	-rm -f ascem/batch/calcite-kinetics-volume-fractions-int.dat
	-rm -f ascem/batch/general-reaction-int.dat
	-rm -f default/multicontinuum/0D_heat_MPH-int.dat
	-rm -f default/multicontinuum/0D_heat_TH-int.dat
	-rm -f default/543/543_hanford_srfcplx_base*.h5
	-rm -f default/anisothermal/thc_1d.h5
	-rm -f default/anisothermal/th_1d.h5
	-rm -f default/condition/dataset/east_surface_cell_centered.h5
	-rm -f default/condition/dataset/east_surface_node_centered.h5
	-rm -f default/condition/dataset/north_surface_cell_centered.h5
	-rm -f default/condition/dataset/north_surface_node_centered.h5
	-rm -f default/condition/dataset/top_surface_cell_centered.h5
	-rm -f default/condition/dataset/top_surface_node_centered.h5
	-rm -f default/condition/dataset/x_line_cell_centered.h5
	-rm -f default/condition/dataset/x_line_node_centered.h5
	-rm -f default/condition/dataset/y_line_cell_centered.h5
	-rm -f default/condition/dataset/y_line_node_centered.h5
	-rm -f default/condition/dataset/z_line_cell_centered.h5
	-rm -f default/condition/dataset/z_line_node_centered.h5
