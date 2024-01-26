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

ifdef RUN_ONLY
	TEST_OPTIONS += --run-only
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
STANDARD_CFG = \
	ascem/batch/batch.cfg \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	co2/sco2.cfg \
	default/543/543.cfg \
	default/anisothermal/anisothermal.cfg \
	default/batch/batch.cfg \
	default/column/column.cfg \
	default/constraint/constraint.cfg \
	default/condition/condition.cfg \
	default/condition/dataset/dataset.cfg \
	default/debug/debug.cfg \
	default/discretization/discretization.cfg \
	default/dispersion/dispersion.cfg \
	default/eos/eos.cfg \
	default/infiltrometer/infiltrometer.cfg \
	default/infiltration1d/infiltration1d.cfg \
	default/input_format/input_format.cfg \
	default/integral_flux/integral_flux.cfg \
	default/materials/anisotropic_tortuosity/anisotropic_tortuosity.cfg \
	default/materials/full_tensor/full_tensor.cfg \
	default/materials/materials.cfg \
	default/materials/mineral_volume/mineral_volume.cfg \
	default/materials/soil_compressibility/soil_compressibility.cfg \
	default/materials/strata_surface/strata_surface.cfg \
	default/multicontinuum/multicontinuum.cfg \
	default/numerical/numerical.cfg \
	default/output/compare/compare.cfg \
	default/output/output.cfg \
	default/restart/restart.cfg \
	default/restart/skip/skip.cfg \
	default/reaction_sandbox/reaction_sandbox.cfg \
	default/regions/cell_mapping/cell_mapping.cfg \
	default/regions/regions.cfg \
	default/salinity/salinity.cfg \
	default/scco2/mphase/mphase.cfg \
	default/srcsink_sandbox/srcsink_sandbox.cfg \
	default/steady/steady.cfg \
	default/tcond/T_dep_tcond.cfg \
	default/timestepping/timestepping.cfg \
	default/velocity/velocity.cfg \
	default/zflow/zflow.cfg \
	general/general.cfg \
	general/hydrate/hydrate.cfg \
	geomech/uniaxial_strain/uniaxial_strain.cfg \
	geomech/uniaxial_strain/uniaxial_strain.cfg \
	geomech/uniaxial_loading/uniaxial_loading.cfg \
	geomech/uniaxial_loading_biot/uniaxial_loading_biot.cfg \
	geomech/terzaghi/terzaghi.cfg \
	geomech/uniaxial_loading/uniaxial_loading.cfg \
	geomech/uniaxial_loading_biot/uniaxial_loading_biot.cfg \
	geomech/terzaghi/terzaghi.cfg \
	geophys/ert/ert.cfg \
	geothermal_hpt/1D_Calcite/calcite_hpt.cfg\
	inversion/inversion.cfg \
	ngee/ngee.cfg \
	pnf/unit_cube/unit_cube.cfg \
	ntrdc/ntrdc.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	shortcourse/regional_doublet/regional_doublet.cfg \
	shortcourse/electron_tower/electron_tower.cfg \
	surfaceflow/inlinesurface/inlinesurface.cfg \
	ufd/ufd.cfg \
	wipp/nuts/nuts.cfg \
	wipp/transport/batch/wipp_batch.cfg \
	wipp/transport/wipp_regular.cfg \
	wipp/wipp_flow/condition/condition.cfg \
	wipp/wipp_flow/gas_generation/gas_generation.cfg \
	wipp/wipp_flow/input_format/input_format.cfg \
	wipp/well_model/well.cfg

ifneq ($(strip $(HYPRE_LIB)),)
STANDARD_CFG += \
	general/cpr_general/cpr_general.cfg \
	wipp/wipp_flow/cpr_wipp/cpr_wipp.cfg \
	wipp/well_model/cpr/cpr_well.cfg
endif

ifdef RUNTIME_ERROR_CHECKING
STANDARD_CFG += \
	runtime/runtime.cfg
endif

#
# domain specific problems
#
GEOCHEMISTRY_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg \
	ascem/batch/batch.cfg \
	default/543/543.cfg \
	default/batch/batch.cfg \
	default/anisothermal/anisothermal.cfg \
	default/constraint/constraint.cfg \
	default/column/column.cfg \
	default/multicontinuum/multicontinuum.cfg \
	geothermal_hpt/1D_Calcite/calcite_hpt.cfg\
	ngee/ngee.cfg \
	shortcourse/1D_Calcite/calcite.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg \
	ufd/ufd.cfg \
	wipp/nuts/nuts.cfg

GEOMECHANICS_CFG = \
	geomech/uniaxial_loading/uniaxial_loading.cfg \
	geomech/uniaxial_loading_biot/uniaxial_loading_biot.cfg \
	geomech/uniaxial_strain/uniaxial_strain.cfg \
	geomech/terzaghi/terzaghi.cfg

FLOW_CFG = \
	default/543/543.cfg \
	default/condition/condition.cfg \
	default/eos/eos.cfg \
	default/infiltration1d/infiltration1d.cfg \
	default/infiltrometer/infiltrometer.cfg \
	default/steady/steady.cfg \
	default/scco2/mphase/mphase.cfg \
	default/zflow/zflow.cfg \
	general/general.cfg \
	general/cpr_general/cpr_general.cfg \
        pnf/unit_cube/unit_cube.cfg \
	ntrdc/ntrdc.cfg \
	shortcourse/1D_variably_saturated_flow/vsat_flow.cfg \
	shortcourse/copper_leaching/cu_leaching.cfg

TRANSPORT_CFG = \
	default/543/543.cfg \
	default/column/column.cfg \
	default/constraint/constraint.cfg \
	default/dispersion/dispersion.cfg \
	default/multicontinuum/multicontinuum.cfg \
	default/velocity/velocity.cfg \
	wipp/transport/batch/wipp_batch.cfg \
	wipp/transport/wipp_regular.cfg \
	wipp/well_model/well.cfg

GEOPHYSICS_CFG = \
	geophys/ert/ert.cfg

INVERSE_CFG = \
	inversion/inversion.cfg

MESH_CFG = \
	default/discretization/discretization.cfg

CHECK_CFG = \
	ascem/1d/1d-calcite/1d-calcite.cfg

TCOND_CFG = \
    default/tcond/T_dep_tcond.cfg

WIPP_CFG = \
	wipp/nuts/nuts.cfg \
	wipp/wipp_flow/gas_generation/gas_generation.cfg \
	wipp/transport/batch/wipp_batch.cfg 

RUNTIME_CFG = \
	runtime/runtime.cfg

test : standard

standard :
ifeq ($(strip $(PARMETIS_LIB)),)
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(STANDARD_CFG)
else
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel unstructured unstructured_parallel \
		--config-files $(STANDARD_CFG)
endif

standard_parallel :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard_parallel \
		--config-files $(STANDARD_CFG)

check : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--suite standard standard_parallel \
		--config-files $(CHECK_CFG)

unstructured :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(STANDARD_CFG) --suite unstructured 

unstructured_parallel :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(STANDARD_CFG) --suite unstructured_parallel

geochemistry :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOCHEMISTRY_CFG) --suite geochemistry 

geomechanics :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOMECHANICS_CFG) --suite geomechanics 
else
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping geomechanics tests."
	@echo "********************************************************"
endif

flow :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(FLOW_CFG) --suite flow

transport :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(TRANSPORT_CFG) --suite transport

geophysics :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(GEOPHYSICS_CFG) --suite geophysics

inverse :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(INVERSE_CFG) --suite inverse

mesh :
ifneq ($(strip $(PARMETIS_LIB)),)
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(MESH_CFG) --suite mesh
else
	@echo "********************************************************"
	@echo "  PFLOTRAN does not appear to be compiled with Parmetis."
	@echo "  Skipping unstructured mesh tests."
	@echo "********************************************************"
endif

ngee-biogeochemistry : 
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		-r ngee --suite biogeochemistry 

test_mimetic :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(MFD_CFG)

tcond_test :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(TCOND_CFG)

wipp_test :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(WIPP_CFG)

runtime_error :
	-$(PYTHON) $(TEST_MANAGER) -e $(PFLOTRAN) $(TEST_OPTIONS) \
		--config-files $(RUNTIME_CFG) --suite runtime

#
# Try to cleanup the regression test output.
#
# NOTE: do NOT add a "*.h5" glob. The .h5 extension is used for both
# input and output files, and a simple glob can not distinguish
# between them. Manually add rm commands for each file that needs to be
# removed.
#
# skip gold files
CLEAN_COMMAND = ! -name '*.gold'
# print full filename followed by null character
CLEAN_COMMAND += -print0
# pipe through xargs (-0 = string is terminated with a null character)
CLEAN_COMMAND += | xargs -0
PLATFORM := $(shell uname)
# This is a workaround for Mac OS which does not accept the '-r' for xargs.
# GNU-based xargs requires '-r' to avoid a failure when the string is empty
# (e.g. no files exists with that suffix).
ifneq ($(PLATFORM),Darwin)
  CLEAN_COMMAND += -r
endif
CLEAN_COMMAND += rm
clean-tests :
	-find . -type f -name '*.testlog' $(CLEAN_COMMAND)
	-find . -type f -name '*.out' $(CLEAN_COMMAND)
	-find . -type f -name '*.tec' $(CLEAN_COMMAND)
	-find . -type f -name '*.pft' $(CLEAN_COMMAND)
	-find . -type f -name '*.ekg' $(CLEAN_COMMAND)
	-find . -type f -name '*-ert-*' $(CLEAN_COMMAND)
	-find . -type f -name '*.mat' $(CLEAN_COMMAND)
	-find . -type f -name '*.bin' $(CLEAN_COMMAND)
	-find . -type f -name '*.rec' $(CLEAN_COMMAND)
	-find . -type f -name '*.xmf' $(CLEAN_COMMAND)
	-find . -type f -name '*.vtk' $(CLEAN_COMMAND)
	-find . -type f -name '*.bio' $(CLEAN_COMMAND)
	-find . -type f -name '*.wf' $(CLEAN_COMMAND)
	-find . -type f -name '*-int.dat' $(CLEAN_COMMAND)
	-find . -type f -name '*-mas.dat' $(CLEAN_COMMAND)
	-find . -type f -name '*.pnl' $(CLEAN_COMMAND)
	-find . -type f -name '*.regression' $(CLEAN_COMMAND)
	-find . -type f -name '*.stdout' $(CLEAN_COMMAND)
	-find . -type f -name '*.old' $(CLEAN_COMMAND)
	-find . -type f -name '*~' $(CLEAN_COMMAND)
	-find . -type f -name 'tmp-restart-*' $(CLEAN_COMMAND)
	-find . -type f -name '*.chk' $(CLEAN_COMMAND)
	-find . -type f -name '*-restart.h5' $(CLEAN_COMMAND)
	-find . -type f -name '*_Jsense*' $(CLEAN_COMMAND)
	-find . -type f -name '*.pyc' $(CLEAN_COMMAND)
	-find . -type f -regex '.*Run[0-9]*.*' $(CLEAN_COMMAND)
	-find . -type d -name '__pycache__' $(CLEAN_COMMAND)
	-rm -f mfd/mfd.h5
	-rm -f default/543/543_flow_and_tracer-np8.h5
	-rm -f default/543/543_flow_and_tracer.h5
	-rm -f default/discretization/mixed_implicit-np4.h5
	-rm -f default/discretization/mixed_implicit.h5
	-rm -f default/eos/*.txt
	-rm -f default/eos/*.dat
	-rm -f general/calcite_flow_and_tran_general.h5
	-rm -f general/ictest1D.h5
	-rm -f general/ictest2D.h5
	-rm -f ngee/sink.h5 
	-rm -f ngee/source.h5
	-rm -f default/543/543_flow_and_tracer_dbase.h5
	-rm -f wipp/bragflo/wipp_flow/gas_generation_2d-np4.h5
	-rm -f wipp/bragflo/wipp_flow/well_production2-ss_mass-0.dat
	-rm -f wipp/bragflo/wipp_flow/gas_generation_2d.h5
	-rm -f wipp/bragflo/wipp_flow/gas_generation_one_cell.h5
	-rm -f wipp/well_model/*.h5
	-rm -f wipp/wipp_flow/gas_generation/*.h5
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
	-rm -f default/salinity/20x20_xz.h5
	-rm -f default/timestepping/20x20_xz.h5
	-rm -f geophys/ert/*.h5
	-rm -f geophys/ert/*.dat
	-rm -f inversion/*.h5
