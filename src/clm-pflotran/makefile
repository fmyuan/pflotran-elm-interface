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
MYFLAGS += ${FC_DEFINE_FLAG}CHUAN_CO2
MYFLAGS += ${FC_DEFINE_FLAG}CLM_PFLOTRAN
MYFLAGS += ${FC_DEFINE_FLAG}STORE_FLOWRATES
MYFLAGS += ${FC_DEFINE_FLAG}GEOMECH
#MYFLAGS += ${FC_DEFINE_FLAG}LEGACY_SATURATION_FUNCTION

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

ifdef scco2
  MYFLAGS += -DCHUAN_CO2
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

# Begin Source Block

util_obj  = \
        ${common_src}appleyard.o \
        ${common_src}generic.o \
	${common_src}input_aux.o \
	${common_src}logging.o \
	${common_src}matrix_block_aux.o \
	${common_src}option.o \
	${common_src}option_flow.o \
	${common_src}option_transport.o \
	${common_src}output_aux.o \
	${common_src}pflotran_constants.o \
	${common_src}pflotran_provenance.o \
	${common_src}string.o \
	${common_src}units.o \
	${common_src}utility.o \
        ${common_src}toil_ims_derivs.o \
	${common_src}variables.o

eos_obj = \
	${common_src}co2_sw.o \
	${common_src}co2_span_wagner_spline.o \
	${common_src}eos.o \
	${common_src}eos_database.o \
	${common_src}eos_gas.o \
	${common_src}eos_oil.o \
	${common_src}gas_eos_mod.o 

mode_aux_obj = \
	${common_src}auxiliary.o \
	${common_src}auxvars_base.o \
	${common_src}auxvars_flow.o \
	${common_src}auxvars_flow_energy.o \
	${common_src}auxvars_toil_ims.o \
	${common_src}auxvars_towg.o \
	${common_src}flash2_aux.o \
	${common_src}general_aux.o \
	${common_src}wipp_flow_aux.o \
	${common_src}pm_base_aux.o \
	${common_src}pm_toil_ims_aux.o \
	${common_src}pm_towg_aux.o \
	${common_src}immis_aux.o \
	${common_src}miscible_aux.o \
	${common_src}mphase_aux.o \
	${common_src}pm_auxiliary.o \
	${common_src}pm_base.o \
	${common_src}pm_base_pointer.o \
	${common_src}pm_general.o \
	${common_src}pm_wipp_flow.o \
	${common_src}pm_bragflo.o \
	${common_src}pm_toil_ims.o \
	${common_src}pm_towg.o \
	${common_src}pm_flash2.o \
	${common_src}pm_immis.o \
	${common_src}pm_mphase.o \
	${common_src}pm_miscible.o \
	${common_src}pm_richards.o \
	${common_src}pm_rt.o \
	${common_src}pm_subsurface_flow.o \
	${common_src}pm_surface_flow.o \
	${common_src}pm_surface_th.o \
	${common_src}pm_surface.o \
	${common_src}pm_th.o \
	${common_src}pm_ufd_decay.o \
	${common_src}pm_waste_form.o \
	${common_src}pm_wipp_srcsink.o \
	${common_src}pm_ufd_biosphere.o \
	${common_src}richards_aux.o \
	${common_src}th_aux.o \
	${common_src}inlinesurface_aux.o

mode_obj = \
	${common_src}flash2.o \
	${common_src}general.o \
	${common_src}general_common.o \
	${common_src}wipp_flow.o \
	${common_src}wipp_flow_common.o \
	${common_src}bragflo.o \
	${common_src}bragflo_common.o \
	${common_src}toil_ims.o \
	${common_src}towg.o \
	${common_src}global.o \
	${common_src}immis.o \
	${common_src}miscible.o \
	${common_src}mphase.o \
	${common_src}reactive_transport.o \
	${common_src}richards.o \
	${common_src}richards_common.o \
	${common_src}th.o \
	${common_src}inlinesurface.o

shared_mode_aux_obj = \
	${common_src}co2eos.o \
	${common_src}co2_span_wagner.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}eos_water.o \
	${common_src}geometry.o \
	${common_src}global_aux.o \
	${common_src}material_aux.o \
	${common_src}reactive_transport_aux.o 

chem_obj = \
	${common_src}transport_constraint.o \
	${common_src}reaction_database.o \
	${common_src}reaction_database_aux.o \
	${common_src}reaction_database_hpt.o \
	${common_src}reaction_gas.o \
	${common_src}reaction_gas_aux.o \
	${common_src}reaction_immobile.o \
	${common_src}reaction_immobile_aux.o \
	${common_src}reaction_microbial.o \
	${common_src}reaction_microbial_aux.o \
	${common_src}reaction_mineral.o \
	${common_src}reaction_mineral_aux.o \
	${common_src}reaction.o \
	${common_src}reaction_aux.o \
	${common_src}reaction_clm.o \
	${common_src}reaction_sandbox.o \
	${common_src}reaction_sandbox_base.o \
	${common_src}reaction_sandbox_clm_cn.o \
	${common_src}reaction_sandbox_example.o \
	${common_src}reaction_sandbox_pnnl_cyber.o \
	${common_src}reaction_sandbox_simple.o \
	${common_src}reaction_sandbox_ufd_wp.o \
	${common_src}reaction_solid_solution.o \
	${common_src}reaction_solid_solution_aux.o \
	${common_src}reaction_surf_complex.o \
	${common_src}reaction_surf_complex_aux.o 

grid_obj = \
	${common_src}connection.o \
	${common_src}communicator_base.o \
	${common_src}communicator_structured.o \
	${common_src}communicator_unstructured.o \
	${common_src}discretization.o \
	${common_src}dm_kludge.o \
	${common_src}grid.o \
	${common_src}grid_structured.o \
	${common_src}grid_unstructured_cell.o \
	${common_src}grid_unstructured_explicit.o \
	${common_src}grid_unstructured_aux.o \
	${common_src}grid_unstructured.o \
	${common_src}grid_unstructured_polyhedra.o

relations_obj = \
	${common_src}saturation.o \
	${common_src}saturation_function.o \
	${common_src}characteristic_curves.o \
	${common_src}characteristic_curves_base.o \
	${common_src}characteristic_curves_common.o \
	${common_src}characteristic_curves_owg.o \
	${common_src}characteristic_curves_wipp.o \
	${common_src}mphase_pckr_mod.o

properties_obj = \
	${common_src}fluid.o \
	${common_src}material.o 

dataset_obj = \
	${common_src}dataset_ascii.o \
	${common_src}dataset_base.o \
	${common_src}dataset_common_hdf5.o \
	${common_src}dataset_global_hdf5.o \
	${common_src}dataset_map_hdf5.o \
	${common_src}dataset.o \
	${common_src}dataset_gridded_hdf5.o \
	${common_src}time_storage.o 

srcsink_obj = \
	${common_src}srcsink_sandbox.o \
	${common_src}srcsink_sandbox_base.o \
	${common_src}srcsink_sandbox_mass_rate.o \
	${common_src}srcsink_sandbox_downreg.o \
	${common_src}srcsink_sandbox_wipp_gas.o \
	${common_src}srcsink_sandbox_wipp_well.o

io_obj = \
	${common_src}hdf5.o \
	${common_src}hdf5_aux.o \
	${common_src}output.o \
	${common_src}output_common.o \
	${common_src}output_ekg.o \
	${common_src}output_hdf5.o \
	${common_src}output_observation.o \
	${common_src}output_surface.o \
	${common_src}output_tecplot.o \
	${common_src}output_vtk.o

misc_obj = \
	${common_src}block_solve.o \
	${common_src}block_tridiag.o \
	${common_src}checkpoint.o \
	${common_src}condition.o \
	${common_src}condition_control.o \
	${common_src}convergence.o \
	${common_src}coupler.o \
	${common_src}data_mediator.o \
	${common_src}data_mediator_dataset.o \
	${common_src}data_mediator_vec.o \
	${common_src}data_mediator_base.o \
	${common_src}debug.o \
	${common_src}e4d_mat_inv.o \
	${common_src}e4d_run.o \
	${common_src}e4d_setup.o \
	${common_src}e4d_vars.o \
	${common_src}factory_hydrogeophysics.o \
	${common_src}factory_pflotran.o \
	${common_src}factory_subsurface.o \
	${common_src}field.o \
	${common_src}hydrostatic.o \
	${common_src}hydrostatic_common.o \
	${common_src}hydrostatic_multi_phase.o \
	${common_src}init_common.o \
	${common_src}init_subsurface.o \
	${common_src}init_subsurface_flow.o \
	${common_src}init_subsurface_transport.o \
	${common_src}integral_flux.o \
	${common_src}lookup_table.o \
	${common_src}matrix_buffer.o \
	${common_src}multisimulation.o \
	${common_src}observation.o \
	${common_src}patch.o \
	${common_src}pmc_auxiliary.o \
	${common_src}pmc_base.o \
	${common_src}pmc_hydrogeophysics.o \
	${common_src}pmc_subsurface.o \
	${common_src}pmc_third_party.o \
	${common_src}realization_base.o \
	${common_src}realization_subsurface.o \
	${common_src}region.o \
	${common_src}regression.o \
	${common_src}secondary_continuum_aux.o\
	${common_src}secondary_continuum.o \
	${common_src}simulation_aux.o \
	${common_src}simulation_base.o \
	${common_src}simulation_hydrogeophysics.o \
	${common_src}simulation_subsurface.o \
	${common_src}spline.o \
	${common_src}strata.o \
	${common_src}timestepper_base.o \
	${common_src}timestepper_BE.o \
	${common_src}transport.o \
	${common_src}waypoint.o \
	${common_src}wipp.o \
	${common_src}wrapper_hydrogeophysics.o

well_obj = \
	${common_src}well_spec_base.o \
	${common_src}well_base.o \
	${common_src}well_flow.o \
	${common_src}well_flow_energy.o \
	${common_src}well_water_injector.o \
	${common_src}well_oil_producer.o \
	${common_src}well_toil_ims.o \
	${common_src}well.o

surface_obj = \
	${common_src}checkpoint_surface.o \
	${common_src}factory_surface.o \
	${common_src}factory_surfsubsurface.o \
	${common_src}init_surface.o \
	${common_src}pmc_surface.o \
	${common_src}simulation_surface.o \
	${common_src}simulation_surfsubsurface.o \
	${common_src}surface_auxiliary.o \
	${common_src}surface_field.o \
	${common_src}surface_flow.o \
	${common_src}surface_global_aux.o \
	${common_src}surface_global.o \
	${common_src}surface_material.o \
	${common_src}realization_surface.o \
	${common_src}surface_th_aux.o \
	${common_src}surface_th.o \
	${common_src}timestepper_surface.o

geomech_obj = \
	${common_src}factory_geomechanics.o \
	${common_src}gauss.o \
	${common_src}geomechanics_auxiliary.o \
	${common_src}geomechanics_realization.o \
	${common_src}geomechanics_discretization.o \
	${common_src}geomech_grid.o \
	${common_src}geomech_grid_aux.o \
	${common_src}geomechanics_condition.o \
	${common_src}geomechanics_coupler.o \
	${common_src}geomechanics_debug.o \
	${common_src}geomechanics_field.o \
	${common_src}geomechanics_force.o \
	${common_src}geomechanics_global.o \
	${common_src}geomechanics_global_aux.o\
	${common_src}geomechanics_material.o \
	${common_src}geomechanics_patch.o \
	${common_src}geomechanics_region.o \
	${common_src}geomechanics_regression.o \
	${common_src}geomechanics_strata.o \
	${common_src}geomechanics_subsurface_properties.o \
	${common_src}output_geomechanics.o \
	${common_src}pm_geomechanics_force.o \
	${common_src}pmc_geomechanics.o \
	${common_src}simulation_geomechanics.o \
	${common_src}shape_function.o \
	${common_src}timestepper_steady.o

solver_obj = \
        ${common_src}preconditioner_cpr.o \
        ${common_src}solver.o \
        ${common_src}solver_cpr.o

clm_interface_obj = \
	${common_src}pflotran_model.o \
	${common_src}clm_pflotran_interface_data.o \
	${common_src}mapping.o

# End Source Block

lib_pflotran_obj = \
	${clm_interface_obj} \
	$(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(geomech_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) $(srcsink_obj) \
	$(io_obj) $(misc_obj) $(surface_obj) $(solver_obj)

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
                   factory_hydrogeophysics.o factory_pflotran.o \
                   factory_subsurface.o factory_surface.o \
                   factory_surfsubsurface.o field.o global_aux.o grid.o \
                   grid_unstructured.o grid_unstructured_aux.o \
                   grid_unstructured_cell.o input_aux.o logging.o mapping.o \
                   material.o multisimulation.o option.o patch.o \
                   pflotran_constants.o realization_base.o \
                   realization_subsurface.o realization_surface.o richards.o \
                   richards_aux.o saturation_function.o simulation_base.o \
                   simulation_subsurface.o simulation_surface.o \
                   simulation_surfsubsurface.o string.o surface_global_aux.o \
                   th.o th_aux.o timestepper_base.o units.o utility.o \
                   variables.o waypoint.o

###############################################################################
# Dependencies
# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
###############################################################################
include ../pflotran/pflotran_dependencies.txt
