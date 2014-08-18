# Makefile for PFLOTRAN coupled with CLM

pflotran_src = ./
common_src   = ./
regression_test_dir = ../../regression_tests
unit_test_dir = ./unittests

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

MYFLAGS = -I.
MYFLAGS += ${FC_DEFINE_FLAG}PROCESS_MODEL
MYFLAGS += ${FC_DEFINE_FLAG}COMPUTE_INTERNAL_MASS_FLUX

ifdef petsc_api_3_3
  MYFLAGS += ${FC_DEFINE_FLAG}HAVE_PETSC_API_3_3
endif

ifdef have_petscpcmg_h
	MYFLAGS += ${FC_DEFINE_FLAG}HAVE_PETSCPCMG_H
endif

ifdef have_hdf5
  MYFLAGS += -I$(HDF5_INCLUDE) -I$(HDF5_LIB) ${FC_DEFINE_FLAG}PETSC_HAVE_HDF5
endif

ifdef store_flowrates
  MYFLAGS += ${FC_DEFINE_FLAG}STORE_FLOWRATES
endif

ifdef clm_pflotran
  MYFLAGS += ${FC_DEFINE_FLAG}CLM_PFLOTRAN
endif

UPDATE_PROVENANCE=0
ifdef provenance
  UPDATE_PROVENANCE=1
endif

CFLAGS   = 
FFLAGS   =
CPPFLAGS = ${MYFLAGS}
FPPFLAGS = ${MYFLAGS}

CLEANFILES       = pflotran libpflotranchem.a make.log pflotran_provenance.F90 \
                   libpflotran.a pflotran_interface

ifdef have_hdf5
LIBS +=  -L${HDF5_LIB} -lhdf5_fortran -lhdf5 -lz 
endif

# Begin Source Block

util_obj  = \
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
	${common_src}variables.o

eos_obj = \
	${common_src}co2eos.o \
	${common_src}co2_sw.o \
	${common_src}co2_sw_rtsafe.o \
	${common_src}co2_span_wagner.o \
	${common_src}co2_span_wagner_spline.o \
	${common_src}eos.o \
	${common_src}eos_gas.o \
	${common_src}eos_water.o \
	${common_src}gas_eos_mod.o 

mode_aux_obj = \
	${common_src}auxiliary.o \
	${common_src}pm_base.o \
	${common_src}pm_richards.o \
	${common_src}pm_rt.o \
	${common_src}pm_subsurface.o \
	${common_src}pm_surface_flow.o \
	${common_src}pm_surface_th.o \
	${common_src}pm_surface.o \
	${common_src}pm_top.o \
	${common_src}pm_th.o \
	${common_src}richards_aux.o \
	${common_src}th_aux.o

mode_obj = \
	${common_src}global.o \
	${common_src}reactive_transport.o \
	${common_src}richards.o \
	${common_src}richards_common.o \
	${common_src}richards_lsm.o \
	${common_src}th.o

shared_mode_aux_obj = \
	${common_src}global_aux.o \
	${common_src}reactive_transport_aux.o 

chem_obj = \
	${common_src}constraint.o \
	${common_src}clm_bgc.o \
	${common_src}database.o \
	${common_src}database_aux.o \
	${common_src}database_hpt.o \
	${common_src}immobile.o \
	${common_src}immobile_aux.o \
	${common_src}mineral.o \
	${common_src}mineral_aux.o \
	${common_src}reaction.o \
	${common_src}reaction_aux.o \
	${common_src}reaction_sandbox.o \
	${common_src}reaction_sandbox_base.o \
	${common_src}reaction_sandbox_clmdec.o \
	${common_src}reaction_sandbox_plantn.o \
	${common_src}reaction_sandbox_nitrification.o \
	${common_src}reaction_sandbox_denitrification.o \
	${common_src}reaction_sandbox_microbial.o \
	${common_src}reaction_sandbox_langmuir.o \
	${common_src}reaction_sandbox_degas.o \
	${common_src}surface_complexation.o \
	${common_src}surface_complexation_aux.o 

grid_obj = \
	${common_src}connection.o \
	${common_src}communicator_base.o \
	${common_src}discretization.o \
	${common_src}dm_kludge.o \
	${common_src}grid.o \
	${common_src}mfd.o \
	${common_src}mfd_aux.o \
	${common_src}unstructured_cell.o \
	${common_src}unstructured_communicator.o \
	${common_src}unstructured_explicit.o \
	${common_src}unstructured_grid_aux.o \
	${common_src}unstructured_grid.o \
	${common_src}unstructured_polyhedra.o \
	${common_src}structured_communicator.o \
	${common_src}structured_grid.o 

relations_obj = \
	${common_src}saturation.o \
	${common_src}saturation_function.o

properties_obj = \
	${common_src}fluid.o \
	${common_src}material_aux.o \
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

io_obj = \
	${common_src}hdf5.o \
	${common_src}hdf5_aux.o \
	${common_src}output.o \
	${common_src}output_common.o \
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
	${common_src}debug.o \
	${common_src}e4d_mat_inv.o \
	${common_src}e4d_run.o \
	${common_src}e4d_setup.o \
	${common_src}e4d_vars.o \
	${common_src}field.o \
	${common_src}geometry.o \
	${common_src}hydrostatic.o \
	${common_src}init.o \
	${common_src}mass_transfer.o \
	${common_src}matrix_buffer.o \
	${common_src}multisimulation.o \
	${common_src}observation.o \
	${common_src}patch.o \
	${common_src}pflotran_factory.o \
	${common_src}pmc_base.o \
	${common_src}pmc_subsurface.o \
	${common_src}realization.o \
	${common_src}realization_base.o \
	${common_src}region.o \
	${common_src}regression.o \
	${common_src}secondary_continuum_aux.o\
	${common_src}secondary_continuum.o \
	${common_src}simulation.o \
	${common_src}simulation_aux.o \
	${common_src}simulation_base.o \
	${common_src}subsurface_factory.o \
	${common_src}subsurface_simulation.o \
	${common_src}solver.o \
	${common_src}spline.o \
	${common_src}strata.o \
	${common_src}timestepper.o \
	${common_src}timestepper_base.o \
	${common_src}timestepper_BE.o \
	${common_src}transport.o \
	${common_src}uniform_velocity.o \
	${common_src}waypoint.o 

surface_obj = \
	${common_src}pmc_surface.o \
	${common_src}timestepper_surface.o \
	${common_src}surface_auxiliary.o \
	${common_src}surface_checkpoint.o \
	${common_src}surface_factory.o \
	${common_src}surfsubsurface_factory.o \
	${common_src}surface_field.o \
	${common_src}surface_flow.o \
	${common_src}surface_init.o \
	${common_src}surface_global_aux.o \
	${common_src}surface_global.o \
	${common_src}surface_material.o \
	${common_src}surface_realization.o \
	${common_src}surface_simulation.o \
	${common_src}surfsubsurface_simulation.o \
	${common_src}surface_th_aux.o \
	${common_src}surface_th.o

# End Source Block

pflotran_obj = $(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) \
	$(io_obj) $(misc_obj) $(surface_obj) ${pflotran_src}pflotran.o

lib_pflotran_obj = $(util_obj) $(eos_obj) $(mode_aux_obj) $(mode_obj) \
	$(shared_mode_aux_obj) $(chem_obj) $(grid_obj) \
	$(relations_obj) $(properties_obj) $(dataset_obj) \
	$(io_obj) $(misc_obj) $(surface_obj)

pflotran : $(pflotran_obj)
	${FLINKER} -o pflotran $(pflotran_obj) ${PETSC_LIB} ${LIBS} 

libpflotran.a : $(lib_pflotran_obj)
	$(AR) $(ARFLAGS) $@ $(lib_pflotran_obj)

pflotran_provenance.F90 : FORCE
ifeq ($(UPDATE_PROVENANCE),1)
	@-rm pflotran_provenance.o pflotran_provenance_module.mod pflotran_provenance.F90
	$(PYTHON) ./pflotran-provenance.py \
		--pflotran-fflags "${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}" \
		--petsc-config "${CONFIGURE_OPTIONS}"
else
	@if test ! -f pflotran_provenance.F90; then cp pflotran_no_provenance.F90 pflotran_provenance.F90; fi
endif

# Should add this to default PETSc targets as well?
%.mod : %.F90
	${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} $<

# NOTE(bja, 2013-07) dependencies for clm interface. Keep these separate so that
# the long list below is easier to update from the main makefile.
clm_pflotran_interface_data.o :
mapping.o : option.o pflotran_constants.o input_aux.o
pflotran_model.o : geomechanics_factory.o hydrogeophysics_factory.o \
             multisimulation.o option.o pflotran_constants.o \
             pflotran_factory.o simulation_base.o subsurface_factory.o \
             surface_factory.o surfsubsurface_factory.o \
             pflotran_provenance.o \
             mapping.o clm_pflotran_interface_data.o

# Dependencies stemming from "use" statements.
# These ensure that the module files are built in the correct order.
auxiliary.o : global_aux.o \
              material_aux.o pflotran_constants.o \
              reactive_transport_aux.o richards_aux.o \
              secondary_continuum_aux.o th_aux.o
block_solve.o :
block_tridiag.o : pflotran_constants.o
clm_bgc.o : clm_pflotran_interface_data.o
checkpoint.o : discretization.o field.o global.o grid.o material.o option.o \
               pflotran_constants.o realization.o string.o variables.o
co2_span_wagner.o : co2_sw_rtsafe.o option.o pflotran_constants.o
co2_span_wagner_spline.o : pflotran_constants.o spline.o
co2_sw.o : co2_span_wagner.o pflotran_constants.o
co2_sw_rtsafe.o : pflotran_constants.o
co2eos.o : pflotran_constants.o
communicator_base.o : dm_kludge.o option.o pflotran_constants.o
condition.o : constraint.o dataset.o dataset_ascii.o dataset_base.o \
              global_aux.o hdf5_aux.o input_aux.o logging.o \
              option.o pflotran_constants.o reaction_aux.o string.o \
              time_storage.o units.o
condition_control.o : condition.o connection.o constraint.o coupler.o \
                      dataset_base.o dataset_common_hdf5.o discretization.o \
                      eos_water.o field.o global_aux.o grid.o \
                      hdf5.o material_aux.o mfd.o option.o patch.o \
                      pflotran_constants.o reaction.o reaction_aux.o \
                      reactive_transport_aux.o realization.o region.o \
                      surface_field.o surface_global_aux.o \
                      surface_realization.o surface_th_aux.o variables.o
connection.o : pflotran_constants.o utility.o
constraint.o : global_aux.o immobile_aux.o input_aux.o logging.o mineral_aux.o \
               option.o pflotran_constants.o reaction_aux.o \
               reactive_transport_aux.o string.o surface_complexation_aux.o \
               units.o
convergence.o : grid.o option.o pflotran_constants.o solver.o
coupler.o : condition.o connection.o grid.o input_aux.o mfd_aux.o option.o \
            pflotran_constants.o region.o string.o structured_grid.o \
            unstructured_explicit.o unstructured_grid_aux.o utility.o
database.o : database_aux.o immobile.o immobile_aux.o input_aux.o \
             mineral.o mineral_aux.o option.o \
             pflotran_constants.o reaction.o reaction_aux.o reaction_sandbox.o \
             string.o surface_complexation_aux.o utility.o
database_aux.o : input_aux.o option.o pflotran_constants.o string.o utility.o
database_hpt.o : database.o database_aux.o input_aux.o mineral_aux.o option.o \
                 pflotran_constants.o reaction.o reaction_aux.o string.o \
                 surface_complexation_aux.o utility.o
dataset.o : dataset_ascii.o dataset_base.o dataset_common_hdf5.o \
            dataset_global_hdf5.o dataset_gridded_hdf5.o dataset_map_hdf5.o \
            dm_kludge.o input_aux.o option.o pflotran_constants.o string.o \
            time_storage.o
dataset_ascii.o : dataset_base.o input_aux.o option.o pflotran_constants.o \
                  string.o time_storage.o units.o utility.o
dataset_base.o : option.o pflotran_constants.o string.o time_storage.o \
                 utility.o
dataset_common_hdf5.o : dataset_base.o hdf5_aux.o input_aux.o logging.o \
                        option.o pflotran_constants.o string.o time_storage.o \
                        units.o
dataset_global_hdf5.o : dataset_base.o dataset_common_hdf5.o dm_kludge.o \
                        logging.o option.o pflotran_constants.o time_storage.o \
                        utility.o
dataset_gridded_hdf5.o : dataset_base.o dataset_common_hdf5.o logging.o \
                         option.o pflotran_constants.o string.o time_storage.o \
                         units.o utility.o
dataset_map_hdf5.o : dataset_base.o dataset_common_hdf5.o input_aux.o \
                     logging.o option.o pflotran_constants.o string.o \
                     time_storage.o units.o utility.o
debug.o : input_aux.o option.o pflotran_constants.o
discretization.o : dm_kludge.o grid.o input_aux.o mfd.o mfd_aux.o option.o \
                   pflotran_constants.o string.o structured_grid.o \
                   unstructured_explicit.o unstructured_grid.o \
                   unstructured_grid_aux.o unstructured_polyhedra.o
dm_kludge.o : pflotran_constants.o unstructured_grid_aux.o
e4d_mat_inv.o :
e4d_run.o : e4d_vars.o
e4d_setup.o : e4d_mat_inv.o e4d_vars.o
e4d_vars.o :
eos.o : eos_gas.o eos_water.o input_aux.o option.o pflotran_constants.o \
        string.o
eos_gas.o : pflotran_constants.o
eos_water.o : pflotran_constants.o
field.o : pflotran_constants.o
fluid.o : input_aux.o option.o pflotran_constants.o string.o
gas_eos_mod.o : pflotran_constants.o
gauss.o : pflotran_constants.o unstructured_cell.o utility.o
geometry.o : input_aux.o option.o pflotran_constants.o
global.o : communicator_base.o connection.o coupler.o field.o global_aux.o \
           grid.o material.o option.o patch.o pflotran_constants.o \
           realization.o realization_base.o variables.o
global_aux.o : option.o pflotran_constants.o utility.o
grid.o : connection.o dm_kludge.o logging.o mfd_aux.o option.o \
         pflotran_constants.o region.o structured_grid.o unstructured_cell.o \
         unstructured_explicit.o unstructured_grid.o unstructured_grid_aux.o \
         unstructured_polyhedra.o utility.o
hdf5.o : discretization.o field.o grid.o hdf5_aux.o logging.o option.o patch.o \
         pflotran_constants.o realization.o realization_base.o region.o \
         unstructured_cell.o utility.o
hdf5_aux.o : logging.o option.o pflotran_constants.o
hydrostatic.o : condition.o connection.o coupler.o dataset_ascii.o \
                dataset_gridded_hdf5.o eos_water.o grid.o \
                option.o pflotran_constants.o region.o structured_grid.o \
                utility.o
immobile.o : immobile_aux.o input_aux.o option.o pflotran_constants.o string.o \
             utility.o
immobile_aux.o : database_aux.o option.o pflotran_constants.o string.o \
                 utility.o
init.o : condition.o condition_control.o connection.o constraint.o \
         convergence.o coupler.o database.o database_hpt.o dataset.o \
         dataset_base.o dataset_common_hdf5.o debug.o discretization.o eos.o \
         eos_water.o field.o fluid.o \
         global.o grid.o hdf5.o input_aux.o \
         logging.o mass_transfer.o material.o material_aux.o mineral.o \
         observation.o option.o output.o output_aux.o \
         output_tecplot.o patch.o pflotran_constants.o reaction.o \
         reaction_aux.o reactive_transport.o realization.o region.o \
         regression.o richards.o saturation_function.o \
         secondary_continuum.o simulation.o solver.o \
         strata.o string.o structured_grid.o surface_field.o \
         surface_flow.o surface_global.o surface_init.o surface_realization.o \
         surface_th.o th.o timestepper.o uniform_velocity.o units.o \
         unstructured_grid.o unstructured_grid_aux.o variables.o waypoint.o
input_aux.o : option.o pflotran_constants.o string.o
logging.o : pflotran_constants.o
mass_transfer.o : dataset_base.o dataset_common_hdf5.o dataset_global_hdf5.o \
                  discretization.o grid.o input_aux.o option.o \
                  pflotran_constants.o string.o
material.o : communicator_base.o dataset_common_hdf5.o global_aux.o \
             input_aux.o material_aux.o option.o pflotran_constants.o \
             saturation_function.o string.o variables.o
material_aux.o : option.o pflotran_constants.o utility.o variables.o
matrix_block_aux.o : option.o pflotran_constants.o utility.o
matrix_buffer.o : grid.o pflotran_constants.o structured_grid.o
mfd.o : connection.o eos_water.o field.o global_aux.o grid.o material_aux.o \
        mfd_aux.o option.o patch.o pflotran_constants.o richards_aux.o \
        unstructured_cell.o unstructured_grid_aux.o
mfd_aux.o : connection.o option.o pflotran_constants.o
mineral.o : database_aux.o global_aux.o input_aux.o material_aux.o \
            mineral_aux.o option.o pflotran_constants.o reaction_aux.o \
            reactive_transport_aux.o string.o units.o \
            utility.o
mineral_aux.o : database_aux.o option.o pflotran_constants.o string.o \
                utility.o
multisimulation.o : input_aux.o option.o pflotran_constants.o simulation.o
observation.o : connection.o input_aux.o option.o pflotran_constants.o \
                region.o string.o
option.o : logging.o option_flow.o option_transport.o pflotran_constants.o
option_flow.o : pflotran_constants.o
option_transport.o : pflotran_constants.o
output.o : connection.o coupler.o debug.o discretization.o field.o grid.o \
           hdf5.o input_aux.o logging.o option.o output_aux.o output_common.o \
           output_hdf5.o output_observation.o output_tecplot.o output_vtk.o \
           patch.o pflotran_constants.o reaction_aux.o realization_base.o \
           variables.o
output_aux.o : input_aux.o option.o pflotran_constants.o string.o variables.o
output_common.o : connection.o coupler.o field.o global_aux.o grid.o hdf5.o \
                  hdf5_aux.o logging.o option.o output_aux.o patch.o \
                  pflotran_constants.o realization_base.o richards_aux.o \
                  unstructured_cell.o unstructured_grid_aux.o variables.o
output_hdf5.o : connection.o coupler.o discretization.o field.o grid.o hdf5.o \
                hdf5_aux.o input_aux.o logging.o option.o output_aux.o \
                output_common.o patch.o pflotran_constants.o \
                pflotran_provenance.o reaction_aux.o realization_base.o \
                string.o unstructured_cell.o unstructured_grid_aux.o \
                variables.o
output_observation.o : connection.o coupler.o discretization.o field.o \
                       global_aux.o grid.o \
                       logging.o material_aux.o \
                       observation.o option.o output_aux.o output_common.o \
                       patch.o pflotran_constants.o reaction_aux.o \
                       reactive_transport.o reactive_transport_aux.o \
                       realization.o realization_base.o region.o richards.o \
                       structured_grid.o th.o utility.o variables.o
output_surface.o : connection.o coupler.o discretization.o field.o grid.o \
                   hdf5.o hdf5_aux.o input_aux.o logging.o option.o \
                   output_aux.o output_common.o output_hdf5.o output_tecplot.o \
                   patch.o pflotran_constants.o reaction_aux.o realization.o \
                   string.o surface_field.o surface_realization.o \
                   unstructured_cell.o unstructured_grid_aux.o utility.o \
                   variables.o
output_tecplot.o : connection.o discretization.o field.o grid.o logging.o \
                   observation.o option.o output_aux.o output_common.o patch.o \
                   pflotran_constants.o reaction_aux.o realization_base.o \
                   region.o secondary_continuum_aux.o string.o \
                   structured_grid.o unstructured_cell.o \
                   unstructured_grid_aux.o utility.o variables.o
output_vtk.o : discretization.o field.o grid.o logging.o option.o output_aux.o \
               output_common.o patch.o pflotran_constants.o reaction_aux.o \
               realization_base.o string.o structured_grid.o variables.o
patch.o : auxiliary.o condition.o connection.o constraint.o coupler.o \
          dataset.o dataset_ascii.o dataset_base.o dataset_common_hdf5.o \
          dataset_gridded_hdf5.o dataset_map_hdf5.o eos_water.o field.o \
          global_aux.o grid.o hydrostatic.o \
          material.o material_aux.o mineral.o mineral_aux.o \
          observation.o option.o output_aux.o \
          pflotran_constants.o reaction.o reaction_aux.o \
          reactive_transport_aux.o region.o richards_aux.o saturation.o \
          saturation_function.o strata.o surface_auxiliary.o \
          surface_complexation_aux.o surface_field.o surface_material.o \
          th_aux.o utility.o variables.o
pflotran.o : \
             multisimulation.o option.o output_aux.o pflotran_constants.o \
             pflotran_factory.o pflotran_provenance.o simulation_base.o \
             subsurface_factory.o surface_factory.o surfsubsurface_factory.o
pflotran_constants.o :
pflotran_factory.o : input_aux.o multisimulation.o option.o \
                     pflotran_constants.o string.o
pflotran_provenance.o :
pm_base.o : option.o output_aux.o pflotran_constants.o realization_base.o
pm_richards.o : pflotran_constants.o pm_base.o pm_subsurface.o richards.o
pm_rt.o : communicator_base.o condition.o discretization.o field.o global.o \
          grid.o mass_transfer.o option.o pflotran_constants.o pm_base.o \
          reaction_aux.o reactive_transport.o realization.o realization_base.o \
          structured_communicator.o unstructured_communicator.o variables.o
pm_subsurface.o : checkpoint.o communicator_base.o condition.o \
                  discretization.o global.o grid.o material.o option.o \
                  pflotran_constants.o pm_base.o realization.o \
                  structured_communicator.o unstructured_communicator.o \
                  variables.o
pm_surface.o : communicator_base.o condition.o discretization.o global.o \
               grid.o option.o pflotran_constants.o pm_base.o \
               surface_checkpoint.o surface_realization.o \
               unstructured_communicator.o
pm_surface_flow.o : condition.o discretization.o pflotran_constants.o \
                    pm_base.o pm_surface.o surface_field.o surface_flow.o
pm_surface_th.o : communicator_base.o condition.o discretization.o grid.o \
                  option.o pflotran_constants.o pm_base.o pm_surface.o \
                  realization.o surface_field.o surface_realization.o \
                  surface_th.o
pm_th.o : communicator_base.o discretization.o global.o grid.o option.o \
          pflotran_constants.o pm_base.o pm_subsurface.o realization.o \
          structured_communicator.o th.o unstructured_communicator.o
pm_top.o : option.o pflotran_constants.o pm_base.o realization.o
pmc_base.o : checkpoint.o logging.o option.o output.o output_aux.o \
             pflotran_constants.o pm_base.o pm_top.o simulation_aux.o \
             timestepper_base.o utility.o waypoint.o
pmc_subsurface.o : connection.o coupler.o discretization.o eos_water.o field.o \
                   grid.o material.o material_aux.o option.o patch.o \
                   pflotran_constants.o pmc_base.o realization.o \
                   realization_base.o string.o variables.o
pmc_surface.o : connection.o coupler.o eos_water.o grid.o option.o output.o \
                output_surface.o patch.o pflotran_constants.o pm_base.o \
                pm_surface_flow.o pmc_base.o realization.o string.o \
                surface_flow.o surface_global_aux.o surface_realization.o \
                surface_th.o surface_th_aux.o timestepper_base.o \
                timestepper_surface.o utility.o
reaction.o : co2_span_wagner.o co2eos.o constraint.o eos_water.o global_aux.o \
             immobile.o immobile_aux.o input_aux.o material_aux.o \
             mineral.o mineral_aux.o option.o \
             pflotran_constants.o reaction_aux.o reaction_sandbox.o \
             reactive_transport_aux.o \
             string.o surface_complexation.o surface_complexation_aux.o \
             units.o utility.o variables.o
reaction_aux.o : database_aux.o immobile_aux.o mineral_aux.o \
                 option.o pflotran_constants.o string.o \
                 surface_complexation_aux.o utility.o
reaction_sandbox.o : global_aux.o input_aux.o material_aux.o option.o \
                     pflotran_constants.o reaction_aux.o \
                     reaction_sandbox_base.o \
                     reactive_transport_aux.o string.o utility.o \
                     reaction_sandbox_clmdec.o \
                     reaction_sandbox_plantn.o \
                     reaction_sandbox_nitrification.o \
                     reaction_sandbox_denitrification.o \
                     reaction_sandbox_microbial.o \
                     reaction_sandbox_langmuir.o \
                     reaction_sandbox_degas.o
reaction_sandbox_base.o : global_aux.o input_aux.o material_aux.o option.o \
                          pflotran_constants.o reaction_aux.o \
                          reactive_transport_aux.o
reaction_sandbox_clmdec.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reaction_sandbox_denitrification.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reaction_sandbox_degas.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o co2eos.o
reaction_sandbox_langmuir.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reaction_sandbox_microbial.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reaction_sandbox_nitrification.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reaction_sandbox_plantn.o : global_aux.o immobile_aux.o input_aux.o \
                            material_aux.o option.o pflotran_constants.o \
                            reaction_aux.o reaction_sandbox_base.o \
                            reactive_transport_aux.o string.o units.o \
                            utility.o clm_bgc.o
reactive_transport.o : condition.o connection.o constraint.o coupler.o debug.o \
                       discretization.o field.o fluid.o global.o global_aux.o \
                       grid.o logging.o mass_transfer.o material.o \
                       material_aux.o option.o output_aux.o output_tecplot.o \
                       patch.o pflotran_constants.o reaction.o reaction_aux.o \
                       reactive_transport_aux.o realization.o region.o \
                       secondary_continuum.o secondary_continuum_aux.o \
                       surface_complexation_aux.o transport.o variables.o
reactive_transport_aux.o : matrix_block_aux.o option.o pflotran_constants.o \
                           reaction_aux.o surface_complexation_aux.o utility.o
realization.o : condition.o connection.o constraint.o coupler.o dataset.o \
                dataset_base.o dataset_common_hdf5.o debug.o discretization.o \
                dm_kludge.o field.o fluid.o grid.o input_aux.o mass_transfer.o \
                material.o material_aux.o mfd.o mfd_aux.o observation.o \
                option.o output_aux.o patch.o pflotran_constants.o reaction.o \
                reaction_aux.o reactive_transport_aux.o realization_base.o \
                region.o saturation_function.o strata.o string.o \
                structured_communicator.o structured_grid.o time_storage.o \
                uniform_velocity.o unstructured_cell.o \
                unstructured_communicator.o unstructured_grid.o \
                unstructured_grid_aux.o variables.o waypoint.o
realization_base.o : communicator_base.o debug.o discretization.o field.o \
                     input_aux.o mass_transfer.o option.o output_aux.o patch.o \
                     pflotran_constants.o reaction_aux.o
region.o : geometry.o input_aux.o logging.o option.o pflotran_constants.o \
           string.o structured_grid.o unstructured_cell.o utility.o
regression.o : discretization.o grid.o input_aux.o option.o output.o \
               output_aux.o output_common.o pflotran_constants.o realization.o \
               string.o timestepper_BE.o utility.o
richards.o : connection.o coupler.o debug.o discretization.o eos_water.o \
             field.o global_aux.o grid.o logging.o mass_transfer.o material.o \
             material_aux.o matrix_buffer.o option.o output_aux.o patch.o \
             pflotran_constants.o realization.o realization_base.o \
             richards_aux.o richards_common.o richards_lsm.o \
             saturation_function.o string.o variables.o
richards_aux.o : eos_water.o global_aux.o material_aux.o matrix_buffer.o \
                 option.o pflotran_constants.o saturation_function.o utility.o
richards_common.o : connection.o eos_water.o global_aux.o material.o \
                    material_aux.o option.o pflotran_constants.o \
                    richards_aux.o saturation_function.o
richards_lsm.o : connection.o coupler.o field.o global_aux.o grid.o logging.o \
                 material.o option.o patch.o pflotran_constants.o \
                 realization.o richards_aux.o richards_common.o \
                 saturation_function.o
saturation.o : condition.o connection.o coupler.o grid.o option.o \
               pflotran_constants.o region.o saturation_function.o
saturation_function.o : input_aux.o option.o pflotran_constants.o string.o \
                        utility.o
secondary_continuum.o : block_solve.o block_tridiag.o condition.o constraint.o \
                        coupler.o eos_water.o global_aux.o grid.o material.o \
                        material_aux.o option.o \
                        pflotran_constants.o reaction.o reaction_aux.o \
                        reactive_transport_aux.o realization.o \
                        secondary_continuum_aux.o string.o utility.o
secondary_continuum_aux.o : option.o pflotran_constants.o \
                            reactive_transport_aux.o utility.o
shape_function.o : gauss.o pflotran_constants.o unstructured_cell.o
simulation.o : option.o \
               pflotran_constants.o reactive_transport.o realization.o \
               regression.o richards.o solver.o surface_realization.o \
               timestepper.o
simulation_aux.o : option.o pflotran_constants.o
simulation_base.o : logging.o option.o output.o output_aux.o \
                    pflotran_constants.o pmc_base.o simulation_aux.o \
                    timestepper_base.o waypoint.o
solver.o : input_aux.o option.o pflotran_constants.o string.o
spline.o : pflotran_constants.o
strata.o : input_aux.o material.o option.o pflotran_constants.o region.o \
           string.o surface_material.o
string.o : pflotran_constants.o
structured_communicator.o : communicator_base.o dm_kludge.o \
                            pflotran_constants.o structured_grid.o
structured_grid.o : connection.o input_aux.o option.o pflotran_constants.o \
                    string.o utility.o
subsurface_factory.o : condition_control.o init.o input_aux.o option.o \
                       output.o output_aux.o pflotran_constants.o pm_base.o \
                       pm_richards.o pm_rt.o pm_subsurface.o \
                       pm_th.o pm_top.o pmc_base.o pmc_subsurface.o \
                       reactive_transport.o realization.o simulation.o \
                       simulation_base.o subsurface_simulation.o timestepper.o \
                       timestepper_BE.o timestepper_base.o
subsurface_simulation.o : logging.o option.o output.o pflotran_constants.o \
                          pmc_base.o pmc_subsurface.o realization.o \
                          regression.o simulation_base.o timestepper_BE.o
surface_auxiliary.o : pflotran_constants.o surface_global_aux.o \
                      surface_th_aux.o
surface_checkpoint.o : discretization.o grid.o option.o output_aux.o \
                       pflotran_constants.o surface_field.o \
                       surface_realization.o
surface_complexation.o : global_aux.o input_aux.o material_aux.o \
                         matrix_block_aux.o option.o pflotran_constants.o \
                         reaction_aux.o reactive_transport_aux.o string.o \
                         surface_complexation_aux.o utility.o
surface_complexation_aux.o : database_aux.o option.o pflotran_constants.o \
                             utility.o
surface_factory.o : condition_control.o init.o input_aux.o logging.o option.o \
                    output.o output_aux.o output_surface.o \
                    pflotran_constants.o pm_base.o pm_surface.o \
                    pm_surface_flow.o pm_surface_th.o pm_top.o pmc_base.o \
                    pmc_surface.o simulation.o simulation_base.o \
                    surface_checkpoint.o surface_realization.o \
                    surface_simulation.o timestepper.o timestepper_base.o \
                    timestepper_surface.o
surface_field.o : pflotran_constants.o
surface_flow.o : connection.o coupler.o debug.o discretization.o dm_kludge.o \
                 eos_water.o global_aux.o grid.o logging.o option.o \
                 output_aux.o patch.o pflotran_constants.o realization.o \
                 realization_base.o string.o surface_field.o \
                 surface_global_aux.o surface_material.o surface_realization.o \
                 variables.o
surface_global.o : connection.o coupler.o discretization.o grid.o option.o \
                   patch.o pflotran_constants.o surface_field.o \
                   surface_global_aux.o surface_realization.o variables.o
surface_global_aux.o : option.o pflotran_constants.o utility.o
surface_init.o : condition.o coupler.o dataset.o dataset_base.o \
                 dataset_common_hdf5.o debug.o discretization.o field.o grid.o \
                 hdf5.o input_aux.o material.o option.o output_aux.o \
                 output_surface.o output_tecplot.o patch.o \
                 pflotran_constants.o region.o solver.o strata.o string.o \
                 structured_grid.o surface_auxiliary.o surface_field.o \
                 surface_material.o surface_realization.o units.o \
                 unstructured_grid.o unstructured_grid_aux.o waypoint.o
surface_material.o : input_aux.o option.o pflotran_constants.o string.o
surface_realization.o : condition.o coupler.o dataset.o dataset_base.o debug.o \
                        discretization.o dm_kludge.o field.o grid.o \
                        input_aux.o option.o output_aux.o patch.o \
                        pflotran_constants.o reaction_aux.o realization.o \
                        realization_base.o region.o strata.o string.o \
                        surface_field.o surface_material.o time_storage.o \
                        unstructured_cell.o unstructured_grid.o \
                        unstructured_grid_aux.o waypoint.o
surface_simulation.o : logging.o option.o output.o pflotran_constants.o \
                       pmc_base.o pmc_surface.o regression.o simulation_base.o \
                       surface_realization.o timestepper_base.o
surface_th.o : connection.o coupler.o debug.o discretization.o dm_kludge.o \
               eos_water.o fluid.o grid.o option.o output_aux.o patch.o \
               pflotran_constants.o realization.o realization_base.o region.o \
               string.o surface_field.o surface_global_aux.o \
               surface_material.o surface_realization.o surface_th_aux.o \
               variables.o
surface_th_aux.o : eos_water.o option.o pflotran_constants.o \
                   surface_global_aux.o
surfsubsurface_factory.o : coupler.o discretization.o dm_kludge.o field.o \
                           grid.o init.o input_aux.o option.o patch.o \
                           pflotran_constants.o pmc_base.o realization.o \
                           region.o simulation.o simulation_aux.o \
                           simulation_base.o string.o subsurface_factory.o \
                           subsurface_simulation.o surface_factory.o \
                           surface_field.o surface_flow.o \
                           surface_realization.o surface_simulation.o \
                           surface_th.o surfsubsurface_simulation.o \
                           timestepper_base.o unstructured_cell.o \
                           unstructured_grid.o unstructured_grid_aux.o
surfsubsurface_simulation.o : logging.o option.o output.o pflotran_constants.o \
                              pmc_base.o pmc_subsurface.o pmc_surface.o \
                              realization.o regression.o simulation_aux.o \
                              simulation_base.o subsurface_simulation.o \
                              surface_realization.o surface_simulation.o \
                              timestepper_base.o
th.o : condition.o connection.o coupler.o debug.o discretization.o eos_water.o \
       field.o fluid.o global_aux.o grid.o material.o material_aux.o option.o \
       output_aux.o patch.o pflotran_constants.o realization.o region.o \
       saturation_function.o secondary_continuum.o secondary_continuum_aux.o \
       string.o th_aux.o utility.o variables.o
th_aux.o : eos_water.o global_aux.o material_aux.o option.o \
           pflotran_constants.o saturation_function.o
time_storage.o : option.o pflotran_constants.o
timestepper.o : checkpoint.o convergence.o discretization.o field.o \
                global.o global_aux.o grid.o input_aux.o \
                logging.o material.o material_aux.o \
                option.o output.o output_aux.o patch.o pflotran_constants.o \
                reaction.o reaction_aux.o reactive_transport.o \
                reactive_transport_aux.o realization.o richards.o \
                secondary_continuum_aux.o solver.o string.o th.o utility.o \
                variables.o waypoint.o
timestepper_BE.o : convergence.o input_aux.o option.o output.o \
                   pflotran_constants.o pm_base.o solver.o string.o \
                   timestepper_base.o utility.o
timestepper_base.o : input_aux.o option.o output.o pflotran_constants.o \
                     pm_base.o string.o waypoint.o
timestepper_surface.o : option.o output.o pflotran_constants.o pm_base.o \
                        pm_surface_flow.o solver.o surface_flow.o \
                        timestepper_base.o waypoint.o
transport.o : connection.o global_aux.o material_aux.o matrix_block_aux.o \
              option.o pflotran_constants.o reactive_transport_aux.o
uniform_velocity.o : input_aux.o logging.o option.o pflotran_constants.o \
                     string.o units.o utility.o
units.o : option.o pflotran_constants.o
unstructured_cell.o : option.o pflotran_constants.o utility.o
unstructured_communicator.o : communicator_base.o dm_kludge.o \
                              pflotran_constants.o unstructured_explicit.o \
                              unstructured_grid.o unstructured_grid_aux.o
unstructured_explicit.o : connection.o geometry.o input_aux.o option.o \
                          pflotran_constants.o string.o \
                          unstructured_grid_aux.o utility.o
unstructured_grid.o : connection.o geometry.o hdf5_aux.o input_aux.o option.o \
                      pflotran_constants.o string.o unstructured_cell.o \
                      unstructured_grid_aux.o utility.o
unstructured_grid_aux.o : geometry.o option.o pflotran_constants.o \
                          unstructured_cell.o utility.o
unstructured_polyhedra.o : connection.o geometry.o input_aux.o option.o \
                           pflotran_constants.o string.o unstructured_cell.o \
                           unstructured_grid.o unstructured_grid_aux.o \
                           utility.o
utility.o : input_aux.o option.o pflotran_constants.o string.o
variables.o :
waypoint.o : option.o output_aux.o pflotran_constants.o utility.o

# developer level regression testing
test : utest rtest

utest : pflotran
	@if [ -d $(unit_test_dir) ]; then \
		$(MAKE) --directory=$(unit_test_dir) test; \
	fi

rtest : pflotran
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) test; \
	fi

clean-tests : FORCE
	@if [ -d $(unit_test_dir) ]; then \
		$(MAKE) --directory=$(unit_test_dir) $(MAKECMDGOALS); \
	fi
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) $(MAKECMDGOALS); \
	fi

# user level tests to verify pflotran is built correctly
check : FORCE
	@if [ -d $(regression_test_dir) ]; then \
		$(MAKE) --directory=$(regression_test_dir) $(MAKECMDGOALS); \
	fi

# null rule to force things to happen
FORCE :

