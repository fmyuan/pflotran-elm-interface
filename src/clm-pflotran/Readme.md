# INSTRUCTIONS #

This is the interface, including minimal PFLOTRAN codes, for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. 

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

***This interface is specifically released for ```PETSC v.3.9x or above```, and ```E3SM v1.1 or above```. ***


## How do I get set up? ##

**(1)** *git clone the repository*, if not yet properly updating codes. (OPTIONAL)
```
git clone https://github.com/fmyuan/pflotran-elm-interface/src/clm-pflotran
```

*NOTE: when git clone, checkout branch 'default-release-v3.9'*


**(2)** *if coupling CLM with PFLOTRAN, need to build a library named as **libpflotran.a**.*
```
cd ./src/clm-pflotran
```

***FIRST***, run makefile to copy or link needed PFLOTRAN source code files (*.F90), if PFLOTRAN codes missing. (OPTIONAL)
```
make copy_common_src

(OR, make link_common_src
to softlink needed PFLOTRAN source codes)

(OR, make clean_common_src
to clean PFLOTRAN source code files, and only leave CLM-PFLOTRAN interface codes, if needed)
```

***SECONDLY***, build the library
```
make PETSC_DIR=$PETSC_DIR smoothing2=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR smoothing2=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0')

```

***FINALLY***, build ELM v1.1 or above  with this library, as usual, BUT must do modifying ELM's makefile or Macro.make as following.




## A FEW Specific Notes on How to modify ELM Macro.make and other machine files##
*I.*  **Macro.make** modified for coupling build -
```
ifeq ($(MODEL), clm) 
  FFLAGS := $(FFLAGS) $(PFLOTRAN_INC)
endif

......

ifeq ($(MODEL), driver) 
  LDFLAGS := $(LDFLAGS) $(PFLOTRAN_LIB)
endif

```

*NOTE*: Modified Macro above requires ***2*** alias $(PFLOTRAN_INC) and $(PFLOTRAN_LIB), set (OR empty) as following.

*II.* **Makefile**
```
# Set PETSc info if it is being used
ifeq ($(strip $(USE_PETSC)), TRUE)
  ifdef PETSC_PATH
    ifndef INC_PETSC
      INC_PETSC:=$(PETSC_PATH)/include
    endif
    ifndef LIB_PETSC
      LIB_PETSC:=$(PETSC_PATH)/lib
    endif
  else
    $(error PETSC_PATH must be defined when USE_PETSC is TRUE)
  endif

  # Get the "PETSC_LIB" list an env var
  include $(PETSC_PATH)/lib/petsc/conf/variables

  # for specific CLM-PFLOTRAN coupling
  ifeq ($(MODEL),clm)
    ifeq ($(strip $(CLM_PFLOTRAN_COUPLED)), TRUE)
      ifneq ($(strip $(CLM_PFLOTRAN_SOURCE_DIR)),)
        CPPDEFS += -DCLM_PFLOTRAN
        ifeq ($(strip $(CLM_PFLOTRAN_COLMODE)), TRUE)
          CPPDEFS += -DCOLUMN_MODE
        endif
      endif
    endif
  endif

endif

```

*III.* **config_compilers.xml** editing ```FFLAGS``` and ```LDFLAGS``` for all or specific compilers. *NOTE*: if this added, No need to modify 'Macro' or 'Macro.make' under case directory.

```
<FFLAGS>
  <!-- -ffree-line-length-none and -ffixed-line-length-none need to be in FFLAGS rather than in FIXEDFLAGS/FREEFLAGS
       so that these are passed to cmake builds (cmake builds don't use FIXEDFLAGS and FREEFLAGS). -->
  <base> -O -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fno-range-check</base>
  <append> -I$(NETCDF_PATH)/include -I$(HDF5_PATH)/include </append>
  <!-- A NOTE here: $(PFLOTRAN_INC) may contain both PETSC and actual PFLOTRAN include dir, or only PETSC include dir, or blanked -->
  <append MODEL="clm"> $(PFLOTRAN_INC) </append>
</FFLAGS>

```
AND,

```
<LDFLAGS>
  <append MODEL="driver"> -L$NETCDF_PATH/lib -Wl,-rpath=$NETCDF_PATH/lib -lnetcdff -lnetcdf </append>
  <append MODEL="driver"> -L$HDF5_PATH/lib -Wl,-rpath=$HDF5_PATH/lib -lhdf5_hl -lhdf5 -lhdf5hl_fortran -lhdf5_fortran </append>
  <!-- A NOTE here: $(PFLOTRAN_LIB) may contain both PETSC libraries and actual PFLOTRAN libray, or only PETSC libraries, or blanked -->
  <append MODEL="driver"> $(PFLOTRAN_LIB) </append>
</LDFLAGS>

```

*IV.* **config_machines.xml** editing for each supported machine (example CADES at ORNL). *NOTE*: after './case.setup', edit 'env_mach_specific.xml' to turn on options.

```
   <!-- for CLM-PFLOTRAN coupling, the PETSC_PATH must be defined specifically upon machines, usually defined in .bashrc -->
   <!-- the following is PETSc v.3.8.x or above -->
   <environment_variables compiler="gnu" mpilib="openmpi">
     <env name="PETSC_PATH">/software/user_tools/current/cades-ccsi/petsc-x/openmpi-1.10-gcc-5.3</env> <!-- PETSc v3.8.x or above -->
   </environment_variables>

   <!-- hack for PFLOTRAN coupling to build ELM model.
        this is a temporary solution, and user must manually edit the following
        in 'env_mach_specific.xml' after case.setup,
        Otherwise, model will build/run as non-PFLOTRAN coupled.
        IF building with PFLOTRAN codes, in addition to the following options,
        you must modify 'CLM_USE_PETSC' in 'env_build.xml' to be 'TRUE'.
   -->
   <environment_variables>
     <env name="CLM_PFLOTRAN_COUPLED">FALSE</env>
     <env name="CLM_PFLOTRAN_COLMODE">FALSE</env>
     <!-- The following pflotran is with PETSc-v3.8.x or above on CADES-->
     <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-interface-v3.8x/src/clm-pflotran</env>

     <!-- by blanking the following 2 names, PETSC libs excluded in e3sm.exe when NOT coupling with PFLOTRAN -->
     <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_DIR}/include</env>
     <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_DIR}/lib -lpetsc -lmetis -lparmetis</env>
   </environment_variables>

```
**NOTE:** You must be sure that $CLM_PFLOTRAN_SOURCE_DIR and $PETSC_PATH are defined in your bash environment setting. Of course the ```libpflotran.a``` are prebuilt as mentioned above. If you DON'T want to include this library in your e3sm.exe (and of course no use of coupled models ), edit 'env_mach_specific.xml' as following:
```
<environment_variables>
  <env name="CLM_PFLOTRAN_COUPLED">FALSE</env>
  <env name="CLM_PFLOTRAN_COLMODE">FALSE</env>
  <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-interface-v3.8x/src/clm-pflotran</env>
<!--
  <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_PATH}/include</env>
  <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_PATH}/lib -lpetsc -lmetis -lparmetis</env>
-->
</environment_variables>
```


***UPDATED: 2019-01-16***
