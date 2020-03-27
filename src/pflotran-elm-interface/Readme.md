# INSTRUCTIONS #

This is the interface, including minimal PFLOTRAN codes, for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. 

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N/REDOX biogeochemistry, i.e. PFLOTRAN.

***This interface is specifically released for ```PETSC v.3.9x or above```, and ```E3SM v1.1 or above```. *** 
*(NOTE:  if would like a test for BGC and TH coupling, branch 'default-release-v3.7' may be what you're look for first)*


## How do I get set up? ##

**(1)** *git clone the repository*, if not yet properly updating codes. (OPTIONAL)
```
git clone https://github.com/fmyuan/pflotran-elm-interface.git
```

*NOTE: when git clone, checkout branch 'default-release-v?.??'*, if need early version of PFLOTRAN with 'v?.??' corresponding to preferred PETSc version. 


**(2)** *if coupling CLM with PFLOTRAN, need to build a library named as **libpflotran.a**.*
```
cd ./src/pflotran-elm-interface
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
make PETSC_DIR=$PETSC_DIR column_mode=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR column_mode=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0' then built codes can be debugged)

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

*II.* **Macro.cmake** (E3SM master since 2019-07)
```
if("${MODEL}" STREQUAL "clm")
set(FFLAGS "${FFLAGS}  $(PFLOTRAN_INC)")
endif()

......

if("${MODEL}" STREQUAL "driver")
......
set(LDFLAGS "${LDFLAGS}  $(PFLOTRAN_LIB)")
endif()

```

*III.* **config_machines.xml** editing ```FFLAGS``` and ```LDFLAGS``` for all or specific compilers. *NOTE*: if this added, No need to modify 'Macro' or 'Macro.make' under case directory.

```
......

<FFLAGS>
  <!-- A NOTE here: $(PFLOTRAN_INC) may contain both PETSC and actual PFLOTRAN include dir, or only PETSC include dir -->
  <append MODEL="clm"> $(PFLOTRAN_INC) </append>
</FFLAGS>

......
<LDFLAGS>
  <!-- A NOTE here: $(PFLOTRAN_LIB) may contain both PETSC libraries and actual PFLOTRAN libray, or only PETSC libraries -->
  <append MODEL="driver"> $(PFLOTRAN_LIB) </append>
</LDFLAGS>


```

*IV.* **config_machines.xml** editing for each supported machine (example CADES at ORNL). *NOTE*: IF NOT, after './case.setup', edit 'env_mach_specific.xml' to turn on options.

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
   -->
   <environment_variables>
     <!-- The following pflotran is with PETSc-v3.8.x or above on CADES-->
     <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-dev/src/pflotran-elm-interface</env>

     <!-- by blanking the following 2 names, PETSC libs excluded in e3sm.exe when NOT coupling with PFLOTRAN -->
     <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_DIR}/include</env>
     <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_DIR}/lib -lpetsc -lmetis -lparmetis</env>
   </environment_variables>

```
**NOTE:** You must be sure that $CLM_PFLOTRAN_SOURCE_DIR and $PETSC_PATH are defined in your bash environment setting. Of course the ```libpflotran.a``` are prebuilt as mentioned above. If you DON'T want to include this library in your e3sm.exe (and of course no use of coupled models ), edit 'env_mach_specific.xml' as following:
```
<environment_variables>
<!--
  <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-dev/src/pflotran-elm-interface</env>
  <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_PATH}/include</env>
  <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_PATH}/lib -lpetsc -lmetis -lparmetis</env>
-->
</environment_variables>
```


***UPDATED: 2019-03-26***
