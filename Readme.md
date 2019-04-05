# README #

This is repository for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. It's keeping updates from PLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The CLM portion is in E3SM repository, branch /fmyuan/lnd/elm-pflotran.

***UPDATE: (2019-04-05) The coupling with E3SM Land Model (ELM) now is working.***

### Branches and Versions ###

* The repository contains master and FIVE(5) branches for different versions of PETSc:

***master*** 

 - the most updated PFLOTRAN codes forked from https://bitbucket.org/pflotran/pflotran
 - only a few minor changes. 
 - NOT suggested to use coupling run, rather as a stand-alone PFLOTRAN. 
 - updated: **2019-04-01**

***(1) pflotran-elm-interface*** 

 - the current development version of clm-pflotran, specifically for testing all thermal-hydrology (TH) and biogeochemistry (C) portion. 
 - **PFLOTRAN**: It's keeping updating with master, and with most recent ***PETSc-dev***

***(2)default-release-v3.7***
 - The current STABLE version of clm-pflotran, specifically for testing all thermal-hydrology (TH) and biogeochemistry (BGC) portion. 
 - **PFLOTRAN**: 6571d9ae3d1cdd699bfc993155c958173e76c20f [6571d9a]. It's older (**Updated: 2017-03-09**) than 'default', and with ***PETSc-dev*** version **3.7.x**.
 - Stable CLM Version: **CLM4_5_35**. (2017-04-25).  *STATUS*: TH & BGC coupling is STABLE. 
 - Target ELM Version: **ELM v1**.  *STATUS*: STABLE for ```BGC```. (2017-06-27); UNDER test for ```TH mode``` 

***(3)default-release-v3.8***
 - for BGC coupled PFLOTRAN, with PETSc 3.8.x. TH coupling is not available.
 - It's working with **ELM v1.1 above**. May not be working with CLM.
 
***(4)default-release-v3.9***
 - for BGC coupled PFLOTRAN, with PETSc 3.9.x. TH coupling is not available.
 - It's working with **ELM v1.1 above**. May not be working with CLM.
 

***(3) default-light*** 

 - **Light-weight version** of PFLOTRAN-elm-interface (Currently updated on Jan-14-2019).  The purpose is to slim the full package of PFLOTRAN-dev for ONLY **THREE (3) MODES:** Richards-MODE, TH-MODE, and CHEMISTRY (i.e. BGC).


### How do I get set up? ###

**(1)** *git clone the repository*
```
git clone https://github.com/fmyuan/pflotran-elm-interface.git
```

**(2)[OPTIONAL]** *git checkout your specific branch, e.g. 'default-release-v3.7'*
```
cd pflotran-interface
git checkout default-release-v3.7
```

**(3a)** *if not coupled with ELM, this repository should be a **stand-alone PFLOTRAN** model.*

- the Source code directory: 
```
cd ./src/pflotran
```

- build it by issuing command:
```
make PETSC_DIR=$PETSC_DIR th_characteristic_curves=TRUE pflotran
(where $PETSC_DIR is your PETSC_DIR directory)
```

- check your repository (regression test):
```
make PETSC_DIR=$PETSC_DIR test
(where $PETSC_DIR is your PETSC_DIR directory)
```

**(3b)** *if coupling ELM with PFLOTRAN, need to build a library named as **libpflotran.a**.*
```
cd ./src/pflotran-elm-interface
```

***FIRST***, run specific make script to softlink PFLOTRAN source code files (*.F90).
```
make link_common_src

(OR, make clean_common_src 
to unlink PFLOTRAN source code files, and only leave PFLOTRAN-ELM interface codes)
```

***SECONDLY***, build the library
```
make PETSC_DIR=$PETSC_DIR th_characteristic_curves=TRUE smoothing2=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR th_characteristic_curves=TRUE smoothing2=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0')

```

***FINALLY***, build CLM (ELM) with this library.

*I.* **Macro (CLM)** or **Macro.make (ELM)** modified, BY any means, for coupling build -
```
ifeq ($(MODEL), clm) 
  ifeq ($(CLM_PFLOTRAN_COLMODE), TRUE) 
    ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
       CPPDEFS +=  -DCOLUMN_MODE 
    endif
  endif

  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     FFLAGS +=  -I$(CLM_PFLOTRAN_SOURCE_DIR)
     CPPDEFS +=  -DCLM_PFLOTRAN 
  endif
endif

......

ifeq ($(MODEL), driver) 
  ifeq ($(CLM_PFLOTRAN_COUPLED), TRUE) 
     LDFLAGS +=  -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB)
  endif
endif

```

*NOTE*: Modified Macro.make above requires ***4*** alias, setted as following.

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

  # Get the "PETSC_LIB" list an env var (UPDATED: 2017-May-19)
  # include $(PETSC_PATH)/conf/variables
  # (1) petsc-git-version: 1a9d3c3c50abf60098813fdf7291fe3540415115
  # (2) in this petsc package, "PETSC_LIB" contains "-L$LIB_PETSC"
  include $(PETSC_PATH)/lib/petsc/conf/variables
  
endif

```

*III.* **config_compilers.xml** editing for all or specific compilers. *NOTE*: if this added, No need to modify 'Macro' or 'Macro.make' under case directory. 

```
<!-- hacking of mach/compiler generated 'Macros.make' for coupling with pflotran -->
<!-- ideally it should go with CLM configuration -->
<compiler>
  <ADD_FFLAGS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE"> -I$(CLM_PFLOTRAN_SOURCE_DIR)</ADD_FFLAGS>
  <ADD_CPPDEFS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE"> -DCLM_PFLOTRAN </ADD_CPPDEFS>
  <ADD_CPPDEFS MODEL="clm" CLM_PFLOTRAN_COUPLED="TRUE" CLM_PFLOTRAN_COLMODE="TRUE"> -DCOLUMN_MODE </ADD_CPPDEFS>
  <ADD_LDFLAGS MODEL="driver" CLM_PFLOTRAN_COUPLED="TRUE"> -L$(CLM_PFLOTRAN_SOURCE_DIR) -lpflotran $(PETSC_LIB)</ADD_LDFLAGS>
</compiler>
<!-- end of hacking 'Macros.make' for coupling with pflotran -->

```


*IV.* **config_machines.xml** editing for each supported machine. *NOTE*: after './case.setup', edit 'env_mach_specific.xml' to turn on options.

```
      <!-- for CLM-PFLOTRAN coupling, the PETSC_PATH must be defined specifically upon machines -->
      <environment_variables>
        <env name="PETSC_PATH" compiler="gnu" mpilib="openmpi">/software/user_tools/current/cades-ccsi/petsc4pf/openmpi-1.10-gcc-5.3</env>      
        <!-- hack for PFLOTRAN coupling (this is a temporary solution, and user must manually edit env_mach_specific.xml after case.setup)-->
        <env name="CLM_PFLOTRAN_COUPLED">FALSE</env>
        <env name="CLM_PFLOTRAN_COLMODE">FALSE</env>
        <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-interface/src/clm-pflotran</env>
      </environment_variables>       

```



*V.* Specifically, for **ELM** build with PFLOTRAN, as external module.
(TODO)
The plan is to copy all source codes ONLY into:
```
      $SRCROOT/components/clm/src/external_models/pflotran-interface/src
```
AND, add the following line in ELM building config file ($SRCROOT/components/clm/bld/configure):
```
       "external_models/sbetr/src/Applications/soil-farm/CENT_ECACNP",
+	     "external_models/pflotran-interface/src/clm-pflotran",
		     "utils", 
		     "cpl" );
```

If setting software environments (specifically PETSc library), the whole ELM building process will automatically setup source codes and dependecies for coupled component, and e3sm.exe usable for either standalone or coupling simulation.



***UPDATED: 2019-04-05***
