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
 

### How do I get set up? ###

**(1)** *git clone the repository*
```
git clone https://github.com/fmyuan/pflotran-elm-interface.git
```

**(2)[OPTIONAL]** *git checkout your specific branch, e.g. 'default-release-v3.7'. The default is 'pflotran-elm-interface'*
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
make PETSC_DIR=$PETSC_DIR pflotran
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

***SECONDLY***, build the library (vertically-only mode ON)
```
make PETSC_DIR=$PETSC_DIR column_mode=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR th_characteristic_curves=TRUE smoothing2=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0')

```

***FINALLY***, build CLM (ELM) with this library.

*I.* **Macro (CLM)** or **Macro.make (ELM)** or **Macro.cmake (ELM master since@2019-07)** modified, BY any means, to include $PFLTRAN_INC and to link $PFLOTRAN_LIB -

```
(Macro/Macro.make)


  ifeq ($(MODEL),clm)
    FFLAGS := $(FFLAGS)  $(PFLOTRAN_INC)
  endif
  
  ......
  
  ifeq ($(MODEL),driver)
    ......
    LDFLAGS := $(LDFLAGS)  $(PFLOTRAN_LIB)
  endif

```
**(ELM master since@2019-07)**
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
*II.* **config_compilers.xml** editing for each supported machine. *NOTE*: after './case.setup', edit **'env_mach_specific.xml'** to modify PETSC_PATH (or PETSC_DIR), CLM_PFLOTRAN_SOURCE_DIR, as **user-defined**.

```
      <!-- for CLM-PFLOTRAN coupling, the PETSC_PATH must be defined specifically upon machines -->
      <environment_variables>
        <env name="PETSC_PATH" compiler="gnu" mpilib="openmpi">/software/user_tools/current/cades-ccsi/petsc4pf/openmpi-1.10-gcc-5.3</env>      
        <!-- hack for PFLOTRAN coupling (this is a temporary solution, and user must manually edit env_mach_specific.xml after case.setup, IF needed)-->
        <env name="CLM_PFLOTRAN_SOURCE_DIR">/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-interface/src/clm-pflotran</env>
        <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_DIR}/include</env>
        <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_DIR}/lib -lpetsc -lmetis -lparmetis</env>
      </environment_variables>       

```



*III.* Specifically, for **ELM** build with PFLOTRAN, as external module.
(e.g. https://github.com/fmyuan/E3SM.git, branch 'elm-pflotran-II')
In this way, ELM source codes will have copied all PFLOTRAN source codes ONLY into:
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



***UPDATED: 2019-08-06***
