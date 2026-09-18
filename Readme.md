# README #

This is repository for coupling PFLOTRAN into CLM/ELM in the DOE-sponsored NGEE-Arctic Project. It keeps updates from PFLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The CLM/ELM portion lives in the E3SM repository (e.g. https://github.com/fmyuan/E3SM.git, branch `elm-pflotran-II`).

***UPDATE: (2019-04-05) The coupling with E3SM Land Model (ELM) now is working.***

### Remotes ###

This repository tracks several upstream/fork remotes:

* `origin` — https://github.com/fmyuan/pflotran-elm-interface.git (this fork; default branch: `pflotran-elm-interface`)
* `bsulman` — https://github.com/bsulman/pflotran-elm-interface.git (collaborator fork, active development)
* `pflotran_bitbucket` — https://bitbucket.org/pflotran/pflotran.git (upstream PFLOTRAN-Dev)

### Branches and Versions ###

The repository carries `master` plus a number of long-lived branches for different coupling targets and PETSc/PFLOTRAN versions. The most actively maintained branches are:

***pflotran-v5.0.0-elm-bgc*** (most recent activity)
 - Latest development branch, coupling against PFLOTRAN v5.0.0-era code for ELM BGC.

***pflotran-elm-bgc***
 - Development branch for ELM BGC coupling.

***pflotran-elm-interface*** (`origin`'s default branch)
 - The long-standing development version of clm-pflotran, for testing thermal-hydrology (TH) and biogeochemistry (C/BGC).
 - Keeps updating with `master` and with recent PETSc.

***bsulman/pflotran-elm-interface***
 - Collaborator (bsulman) fork/branch of the interface, used for active co-development.

***master***
 - The most updated PFLOTRAN codes forked from https://bitbucket.org/pflotran/pflotran, with only a few minor changes.
 - NOT suggested for coupled runs — use as a stand-alone PFLOTRAN instead.

Older, less actively maintained release branches (kept for reference/reproducibility of past studies):
`default-release-v2020`, `default-release-v2021`, `default-release-v2021a`, `default-release-v3.7`, `default-release-v3.8`, `default-release-v3.9`, `default-release-v3.10`, `default-release-v3.11`, `default-release-v3.12`, `simpleTH`, `THonly`, `Water3PhaseFlow`.

Check `git log -1 <branch>` for a given branch's last update date, and its own `Readme.md`/`readme.rst` (where present) for version-specific notes, before using it.


### How do I get set up? ###

**(1)** *git clone the repository*
```
git clone https://github.com/fmyuan/pflotran-elm-interface.git
```

**(2)[OPTIONAL]** *git checkout your specific branch, e.g. 'pflotran-elm-bgc'. The default is 'pflotran-elm-interface'*
```
cd pflotran-elm-interface
git checkout pflotran-elm-bgc
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



***UPDATED: 2026-09-18***
