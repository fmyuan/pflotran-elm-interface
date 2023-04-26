module slatec_pchip_module
#include "petsc/finclude/petscsys.h"
implicit none

! **************************************************************************** !
!
! This module contains selected subroutines from the PCHIP package of the SLATEC
! public domain library. The code has been minimally modified to conform to 
! PFLOTRAN conventions. These modifications include:
! 1. Using variable type definitions in petsc/finclude/petscsys.h instead
!    native Fortran types.
! 2. Replacing fixed format continuation characters with free-format ampersands.
! 3. Removing declaration of local functions 
! 4. Replace DATA/SAVE statements with PARAMETER for constants
!
! XERMSG functionality has not been replicated, but could be integrated with
! existing PFLOTRAN error message routines.
!
! Imported SLATEC/PCHIP decks
! 1. PCHDOC
! 2. PCHIM
! 3. PCHST
! 4. CHFDV
! 5. PCHFD

! Importer: Matthew Paul
! Date:     04/25/2023
!
! **************************************************************************** !

private XERMSG
public PCHDOC
public PCHIM
public PCHST
public CHFDV
public PCHFD

contains

subroutine XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
  character(*) :: LIBRAR, SUBROU, MESSG
  PetscInt :: NERR, LEVEL
end subroutine XERMSG

!DECK PCHDOC
      SUBROUTINE PCHDOC
!***BEGIN PROLOGUE  PCHDOC
!***PURPOSE  Documentation for PCHIP, a Fortran package for piecewise
!            cubic Hermite interpolation of data.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A, Z
!***TYPE      ALL (PCHDOC-A)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, DOCUMENTATION,
!             MONOTONE INTERPOLATION, PCHIP,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!            PCHIP:  Piecewise Cubic Hermite Interpolation Package
!
!      This document describes the contents of PCHIP, which is a
!   Fortran package for piecewise cubic Hermite interpolation of data.
!   It features software to produce a monotone and "visually pleasing"
!   interpolant to monotone data.  As is demonstrated in Reference 4,
!   such an interpolant may be more reasonable than a cubic spline if
!   the data contains both "steep" and "flat" sections.  Interpola-
!   tion of cumulative probability distribution functions is another
!   application.  (See References 2-4 for examples.)
!
!
!      All piecewise cubic functions in PCHIP are represented in
!   cubic Hermite form; that is, f(x) is determined by its values
!   F(I) and derivatives D(I) at the breakpoints X(I), I=1(1)N.
!   Throughout the package a PCH function is represented by the
!   five variables  N, X, F, D, INCFD:
!     N     - number of data points;
!     X     - abscissa values for the data points;
!     F     - ordinates (function values) for the data points;
!     D     - slopes (derivative values) at the data points;
!     INCFD - increment between successive elements in the F- and
!             D-arrays (more on this later).
!   These appear together and in the same order in all calls.
!
!      The double precision equivalents of the PCHIP routines are
!   obtained from the single precision names by prefixing the
!   single precision names with a D.  For example, the double
!   precision equivalent of PCHIM is DPCHIM.
!
!      The contents of the package are as follows:
!
!   1. Determine Derivative Values.
!
!      NOTE:  These routines provide alternate ways of determining D
!             if these values are not already known.
!
!         PCHIM -- Piecewise Cubic Hermite Interpolation to Monotone
!               data.
!               Used if the data are monotonic or if the user wants
!               to guarantee that the interpolant stays within the
!               limits of the data.  (See Reference 3.)
!
!         PCHIC -- Piecewise Cubic Hermite Interpolation Coefficients.
!               Used if neither of the above conditions holds, or if
!               the user wishes control over boundary derivatives.
!               Will generally reproduce monotonicity on subintervals
!               over which the data are monotonic.
!
!         PCHSP -- Piecewise Cubic Hermite Spline.
!               Produces a cubic spline interpolator in cubic Hermite
!               form.  Provided primarily for easy comparison of the
!               spline with other piecewise cubic interpolants.  (A
!               modified version of de Boor's CUBSPL, Reference 1.)
!
!   2. Evaluate, Differentiate, or Integrate Resulting PCH Function.
!
!      NOTE:  If derivative values are available from some other
!             source, these routines can be used without calling
!             any of the previous routines.
!
!         CHFEV -- Cubic Hermite Function EValuator.
!               Evaluates a single cubic Hermite function at an array
!               of points.  Used when the interval is known, as in
!               graphing applications.  Called by PCHFE.
!
!         PCHFE -- Piecewise Cubic Hermite Function Evaluator.
!               Used when the interval is unknown or the evaluation
!               array spans more than one data interval.
!
!         CHFDV -- Cubic Hermite Function and Derivative Evaluator.
!               Evaluates a single cubic Hermite function and its
!               first derivative at an array of points.  Used when
!               the interval is known, as in graphing applications.
!               Called by PCHFD.
!
!         PCHFD -- Piecewise Cubic Hermite Function and Derivative
!               Evaluator.
!               Used when the interval is unknown or the evaluation
!               array spans more than one data interval.
!
!         PCHID -- Piecewise Cubic Hermite Integrator, Data Limits.
!               Computes the definite integral of a piecewise cubic
!               Hermite function when the integration limits are data
!               points.
!
!         PCHIA -- Piecewise Cubic Hermite Integrator, Arbitrary Limits.
!               Computes the definite integral of a piecewise cubic
!               Hermite function over an arbitrary finite interval.
!
!   3. Utility routines.
!
!         PCHBS -- Piecewise Cubic Hermite to B-Spline converter.
!               Converts a PCH function to B-representation, so that
!               it can be used with other elements of the B-spline
!               package (see BSPDOC).
!
!         PCHCM -- Piecewise Cubic Hermite, Check Monotonicity of.
!               Checks the monotonicity of an arbitrary PCH function.
!               Might be used with PCHSP to build a polyalgorithm for
!               piecewise C-2 interpolation.
!
!   4. Internal routines.
!
!         CHFIE -- Cubic Hermite Function Integral Evaluator.
!               (Real function called by PCHIA.)
!
!         CHFCM -- Cubic Hermite Function, Check Monotonicity of.
!               (Integer function called by PCHCM.)
!
!         PCHCE -- PCHIC End Derivative Setter.
!               (Called by PCHIC.)
!
!         PCHCI -- PCHIC Initial Derivative Setter.
!               (Called by PCHIC.)
!
!         PCHCS -- PCHIC Monotonicity Switch Derivative Setter.
!               (Called by PCHIC.)
!
!         PCHDF -- PCHIP Finite Difference Formula.
!               (Real function called by PCHCE and PCHSP.)
!
!         PCHST -- PCHIP Sign Testing Routine.
!               (Real function called by various PCHIP routines.)
!
!         PCHSW -- PCHCS Switch Excursion Adjuster.
!               (Called by PCHCS.)
!
!   The calling sequences for these routines are described in the
!   prologues of the respective routines.
!
!
!      INCFD, the increment between successive elements in the F-
!   and D-arrays is included in the representation of a PCH function
!   in this package to facilitate two-dimensional applications.  For
!   "normal" usage INCFD=1, and F and D are one-dimensional arrays.
!   one would call PCHxx (where "xx" is "IM", "IC", or "SP") with
!
!              N, X, F, D, 1  .
!
!   Suppose, however, that one has data on a rectangular mesh,
!
!         F2D(I,J) = value at (X(I), Y(J)),  I=1(1)NX,
!                                            J=1(1)NY.
!   Assume the following dimensions:
!
!         REAL  X(NXMAX), Y(NYMAX)
!         REAL  F2D(NXMAX,NYMAX), FX(NXMAX,NYMAX), FY(NXMAX,NYMAX)
!
!   where  2.LE.NX.LE.NXMAX AND 2.LE.NY.LE.NYMAX .  To interpolate
!   in X along the line  Y = Y(J), call PCHxx with
!
!              NX, X, F2D(1,J), FX(1,J), 1  .
!
!   To interpolate along the line X = X(I), call PCHxx with
!
!              NY, Y, F2D(I,1), FY(I,1), NXMAX  .
!
!   (This example assumes the usual columnwise storage of 2-D arrays
!    in Fortran.)
!
!***REFERENCES  1. Carl de Boor, A Practical Guide to Splines, Springer-
!                 Verlag, New York, 1978 (esp. Chapter IV, pp.49-62).
!               2. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
!                 Package, Report UCRL-87285, Lawrence Livermore Natio-
!                 nal Laboratory, July 1982.  [Poster presented at the
!                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
!               3. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               4. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811106  DATE WRITTEN
!   870930  Updated Reference 3.
!   890414  Changed PCHMC and CHFMC to PCHCM and CHFCM, respectively,
!           and augmented description of PCHCM.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   910826  1. Revised purpose, clarified role of argument INCFD,
!              corrected error in example, and removed redundant
!              reference list.
!           2. Added description of PCHBS.  (FNF)
!   920429  Revised format and order of references.  (WRB,FNF)
!   930505  Changed CHFIV to CHFIE.  (FNF)
!***END PROLOGUE  PCHDOC
!-----------------------------------------------------------------------
!     THIS IS A DUMMY SUBROUTINE, AND SHOULD NEVER BE CALLED.
!
!***FIRST EXECUTABLE STATEMENT  PCHDOC
      RETURN
!------------- LAST LINE OF PCHDOC FOLLOWS -----------------------------
      END

!DECK PCHIM
      SUBROUTINE PCHIM (N, X, F, D, INCFD, IERR)
!***BEGIN PROLOGUE  PCHIM
!***PURPOSE  Set derivatives needed to determine a monotone piecewise
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See PCHIC if user control is
!            desired over boundary or switch conditions.)
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      SINGLE PRECISION (PCHIM-S, DPCHIM-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubic
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See PCHIC if user control of boundary con-
!     ditions is desired.)
!
!     If the data are only piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See PCHIC if user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by PCHFE or PCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  PCHIM (N, X, F, D, INCFD, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of dependent variable values to be inter-
!           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
!           PCHIM is designed for monotonic data, but it will work for
!           any F-array.  It will force extrema at points where mono-
!           tonicity switches direction.  If some other treatment of
!           switch points is desired, PCHIC should be used instead.
!                                     -----
!     D -- (output) real array of derivative values at the data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               2. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  PCHST, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820201  1. Introduced  PCHST  to reduce possible over/under-
!             flow problems.
!           2. Rearranged derivative formula for same reason.
!   820602  1. Modified end conditions to be continuous functions
!             of data when monotonicity switches in next interval.
!           2. Modified formulas so end conditions are less prone
!             of over/underflow problems.
!   820803  Minor cosmetic changes for release 1.
!   870813  Updated Reference 1.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  PCHIM
!  Programming notes:
!
!     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. To produce a double precision version, simply:
!        a. Change PCHIM to DPCHIM wherever it occurs,
!        b. Change PCHST to DPCHST wherever it occurs,
!        c. Change all references to the Fortran intrinsics to their
!           double precision equivalents,
!        d. Change the real declarations to double precision, and
!        e. Change the constants ZERO and THREE to double precision.
!
!  DECLARE ARGUMENTS.
!
      PetscInt :: N, INCFD, IERR
      PetscReal :: X(*), F(INCFD,*), D(INCFD,*)
!
!  DECLARE LOCAL VARIABLES.
!
      PetscInt :: I, NLESS1
      PetscReal :: DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE, &
            H1, H2, HSUM, HSUMT3, W1, W2
      PetscReal, parameter :: ZERO = 0d0, THREE = 3d0
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
!
!  NORMAL CASE  (N .GE. 3).
!
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
      ENDIF
!
!  LOOP THROUGH INTERIOR POINTS.
!
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
!
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
         D(1,I) = ZERO
         IF ( PCHST(DEL1,DEL2) )  42, 41, 45
!
!        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( PCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
!
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = MAX( ABS(DEL1), ABS(DEL2) )
         DMIN = MIN( ABS(DEL1), ABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
!
   50 CONTINUE
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
      ENDIF
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'PCHIM', &
         'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'PCHIM', 'INCREMENT LESS THAN ONE', IERR, &
         1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'PCHIM', 'X-ARRAY NOT STRICTLY INCREASING' &
         , IERR, 1)
      RETURN
!------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END

!DECK PCHST
      REAL FUNCTION PCHST (ARG1, ARG2)
!***BEGIN PROLOGUE  PCHST
!***SUBSIDIARY
!***PURPOSE  PCHIP Sign-Testing Routine
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHST-S, DPCHST-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         PCHST:  PCHIP Sign-Testing Routine.
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
!     The object is to do this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!  Fortran intrinsics used:  SIGN.
!
!***SEE ALSO  PCHCE, PCHCI, PCHCS, PCHIM
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHST
!
!**End
!
!  DECLARE ARGUMENTS.
!
      PetscReal :: ARG1, ARG2
!
!  DECLARE LOCAL VARIABLES.
!
      PetscReal, parameter :: ONE = 1d0, ZERO = 0d0
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
!
      RETURN
!------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END

!DECK CHFDV
      SUBROUTINE CHFDV (X1, X2, F1, F2, D1, D2, NE, XE, FE, DE, NEXT, IERR)
!***BEGIN PROLOGUE  CHFDV
!***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
!            first derivative at an array of points.  While designed for
!            use by PCHFD, it may be useful directly as an evaluator
!            for a piecewise cubic Hermite function in applications,
!            such as graphing, where the interval is known in advance.
!            If only function values are required, use CHFEV instead.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H1
!***TYPE      SINGLE PRECISION (CHFDV-S, DCHFDV-D)
!***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION, PCHIP
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!        CHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use CHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
!
!        CALL  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1.EQ.X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real array of values of the cubic function defined
!           by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1.EQ.X2 .
!                (Output arrays have not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811019  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CHFDV
!  Programming notes:
!
!     To produce a double precision version, simply:
!        a. Change CHFDV to DCHFDV wherever it occurs,
!        b. Change the real declaration to double precision, and
!        c. Change the constant ZERO to double precision.
!
!  DECLARE ARGUMENTS.
!
      PetscInt :: NE, NEXT(2), IERR
      PetscReal :: X1, X2, F1, F2, D1, D2, XE(*), FE(*), DE(*)
!
!  DECLARE LOCAL VARIABLES.
!
      PetscInt :: I
      PetscReal :: C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA
      PetscReal, parameter :: ZERO = 0d0
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  CHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
!
!  INITIALIZE.
!
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = MIN(ZERO, H)
      XMA = MAX(ZERO, H)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
!                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
!
!  EVALUATION LOOP.
!
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
!          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
!
!  NORMAL RETURN.
!
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'CHFDV', &
         'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
      RETURN
!
 5002 CONTINUE
!     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'CHFDV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)
      RETURN
!------------- LAST LINE OF CHFDV FOLLOWS ------------------------------
      END

!DECK PCHFD
      SUBROUTINE PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!***BEGIN PROLOGUE  PCHFD
!***PURPOSE  Evaluate a piecewise cubic Hermite function and its first
!            derivative at an array of points.  May be used by itself
!            for Hermite interpolation, or as an evaluator for PCHIM
!            or PCHIC.  If only function values are required, use
!            PCHFE instead.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H1
!***TYPE      SINGLE PRECISION (PCHFD-S, DPCHFD-D)
!***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
!             HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC EVALUATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          PCHFD:  Piecewise Cubic Hermite Function and Derivative
!                  evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!     gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use PCHFE, instead.
!
!     To provide compatibility with PCHIM and PCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
!        LOGICAL  SKIP
!
!        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in PCHIM or PCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real array of values of the cubic Hermite function
!           defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine CHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CHFDV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811020  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870707  Minor cosmetic changes to prologue.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  PCHFD
!  Programming notes:
!
!     1. To produce a double precision version, simply:
!        a. Change PCHFD to DPCHFD, and CHFDV to DCHFDV, wherever they
!           occur,
!        b. Change the real declaration to double precision,
!
!     2. Most of the coding between the call to CHFDV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. CHFDV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of PCHFD that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points .LT.
!        X(N-1), followed by points .GT.X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.
!
!  DECLARE ARGUMENTS.
!
      PetscInt :: N, INCFD, NE, IERR
      PetscReal :: X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*), DE(*)
      PetscBool :: SKIP
!
!  DECLARE LOCAL VARIABLES.
!
      PetscInt :: I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  PCHFD
      IF (SKIP)  GO TO 5
!
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
         IF (JFIRST .GT. NE)  GO TO 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
!
   40    CONTINUE
         NJ = J - JFIRST
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
         IF (NJ .EQ. 0)  GO TO 50
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
        CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
                    NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
!       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
!
         IF (NEXT(2) .EQ. 0)  GO TO 42
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
            IF (IR .LT. N)  GO TO 41
!           IF (IR .EQ. N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
!           ENDIF
!        ENDIF
   42    CONTINUE
!
         IF (NEXT(1) .EQ. 0)  GO TO 49
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
            IF (IR .GT. 2)  GO TO 43
!           IF (IR .EQ. 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN CHFDV.
               GO TO 5005
!
   45          CONTINUE
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!              AT THIS POINT, EITHER  XE(J) .LT. X(1)
!                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
               IR = MAX(1, I-1)
!           ENDIF
!        ENDIF
   49    CONTINUE
!
         JFIRST = J
!
!     END OF IR-LOOP.
!
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
!
!  NORMAL RETURN.
!
 5000 CONTINUE
      RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'PCHFD', &
         'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
      RETURN
!
 5002 CONTINUE
!     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'PCHFD', 'INCREMENT LESS THAN ONE', IERR, &
         1)
      RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERMSG ('SLATEC', 'PCHFD', 'X-ARRAY NOT STRICTLY INCREASING' &
         , IERR, 1)
      RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
      IERR = -4
      CALL XERMSG ('SLATEC', 'PCHFD', &
         'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
      RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM CHFDV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      CALL XERMSG ('SLATEC', 'PCHFD', &
         'ERROR RETURN FROM CHFDV -- FATAL', IERR, 2)
      RETURN
!------------- LAST LINE OF PCHFD FOLLOWS ------------------------------
      END

end module slatec_pchip_module
