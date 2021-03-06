#include "reflection.h"
module Test_UnixProcess_mod
   use TestSuite_mod
   use Assert_mod
   use Exception_mod
   use UnixProcess_mod
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('StringConversionUtilities')
#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

      ADD(testIsActive)
      ADD(testGetLine)
      ADD(testGetLine2)

   end function suite

   !------
   ! A bit self-referential, but at least it serves to drive
   ! development.Start a background command that persists.
   !------
   subroutine testIsActive()
      type (UnixProcess) :: process

      process = UnixProcess('sleep 10', runInBackground=.true.)
      call assertTrue(process%isActive(),'hmm')
      if (anyExceptions()) return

      call process%terminate()
      call assertFalse(process%isActive(),'huh')
      
   end subroutine testIsActive

   subroutine testGetLine()
      type (UnixProcess) :: process

      character(len=:), allocatable :: line

      process = UnixProcess('echo hello')
      line = process%getLine()
      call assertEqual('hello', trim(line))
      
   end subroutine testGetLine

   ! In this test, getLine is called twice because we first
   ! need to get the pid for the background process.
   subroutine testGetLine2()
      type (UnixProcess) :: process

      character(len=:), allocatable :: line

      process = UnixProcess('echo hello', runInBackground=.true.)
      line = process%getLine()

      call assertEqual('hello', trim(line))
      
   end subroutine testGetLine2

end module Test_UnixProcess_mod
