module Shape_class
  
  implicit none
  
  private

  type, public :: shape
  
  ! empty instance variables
  
  contains
  
    procedure, public :: CalculateArea => CalculateShapeArea
    procedure, public :: CalculatePerimeter => CalculateShapePerimeter
    procedure, public :: PrintDescription => PrintShapeDescription

  end type shape
  
contains

!******************************************************************************!

real function CalculateShapeArea(this)
  
  class(shape) :: this
  
  CalculateShapeArea = 0.
  
end function CalculateShapeArea

!******************************************************************************!

real function CalculateShapePerimeter(this)
  
  class(shape) :: this
  
  CalculateShapePerimeter = 0.
  
end function CalculateShapePerimeter

!******************************************************************************!

character(len=50) function PrintShapeDescription(this)
  
  class(shape) :: this
  
  PrintShapeDescription = ''
  
end function PrintShapeDescription

end module Shape_class

!******************************************************************************!
!******************************************************************************!

module Circle_class

  use Shape_class

  implicit none
  
  private

  type, public, extends(shape) :: circle_type
  
    real :: r = 0.
  
  contains
  
    procedure, public :: CalculateArea => CalculateCircleArea
    procedure, public :: CalculatePerimeter => CalculateCirclePerimeter
    procedure, public :: PrintDescription => PrintCircleDescription
  
  end type circle_type
  
  real, parameter :: PI = 3.141593

  public :: CreateCircle
  
contains

!******************************************************************************!

function CreateCircle(r)

  real :: r
  
  class(circle_type), pointer :: CreateCircle

  allocate(CreateCircle)
  CreateCircle%r = r

end function CreateCircle

!******************************************************************************!

real function CalculateCircleArea(this)

  class(circle_type) :: this
  
  CalculateCircleArea = PI * this%r**2.

end function CalculateCircleArea

!******************************************************************************!

real function CalculateCirclePerimeter(this)

  class(circle_type) :: this
  
  CalculateCirclePerimeter = 2. * PI * this%r

end function CalculateCirclePerimeter

!******************************************************************************!

character(len=50) function PrintCircleDescription(this)
  
  class(circle_type) :: this
  
  write(PrintCircleDescription,'(a,f6.2)') 'Circle of radius ', this%r
  
end function PrintCircleDescription

end module Circle_class

!******************************************************************************!
!******************************************************************************!

module Rectangle_class

  use Shape_class

  implicit none
  
  private

  type, public, extends(shape) :: rectangle_type
  
    real :: l = 0.
    real :: w = 0.
  
  contains
  
    procedure, public :: CalculateArea => CalculateRectangleArea
    procedure, public :: CalculatePerimeter => CalculateRectanglePerimeter
    procedure, public :: PrintDescription => PrintRectangleDescription
  
  end type rectangle_type

  public :: CreateRectangle

contains

!******************************************************************************!

function CreateRectangle(l,w)

  real :: l
  real :: w

  class(rectangle_type), pointer :: CreateRectangle
  
  allocate(CreateRectangle)
  CreateRectangle%l = l
  CreateRectangle%w = w

end function CreateRectangle

!******************************************************************************!

real function CalculateRectangleArea(this)

  class(rectangle_type) :: this
  
  CalculateRectangleArea = this%l * this%w

end function CalculateRectangleArea

!******************************************************************************!

real function CalculateRectanglePerimeter(this)

  class(rectangle_type) :: this
  
  CalculateRectanglePerimeter = 2. * this%l + 2. * this%w

end function CalculateRectanglePerimeter

!******************************************************************************!

character(len=50) function PrintRectangleDescription(this)
  
  class(rectangle_type) :: this
  
  write(PrintRectangleDescription,'(a,f6.2,a,f6.2)') 'Rectangle of length ',  &
    this%l, ' and width ', this%w
  
end function PrintRectangleDescription

end module Rectangle_class

!******************************************************************************!
!******************************************************************************!

module Square_class

  use Rectangle_class

  implicit none
  
  private

  type, public, extends(rectangle_type) :: square_type
  
  ! no new variables
  
  contains
  
    procedure, public :: PrintDescription => PrintSquareDescription
  
  end type square_type

  public :: CreateSquare
  
contains

!******************************************************************************!

function CreateSquare(l)

  real :: l

  type(square_type), pointer :: CreateSquare
  
  allocate(CreateSquare)
  CreateSquare%l = l
  CreateSquare%w = l

end function CreateSquare

!******************************************************************************!

character(len=50) function PrintSquareDescription(this)
  
  class(square_type) :: this
  
  write(PrintSquareDescription,'(a,f6.2)') 'Square of length ', this%l
  
end function PrintSquareDescription

end module Square_class

!******************************************************************************!
!******************************************************************************!
!******************************************************************************!

program test_shape

  use Circle_class
  use Rectangle_class
  use Shape_class
  use Square_class
  
  implicit none
  
  type(circle_type), pointer :: circle
  type(square_type), pointer :: square
  type(rectangle_type), pointer :: rectangle
  
  type :: shape_ptr
    class(shape), pointer :: p
  end type shape_ptr
  
  integer :: i
  type(shape_ptr) :: shapes(3)
  
  write(*,'(a)') 'Beginning of Fortran 2003 Shape test...'
  
  circle => CreateCircle(2.)
  square => CreateSquare(3.)
  rectangle => CreateRectangle(4.,1.)
  
  shapes(1)%p => circle
  shapes(2)%p => square
  shapes(3)%p => rectangle

  do i = 1, size(shapes)
    write(*,'(/a)') shapes(i)%p%PrintDescription()
    write(*,'(a,f8.4)') 'Area      = ', shapes(i)%p%CalculateArea()
    write(*,'(a,f8.4)') 'Perimeter = ', shapes(i)%p%CalculatePerimeter()
  enddo

  do i = 1, size(shapes)
    deallocate(shapes(i)%p)
    nullify(shapes(i)%p)
  enddo
  nullify(circle)
  nullify(square)
  nullify(rectangle)
  
  write(*,'(/a)') 'Successful completion of Fortran 2003 Shape test!'

end program test_shape
