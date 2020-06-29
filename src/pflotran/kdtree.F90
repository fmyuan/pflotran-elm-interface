! Author: Matthew Kennel, Institute for Nonlinear Science (2004)
! Software: Kdtree2
! Access Date: 4/14/20
! Repository: jmhodges/kdtree2
! Citation: https://arxiv.org/abs/physics/0408067
! License: Academic Free Licence (AFL)


module Kdtree_module

#include "petsc/finclude/petscsys.h"
 
  use petscsnes

  implicit none

  PetscInt :: bucket_size = 12   

  type, public :: kdtree_result
    PetscReal :: dis
    PetscInt :: idx 
  end type kdtree_result

  ! The priority queue consists of elements
  ! priority(1:heap_size), with associated payload(:).
  !
  ! There are heap_size active elements.
  ! Assumes the allocation is always sufficient.  Will NOT increase it
  ! to match.
  type, public :: pq
    PetscInt :: heap_size != 0
    type(kdtree_result), pointer :: elems(:)
   end type pq
   
  type, public :: interval
    PetscReal :: lower, upper 
  end type interval

  ! an internal tree node
  type, public :: tree_node
    ! the dimension to cut
    PetscInt :: cut_dim
    ! where to cut the dimension
    PetscReal :: cut_val
    ! improved cutoffs knowing the spread in child boxes.
    PetscReal :: cut_val_left, cut_val_right
    PetscInt :: l, u
    type(tree_node), pointer :: left, right
    type(interval), pointer :: box(:)! => null()
  end type tree_node

  ! Global information about the tree, one per tree
  type, public :: kdtree
    ! dimensionality and total # of points
    PetscInt :: dimen, n 
    ! pointer to the actual data array
    PetscReal, pointer :: the_data(:,:)
    ! permuted index into the data, so that indexes[l..u] of some
    ! bucket represent the indexes of the actual points in that
    ! bucket.
    PetscInt, pointer :: ind(:)
    PetscBool :: sort 
    PetscBool :: rearrange 
    ! if (rearrange .eqv. .true.) then rearranged_data has been
    ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
    ! permitting search to use more cache-friendly rearranged_data, at
    ! some initial computation and storage cost.
    PetscReal, pointer :: rearranged_data(:,:) 
    ! root pointer of the tree
    type (tree_node), pointer :: root
  end type kdtree

  type, public :: tree_search_record    
    ! One of these is created for each search.
    ! Many fields are copied from the tree structure, in order to
    ! speed up the search.
    PetscInt :: dimen, nn, nfound
    PetscReal :: ballsize
    PetscInt :: centeridx, correltime
    ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
    PetscInt :: nalloc
    PetscBool :: rearrange, overflow
    PetscReal, pointer :: qv(:)
    type(kdtree_result), pointer :: results(:)
    type(pq) :: pq
    PetscReal, pointer :: data(:,:)
    PetscInt, pointer :: ind(:)
  end type tree_search_record

  ! A GLOBAL VARIABLE for search
  type(tree_search_record),target :: sr   

contains

! ************************************************************************** !

function KdtreeCreate()

  implicit none

  type(kdtree), pointer :: KdtreeCreate
  
  type(kdtree), pointer :: aux

  allocate(aux)
  aux%dimen = 0
  aux%n = 0
  aux%sort = PETSC_FALSE
  aux%rearrange = PETSC_FALSE
  nullify(aux%the_data)
  nullify(aux%ind)
  nullify(aux%rearranged_data)
  nullify(aux%root)

  KdtreeCreate => aux
  
end function KdtreeCreate

! ************************************************************************** !

subroutine KdtreeConstruct(kdtree_aux,input_data,dim,sort,rearrange)  
  
  implicit none  

  type(kdtree), pointer :: kdtree_aux !mr
  PetscInt, optional :: dim
  PetscBool, optional :: sort, rearrange
  PetscReal, target :: input_data(:,:)
  PetscInt :: i

 ! allocate(mr)

  kdtree_aux%the_data => input_data
    
  if (present(dim)) then
    kdtree_aux%dimen = dim
  else
    kdtree_aux%dimen = size(input_data,1)
  endif
  kdtree_aux%n = size(input_data,2)

  call BuildTree(kdtree_aux)

  if (present(sort)) then
    kdtree_aux%sort = sort
  else
    kdtree_aux%sort = PETSC_FALSE
  endif

  if (present(rearrange)) then
    kdtree_aux%rearrange = rearrange
  else
    kdtree_aux%rearrange = PETSC_TRUE
  endif

  if (kdtree_aux%rearrange) then
    allocate(kdtree_aux%rearranged_data(kdtree_aux%dimen,kdtree_aux%n))
    do i=1,kdtree_aux%n
      kdtree_aux%rearranged_data(:,i) = kdtree_aux%the_data(:, &
      kdtree_aux%ind(i))
    enddo
  else
    nullify(kdtree_aux%rearranged_data)
  endif

end subroutine KdtreeConstruct

! ************************************************************************** !

subroutine BuildTree(tp)

  implicit none
  
  type(kdtree), pointer :: tp
  PetscInt :: j
  type(tree_node), pointer :: dummy => null()

 ! nullify(dummy%box)
  allocate (tp%ind(tp%n))
  forall (j=1:tp%n)
    tp%ind(j) = j
  end forall

  call BuildTreeForRange(tp,1,tp%n,dummy,tp%root)
  
end subroutine BuildTree

! ************************************************************************** !

recursive subroutine BuildTreeForRange(tp,l,u,parent,res)
  
  implicit none
  
  ! Function Return Cut_value 
  type(tree_node), pointer :: res

  ! Structure Arguments 
  type(kdtree), pointer :: tp
  type(tree_node),pointer :: parent

  ! Scalar Arguments

  PetscInt :: l, u

  ! Local Scalars
  PetscInt :: i, c, m, dimen
  PetscBool :: recompute

  PetscReal :: average

  ! first compute min and max
  dimen = tp%dimen
  allocate(res)
  allocate(res%box(dimen))

  ! First, compute an APPROXIMATE bounding box of all points &
  ! associated with this node.
  if (u < l) then
    ! no points in this box
    nullify(res)
    return
  endif

  if ((u - l) <= bucket_size) then
         
    ! always compute true bounding box for terminal nodes.        
    do i=1,dimen
      call SpreadInCoordinate(tp,i,l,u,res%box(i))
    enddo

    res%cut_dim = 0
    res%cut_val = 0.0
    res%l = l
    res%u = u
    res%left => null()
    res%right => null()
  else
         
    ! modify approximate bounding box.  This will be an
    ! overestimate of the true bounding box, as we are only recomputing
    ! the bounding box for the dimension that the parent split on.
    !
    ! Going to a true bounding box computation would significantly
    ! increase the time necessary to build the tree, and usually
    ! has only a very small difference.  This box is not used
    ! for searching but only for deciding which coordinate to split on.
         
    do i=1,dimen
      recompute= PETSC_TRUE !.true.
      if (associated(parent)) then
        if (i .ne. parent%cut_dim) then
          recompute= PETSC_FALSE !.false.
        endif
      endif
      if (recompute) then
        call SpreadInCoordinate(tp,i,l,u,res%box(i))
      else
        res%box(i) = parent%box(i)
      endif
    enddo
  
    c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         
    ! c is the identity of which coordinate has the greatest spread.
    if (PETSC_FALSE) then
      ! select exact median to have fully balanced tree.      
      m = (l+u)/2
      call SelectOnCoordinate(tp%the_data,tp%ind,c,m,l,u)
    else
            
      ! select point halfway between min and max, as per A. Moore,
      ! who says this helps in some degenerate cases, or
      ! actual arithmetic average.

      if (PETSC_TRUE) then
        ! actually compute average
        average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1)
      else
        average = (res%box(c)%upper + res%box(c)%lower)/2.0
      endif

      res%cut_val = average
      m = SelectOnCoordinateValue(tp%the_data,tp%ind,c,average,l,u)
    endif

    ! moves indexes around
    res%cut_dim = c
    res%l = l
    res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

    call BuildTreeForRange(tp,l,m,res,res%left)
    call BuildTreeForRange(tp,m+1,u,res,res%right)

    if (associated(res%right) .eqv. PETSC_FALSE) then
      res%box = res%left%box
      res%cut_val_left = res%left%box(c)%upper
      res%cut_val = res%cut_val_left
    elseif (associated(res%left) .eqv. PETSC_FALSE) then
      res%box = res%right%box
      res%cut_val_right = res%right%box(c)%lower
      res%cut_val = res%cut_val_right
    else
      res%cut_val_right = res%right%box(c)%lower
      res%cut_val_left = res%left%box(c)%upper
      res%cut_val = (res%cut_val_left + res%cut_val_right)/2


      ! now remake the true bounding box for self.
      ! Since we are taking unions (in effect) of a tree structure,
      ! this is much faster than doing an exhaustive
      ! search over all points
      res%box%upper = max(res%left%box%upper,res%right%box%upper)
      res%box%lower = min(res%left%box%lower,res%right%box%lower)
     endif
  endif
  
end subroutine BuildTreeForRange

! ************************************************************************** !

function SelectOnCoordinateValue(v,ind,c,alpha,li,ui)

  ! Move elts of ind around between l and u, so that all points
  ! <= than alpha (in c cooordinate) are first, and then
  ! all points > alpha are second.

  !
  ! Algorithm (matt kennel).
  !
  ! Consider the list as having three parts: on the left,
  ! the points known to be <= alpha.  On the right, the points
  ! known to be > alpha, and in the middle, the currently unknown
  ! points.   The algorithm is to scan the unknown points, starting
  ! from the left, and swapping them so that they are added to
  ! the left stack or the right stack, as appropriate.
  !
  ! The algorithm finishes when the unknown stack is empty.

  implicit none

  PetscInt :: c, li, ui
  PetscReal :: alpha

  PetscReal :: v(1:,1:)

  PetscInt :: ind(1:), tmp, lb, rb

  PetscInt :: SelectOnCoordinateValue

  lb = li; rb = ui

  do while (lb < rb)
    if ( v(c,ind(lb)) <= alpha ) then
      lb = lb + 1
    else
      tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
      rb = rb - 1
    endif
  enddo

  if (v(c,ind(lb)) <= alpha) then
    SelectOnCoordinateValue = lb
  else
    SelectOnCoordinateValue = lb - 1
  endif

end function SelectOnCoordinateValue

! ************************************************************************** !

subroutine SelectOnCoordinate(v,ind,c,k,li,ui)
  ! Move elts of ind around between l and u, so that the kth
  ! element
  ! is >= those below, <= those above, in the coordinate c.

  implicit none

  PetscInt :: c, k, li, ui
  PetscInt :: i, l, m, s, t, u, ind(:)
  PetscReal ::  v(:,:)

  l = li
  u = ui
  do while (l < u)
    t = ind(l)
    m = l
    do i = l + 1, u
      if (v(c,ind(i)) < v(c,t)) then
        m = m + 1
        s = ind(m)
        ind(m) = ind(i)
        ind(i) = s
      endif
    enddo
    s = ind(l)
    ind(l) = ind(m)
    ind(m) = s
    if (m <= k) l = m + 1
    if (m >= k) u = m - 1
  enddo
 
end subroutine SelectOnCoordinate

! ************************************************************************** !

subroutine SpreadInCoordinate(tp,c,l,u,interv)
  ! the spread in coordinate 'c', between l and u.
  !
  ! Return lower bound in 'smin', and upper in 'smax',

  implicit none
  
  type(kdtree), pointer :: tp
  type(interval), intent(out) :: interv

  PetscInt :: c, l, u

  PetscReal :: last, lmax, lmin, t, smin, smax
  PetscInt :: i, ulocal

  PetscReal, pointer :: v(:,:)
  PetscInt, pointer :: ind(:)

  v => tp%the_data(1:,1:)
  ind => tp%ind(1:)
  smin = v(c,ind(l))
  smax = smin

  ulocal = u

  do i = l + 2, ulocal, 2
    lmin = v(c,ind(i-1))
    lmax = v(c,ind(i))
    if (lmin > lmax) then
      t = lmin
      lmin = lmax
      lmax = t
    endif
    if (smin > lmin) smin = lmin
    if (smax < lmax) smax = lmax
  enddo
 
  if (i == ulocal + 1) then
    last = v(c,ind(ulocal))
    if (smin > last) smin = last
    if (smax < last) smax = last
  endif

   interv%lower = smin
   interv%upper = smax

end subroutine SpreadInCoordinate

! ************************************************************************** !

subroutine KdtreeDestroy(tp)

  ! Deallocates all memory for the tree, except input data matrix

  implicit none

  type(kdtree), pointer :: tp

  call KdtreeDestroyNode(tp%root)

  deallocate (tp%ind)
  nullify (tp%ind)

  if (tp%rearrange) then
    deallocate(tp%rearranged_data)
    nullify(tp%rearranged_data)
  endif

  deallocate(tp)
 end subroutine KdtreeDestroy

recursive subroutine KdtreeDestroyNode(np)

  implicit none

  type (tree_node), pointer :: np

  if (associated(np%left)) then
    call KdtreeDestroyNode(np%left)
    nullify (np%left)
  endif
  if (associated(np%right)) then
    call KdtreeDestroyNode(np%right)
    nullify (np%right)
  endif
  if (associated(np%box)) deallocate(np%box)
  deallocate(np)
  
end subroutine KdtreeDestroyNode

! ************************************************************************** !

subroutine KdtreeNNearest(tp,qv,nn,results)
  
  ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
  ! returning their indexes and distances in 'indexes' and 'distances'
  ! arrays already allocated passed to this subroutine.

  implicit none
  
  type(kdtree), pointer :: tp
  PetscReal, target :: qv(:)
  PetscInt :: nn
  type(kdtree_result), target :: results(:)
   
  sr%ballsize = huge(1.0)
  sr%qv => qv
  sr%nn = nn
  sr%nfound = 0
  sr%centeridx = -1
  sr%correltime = 0
  sr%overflow = PETSC_FALSE

  sr%results => results

  sr%nalloc = nn 

  sr%ind => tp%ind
  sr%rearrange = tp%rearrange
  if (tp%rearrange) then
    sr%Data => tp%rearranged_data
  else
    sr%Data => tp%the_data
  endif
  sr%dimen = tp%dimen

  sr%pq = KdtreeCreatePQ(results)

  call KdtreeSearch(tp%root)   

  if (tp%sort) then
    call KdtreeSortResults(nn, results)
  endif

  
end subroutine KdtreeNNearest

! ************************************************************************** !

recursive subroutine KdtreeSearch(node)
  !
  ! This is the innermost core routine of the kd-tree search.  Along
  ! with "process_terminal_node", it is the performance bottleneck.
  !
  ! This version uses a logically complete secondary search of
  ! "box in bounds", whether the sear
  !

  implicit none
  
  type(tree_node), pointer :: node
  type(tree_node), pointer :: ncloser, nfarther

  PetscInt :: cut_dim, i
  PetscReal :: qval, dis, ballsize
  PetscReal, pointer :: qv(:)

  type(interval), pointer :: box(:)

  if ((associated(node%left) .and. associated(node%right)) .eqv. PETSC_FALSE) then
    ! we are on a terminal node
    if (sr%nn == 0) then
      call ProcessTerminalNodeFixedBall(node)
    else
      call ProcessTerminalNode(node)
    endif
  else
    ! we are not on a terminal node
    qv => sr%qv(1:)
    cut_dim = node%cut_dim
    qval = qv(cut_dim)

    if (qval < node%cut_val) then
      ncloser => node%left
      nfarther => node%right
      dis = (node%cut_val_right - qval)**2
    else
      ncloser => node%right
      nfarther => node%left
      dis = (node%cut_val_left - qval)**2
    endif

    if (associated(ncloser)) call KdtreeSearch(ncloser)

    ! we may need to search the second node.
    if (associated(nfarther)) then
      ballsize = sr%ballsize
      if (dis <= ballsize) then
             
        ! we do this separately as going on the first cut dimen is often
        ! a good idea.
        ! note that if extra**2 < sr%ballsize, then the next
        ! check will also be false.             
        box => node%box(1:)
        do i=1,sr%dimen
          if (i /= cut_dim) then
            dis = dis + KdtreeDisFromBnd(qv(i),box(i)%lower,box(i)%upper)
            if (dis > ballsize) then
              return
            endif
          endif
        enddo

        call KdtreeSearch(nfarther)
      endif
    endif
  endif
 
end subroutine KdtreeSearch

! ************************************************************************** !

function KdtreeDisFromBnd(x,amin,amax) result (res)

  implicit none
  
  PetscReal :: x, amin, amax
  PetscReal :: res
  
  if (x > amax) then
    res = (x - amax)**2;
    return
  else
    if (x < amin) then
      res = (amin - x)**2;
      return
    else
      res = 0.0
      return
    endif
  endif
  return
  
end function KdtreeDisFromBnd

! ************************************************************************** !

subroutine ProcessTerminalNode(node)
    
  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure.

  implicit none
  
  type (tree_node), pointer :: node
  !
  PetscReal, pointer :: qv(:), data(:,:)
  PetscInt, pointer :: ind(:)

  PetscInt :: dimen, i , indexofi, k, centeridx, correltime
  PetscReal :: ballsize, sd, newpri
  PetscBool :: rearrange
  type(pq), pointer :: pqp
    
  qv => sr%qv(1:)
  pqp => sr%pq
  dimen = sr%dimen
  ballsize = sr%ballsize
  rearrange = sr%rearrange
  ind => sr%ind(1:)
  data => sr%Data(1:,1:)
  centeridx = sr%centeridx
  correltime = sr%correltime

  do i = node%l, node%u
    if (rearrange) then
      sd = 0.0
      do k = 1,dimen
        sd = sd + (data(k,i) - qv(k))**2
        if (sd>ballsize) exit
      enddo
      if (sd>ballsize) cycle
      indexofi = ind(i)  
    else
      indexofi = ind(i)
      sd = 0.0
      do k = 1,dimen
        sd = sd + (data(k,indexofi) - qv(k))**2
        if (sd > ballsize) exit 
      enddo
      if (sd > ballsize) cycle
    endif

    if (centeridx > 0) then 
      if (abs(indexofi - centeridx) < correltime) cycle
    endif

    if (sr%nfound < sr%nn) then
          
      sr%nfound = sr%nfound + 1
      newpri = KdtreePQInsert(pqp,sd,indexofi)
      if (sr%nfound == sr%nn) ballsize = newpri
    else
          
      ballsize = KdtreePQReplaceMax(pqp,sd,indexofi)
    endif
  enddo 
    
  sr%ballsize = ballsize

end subroutine ProcessTerminalNode

! ************************************************************************** !

subroutine ProcessTerminalNodeFixedBall(node)
    
  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure, i.e.
  ! save all within a fixed ball.

  implicit none
  
  type(tree_node), pointer :: node

  PetscReal, pointer :: qv(:), data(:,:)
  PetscInt, pointer :: ind(:)

  PetscInt :: nfound, dimen, i, indexofi, k
  PetscInt :: centeridx, correltime, nn
  PetscReal :: ballsize, sd
  PetscBool :: rearrange
    
  ! copy values from sr to local variables
    
  qv => sr%qv(1:)
  dimen = sr%dimen
  ballsize = sr%ballsize
  rearrange = sr%rearrange
  ind => sr%ind(1:)
  data => sr%Data(1:,1:)
  centeridx = sr%centeridx
  correltime = sr%correltime
  nn = sr%nn 
  nfound = sr%nfound

  ! search through terminal bucket.
  do i = node%l, node%u

    if (rearrange) then
      sd = 0.0
      do k = 1,dimen
        sd = sd + (data(k,i) - qv(k))**2
        if (sd > ballsize) exit 
      enddo
      if (sd > ballsize) cycle
      indexofi = ind(i) 
    else
      indexofi = ind(i)
      sd = 0.0
      do k = 1,dimen
        sd = sd + (data(k,indexofi) - qv(k))**2
        if (sd>ballsize) exit 
      enddo
      if (sd > ballsize) cycle
    endif

    if (centeridx > 0) then 
      if (abs(indexofi - centeridx) < correltime) cycle 
    endif

    nfound = nfound+1
    if (nfound > sr%nalloc) then
      sr%overflow = PETSC_TRUE
    else
      sr%results(nfound)%dis = sd
      sr%results(nfound)%idx = indexofi
    endif
  enddo
     
  sr%nfound = nfound
  
end subroutine ProcessTerminalNodeFixedBall

! ************************************************************************** !

subroutine KdtreeSortResults(nfound,results)

  !  Use after search to sort results(1:nfound) in order of increasing
  !  distance.

  implicit none

  PetscInt :: nfound
  type(kdtree_result), target :: results(:)

  if (nfound > 1) call KdtreeHeapsortStruct(results,nfound)

end subroutine KdtreeSortResults

! ************************************************************************** !

subroutine KdtreeHeapsortStruct(a,n)
  
  ! Sort a(1:n) in ascending order

  implicit none
               
  PetscInt :: n
  type(kdtree_result),intent(inout) :: a(:)
  type(kdtree_result) :: value

  PetscInt :: i, j, ileft, iright

  ileft = n/2 + 1
  iright = n

  if(n == 1) return

  do
    if(ileft > 1)then
      ileft = ileft - 1
      value = a(ileft)
    else
      value = a(iright)
      a(iright) = a(1)
      iright = iright - 1
      if (iright == 1) then
        a(1) = value
        return
      endif
    endif
    i = ileft
    j = 2 * ileft
    do while (j <= iright)
      if(j < iright) then
        if(a(j)%dis < a(j+1)%dis) j = j + 1
      endif
      if(value%dis < a(j)%dis) then
        a(i) = a(j);
        i = j
        j = j + j
      else
        j = iright + 1
      endif
    enddo
    a(i) = value
  enddo
 
end subroutine KdtreeHeapsortStruct

!************************************************************************** !

function KdtreeCreatePQ(results_in)
  !
  ! Create a priority queue from ALREADY allocated
  ! array pointers for storage.  NOTE! It will NOT
  ! add any alements to the heap, i.e. any existing
  ! data in the input arrays will NOT be used and may
  ! be overwritten.
  !

  implicit none
  
  type(kdtree_result), target:: results_in(:)
  type(pq) :: KdtreeCreatePQ
  PetscInt :: nalloc

  nalloc = size(results_in,1)
  KdtreeCreatePQ%elems => results_in
  KdtreeCreatePQ%heap_size = 0
  
end function KdtreeCreatePQ

!************************************************************************** !

subroutine KdtreeHeapify(a,i_in)
  !
  ! take a heap rooted at 'i' and force it to be in the
  ! heap canonical form.   This is performance critical
  ! and has been tweaked a little to reflect this.
  !
  
  type(pq), pointer :: a
  PetscInt :: i_in
  PetscInt :: i, l, r, largest
  PetscReal :: pri_i, pri_l, pri_r, pri_largest
  type(kdtree_result) :: temp

  i = i_in

  do while (largest /= i)
    l = 2 * i
    r = l + 1

    if (l > a%heap_size) then
      exit
    else
      pri_i = a%elems(i)%dis
      pri_l = a%elems(l)%dis
      if (pri_l > pri_i) then
        largest = l
        pri_largest = pri_l
      else
        largest = i
        pri_largest = pri_i
      endif

      if (r < a%heap_size) then
        pri_r = a%elems(r)%dis
        if (pri_r > pri_largest) then
          largest = r
        endif
      endif
    endif

    if (largest /= i) then

      temp = a%elems(i)
      a%elems(i) = a%elems(largest)
      a%elems(largest) = temp

      i = largest   
    endif
  enddo 
  
end subroutine KdtreeHeapify

!************************************************************************** !

subroutine KdtreePQExtractMax(a,e)
  !
  ! return the priority and payload of maximum priority
  ! element, and remove it from the queue.
  ! (equivalent to 'pop()' on a stack)
  !

  implicit none
  
  type(pq),pointer :: a
  type(kdtree_result), intent(out) :: e

  e = a%elems(1)

  a%elems(1) = a%elems(a%heap_size)
  a%heap_size = a%heap_size-1
  call KdtreeHeapify(a,1)

  return

end subroutine  KdtreePQExtractMax

!************************************************************************** !

function KdtreePQInsert(a,dis,idx)
  !
  ! Insert a new element and return the new maximum priority,
  ! which may or may not be the same as the old maximum priority.
  !
  implicit none
  
  type(pq), pointer :: a
  PetscReal :: dis
  PetscInt :: idx
  PetscInt :: i, isparent
  PetscReal:: parentdis
  PetscReal :: KdtreePQInsert

  a%heap_size = a%heap_size + 1
  i = a%heap_size

  do while (i > 1)
    isparent = int(i/2)
    parentdis = a%elems(isparent)%dis
    if (dis > parentdis) then
      a%elems(i)%dis = parentdis
      a%elems(i)%idx = a%elems(isparent)%idx
      i = isparent
    else
      exit
    endif
  enddo

  a%elems(i)%dis = dis
  a%elems(i)%idx = idx

  KdtreePQInsert = a%elems(1)%dis
  return

end function KdtreePQInsert

!************************************************************************** !

function KdtreePQReplaceMax(a,dis,idx)
  !
  ! Replace the extant maximum priority element
  ! in the PQ with (dis,idx).  Return
  ! the new maximum priority, which may be larger
  ! or smaller than the old one.
  !

  implicit none
  
  type(pq), pointer :: a
  PetscReal :: dis
  PetscInt :: idx
  PetscInt :: parent, child, N
  PetscReal :: prichild, prichildp1

  type(kdtree_result) :: etmp
  PetscReal :: KdtreePQReplaceMax


  N = a%heap_size
  if (N >= 1) then
    parent = 1
    child = 2

    do while (child <= N)
      prichild = a%elems(child)%dis

      if (child < N) then
        prichildp1 = a%elems(child + 1)%dis
        if (prichild < prichildp1) then
          child = child + 1
          prichild = prichildp1
        endif
      endif

      if (dis >= prichild) then
        exit 
      else
        a%elems(parent) = a%elems(child)
        parent = child
        child = 2 * parent
      endif
    enddo 
    a%elems(parent)%dis = dis
    a%elems(parent)%idx = idx
    KdtreePQReplaceMax = a%elems(1)%dis
  else
    a%elems(1)%dis = dis
    a%elems(1)%idx = idx
    KdtreePQReplaceMax = dis
  endif

  return
  
end function KdtreePQReplaceMax

end module Kdtree_module
