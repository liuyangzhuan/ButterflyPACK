! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so.

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).

#include "ButterflyPACK_config.fi"
module BPACK_linkedlist

!**** declarations for list type with unlimited polymorphism
type :: nod
	type(nod), pointer :: next => null()
	class(*), allocatable :: item
	contains
#ifdef HAVE_FINAL
	final :: nod_finalizer
#endif
end type nod

type :: list
	integer :: num_nods = 0
	integer :: idx=0
	type(nod), pointer :: head => null()
	type(nod), pointer :: tail => null()
	procedure(nod_score),nopass,pointer :: FuncScore=>null()
	contains
#ifdef HAVE_FINAL
	final :: list_finalizer
#endif
	procedure :: len => list_length
end type list

! interfaces:
interface list
	module procedure :: list_constructor
end interface list

interface nod
	module procedure :: nod_constructor
end interface nod

interface append
	module procedure :: list_append_item
end interface append

! interface get_item
	! module procedure :: list_get_item_character
	! module procedure :: nod_get_item_character
! end interface get_item

interface get_nod
	module procedure :: list_get_nod
end interface get_nod

interface assignment(=)
	module procedure :: nod_assign_nod_to_nod
end interface assignment(=)

abstract interface
	function nod_score(this) result(score)
		import :: nod
		implicit none
		type(nod)::this
		real(kind=8)::score
	end function nod_score
end interface
contains


!===============================================================================
!  list_append_item:
!
!    Finalizes the components of the given list.
!
subroutine list_append_item( this, item )
type(list), intent(inout) :: this
class(*), intent(in) :: item
if(associated(this%tail)) then
  allocate(this%tail%next)
  this%tail%next=nod(item)
  this%tail => this%tail%next
else
  allocate(this%head)
  this%head=nod(item)
  this%tail => this%head
end if
this%num_nods = this%num_nods + 1
end subroutine list_append_item
!===============================================================================


!===============================================================================
!  list_remove_headitem:
!
!    Remove the head component of the given list.
!
subroutine list_remove_headitem( this )
type(list), intent(inout) :: this
type(nod),pointer:: cur
if(associated(this%head))then
	cur=>this%head%next
	call nod_finalizer(this%head)
	deallocate(this%head)
	this%head=>cur
	this%num_nods = this%num_nods - 1
endif
if(this%num_nods==0)this%tail=>null()
end subroutine list_remove_headitem
!===============================================================================


!===============================================================================
!  list_constructor:
!
!    Returns an uninitialized list.
!
function list_constructor( ) result( val )
type(list) :: val
val%num_nods=0
end function list_constructor
!===============================================================================

!===============================================================================
!  list_finalizer:
!
!    Finalizes the components of the given list.
!
subroutine list_finalizer( this )
type(list), intent(inout) :: this

do while(this%num_nods>0)
	call list_remove_headitem(this)
enddo
end subroutine list_finalizer
!===============================================================================

! !===============================================================================
! !  list_get_item_character:
! !
! !    Gets the i-th nod from this list and returns item value if it is a
! !    character.
! !
! !      STAT   ERRMSG
! !        -1   item found but not of type character
! !        -2   nod found but item not allocated
! !        -3   nod not found (inod exceeds list bounds)
! !
! subroutine list_get_item_character( this, inod, chVal, stat, errmsg )
! type(list), intent(in) :: this
! integer, intent(in) :: inod
! character(*), intent(out) :: chVal
! integer, intent(out), optional :: stat
! character(*), intent(out), optional :: errmsg
! ! local variables:
! type(nod) :: nVal
! integer    :: istat

! call get_nod(this, inod, nVal, stat=istat)
! if (istat == 0) then
  ! call get_item(nVal, chVal, stat=istat)
! else
  ! istat = -3
! end if

! if (present(stat)) stat = istat

! if (present(errmsg)) then
  ! select case (istat)
  ! case (-1)
	! errmsg = 'item found but not of type character'
  ! case (-2)
	! errmsg = 'nod found but item not allocated'
  ! case (-3)
	! errmsg = 'nod not found (inod exceeds list bounds)'
  ! case default
	! errmsg = ''
  ! end select
! end if
! end subroutine list_get_item_character
! !===============================================================================

!===============================================================================
!  list_get_nod:
!
!    Gets the i-th nod from this list.
!
!      STAT   ERRMSG
!        -1   nod not found (inod exceeds list bounds)
!
subroutine list_get_nod( this, inod, nVal, stat, errmsg )
type(list), intent(in) :: this
integer, intent(in) :: inod
type(nod), intent(out) :: nVal
integer, intent(out), optional :: stat
character(*), intent(out), optional :: errmsg
! local variables:
integer :: i, istat, list_len
type(nod), pointer :: current_nod


if (inod < 1) then
  istat = -1
else if (inod > this%len()) then
  istat = -1
else
  istat = 0

  current_nod => this%head
  do i = 2, inod
	current_nod => current_nod%next
  end do
  nVal = current_nod
  current_nod => null()
end if

if (present(stat)) stat = istat

if (present(errmsg)) then
  select case (istat)
  case (-1)
	errmsg = 'nod not found (inod exceeds list bounds)'
  case default
	errmsg = ''
  end select
end if
end subroutine list_get_nod
!===============================================================================



!===============================================================================
!  list_print_nod_score:
!
!    Gets the i-th nod from this list.
!
!      STAT   ERRMSG
!        -1   nod not found (inod exceeds list bounds)
!
subroutine list_print_scores( this,FuncScore )
type(list), intent(in) :: this
! local variables:
integer :: i, istat, list_len
type(nod), pointer :: current_nod
procedure(nod_score)::FuncScore
  current_nod => this%head
  do i = 1, this%num_nods
	write(*,*)i,FuncScore(current_nod)
	current_nod => current_nod%next
  end do

end subroutine list_print_scores
!===============================================================================



!===============================================================================
!  list_length:
!
!    Returns the number of nods in the given list
!
function list_length( self ) result( val )
class(list), intent(in) :: self
integer :: val

val = self%num_nods
end function list_length
!===============================================================================

!===============================================================================
!  nod_assign_nod_to_nod:
!
!    Returns .TRUE. if the given nod has an item of type complex with value
!    equal to the given complex value.
!
subroutine nod_assign_nod_to_nod( LHS, RHS )
type(nod), intent(inout) :: LHS
type(nod), intent(in) :: RHS
type(nod),pointer:: cur
class(*),pointer::ptrl,ptrr,ptr
integer ii

if (allocated(LHS%item))then
	select TYPE(ptrl=>LHS%item)
		type is (list)
			call list_finalizer(ptrl)
		class default ! intrinsic types and derived types with scalar,allocatable, or pointers
	end select
	deallocate(LHS%item)
endif

if(allocated(RHS%item))then
	allocate(LHS%item, source=RHS%item)
	select TYPE(ptrr=>RHS%item)
		type is (list)
			select TYPE(ptrl=>LHS%item)
			type is (list)
				ptrl%num_nods=0
				ptrl%head=>null()
				ptrl%tail=>null()
				ptrl%idx=ptrr%idx
				cur=>ptrr%head
				do ii=1,ptrr%num_nods
					call list_append_item(ptrl, cur%item)
					cur=>cur%next
				enddo
			class default ! intrinsic types and derived types with scalar,allocatable, or pointers
				stop
			end select
		class default ! intrinsic types and derived types with scalar,allocatable, or pointers
			! LHS%item=RHS%item
	end select
endif

end subroutine nod_assign_nod_to_nod
!===============================================================================

!===============================================================================
!  nod_constructor:
!
!    Returns a nod constructed from the given item.
!
function nod_constructor( item ) result( val )
class(*), intent(in), optional :: item
type(nod) :: val
type(nod),pointer:: cur
integer ii
class(*),pointer:: ptr
if (present(item))then
	allocate(val%item, source=item)
	select TYPE(item)
		type is (list)
			select TYPE(ptr=>val%item)
				type is (list)
					ptr%num_nods=0
					ptr%head=>null()
					ptr%tail=>null()
					ptr%idx=item%idx
					cur=>item%head
					do ii=1,item%num_nods
						call list_append_item(ptr, cur%item)
						cur=>cur%next
					enddo
				class default ! intrinsic types and derived types with scalar,allocatable, or pointers
					stop
			end select
		class default ! intrinsic types and derived types with scalar,allocatable, or pointers
			! val%item=item
	end select
endif
end function nod_constructor
!===============================================================================

!===============================================================================
!  nod_finalizer:
!
!    Finalizes the components of the given nod.
!
subroutine nod_finalizer( this )
type(nod), intent(inout) :: this
class(*),pointer::ptr
if (associated(this%next)) nullify(this%next)
if (allocated(this%item))then
	select TYPE(ptr=>this%item)
		type is (list)
			call list_finalizer(ptr)
		class default ! intrinsic types and derived types with scalar,allocatable, or pointers
	end select
	deallocate(this%item)
endif
end subroutine nod_finalizer
!===============================================================================

! !===============================================================================
! !  nod_get_item_character:
! !
! !    Returns .TRUE. if the given nod has an item of type character with value
! !    equal to the given character value.
! !
! !      STAT   ERRMSG
! !        -1   nod item not of type character
! !        -2   nod item not allocated
! !
! subroutine nod_get_item_character( this, sVal, stat, errmsg )
! type(nod), intent(in) :: this
! character(*), intent(out) :: sVal
! integer, intent(out), optional :: stat
! character(*), intent(out), optional :: errmsg
! ! local variables:
! integer :: istat

! if (allocated(this%item)) then
  ! select type (item => this%item)
  ! type is (character(*))
	! sVal = item
	! istat = 0
  ! class default
	! istat = -1
  ! end select
! else
  ! istat = -2
! end if

! if (present(stat)) stat = istat

! if (present(errmsg)) then
  ! select case (istat)
  ! case (-1)
	! errmsg = 'nod item is not of type character'
  ! case (-2)
	! errmsg = 'nod item is not allocated'
  ! case default
	! errmsg = ''
  ! end select
! end if
! end subroutine nod_get_item_character
! !===============================================================================



function nod_score_integer(this) result(score)
	implicit none
	type(nod)::this
	real(kind=8)::score
	class(*),pointer::ptr

	select TYPE(ptr=>this%item)
		type is (integer)
			score=dble(ptr)
		class default
			write(*,*)'unexpected item type in nod_score_integer'
			stop
	end select
end function nod_score_integer

function nod_score_dble(this) result(score)
	implicit none
	type(nod)::this
	real(kind=8)::score
	class(*),pointer::ptr

	select TYPE(ptr=>this%item)
		type is (real(kind=8))
			score=dble(ptr)
		class default
			write(*,*)'unexpected item type in nod_score_dble'
			stop
	end select
end function nod_score_dble





!===============================================================================
!  MergeSort:
!
!    mergesort on a linked list
!
recursive subroutine MergeSort(headRef,FuncScore)
	implicit none
	type(nod),pointer::headRef,head,a,b,tmp
	procedure(nod_score)::FuncScore
	integer ii

	head => headRef
	if(associated(head))then
	if(associated(head%next))then
		call FrontBackSplit(head,a,b)
		call MergeSort(a,FuncScore)
		call MergeSort(b,FuncScore)
		call SortedMerge(a,b,headRef,FuncScore)
	endif
	endif
end subroutine MergeSort

recursive subroutine SortedMerge(a,b,result,FuncScore)
	implicit none
	type(nod),pointer::result,a,b
	procedure(nod_score)::FuncScore
	result=>null()
	if(.not. associated(a))then
		result=>b
	else if(.not. associated(b))then
		result=>a
	else
		if(FuncScore(a)<=FuncScore(b))then
			result=>a
			a=>a%next
			call SortedMerge(a,b,result%next,FuncScore)
		else
			result=>b
			b=>b%next
			call SortedMerge(a,b,result%next,FuncScore)
		endif
	endif
end subroutine SortedMerge

subroutine FrontBackSplit(source,frontRef,backRef)
	implicit none
	type(nod),pointer::source,frontRef,backRef,fast,slow
	slow=>source
	fast=>source%next
	do while(associated(fast))
		fast =>fast%next
		if(associated(fast))then
			slow =>slow%next
			fast =>fast%next
		endif
	enddo
	frontRef=>source
	backRef=>slow%next
	slow%next=>null()
end subroutine FrontBackSplit



end module BPACK_linkedlist
