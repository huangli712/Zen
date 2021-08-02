!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : linkedlist
!!! source  : m_linkedlist.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           07/29/2021 by li huang (last modified)
!!! purpose : this purpose of this module is to implement a typical and
!!!           useful data structure --- linked list. it is a generic
!!!           linked list, capable of storing arbitrary data.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module linkedlist
     implicit none

!!========================================================================
!!>>> declare global data types                                        <<<
!!========================================================================

!! module variables
     ! a public variable used as a mold for transfer() subroutine.
     integer, allocatable :: list_d(:)

!! module structs
     ! node for the linked list. it contains a pointer pointing to next
     ! node, and the integer pointer array is used to store data.
     type list_t
         private
         !
         integer, pointer :: data(:) => null()
         type (list_t), pointer :: next => null()
     end type list_t

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: list_d
     public :: list_t

     public :: list_init
     public :: list_free

     public :: list_insert
     public :: list_put
     public :: list_get
     public :: list_next
     public :: list_count

  contains ! encapsulated functionality

!!
!! @sub list_init
!!
!! initialize a head node [self] for a list and optionally store the
!! provided data.
!!
  subroutine list_init(self, data)
     implicit none

!! external arguments
     ! pointer to new linked list
     type (list_t), pointer :: self

     ! the data for the first element
     integer, intent(in), optional :: data(:)

!! [body

     ! allocate memory for linked list
     allocate(self)
     nullify(self%next)

     ! check whether we should make an empty node
     if ( present(data) ) then
         allocate( self%data( size(data) ) )
         self%data = data
     else
         nullify(self%data)
     endif ! back if ( present(data) ) block

!! body]

     return
  end subroutine list_init

!!
!! @sub list_free
!!
!! free the entire list and all data, beginning at node [self].
!!
  subroutine list_free(self)
     implicit none

!! external arguments
     ! pointer to the list to be destroyed
     type (list_t), pointer :: self

!! local variables
     ! pointer to the current node
     type (list_t), pointer :: curr

     ! pointer to the next node
     type (list_t), pointer :: next

!! [body

     ! go through the whole linked list
     curr => self
     !
     do while ( associated(curr) )
         ! get next node
         next => curr%next
         !
         ! release memory for the internal data
         if ( associated(curr%data) ) then
             deallocate(curr%data)
             nullify(curr%data)
         endif ! back if ( associated(curr%data) ) block
         !
         ! release memory for the node itself
         deallocate(curr)
         nullify(curr)
         !
         ! point to next node
         curr => next
     enddo ! over do while loop

!! body]

     return
  end subroutine list_free

!!
!! @sub list_insert
!!
!! insert a node containing data (optional) after node [self].
!!
  subroutine list_insert(self, data)
     implicit none

!! external arguments
     ! element in the linked list after which
     ! the new element should be inserted.
     type (list_t), pointer :: self

     ! the data for the new element
     integer, intent(in), optional :: data(:)

!! local variables
     ! pointer to new node
     type (list_t), pointer :: next

!! [body

     ! allocate memory for new node
     allocate(next)

     ! whether we should build an empty node
     if ( present(data) ) then
         allocate( next%data( size(data) ) )
         next%data = data
     else
         nullify(next%data)
     endif ! back if ( present(data) ) block

     ! update the linked list
     next%next => self%next
     self%next => next

!! body]

     return
  end subroutine list_insert

!!
!! @sub list_put
!!
!! store the encoded data in list node [self].
!!
  subroutine list_put(self, data)
     implicit none

!! external arguments
     ! element in the linked list
     type (list_t), pointer :: self

     ! the data to be stored
     integer, intent(in) :: data(:)

!! [body

     ! release old memory at first
     if ( associated(self%data) ) then
         deallocate(self%data)
         nullify(self%data)
     endif ! back if ( associated(self%data) ) block

     ! allocate new memory
     allocate( self%data( size(data) ) )

     ! save the data
     self%data = data

!! body]

     return
  end subroutine list_put

!!
!! @fun list_get
!!
!! return the data stored in the node [self].
!!
  function list_get(self) result(data)
     implicit none

!! external arguments
     ! node in the linked list
     type (list_t), pointer :: self

     ! function value, the node's data
     integer, pointer :: data(:)

!! [body

     data => self%data

!! body]

     return
  end function list_get

!!
!! @fun list_next
!!
!! return the next node after node [self].
!!
  function list_next(self) result(next)
     implicit none

!! external arguments
     ! pointer to the list
     type (list_t), pointer :: self

     ! function value, pointer to the next node
     type (list_t), pointer :: next

!! [body

     next => self%next

!! body]

     return
  end function list_next

!!
!! @fun list_count
!!
!! count the number of nodes in the list [self]. in fact, this function
!! can be used to return the number of nodes after a given node.
!!
  function list_count(self) result(counter)
     implicit none

!! external arguments
     ! pointer to the list
     type (list_t), pointer :: self

     ! function value
     integer :: counter

!! local variables
     ! pointer to current node
     type (list_t), pointer :: curr

!! [body

     if ( associated(self) ) then
         counter = 1
         curr => self
         do while ( associated(curr%next) )
             counter = counter + 1
             curr => curr%next
         enddo ! over do while loop
     else
         counter = 0
     endif ! back if ( associated(self) ) block

!! body]

     return
  end function list_count

  end module linkedlist
