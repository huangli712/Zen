!!! project : flink @ sakura
!!! program : sparse
!!! source  : m_sparse.f90
!!! type    : module
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 02/01/2010 by li huang (created)
!!!           01/16/2025 by li huang (last modified)
!!! purpose : the purpose of this module is to implement important sparse
!!!           matrix/vector operations, including matrix multiplication,
!!!           format conversion, etc. the internal format of sparse matrix
!!!           used in this module is CSR (compressed sparse row) format.
!!! status  : unstable
!!! comment : only support real(dp) and complex(dp) data types
!!!-----------------------------------------------------------------------

  module sparse
     implicit none

!!========================================================================
!!>>> declare local parameters                                         <<<
!!========================================================================

!! module parameters
     ! dp: double precision for real and complex number
     integer, private, parameter :: dp    = kind(1.0d0)

     ! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare global structures                                        <<<
!!========================================================================

     ! sparse_t: abstract struct for sparse matrix
     type, private, abstract :: sparse_t

         ! nrows: number of rows
         integer :: nrows = 0

         ! ncols: number of columns
         integer :: ncols = 0

         ! nnz: number of non-zero values
         integer :: nnz   = 0

     end type sparse_t

     ! csr_t: sparse matrix in compressed sparse row format
     type, private, extends(sparse_t) :: csr_t

         ! rowptr: matrix row pointer
         integer, allocatable :: rowptr(:)

         ! colptr: matrix column pointer
         integer, allocatable :: colptr(:)

     end type csr_t

     ! csr_d: compressed sparse row format, real(dp) version
     type, public, extends(csr_t) :: csr_d
         !
         real(dp), allocatable :: V(:)
         !
     end type csr_d

     ! csr_z: compressed sparse row format, complex(dp) version
     type, public, extends(csr_t) :: csr_z
         !
         complex(dp), allocatable :: V(:)
         !
     end type csr_z

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! CSR -> allocation
     private :: csr_alloc_d ! real(dp) version
     private :: csr_alloc_z ! complex(dp) version
     !
     private :: csr_alloc_d_t
     private :: csr_alloc_z_t

     ! CSR -> deallocation
     private :: csr_free_d  ! real(dp) version
     private :: csr_free_z  ! complex(dp) version
     !
     private :: csr_free_d_t
     private :: csr_free_z_t

     ! CSR -> print
     private :: csr_print_d ! real(dp) version
     private :: csr_print_z ! complex(dp) version
     !
     private :: csr_print_d_t
     private :: csr_print_z_t

     ! CSR -> DNS (conversion)
     private :: csr_dns_d   ! real(dp) version
     private :: csr_dns_z   ! complex(dp) version
     !
     private :: csr_dns_d_t
     private :: csr_dns_z_t

     ! DNS -> CSR (conversion)
     private :: dns_csr_d   ! real(dp) version
     private :: dns_csr_z   ! complex(dp) version
     !
     private :: dns_csr_d_t
     private :: dns_csr_z_t

     ! CSR -> CSR (copy)
     private :: csr_csr_d   ! real(dp) version
     private :: csr_csr_z   ! complex(dp) version
     !
     private :: csr_csr_d_t
     private :: csr_csr_z_t

     ! CSR -> getter
     private :: get_csr_d   ! real(dp) version
     private :: get_csr_z   ! complex(dp) version
     !
     private :: get_csr_d_t
     private :: get_csr_z_t

     ! CSR X VEC
     private :: csr_mv_d    ! real(dp) version
     private :: csr_mv_z    ! complex(dp) version
     !
     private :: csr_mv_d_t
     private :: csr_mv_z_t

     ! CSR X CSR
     private :: csr_mm_d    ! real(dp) version
     private :: csr_mm_z    ! complex(dp) version
     !
     private :: csr_mm_d_t
     private :: csr_mm_z_t

     ! CSR X DIA
     private :: csr_md_d    ! real(dp) version
     private :: csr_md_z    ! complex(dp) version
     !
     private :: csr_md_d_t
     private :: csr_md_z_t

     ! DIA X CSR
     private :: csr_dm_d    ! real(dp) version
     private :: csr_dm_z    ! complex(dp) version
     !
     private :: csr_dm_d_t
     private :: csr_dm_z_t

     ! CSR + CSR
     private :: csr_plus_d  ! real(dp) version
     private :: csr_plus_z  ! complex(dp) version
     !
     private :: csr_plus_d_t
     private :: csr_plus_z_t

!!========================================================================
!!>>> declare interface and module procedure                           <<<
!!========================================================================

     public :: csr_alloc
     interface csr_alloc
         module procedure csr_alloc_d
         module procedure csr_alloc_z
         module procedure csr_alloc_d_t
         module procedure csr_alloc_z_t
     end interface csr_alloc

     public :: csr_free
     interface csr_free
         module procedure csr_free_d
         module procedure csr_free_z
         module procedure csr_free_d_t
         module procedure csr_free_z_t
     end interface csr_free

     public :: csr_print
     interface csr_print
         module procedure csr_print_d
         module procedure csr_print_z
         module procedure csr_print_d_t
         module procedure csr_print_z_t
     end interface csr_print

     public :: csr_dns
     interface csr_dns
         module procedure csr_dns_d
         module procedure csr_dns_z
         module procedure csr_dns_d_t
         module procedure csr_dns_z_t
     end interface csr_dns

     public :: dns_csr
     interface dns_csr
         module procedure dns_csr_d
         module procedure dns_csr_z
         module procedure dns_csr_d_t
         module procedure dns_csr_z_t
     end interface dns_csr

     public :: csr_csr
     interface csr_csr
         module procedure csr_csr_d
         module procedure csr_csr_z
         module procedure csr_csr_d_t
         module procedure csr_csr_z_t
     end interface csr_csr

     public :: get_csr
     interface get_csr
         module procedure get_csr_d
         module procedure get_csr_z
         module procedure get_csr_d_t
         module procedure get_csr_z_t
     end interface get_csr

     public :: csr_mv
     interface csr_mv
         module procedure csr_mv_d
         module procedure csr_mv_z
         module procedure csr_mv_d_t
         module procedure csr_mv_z_t
     end interface csr_mv

     public :: csr_mm
     interface csr_mm
         module procedure csr_mm_d
         module procedure csr_mm_z
         module procedure csr_mm_d_t
         module procedure csr_mm_z_t
     end interface csr_mm

     public :: csr_md
     interface csr_md
         module procedure csr_md_d
         module procedure csr_md_z
         module procedure csr_md_d_t
         module procedure csr_md_z_t
     end interface csr_md

     public :: csr_dm
     interface csr_dm
         module procedure csr_dm_d
         module procedure csr_dm_z
         module procedure csr_dm_d_t
         module procedure csr_dm_z_t
     end interface csr_dm

     public :: csr_plus
     interface csr_plus
         module procedure csr_plus_d
         module procedure csr_plus_z
         module procedure csr_plus_d_t
         module procedure csr_plus_z_t
     end interface csr_plus

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory                                                  <<<
!!========================================================================

!!
!! @sub csr_alloc_d
!!
!! allocates memory for a row-stored sparse matrix.
!!
  subroutine csr_alloc_d(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, allocatable, intent(inout) :: ia(:)
     integer, allocatable, intent(inout) :: ja(:)
     real(dp), allocatable, intent(inout) :: a(:)

!! local variables
     ! status flag
     integer :: istat

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! allocate memory
     allocate(ia(nrows+1), stat = istat)
     allocate(ja(nnz), stat = istat)
     allocate(a(nnz), stat = istat)
     !
     if ( istat /= 0 ) then
         write(mystd,'(a)') 'sparse: can not allocate enough memory in csr_alloc_d'
         STOP
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     ia = 0
     ja = 0
     a = 0.0_dp

!! body]

     return
  end subroutine csr_alloc_d

!!
!! @sub csr_alloc_z
!!
!! allocates memory for a row-stored sparse matrix.
!!
  subroutine csr_alloc_z(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, allocatable, intent(inout) :: ia(:)
     integer, allocatable, intent(inout) :: ja(:)
     complex(dp), allocatable, intent(inout) :: a(:)

!! local variables
     ! status flag
     integer :: istat

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! allocate memory
     allocate(ia(nrows+1), stat = istat)
     allocate(ja(nnz), stat = istat)
     allocate(a(nnz), stat = istat)
     !
     if ( istat /= 0 ) then
         write(mystd,'(a)') 'sparse: can not allocate enough memory in csr_alloc_z'
         STOP
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     ia = 0
     ja = 0
     a = dcmplx(0.0_dp, 0.0_dp)

!! body]

     return
  end subroutine csr_alloc_z

!!
!! @sub csr_alloc_d_t
!!
!! allocates memory for a row-stored sparse matrix.
!!
  subroutine csr_alloc_d_t(nrows, ncols, nnz, csr)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! csr, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csr

!! local variables
     ! status flag
     integer :: istat

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! allocate memory
     allocate(csr%rowptr(nrows+1), stat = istat)
     allocate(csr%colptr(nnz), stat = istat)
     allocate(csr%V(nnz), stat = istat)
     !
     if ( istat /= 0 ) then
         write(mystd,'(a)') 'sparse: can not allocate enough memory in csr_alloc_d_t'
         STOP
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     csr%rowptr = 0
     csr%colptr = 0
     csr%V = 0.0_dp
     !
     csr%nrows = nrows
     csr%ncols = ncols
     csr%nnz   = nnz

!! body]

     return
  end subroutine csr_alloc_d_t

!!
!! @sub csr_alloc_z_t
!!
!! allocates memory for a row-stored sparse matrix.
!!
  subroutine csr_alloc_z_t(nrows, ncols, nnz, csr)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! csr, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csr

!! local variables
     ! status flag
     integer :: istat

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! allocate memory
     allocate(csr%rowptr(nrows+1), stat = istat)
     allocate(csr%colptr(nnz), stat = istat)
     allocate(csr%V(nnz), stat = istat)
     !
     if ( istat /= 0 ) then
         write(mystd,'(a)') 'sparse: can not allocate enough memory in csr_alloc_z_t'
         STOP
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     csr%rowptr = 0
     csr%colptr = 0
     csr%V = dcmplx(0.0_dp, 0.0_dp)
     !
     csr%nrows = nrows
     csr%ncols = ncols
     csr%nnz   = nnz

!! body]

     return
  end subroutine csr_alloc_z_t

!!========================================================================
!!>>> deallocate memory                                                <<<
!!========================================================================

!!
!! @sub csr_free_d
!!
!! deallocates memory for a row-stored sparse matrix.
!!
  subroutine csr_free_d(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(inout) :: nrows

     ! column dimension of dense matrix
     integer, intent(inout) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(inout) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, allocatable, intent(inout) :: ia(:)
     integer, allocatable, intent(inout) :: ja(:)
     real(dp), allocatable, intent(inout) :: a(:)

!! [body

     if ( allocated(ia) ) deallocate(ia)
     if ( allocated(ja) ) deallocate(ja)
     if ( allocated( a) ) deallocate( a)
     !
     nrows = 0
     ncols = 0
     nnz   = 0

!! body]

     return
  end subroutine csr_free_d

!!
!! @sub csr_free_z
!!
!! deallocates memory for a row-stored sparse matrix.
!!
  subroutine csr_free_z(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(inout) :: nrows

     ! column dimension of dense matrix
     integer, intent(inout) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(inout) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, allocatable, intent(inout) :: ia(:)
     integer, allocatable, intent(inout) :: ja(:)
     complex(dp), allocatable, intent(inout) :: a(:)

!! [body

     if ( allocated(ia) ) deallocate(ia)
     if ( allocated(ja) ) deallocate(ja)
     if ( allocated( a) ) deallocate( a)
     !
     nrows = 0
     ncols = 0
     nnz   = 0

!! body]

     return
  end subroutine csr_free_z

!!
!! @sub csr_free_d_t
!!
!! deallocates memory for a row-stored sparse matrix.
!!
  subroutine csr_free_d_t(csr)
     implicit none

!! external arguments
     ! csr, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csr

!! [body

     if ( allocated(csr%rowptr) ) deallocate(csr%rowptr)
     if ( allocated(csr%colptr) ) deallocate(csr%colptr)
     if ( allocated(csr%V     ) ) deallocate(csr%V     )
     !
     csr%nrows = 0
     csr%ncols = 0
     csr%nnz   = 0

!! body]

     return
  end subroutine csr_free_d_t

!!
!! @sub csr_free_z_t
!!
!! deallocates memory for a row-stored sparse matrix.
!!
  subroutine csr_free_z_t(csr)
     implicit none

!! external arguments
     ! csr, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csr

!! [body

     if ( allocated(csr%rowptr) ) deallocate(csr%rowptr)
     if ( allocated(csr%colptr) ) deallocate(csr%colptr)
     if ( allocated(csr%V     ) ) deallocate(csr%V     )
     !
     csr%nrows = 0
     csr%ncols = 0
     csr%nnz   = 0

!! body]

     return
  end subroutine csr_free_z_t

!!========================================================================
!!>>> input / output                                                   <<<
!!========================================================================

!!
!! @sub csr_print_d
!!
!! prints a row-stored sparse matrix.
!!
  subroutine csr_print_d(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     write(mystd,'(a,i4)') 'number of rows: ', nrows
     write(mystd,'(a,i4)') 'number of columns: ', ncols
     write(mystd,'(a,i4)') 'maximum number of nonzero elements: ', nnz
     !
     write(mystd,'(a)') 'ia:'
     do i=1,nrows+1
         write(mystd,'(2i4)') i, ia(i)
     enddo ! over i={1,nrows+1} loop
     !
     write(mystd,'(a)') 'ja:'
     do i=1,nnz
         write(mystd,'(2i4)') i, ja(i)
     enddo ! over i={1,nnz} loop
     !
     write(mystd,'(a)') 'a:'
     do i=1,nnz
         write(mystd,'(i4,f16.8)') i, a(i)
     enddo ! over i={1,nnz} loop

!! body]

     return
  end subroutine csr_print_d

!!
!! @sub csr_print_z
!!
!! prints a row-stored sparse matrix.
!!
  subroutine csr_print_z(nrows, ncols, nnz, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     write(mystd,'(a,i4)') 'number of rows: ', nrows
     write(mystd,'(a,i4)') 'number of columns: ', ncols
     write(mystd,'(a,i4)') 'maximum number of nonzero elements: ', nnz
     !
     write(mystd,'(a)') 'ia:'
     do i=1,nrows+1
         write(mystd,'(2i4)') i, ia(i)
     enddo ! over i={1,nrows+1} loop
     !
     write(mystd,'(a)') 'ja:'
     do i=1,nnz
         write(mystd,'(2i4)') i, ja(i)
     enddo ! over i={1,nnz} loop
     !
     write(mystd,'(a)') 'a:'
     do i=1,nnz
         write(mystd,'(i4,2f16.8)') i, a(i)
     enddo ! over i={1,nnz} loop

!! body]

     return
  end subroutine csr_print_z

!!
!! @sub csr_print_d_t
!!
!! prints a row-stored sparse matrix.
!!
  subroutine csr_print_d_t(csr)
     implicit none

!! external arguments
     ! csr, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csr

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( csr%nrows <= 0 .or. csr%ncols <= 0 .or. csr%nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     write(mystd,'(a,i4)') 'number of rows: ', csr%nrows
     write(mystd,'(a,i4)') 'number of columns: ', csr%ncols
     write(mystd,'(a,i4)') 'maximum number of nonzero elements: ', csr%nnz
     !
     write(mystd,'(a)') 'ia:'
     do i=1,csr%nrows+1
         write(mystd,'(2i4)') i, csr%rowptr(i)
     enddo ! over i={1,csr%nrows+1} loop
     !
     write(mystd,'(a)') 'ja:'
     do i=1,csr%nnz
         write(mystd,'(2i4)') i, csr%colptr(i)
     enddo ! over i={1,csr%nnz} loop
     !
     write(mystd,'(a)') 'a:'
     do i=1,csr%nnz
         write(mystd,'(i4,f16.8)') i, csr%V(i)
     enddo ! over i={1,csr%nnz} loop

!! body]

     return
  end subroutine csr_print_d_t

!!
!! @sub csr_print_z_t
!!
!! prints a row-stored sparse matrix.
!!
  subroutine csr_print_z_t(csr)
     implicit none

!! external arguments
     ! csr, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csr

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( csr%nrows <= 0 .or. csr%ncols <= 0 .or. csr%nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     write(mystd,'(a,i4)') 'number of rows: ', csr%nrows
     write(mystd,'(a,i4)') 'number of columns: ', csr%ncols
     write(mystd,'(a,i4)') 'maximum number of nonzero elements: ', csr%nnz
     !
     write(mystd,'(a)') 'ia:'
     do i=1,csr%nrows+1
         write(mystd,'(2i4)') i, csr%rowptr(i)
     enddo ! over i={1,csr%nrows+1} loop
     !
     write(mystd,'(a)') 'ja:'
     do i=1,csr%nnz
         write(mystd,'(2i4)') i, csr%colptr(i)
     enddo ! over i={1,csr%nnz} loop
     !
     write(mystd,'(a)') 'a:'
     do i=1,csr%nnz
         write(mystd,'(i4,2f16.8)') i, csr%V(i)
     enddo ! over i={1,csr%nnz} loop

!! body]

     return
  end subroutine csr_print_z_t

!!========================================================================
!!>>> convert sparse matrix into dense matrix                          <<<
!!========================================================================

!!
!! @sub csr_dns_d
!!
!! converts a row-stored sparse matrix into a densely stored one.
!!
  subroutine csr_dns_d(nrows, ncols, nnz, ia, ja, a, dns)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! array where to store dense matrix
     real(dp), intent(out) :: dns(nrows,ncols)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init dns matrix
     dns = 0.0_dp

     ! convert sparse matrix to dense matrix
     do i=1,nrows
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncols ) then
                 write(mystd,'(a)') 'sparse: error in csr_dns_d'
                 STOP
             endif ! back if ( j > ncols ) block
             dns(i,j) = a(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_dns_d

!!
!! @sub csr_dns_z
!!
!! converts a row-stored sparse matrix into a densely stored one.
!!
  subroutine csr_dns_z(nrows, ncols, nnz, ia, ja, a, dns)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! array where to store dense matrix
     complex(dp), intent(out) :: dns(nrows,ncols)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init dns matrix
     dns = dcmplx(0.0_dp, 0.0_dp)

     ! convert sparse matrix to dense matrix
     do i=1,nrows
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncols ) then
                 write(mystd,'(a)') 'sparse: error in csr_dns_z'
                 STOP
             endif ! back if ( j > ncols ) block
             dns(i,j) = a(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_dns_z

!!
!! @sub csr_dns_d_t
!!
!! converts a row-stored sparse matrix into a densely stored one.
!!
  subroutine csr_dns_d_t(nrows, ncols, csr, dns)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! csr, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csr

     ! array where to store dense matrix
     real(dp), intent(out) :: dns(nrows,ncols)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init dns matrix
     dns = 0.0_dp

     ! convert sparse matrix to dense matrix
     do i=1,csr%nrows
         do k=csr%rowptr(i),csr%rowptr(i+1)-1
             j = csr%colptr(k)
             if ( j > csr%ncols ) then
                 write(mystd,'(a)') 'sparse: error in csr_dns_d_t'
                 STOP
             endif ! back if ( j > csr%ncols ) block
             dns(i,j) = csr%V(k)
         enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop
     enddo ! over i={1,csr%nrows} loop

!! body]

     return
  end subroutine csr_dns_d_t

!!
!! @sub csr_dns_z_t
!!
!! converts a row-stored sparse matrix into a densely stored one.
!!
  subroutine csr_dns_z_t(nrows, ncols, csr, dns)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! csr, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csr

     ! array where to store dense matrix
     complex(dp), intent(out) :: dns(nrows,ncols)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init dns matrix
     dns = dcmplx(0.0_dp, 0.0_dp)

     ! convert sparse matrix to dense matrix
     do i=1,csr%nrows
         do k=csr%rowptr(i),csr%rowptr(i+1)-1
             j = csr%colptr(k)
             if ( j > csr%ncols ) then
                 write(mystd,'(a)') 'sparse: error in csr_dns_z_t'
                 STOP
             endif ! back if ( j > csr%ncols ) block
             dns(i,j) = csr%V(k)
         enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop
     enddo ! over i={1,csr%nrows} loop

!! body]

     return
  end subroutine csr_dns_z_t

!!========================================================================
!!>>> convert dense matrix into sparse matrix                          <<<
!!========================================================================

!!
!! @sub dns_csr_d
!!
!! converts a densely stored matrix into a row orientied compactly
!! sparse matrix.
!!
  subroutine dns_csr_d(nrows, ncols, nnz, dns, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! input densely stored matrix
     real(dp), intent(in) :: dns(nrows,ncols)

     ! ia, ja, a, output matrix in compressed sparse row format
     integer, intent(out) :: ia(nrows+1)
     integer, intent(out) :: ja(nnz)
     real(dp), intent(out) :: a(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix
     a = 0.0_dp
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrows
         do j=1,ncols
             if ( dns(i,j) == 0.0_dp ) CYCLE
             ja(k) = j
             a(k) = dns(i,j)
             k = k + 1
             if ( k > nnz ) then
                 write(mystd,'(a)') 'sparse: error in dns_csr_d'
                 STOP
             endif ! back if ( k > nnz ) block
         enddo ! over j={1,ncols} loop
         ia(i+1) = k
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine dns_csr_d

!!
!! @sub dns_csr_z
!!
!! converts a densely stored matrix into a row orientied compactly
!! sparse matrix.
!!
  subroutine dns_csr_z(nrows, ncols, nnz, dns, ia, ja, a)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! input densely stored matrix
     complex(dp), intent(in) :: dns(nrows,ncols)

     ! ia, ja, a, output matrix in compressed sparse row format
     integer, intent(out) :: ia(nrows+1)
     integer, intent(out) :: ja(nnz)
     complex(dp), intent(out) :: a(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix
     a = dcmplx(0.0_dp, 0.0_dp)
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrows
         do j=1,ncols
             if ( real( dns(i,j) ) == 0.0_dp .and. &
                & aimag( dns(i,j) ) == 0.0_dp ) CYCLE
             ja(k) = j
             a(k) = dns(i,j)
             k = k + 1
             if ( k > nnz ) then
                 write(mystd,'(a)') 'sparse: error in dns_csr_z'
                 STOP
             endif ! back if ( k > nnz ) block
         enddo ! over j={1,ncols} loop
         ia(i+1) = k
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine dns_csr_z

!!
!! @sub dns_csr_d_t
!!
!! converts a densely stored matrix into a row orientied compactly
!! sparse matrix.
!!
  subroutine dns_csr_d_t(nrows, ncols, dns, csr)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! input densely stored matrix
     real(dp), intent(in) :: dns(nrows,ncols)

     ! csr, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csr

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix
     csr%V = 0.0_dp
     csr%rowptr = 0
     csr%colptr = 0

     k = 1
     csr%rowptr(1) = 1
     do i=1,nrows
         do j=1,ncols
             if ( dns(i,j) == 0.0_dp ) CYCLE
             csr%colptr(k) = j
             csr%V(k) = dns(i,j)
             k = k + 1
             if ( k > csr%nnz ) then
                 write(mystd,'(a)') 'sparse: error in dns_csr_d_t'
                 STOP
             endif ! back if ( k > csr%nnz ) block
         enddo ! over j={1,ncols} loop
         csr%rowptr(i+1) = k
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine dns_csr_d_t

!!
!! @sub dns_csr_z_t
!!
!! converts a densely stored matrix into a row orientied compactly
!! sparse matrix.
!!
  subroutine dns_csr_z_t(nrows, ncols, dns, csr)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! input densely stored matrix
     complex(dp), intent(in) :: dns(nrows,ncols)

     ! csr, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csr

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix
     csr%V = dcmplx(0.0_dp, 0.0_dp)
     csr%rowptr = 0
     csr%colptr = 0

     k = 1
     csr%rowptr(1) = 1
     do i=1,nrows
         do j=1,ncols
             if ( real( dns(i,j) ) == 0.0_dp .and. &
                & aimag( dns(i,j) ) == 0.0_dp ) CYCLE
             csr%colptr(k) = j
             csr%V(k) = dns(i,j)
             k = k + 1
             if ( k > csr%nnz ) then
                 write(mystd,'(a)') 'sparse: error in dns_csr_z_t'
                 STOP
             endif ! back if ( k > csr%nnz ) block
         enddo ! over j={1,ncols} loop
         csr%rowptr(i+1) = k
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine dns_csr_z_t

!!========================================================================
!!>>> copy sparse matrix to sparse matrix                              <<<
!!========================================================================

!!
!! @sub csr_csr_d
!!
!! copy data between two row orientied compactly sparse matrices.
!!
  subroutine csr_csr_d(nrows, nnz, ia, ja, a, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     real(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do i=ia(1),ia(nrows+1)-1
         jb(i) = ja(i)
     enddo ! over i={ia(1),ia(nrows+1)-1} loop

     do i=ia(1),ia(nrows+1)-1
         b(i) = a(i)
     enddo ! over i={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_csr_d

!!
!! @sub csr_csr_z
!!
!! copy data between two row orientied compactly sparse matrices.
!!
  subroutine csr_csr_z(nrows, nnz, ia, ja, a, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     complex(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do i=ia(1),ia(nrows+1)-1
         jb(i) = ja(i)
     enddo ! over i={ia(1),ia(nrows+1)-1} loop

     do i=ia(1),ia(nrows+1)-1
         b(i) = a(i)
     enddo ! over i={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_csr_z

!!
!! @sub csr_csr_d_t
!!
!! copy data between two row orientied compactly sparse matrices.
!!
  subroutine csr_csr_d_t(csra, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csra

     ! csrb, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do i=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(i) = csra%colptr(i)
     enddo ! over i={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

     do i=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%V(i) = csra%V(i)
     enddo ! over i={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_csr_d_t

!!
!! @sub csr_csr_z_t
!!
!! copy data between two row orientied compactly sparse matrices.
!!
  subroutine csr_csr_z_t(csra, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csra

     ! csrb, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do i=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(i) = csra%colptr(i)
     enddo ! over i={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

     do i=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%V(i) = csra%V(i)
     enddo ! over i={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_csr_z_t

!!========================================================================
!!>>> visit elements in sparse matrix                                  <<<
!!========================================================================

!!
!! @fun get_csr_d
!!
!! this function returns the element a(i,j) of matrix a, which is stored
!! in compressed sparse row format.
!!
  real(dp) &
  function get_csr_d(i, j, nrows, nnz, ia, ja, a) result(elm)
     implicit none

!! external arguments
     ! the row index of the element sought
     integer, intent(in) :: i

     ! the column index of the element sought
     integer, intent(in) :: j

     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

!! local variables
     ! loop index
     integer :: k

     ! memory address of a(i,j)
     integer :: addr

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! initialization
     addr = 0
     elm = 0.0_dp

     ! scan the row - exit as soon as a(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

     ! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = a(addr)
     endif ! back if ( addr /= 0 ) block

!! body]

     return
  end function get_csr_d

!!
!! @fun get_csr_z
!!
!! this function returns the element a(i,j) of matrix a, which is stored
!! in compressed sparse row format.
!!
  complex(dp) &
  function get_csr_z(i, j, nrows, nnz, ia, ja, a) result(elm)
     implicit none

!! external arguments
     ! the row index of the element sought
     integer, intent(in) :: i

     ! the column index of the element sought
     integer, intent(in) :: j

     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

!! local variables
     ! loop index
     integer :: k

     ! memory address of a(i,j)
     integer :: addr

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! initialization
     addr = 0
     elm = dcmplx(0.0_dp, 0.0_dp)

     ! scan the row - exit as soon as a(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

     ! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = a(addr)
     endif ! back if ( addr /= 0 ) block

!! body]

     return
  end function get_csr_z

!!
!! @fun get_csr_d_t
!!
!! this function returns the element a(i,j) of matrix a, which is stored
!! in compressed sparse row format.
!!
  real(dp) &
  function get_csr_d_t(i, j, csr) result(elm)
     implicit none

!! external arguments
     ! the row index of the element sought
     integer, intent(in) :: i

     ! the column index of the element sought
     integer, intent(in) :: j

     ! csr, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csr

!! local variables
     ! loop index
     integer :: k

     ! memory address of a(i,j)
     integer :: addr

!! [body

     ! initialization
     addr = 0
     elm = 0.0_dp

     ! scan the row - exit as soon as a(i,j) is found
     do k=csr%rowptr(i),csr%rowptr(i+1)-1
         if ( csr%colptr(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( csr%colptr(k) == j ) block
     enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop

     ! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = csr%V(addr)
     endif ! back if ( addr /= 0 ) block

!! body]

     return
  end function get_csr_d_t

!!
!! @fun get_csr_z_t
!!
!! this function returns the element a(i,j) of matrix a, which is stored
!! in compressed sparse row format.
!!
  complex(dp) &
  function get_csr_z_t(i, j, csr) result(elm)
     implicit none

!! external arguments
     ! the row index of the element sought
     integer, intent(in) :: i

     ! the column index of the element sought
     integer, intent(in) :: j

     ! csr, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csr

!! local variables
     ! loop index
     integer :: k

     ! memory address of a(i,j)
     integer :: addr

!! [body

     ! initialization
     addr = 0
     elm = dcmplx(0.0_dp, 0.0_dp)

     ! scan the row - exit as soon as a(i,j) is found
     do k=csr%rowptr(i),csr%rowptr(i+1)-1
         if ( csr%colptr(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( csr%colptr(k) == j ) block
     enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop

     ! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = csr%V(addr)
     endif ! back if ( addr /= 0 ) block

!! body]

     return
  end function get_csr_z_t

!!========================================================================
!!>>> sparse matrix-vector multiplication                              <<<
!!========================================================================

!!
!! @sub csr_mv_d
!!
!! multiplies a matrix by a vector using the dot product form.
!!
  subroutine csr_mv_d(nrows, ncols, nnz, ia, ja, a, x, y)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! vector, length equal to the column dimension of the dense matrix
     real(dp), intent(in) :: x(ncols)

     ! vector, real array of length nrows, containing the product y = A . x
     real(dp), intent(out) :: y(nrows)

!! local variables
     ! loop index
     integer :: i
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! zero out output vector
     y = 0.0_dp

     ! compute the inner product of row i with vector x
     do i=1,nrows
         do k=ia(i),ia(i+1)-1
             y(i) = y(i) + a(k) * x( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_mv_d

!!
!! @sub csr_mv_z
!!
!! multiplies a matrix by a vector using the dot product form.
!!
  subroutine csr_mv_z(nrows, ncols, nnz, ia, ja, a, x, y)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! vector, length equal to the column dimension of the dense matrix
     complex(dp), intent(in) :: x(ncols)

     ! vector, complex(dp) array of length nrows, containing the product y = A . x
     complex(dp), intent(out) :: y(nrows)

!! local variables
     ! loop index
     integer :: i
     integer :: k

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! zero out output vector
     y = dcmplx(0.0_dp, 0.0_dp)

     ! compute the inner product of row i with vector sx
     do i=1,nrows
         do k=ia(i),ia(i+1)-1
             y(i) = y(i) + a(k) * x( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_mv_z

!!
!! @sub csr_mv_d_t
!!
!! multiplies a matrix by a vector using the dot product form.
!!
  subroutine csr_mv_d_t(nrows, ncols, csr, x, y)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! csr, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csr

     ! vector, length equal to the column dimension of the dense matrix
     real(dp), intent(in) :: x(ncols)

     ! vector, real array of length nrows, containing the product y = A . x
     real(dp), intent(out) :: y(nrows)

!! local variables
     ! loop index
     integer :: i
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! zero out output vector
     y = 0.0_dp

     ! compute the inner product of row i with vector x
     do i=1,nrows
         do k=csr%rowptr(i),csr%rowptr(i+1)-1
             y(i) = y(i) + csr%V(k) * x( csr%colptr(k) )
         enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_mv_d_t

!!
!! @sub csr_mv_z_t
!!
!! multiplies a matrix by a vector using the dot product form.
!!
  subroutine csr_mv_z_t(nrows, ncols, csr, x, y)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! column dimension of dense matrix
     integer, intent(in) :: ncols

     ! csr, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csr

     ! vector, length equal to the column dimension of the dense matrix
     complex(dp), intent(in) :: x(ncols)

     ! vector, complex(dp) array of length nrows, containing the product y = A . x
     complex(dp), intent(out) :: y(nrows)

!! local variables
     ! loop index
     integer :: i
     integer :: k

!! [body

     ! check dimensions
     if ( csr%nrows /= nrows .or. csr%ncols /= ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! zero out output vector
     y = dcmplx(0.0_dp, 0.0_dp)

     ! compute the inner product of row i with vector sx
     do i=1,nrows
         do k=csr%rowptr(i),csr%rowptr(i+1)-1
             y(i) = y(i) + csr%V(k) * x( csr%colptr(k) )
         enddo ! over k={csr%rowptr(i),csr%rowptr(i+1)-1} loop
     enddo ! over i={1,nrows} loop

!! body]

     return
  end subroutine csr_mv_z_t

!!========================================================================
!!>>> sparse matrix-matrix multiplication                              <<<
!!========================================================================

!!
!! @sub csr_mm_d
!!
!! performs the matrix by matrix product C = A * B.
!!
  subroutine csr_mm_d(nrows, ndims, ncols, nnz, ia, ja, a, ib, jb, b, ic, jc, c)
     implicit none

!! external arguments
     ! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in) :: nrows

     ! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in) :: ndims

     ! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in) :: ncols

     ! the length of the arrays c and jc.
     !
     ! this subroutine will stop if the result matrix C has a number of
     ! elements that exceeds nnz.
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     ! matrix A -> shape(nrows,ndims)
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! ib, jb, b, input matrix in compressed sparse row format
     ! matrix B -> shape(ndims,ncols)
     integer, intent(in) :: ib(ndims+1)
     integer, intent(in) :: jb(nnz)
     real(dp), intent(in) :: b(nnz)

     ! ic, jc, c, output matrix in compressed sparse row format
     ! matrix C -> shape(nrows,ncols)
     integer, intent(out) :: ic(nrows+1)
     integer, intent(out) :: jc(nnz)
     real(dp), intent(out) :: c(nnz)

!! local variables
     ! loop index
     integer :: i, j, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in
     ! matrix B, which is an array that has nonzero value if the column
     ! index already exist, in which case the value is the index of
     ! that column.
     integer :: iw(ncols)

     ! dummy real(dp) variables, used to improve the ratio of floating
     ! point operations to memory accesses.
     real(dp) :: atmp, btmp

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ndims <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     ic(1) = 1

     q = 0
     do i=1,nrows
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = a(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = b(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     c(q) = atmp * btmp
                 else
                     c(p) = c(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrows} loop

     ! check the number of nonzero elements
     if ( q > nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_mm_d'
         STOP
     endif ! back if ( q > nnz ) block

!! body]

     return
  end subroutine csr_mm_d

!!
!! @sub csr_mm_z
!!
!! performs the matrix by matrix product C = A * B.
!!
  subroutine csr_mm_z(nrows, ndims, ncols, nnz, ia, ja, a, ib, jb, b, ic, jc, c)
     implicit none

!! external arguments
     ! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in) :: nrows

     ! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in) :: ndims

     ! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in) :: ncols

     ! the length of the arrays c and jc.
     !
     ! this subroutine will stop if the result matrix C has a number of
     ! elements that exceeds nnz.
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     ! matrix A -> shape(nrows,ndims)
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! ib, jb, b, input matrix in compressed sparse row format
     ! matrix B -> shape(ndims,ncols)
     integer, intent(in) :: ib(ndims+1)
     integer, intent(in) :: jb(nnz)
     complex(dp), intent(in) :: b(nnz)

     ! ic, jc, c, output matrix in compressed sparse row format
     ! matrix C -> shape(nrows,ncols)
     integer, intent(out) :: ic(nrows+1)
     integer, intent(out) :: jc(nnz)
     complex(dp), intent(out) :: c(nnz)

!! local variables
     ! loop index
     integer :: i, j, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in
     ! matrix B, which is an array that has nonzero value if the column
     ! index already exist, in which case the value is the index of
     ! that column.
     integer :: iw(ncols)

     ! dummy complex(dp) variables, used to improve the ratio of
     ! floating point operations to memory accesses.
     complex(dp) :: atmp, btmp

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ndims <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     ic(1) = 1

     q = 0
     do i=1,nrows
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = a(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = b(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     c(q) = atmp * btmp
                 else
                     c(p) = c(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrows} loop

     ! check the number of nonzero elements
     if ( q > nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_mm_z'
         STOP
     endif ! back if ( q > nnz ) block

!! body]

     return
  end subroutine csr_mm_z

!!
!! @sub csr_mm_d_t
!!
!! performs the matrix by matrix product C = A * B.
!!
  subroutine csr_mm_d_t(csra, csrb, csrc)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csra

     ! csrb, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csrb

     ! csrc, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csrc

!! local variables
     ! loop index
     integer :: i, j, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in
     ! matrix B, which is an array that has nonzero value if the column
     ! index already exist, in which case the value is the index of
     ! that column.
     integer :: iw(csrb%ncols)

     ! dummy real(dp) variables, used to improve the ratio of
     ! floating point operations to memory accesses.
     real(dp) :: atmp, btmp

!! [body

     ! check dimensions
     if ( csra%nrows /= csrc%nrows .or. &
          csra%ncols /= csrb%nrows .or. &
          csrb%ncols /= csrc%ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block
     !
     if ( csra%nnz /= csrb%nnz .or. &
          csra%nnz /= csrc%nnz .or. &
          csrb%nnz /= csrc%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     csrc%rowptr(1) = 1

     q = 0
     do i=1,csrc%nrows
         do ka=csra%rowptr(i),csra%rowptr(i+1)-1
             j = csra%colptr(ka)
             atmp = csra%V(ka)
             do kb=csrb%rowptr(j),csrb%rowptr(j+1)-1
                 k = csrb%colptr(kb)
                 btmp = csrb%V(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     csrc%colptr(q) = k
                     csrc%V(q) = atmp * btmp
                 else
                     csrc%V(p) = csrc%V(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={csrb%rowptr(j),csrb%rowptr(j+1)-1} loop
         enddo ! over ka={csra%rowptr(i),csra%rowptr(i+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=csrc%rowptr(i),q
             iw( csrc%colptr( k ) ) = 0
         enddo ! over k={csrc%rowptr(i),q} loop
         csrc%rowptr(i+1) = q + 1
     enddo ! over i={1,csrc%nrows} loop

     ! check the number of nonzero elements
     if ( q > csrc%nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_mm_d_t'
         STOP
     endif ! back if ( q > csrc%nnz ) block

!! body]

     return
  end subroutine csr_mm_d_t

!!
!! @sub csr_mm_z_t
!!
!! performs the matrix by matrix product C = A * B.
!!
  subroutine csr_mm_z_t(csra, csrb, csrc)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csra

     ! csrb, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csrb

     ! csrc, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csrc

!! local variables
     ! loop index
     integer :: i, j, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in
     ! matrix B, which is an array that has nonzero value if the column
     ! index already exist, in which case the value is the index of
     ! that column.
     integer :: iw(csrb%ncols)

     ! dummy complex(dp) variables, used to improve the ratio of
     ! floating point operations to memory accesses.
     complex(dp) :: atmp, btmp

!! [body

     ! check dimensions
     if ( csra%nrows /= csrc%nrows .or. &
          csra%ncols /= csrb%nrows .or. &
          csrb%ncols /= csrc%ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block
     !
     if ( csra%nnz /= csrb%nnz .or. &
          csra%nnz /= csrc%nnz .or. &
          csrb%nnz /= csrc%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     csrc%rowptr(1) = 1

     q = 0
     do i=1,csrc%nrows
         do ka=csra%rowptr(i),csra%rowptr(i+1)-1
             j = csra%colptr(ka)
             atmp = csra%V(ka)
             do kb=csrb%rowptr(j),csrb%rowptr(j+1)-1
                 k = csrb%colptr(kb)
                 btmp = csrb%V(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     csrc%colptr(q) = k
                     csrc%V(q) = atmp * btmp
                 else
                     csrc%V(p) = csrc%V(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={csrb%rowptr(j),csrb%rowptr(j+1)-1} loop
         enddo ! over ka={csra%rowptr(i),csra%rowptr(i+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=csrc%rowptr(i),q
             iw( csrc%colptr( k ) ) = 0
         enddo ! over k={csrc%rowptr(i),q} loop
         csrc%rowptr(i+1) = q + 1
     enddo ! over i={1,csrc%nrows} loop

     ! check the number of nonzero elements
     if ( q > csrc%nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_mm_z_t'
         STOP
     endif ! back if ( q > csrc%nnz ) block

!! body]

     return
  end subroutine csr_mm_z_t

!!========================================================================
!!>>> sparse matrix-diagonal matrix multiplication                     <<<
!!========================================================================

!!
!! @sub csr_md_d
!!
!! performs the matrix by matrix product B = A * Diag.
!!
  subroutine csr_md_d(nrows, nnz, ia, ja, a, diag, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! diagonal matrix stored as a vector diag
     real(dp), intent(in) :: diag(nrows)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     real(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     b = 0.0_dp
     ib = 0
     jb = 0

     ! scale each element
     do i=1,nrows
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrows} loop

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do k=ia(1),ia(nrows+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_md_d

!!
!! @sub csr_md_z
!!
!! performs the matrix by matrix product B = A * Diag.
!!
  subroutine csr_md_z(nrows, nnz, ia, ja, a, diag, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! diagonal matrix stored as a vector diag
     complex(dp), intent(in) :: diag(nrows)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     complex(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     b = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

     ! scale each element
     do i=1,nrows
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrows} loop

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do k=ia(1),ia(nrows+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_md_z

!!
!! @sub csr_md_d_t
!!
!! performs the matrix by matrix product B = A * Diag.
!!
  subroutine csr_md_d_t(csra, diag, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csra

     ! diagonal matrix stored as a vector diag
     real(dp), intent(in) :: diag(csra%nrows)

     ! csrb, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     csrb%V = 0.0_dp
     csrb%rowptr = 0
     csrb%colptr = 0

     ! scale each element
     do i=1,csra%nrows
         k1 = csra%rowptr(i)
         k2 = csra%rowptr(i+1) - 1
         do k=k1,k2
             csrb%V(k) = csra%V(k) * diag( csra%colptr(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,csra%nrows} loop

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do k=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(k) = csra%colptr(k)
     enddo ! over k={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_md_d_t

!!
!! @sub csr_md_z_t
!!
!! performs the matrix by matrix product B = A * Diag.
!!
  subroutine csr_md_z_t(csra, diag, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csra

     ! diagonal matrix stored as a vector diag
     complex(dp), intent(in) :: diag(csra%nrows)

     ! csrb, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     csrb%V = dcmplx(0.0_dp, 0.0_dp)
     csrb%rowptr = 0
     csrb%colptr = 0

     ! scale each element
     do i=1,csra%nrows
         k1 = csra%rowptr(i)
         k2 = csra%rowptr(i+1) - 1
         do k=k1,k2
             csrb%V(k) = csra%V(k) * diag( csra%colptr(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,csra%nrows} loop

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do k=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(k) = csra%colptr(k)
     enddo ! over k={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_md_z_t

!!========================================================================
!!>>> diagonal matrix-sparse matrix multiplication                     <<<
!!========================================================================

!!
!! @sub csr_dm_d
!!
!! performs the matrix by matrix product B = Diag * A.
!!
  subroutine csr_dm_d(nrows, nnz, diag, ia, ja, a, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! diagonal matrix stored as a vector diag
     real(dp), intent(in) :: diag(nrows)

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     real(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     b = 0.0_dp
     ib = 0
     jb = 0

     ! normalize each row
     do i=1,nrows
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrows} loop

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do k=ia(1),ia(nrows+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_dm_d

!!
!! @sub csr_dm_z
!!
!! performs the matrix by matrix product B = Diag * A.
!!
  subroutine csr_dm_z(nrows, nnz, diag, ia, ja, a, ib, jb, b)
     implicit none

!! external arguments
     ! row dimension of dense matrix
     integer, intent(in) :: nrows

     ! maximum number of nonzero elements allowed
     integer, intent(in) :: nnz

     ! diagonal matrix stored as a vector diag
     complex(dp), intent(in) :: diag(nrows)

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! ib, jb, b, output matrix in compressed sparse row format
     integer, intent(out) :: ib(nrows+1)
     integer, intent(out) :: jb(nnz)
     complex(dp), intent(out) :: b(nnz)

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     b = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

     ! normalize each row
     do i=1,nrows
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrows} loop

     do i=1,nrows+1
         ib(i) = ia(i)
     enddo ! over i={1,nrows+1} loop

     do k=ia(1),ia(nrows+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrows+1)-1} loop

!! body]

     return
  end subroutine csr_dm_z

!!
!! @sub csr_dm_d_t
!!
!! performs the matrix by matrix product B = Diag * A.
!!
  subroutine csr_dm_d_t(diag, csra, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csra

     ! diagonal matrix stored as a vector diag
     real(dp), intent(in) :: diag(csra%nrows)

     ! csrb, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     csrb%V = 0.0_dp
     csrb%rowptr = 0
     csrb%colptr = 0

     ! normalize each row
     do i=1,csra%nrows
         k1 = csra%rowptr(i)
         k2 = csra%rowptr(i+1) - 1
         do k=k1,k2
             csrb%V(k) = csra%V(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,csra%nrows} loop

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do k=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(k) = csra%colptr(k)
     enddo ! over k={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_dm_d_t

!!
!! @sub csr_dm_z_t
!!
!! performs the matrix by matrix product B = Diag * A.
!!
  subroutine csr_dm_z_t(diag, csra, csrb)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csra

     ! diagonal matrix stored as a vector diag
     complex(dp), intent(in) :: diag(csra%nrows)

     ! csrb, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csrb

!! local variables
     ! loop index
     integer :: i
     integer :: k

     ! loop index
     integer :: k1
     integer :: k2

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%nnz   /= csrb%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init sparse matrix B
     csrb%V = dcmplx(0.0_dp, 0.0_dp)
     csrb%rowptr = 0
     csrb%colptr = 0

     ! normalize each row
     do i=1,csra%nrows
         k1 = csra%rowptr(i)
         k2 = csra%rowptr(i+1) - 1
         do k=k1,k2
             csrb%V(k) = csra%V(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,csra%nrows} loop

     do i=1,csra%nrows+1
         csrb%rowptr(i) = csra%rowptr(i)
     enddo ! over i={1,csra%nrows+1} loop

     do k=csra%rowptr(1),csra%rowptr(csra%nrows+1)-1
         csrb%colptr(k) = csra%colptr(k)
     enddo ! over k={csra%rowptr(1),csra%rowptr(csra%nrows+1)-1} loop

!! body]

     return
  end subroutine csr_dm_z_t

!!========================================================================
!!>>> sparse matrix-matrix addition                                    <<<
!!========================================================================

!!
!! @sub csr_plus_d
!!
!! performs the CSR matrix sum C = A + B.
!!
  subroutine csr_plus_d(nrows, ncols, nnz, ia, ja, a, ib, jb, b, ic, jc, c)
     implicit none

!! external arguments
     ! row dimension of A and B.
     integer, intent(in) :: nrows

     ! column dimension of A and B.
     integer, intent(in) :: ncols

     ! the length of the arrays c and jc.
     !
     ! this subroutine will stop if the result matrix C has a number of
     ! elements that exceeds nnz.
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     real(dp), intent(in) :: a(nnz)

     ! ib, jb, b, input matrix in compressed sparse row format
     integer, intent(in) :: ib(nrows+1)
     integer, intent(in) :: jb(nnz)
     real(dp), intent(in) :: b(nnz)

     ! ic, jc, c, output matrix in compressed sparse row format
     integer, intent(out) :: ic(nrows+1)
     integer, intent(out) :: jc(nnz)
     real(dp), intent(out) :: c(nnz)

!! local variables
     ! loop index
     integer :: i, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in A
     integer :: iw(ncols)

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     ic(1) = 1

     q = 0
     do i=1,nrows
         do ka=ia(i),ia(i+1)-1
             q = q + 1
             k = ja(ka)
             iw(k) = q
             jc(q) = k
             c(q) = a(ka)
         enddo ! over ka={ia(i),ia(i+1)-1} loop
         !
         do kb=ib(i),ib(i+1)-1
             k = jb(kb)

             p = iw(k)
             if ( p == 0 ) then
                 q = q + 1
                 iw(k) = q
                 jc(q) = k
                 c(q) = b(kb)
             else
                 c(p) = c(p) + b(kb)
             endif ! back if ( p == 0 ) block
         enddo ! over kb={ib(j),ib(j+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc(k) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrows} loop

     ! check the number of nonzero elements
     if ( q > nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_plus_d'
         STOP
     endif ! back if ( q > nnz ) block

!! body]

     return
  end subroutine csr_plus_d

!!
!! @sub csr_plus_z
!!
!! performs the CSR matrix sum C = A + B.
!!
  subroutine csr_plus_z(nrows, ncols, nnz, ia, ja, a, ib, jb, b, ic, jc, c)
     implicit none

!! external arguments
     ! row dimension of A and B.
     integer, intent(in) :: nrows

     ! column dimension of A and B.
     integer, intent(in) :: ncols

     ! the length of the arrays c and jc.
     !
     ! this subroutine will stop if the result matrix C has a number of
     ! elements that exceeds nnz.
     integer, intent(in) :: nnz

     ! ia, ja, a, input matrix in compressed sparse row format
     integer, intent(in) :: ia(nrows+1)
     integer, intent(in) :: ja(nnz)
     complex(dp), intent(in) :: a(nnz)

     ! ib, jb, b, input matrix in compressed sparse row format
     integer, intent(in) :: ib(nrows+1)
     integer, intent(in) :: jb(nnz)
     complex(dp), intent(in) :: b(nnz)

     ! ic, jc, c, output matrix in compressed sparse row format
     integer, intent(out) :: ic(nrows+1)
     integer, intent(out) :: jc(nnz)
     complex(dp), intent(out) :: c(nnz)

!! local variables
     ! loop index
     integer :: i, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in A
     integer :: iw(ncols)

!! [body

     ! check dimensions
     if ( nrows <= 0 .or. ncols <= 0 .or. nnz <= 0 ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     ic(1) = 1

     q = 0
     do i=1,nrows
         do ka=ia(i),ia(i+1)-1
             q = q + 1
             k = ja(ka)
             iw(k) = q
             jc(q) = k
             c(q) = a(ka)
         enddo ! over ka={ia(i),ia(i+1)-1} loop
         !
         do kb=ib(i),ib(i+1)-1
             k = jb(kb)

             p = iw(k)
             if ( p == 0 ) then
                 q = q + 1
                 iw(k) = q
                 jc(q) = k
                 c(q) = b(kb)
             else
                 c(p) = c(p) + b(kb)
             endif ! back if ( p == 0 ) block
         enddo ! over kb={ib(j),ib(j+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc(k) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrows} loop

     ! check the number of nonzero elements
     if ( q > nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_plus_z'
         STOP
     endif ! back if ( q > nnz ) block

!! body]

     return
  end subroutine csr_plus_z

!!
!! @sub csr_plus_d_t
!!
!! performs the CSR matrix sum C = A + B.
!!
  subroutine csr_plus_d_t(csra, csrb, csrc)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csra

     ! csrb, an input matrix in compressed sparse row format
     type (csr_d), intent(in) :: csrb

     ! csrc, an output matrix in compressed sparse row format
     type (csr_d), intent(inout) :: csrc

!! local variables
     ! loop index
     integer :: i, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in A
     integer :: iw(csra%ncols)

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%nrows /= csrc%nrows .or. &
          csrb%nrows /= csrc%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%ncols /= csrc%ncols .or. &
          csrb%ncols /= csrc%ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block
     !
     if ( csra%nnz /= csrb%nnz .or. &
          csra%nnz /= csrc%nnz .or. &
          csrb%nnz /= csrc%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     csrc%rowptr(1) = 1

     q = 0
     do i=1,csra%nrows
         do ka=csra%rowptr(i),csra%rowptr(i+1)-1
             q = q + 1
             k = csra%colptr(ka)
             iw(k) = q
             csrc%colptr(q) = k
             csrc%V(q) = csra%V(ka)
         enddo ! over ka={csra%rowptr(i),csra%rowptr(i+1)-1} loop
         !
         do kb=csrb%rowptr(i),csrb%rowptr(i+1)-1
             k = csrb%colptr(kb)

             p = iw(k)
             if ( p == 0 ) then
                 q = q + 1
                 iw(k) = q
                 csrc%colptr(q) = k
                 csrc%V(q) = csrb%V(kb)
             else
                 csrc%V(p) = csrc%V(p) + csrb%V(kb)
             endif ! back if ( p == 0 ) block
         enddo ! over kb={csrb%rowptr(j),csrb%rowptr(j+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=csrc%rowptr(i),q
             iw( csrc%colptr(k) ) = 0
         enddo ! over k={csrc%rowptr(i),q} loop
         csrc%rowptr(i+1) = q + 1
     enddo ! over i={1,csra%nrows} loop

     ! check the number of nonzero elements
     if ( q > csrc%nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_plus_d_t'
         STOP
     endif ! back if ( q > csrc%nnz ) block

!! body]

     return
  end subroutine csr_plus_d_t

!!
!! @sub csr_plus_z_t
!!
!! performs the CSR matrix sum C = A + B.
!!
  subroutine csr_plus_z_t(csra, csrb, csrc)
     implicit none

!! external arguments
     ! csra, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csra

     ! csrb, an input matrix in compressed sparse row format
     type (csr_z), intent(in) :: csrb

     ! csrc, an output matrix in compressed sparse row format
     type (csr_z), intent(inout) :: csrc

!! local variables
     ! loop index
     integer :: i, k

     ! loop index
     integer :: ka, kb

     ! dummy integer variables
     integer :: p, q

     ! integer work array of length equal to the number of columns in A
     integer :: iw(csra%ncols)

!! [body

     ! check dimensions
     if ( csra%nrows /= csrb%nrows .or. &
          csra%nrows /= csrc%nrows .or. &
          csrb%nrows /= csrc%nrows .or. &
          csra%ncols /= csrb%ncols .or. &
          csra%ncols /= csrc%ncols .or. &
          csrb%ncols /= csrc%ncols ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block
     !
     if ( csra%nnz /= csrb%nnz .or. &
          csra%nnz /= csrc%nnz .or. &
          csrb%nnz /= csrc%nnz ) then
         write(mystd,'(a)') 'sparse: wrong dimensions for sparse matrix'
         STOP
     endif ! back if block

     ! init work array
     iw = 0

     ! init sparse matrix C
     csrc%rowptr(1) = 1

     q = 0
     do i=1,csra%nrows
         do ka=csra%rowptr(i),csra%rowptr(i+1)-1
             q = q + 1
             k = csra%colptr(ka)
             iw(k) = q
             csrc%colptr(q) = k
             csrc%V(q) = csra%V(ka)
         enddo ! over ka={csra%rowptr(i),csra%rowptr(i+1)-1} loop
         !
         do kb=csrb%rowptr(i),csrb%rowptr(i+1)-1
             k = csrb%colptr(kb)

             p = iw(k)
             if ( p == 0 ) then
                 q = q + 1
                 iw(k) = q
                 csrc%colptr(q) = k
                 csrc%V(q) = csrb%V(kb)
             else
                 csrc%V(p) = csrc%V(p) + csrb%V(kb)
             endif ! back if ( p == 0 ) block
         enddo ! over kb={csrb%rowptr(j),csrb%rowptr(j+1)-1} loop

         ! done this row i, so set work array to zero again
         do k=csrc%rowptr(i),q
             iw( csrc%colptr(k) ) = 0
         enddo ! over k={csrc%rowptr(i),q} loop
         csrc%rowptr(i+1) = q + 1
     enddo ! over i={1,csra%nrows} loop

     ! check the number of nonzero elements
     if ( q > csrc%nnz ) then
         write(mystd,'(a)') 'sparse: error in csr_plus_z_t'
         STOP
     endif ! back if ( q > csrc%nnz ) block

!! body]

     return
  end subroutine csr_plus_z_t

  end module sparse
