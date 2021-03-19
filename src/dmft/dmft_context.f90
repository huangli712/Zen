!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_group
!!!           dmft_window
!!!           dmft_lattice
!!!           dmft_bzone
!!!           dmft_tetra
!!!           dmft_eigen
!!!           dmft_projs
!!!           dmft_fmesh
!!!           dmft_sigma
!!!           dmft_green
!!!           dmft_weiss
!!!           context    module
!!! source  : dmft_context.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           03/19/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @mod dmft_group
!!
  module dmft_group
     implicit none

     integer, public, save :: max_ndim
     character(len=4), public, save, allocatable :: shell(:)
     integer, public, save, allocatable :: site(:)
     integer, public, save, allocatable :: l(:)
     logical, public, save, allocatable :: corr(:)
     integer, public, save, allocatable :: ndim(:)

  end module dmft_group

!!
!! @mod dmft_window
!!
  module dmft_window
     implicit none

     integer, public, save :: max_nbnd
     integer, public, save, allocatable :: bmin(:)
     integer, public, save, allocatable :: bmax(:)
     integer, public, save, allocatable :: nbnd(:)
     integer, public, save, allocatable :: kwin(:,:,:,:)

  end module dmft_window

!!
!! @mod dmft_lattice
!!
  module dmft_lattice
     use constants, only : dp

     implicit none

     character(len=2), public, save, allocatable :: sorts(:)
     character(len=2), public, save, allocatable :: atoms(:)
     integer,  public, save, allocatable :: sortn(:)
     real(dp), public, save, allocatable :: lvect(:,:)
     real(dp), public, save, allocatable :: coord(:,:)

  end module dmft_lattice

!!
!! @mod dmft_bzone
!!
  module dmft_bzone
     use constants, only : dp

     implicit none

     real(dp), public, save, allocatable :: kmesh(:,:)
     real(dp), public, save, allocatable :: weight(:)

  end module dmft_bzone

!!
!! @mod dmft_tetra
!!
  module dmft_tetra
     implicit none

     integer, public, save, allocatable :: tetra(:,:)

  end module dmft_tetra

!!
!! @mod dmft_eigen
!!
  module dmft_eigen
     use constants, only : dp

     implicit none

     real(dp), public, save, allocatable :: enk(:,:,:)
     real(dp), public, save, allocatable :: occupy(:,:,:)

  end module dmft_eigen

!!
!! @mod dmft_projs
!!
  module dmft_projs
     use constants, only : dp

     implicit none

     complex(dp), public, save, allocatable :: psichi(:,:,:,:,:)

  end module dmft_projs

!!
!! @mod dmft_fmesh
!!
  module dmft_fmesh
     use constants, only : dp

     implicit none

     real(dp), public, save, allocatable :: imesh(:)
     real(dp), public, save, allocatable :: rmesh(:)

  end module dmft_fmesh

!!
!! @mod dmft_sigma
!!
  module dmft_sigma
     implicit none

     complex(dp), public, save, allocatable :: sigdc
     complex(dp), public, save, allocatable :: sig_l
     complex(dp), public, save, allocatable :: sig_k

  end module dmft_sigma

!!
!! @mod dmft_green
!!
  module dmft_green
     implicit none

  end module dmft_green

!!
!! @mod dmft_weiss
!!
  module dmft_weiss
     implicit none

  end module dmft_weiss

  module context
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : ngrp, nwnd
     use control, only : nband, nkpt, nspin
     use control, only : natom, nsort
     use control, only : mfreq, nfreq
     use control, only : ntet

     use dmft_group
     use dmft_window
     use dmft_lattice
     use dmft_bzone
     use dmft_tetra
     use dmft_eigen
     use dmft_projs
     use dmft_fmesh
     use dmft_sigma

     implicit none

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! declaration of module procedures: allocate memory
     public :: cat_alloc_group
     public :: cat_alloc_window
     public :: cat_alloc_lattice
     public :: cat_alloc_bzone
     public :: cat_alloc_tetra
     public :: cat_alloc_eigen
     public :: cat_alloc_projs
     public :: cat_alloc_fmesh
     public :: cat_alloc_sigma

! declaration of module procedures: deallocate memory
     public :: cat_free_group
     public :: cat_free_window
     public :: cat_free_lattice
     public :: cat_free_bzone
     public :: cat_free_tetra
     public :: cat_free_eigen
     public :: cat_free_projs
     public :: cat_free_fmesh
     public :: cat_free_sigma

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_group
!!
  subroutine cat_alloc_group()
     implicit none

! allocate memory
     allocate(site(ngrp), stat = istat)
     allocate(l(ngrp), stat = istat)
     allocate(corr(ngrp), stat = istat)
     allocate(shell(ngrp), stat = istat)
     allocate(ndim(ngrp), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_group','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     site = 0
     l = 0
     corr = .false.
     shell = 's'
     ndim = 0

     return
  end subroutine cat_alloc_group

!!
!! @sub cat_alloc_window
!!
  subroutine cat_alloc_window()
     implicit none

! allocate memory
     allocate(bmin(nwnd), stat = istat)
     allocate(bmax(nwnd), stat = istat)
     allocate(nbnd(nwnd), stat = istat)
     allocate(kwin(nkpt,nspin,2,nwnd), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_window','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     bmin = 0
     bmax = 0
     nbnd = 0
     kwin = 0

     return
  end subroutine cat_alloc_window

!!
!! @sub cat_alloc_lattice
!!
  subroutine cat_alloc_lattice()
     implicit none

! allocate memory
     allocate(sorts(nsort), stat = istat)
     allocate(atoms(natom), stat = istat)
     allocate(sortn(nsort), stat = istat)
     allocate(lvect(3,3), stat = istat)
     allocate(coord(natom,3), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_lattice','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sorts = 'X'
     atoms = 'X'
     sortn = 0
     lvect = zero
     coord = zero

     return
  end subroutine cat_alloc_lattice

!!
!! @sub cat_alloc_bzone
!!
  subroutine cat_alloc_bzone()
     implicit none

! allocate memory
     allocate(kmesh(nkpt,3), stat = istat)
     allocate(weight(nkpt), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_bzone','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     kmesh = zero
     weight = zero

     return
  end subroutine cat_alloc_bzone

!!
!! @sub cat_alloc_tetra
!!
  subroutine cat_alloc_tetra()
     implicit none

! allocate memory
     allocate(tetra(ntet,5), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_tetra','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tetra = 0

     return
  end subroutine cat_alloc_tetra

!!
!! @sub cat_alloc_eigen
!!
  subroutine cat_alloc_eigen()
     implicit none

! allocate memory
     allocate(enk(nband,nkpt,nspin), stat = istat)
     allocate(occupy(nband,nkpt,nspin), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_eigen','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     enk = zero
     occupy = zero

     return
  end subroutine cat_alloc_eigen

!!
!! @sub cat_alloc_projs
!!
  subroutine cat_alloc_projs()
     implicit none

! allocate memory
     allocate(psichi(max_ndim,max_nbnd,nkpt,nspin,ngrp), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_projs','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     psichi = czero

     return
  end subroutine cat_alloc_projs

!!
!! @sub cat_alloc_fmesh
!!
  subroutine cat_alloc_fmesh()
     implicit none

! allocate memory
     allocate(imesh(mfreq), stat = istat)
     allocate(rmesh(nfreq), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fmesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     imesh = zero
     rmesh = zero

     return
  end subroutine cat_alloc_fmesh

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_group
!!
  subroutine cat_free_group()
     implicit none

     return
  end subroutine cat_free_group

!!
!! @sub cat_free_window
!!
  subroutine cat_free_window()
     implicit none

     return
  end subroutine cat_free_window

!!
!! @sub cat_free_lattice
!!
  subroutine cat_free_lattice()
     implicit none

     return
  end subroutine cat_free_lattice

!!
!! @sub cat_free_bzone
!!
  subroutine cat_free_bzone()
     implicit none

     return
  end subroutine cat_free_bzone

!!
!! @sub cat_free_tetra
!!
  subroutine cat_free_tetra()
     implicit none

     return
  end subroutine cat_free_tetra

!!
!! @sub cat_free_eigen
!!
  subroutine cat_free_eigen()
     implicit none

     return
  end subroutine cat_free_eigen

!!
!! @sub cat_free_projs
!!
  subroutine cat_free_projs()
     implicit none

     return
  end subroutine cat_free_projs

!!
!! @sub cat_free_fmesh
!!
  subroutine cat_free_fmesh()
     implicit none

     return
  end subroutine cat_free_fmesh

!!
!! @sub cat_free_sigma
!!
  subroutine cat_free_sigma()
     implicit none

     return
  end subroutine cat_free_sigma

  end module context
