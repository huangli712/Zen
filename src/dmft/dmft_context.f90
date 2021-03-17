!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_group
!!!           dmft_window
!!!           dmft_lattice
!!!           dmft_kmesh
!!!           dmft_fmesh
!!!           dmft_tetra
!!!           dmft_eigen
!!!           dmft_projs
!!!           dmft_sigma
!!!           dmft_green
!!!           dmft_weiss
!!!           context    module
!!! source  : dmft_context.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           03/17/2021 by li huang (last modified)
!!! purpose : 
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @mod dmft_group
!!
  module dmft_group
     implicit none

     integer, public, save, allocatable :: site(:,:)
     integer, public, save, allocatable :: l(:,:)
     logical, public, save, allocatable :: corr(:,:)
     character(len=4), public, save, allocatable :: shell(:,:)
     integer, public, save, allocatable :: ndim(:,:)

  end module dmft_group

!!
!! @mod dmft_window
!!
  module dmft_window
     implicit none

     integer, public, save, allocatable :: bmin(:,:)
     integer, public, save, allocatable :: bmax(:,:)
     integer, public, save, allocatable :: nbnd(:,:)
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
!! @mod dmft_kmesh
!!
  module dmft_kmesh
     use constants, only : dp

     implicit none

     real(dp), public, save, allocatable :: kpts(:,:)
     real(dp), public, save, allocatable :: wghts(:)

  end module dmft_kmesh

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
  end module dmft_eigen

!!
!! @mod dmft_projs
!!
  module dmft_projs
     use constants, only : dp

     implicit none

     complex(dp), public, save, allocatable :: psichi(:)

  end module dmft_projs

!!
!! @mod dmft_sigma
!!
  module dmft_sigma
     implicit none
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

     use dmft_lattice

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
     public :: cat_alloc_lattice

! declaration of module procedures: deallocate memory
     public :: cat_free_lattice

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

  subroutine cat_alloc_lattice()
     implicit none

     return
  end subroutine cat_alloc_lattice

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

  subroutine cat_free_lattice()
     implicit none

     return
  end subroutine cat_free_lattice

  end module context
