!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : context    module
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
!! @mod dmft_maps
!!
  module dmft_maps
  end module dmft_maps

!!
!! @mod dmft_atoms
!!
  module dmft_atoms
  end module dmft_atoms

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
  end module dmft_sigma

!!
!! @mod dmft_green
!!
  module dmft_green
  end module dmft_green

  module context
     use constants, only : dp
     use constants, only : zero, czero

     implicit none

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

  end module context
