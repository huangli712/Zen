!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : control    module
!!! source  : dmft_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           04/02/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! code name of the current dynamical mean-field theory engine
!!
     character(len = 09), public, save :: cname = 'JACARANDA'

!!========================================================================
!!>>> integer variables: from dmft.in                                  <<<
!!========================================================================

!!
!! @var task
!!
!! control flag, determine the running mode of the code
!!
!! if task == 0:
!!     search the fermi level only
!!
!! if task == 1:
!!     calculate the local hybridization function. of course, the fermi
!!     level should be updated as well
!!
!! if task == 2:
!!     calculate charge correction due to the electronic correlation
!!
     integer, public, save :: task   = 1

!!
!! @var axis
!!
!! control flag, determine the axis for brillouin zone integration
!!
!! if axis == 1:
!!     imaginary axis
!!
!! if axis == 2:
!!     real axis
!!
     integer, public, save :: axis   = 1

!!========================================================================
!!>>> logical variables: from dmft.in                                  <<<
!!========================================================================

!!
!! @var lfermi
!!
!! control flag, determine whether the fermi level should be updated
!!
!! if lfermi == .true.
!!     search the fermi level
!!
!! if lfermi == .false.
!!     fix the fermi leve. in other words, the dft fermi level is used
!!
     logical, public, save :: lfermi = .true.

!!
!! @var ltetra
!!
!! control flag, determine whether the analytical tetrahedron method is
!! used to do the brillouin zone integration
!!
!! if ltetra == .true.
!!     use the analytical tetrahedron algorithm
!!
!! if ltetra == .false.
!!     use the direct algorithm
!!
     logical, public, save :: ltetra = .true.

!!========================================================================
!!>>> real variables: from dmft.in                                     <<<
!!========================================================================

!!
!! @var beta
!!
!! inverse temperature, \beta = 1 / T
!!
     real(dp), public, save :: beta  = 8.00_dp

!!
!! @var mc
!!
!! convergence criterion for fermi level search
!!
     real(dp), public, save :: mc    = 0.0001_dp

!!========================================================================
!!>>> integer variables: from params.ir                                <<<
!!========================================================================

!!========================================================================
!!>>> real variables: from params.ir                                   <<<
!!========================================================================








     integer, public, save :: nsort  = 3
     integer, public, save :: natom  = 5
     integer, public, save :: nband  = 30
     integer, public, save :: nkpt   = 729
     integer, public, save :: nspin  = 1
     integer, public, save :: ntet   = 4374
     integer, public, save :: nsite  = 1
     integer, public, save :: ngrp   = 1
     integer, public, save :: nwnd   = 1
     integer, public, save :: mfreq  = 8193
     integer, public, save :: nfreq  = 513


     real(dp), public, save :: scal  = 4.00_dp
     real(dp), public, save :: volt  = 1.00_dp
     real(dp), public, save :: fermi = 0.00_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

!!
!! @var nprocs
!!
!! number of processors: default value 1
!!
     integer, public, save :: nprocs = 1

!!
!! @var myid
!!
!! the id of current process: default value 0
!!
     integer, public, save :: myid   = 0

!!
!! @var master
!!
!! denote as the controller process: default value 0
!!
     integer, public, save :: master = 0

!!
!! @var cid
!!
!! the id of current process in cartesian topology (cid == myid)
!!
     integer, public, save :: cid    = 0

!!
!! @var cx
!!
!! the x coordinates of current process in cartesian topology
!!
     integer, public, save :: cx     = 0

!!
!! @var cy
!!
!! the y coordinates of current process in cartesian topology
!!
     integer, public, save :: cy     = 0

  end module control
