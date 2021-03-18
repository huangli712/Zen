!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : control    module
!!! source  : dmft_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           03/18/2021 by li huang (last modified)
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
!!>>> integer variables                                                <<<
!!========================================================================

     integer, public, save :: task   = 1
     integer, public, save :: axis   = 1

     integer, public, save :: nsort  = 3
     integer, public, save :: natom  = 5
     integer, public, save :: nband  = 30
     integer, public, save :: nkpt   = 729
     integer, public, save :: nspin  = 1
     integer, public, save :: ntet   = 4374
     integer, public, save :: ngrp   = 1
     integer, public, save :: nwnd   = 1
     integer, public, save :: mfreq  = 8193
     integer, public, save :: nfreq  = 513

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

     real(dp), public, save :: beta  = 8.00_dp
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
