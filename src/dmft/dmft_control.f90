!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : control    module
!!! source  : dmft_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           02/26/2021 by li huang (last modified)
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

     integer, public, save :: nband = 1

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
