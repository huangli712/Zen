!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : control    module
!!!           version    module
!!! source  : dmft_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           05/09/2021 by li huang (last modified)
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

!!
!! @var model
!!
!! name of strongly correlated materials or models
!!
     character(len = 99), public, save :: model = 'SrVO3'

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
!!     level might be updated as well (depends on lfermi)
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
!!     fix the fermi level. in other words, the dft fermi level is used
!!
     logical, public, save :: lfermi = .true.

!!
!! @var ltetra
!!
!! control flag, determine whether the analytical tetrahedron method is
!! used to perform the brillouin zone integration
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

!!
!! @var nsort
!!
!! number of atomic sorts in the model
!!
     integer, public, save :: nsort  = 3

!!
!! @var natom
!!
!! number of atoms in the model
!!
     integer, public, save :: natom  = 5

!!
!! @var nband
!!
!! number of dft bands
!!
     integer, public, save :: nband  = 30

!!
!! @var nkpt
!!
!! number of k-mesh points
!!
     integer, public, save :: nkpt   = 729

!!
!! @var nspin
!!
!! number of spin orientations
!!
     integer, public, save :: nspin  = 1

!!
!! @var ntet
!!
!! number of tetrahedra. note that ntet = 1 means that the tetrahedron
!! data are absent and ltetra must be .false.
!!
     integer, public, save :: ntet   = 1

!!
!! @var ngrp
!!
!! number of groups of projectors, which are used to create the Hilbert
!! subspace for correlated orbitals
!!
     integer, public, save :: ngrp   = 1

!!
!! @var nwnd
!!
!! number of energy windows or band windows, which are used to restrict
!! the correlated subspace
!!
     integer, public, save :: nwnd   = 1

!!
!! @var nsite
!!
!! number of correlated electron problems, i.e, number of impurity sites
!! in which the correlated effect is considered
!!
     integer, public, save :: nsite  = 1

!!
!! @var nmesh
!!
!! number of frequency points
!!
     integer, public, save :: nmesh  = 8193

!!========================================================================
!!>>> real variables: from params.ir                                   <<<
!!========================================================================

!!
!! @var scale
!!
!! an universal scaling factor for the lattice constants
!!
     real(dp), public, save :: scale = 4.00_dp

!!
!! @var fermi
!!
!! default fermi level, which is usually taken from the dft calculations
!!
     real(dp), public, save :: fermi = 0.00_dp

!!
!! @var volt
!!
!! volume of a tetrahedron, which is used to renormalize the integration
!! by using the analytical tetrahedron algorithm
!!
     real(dp), public, save :: volt  = 1.00_dp

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



  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL = 'v0.1.2 @ 2021.05.09D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.1.2'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2021.05.09'

!!
!! @var V_STAT
!!
!! version string, only status info., D means devel, T testing, R released.
!!
     character(len=01), public, parameter :: V_STAT = 'D'

!!
!! @var V_AUTH
!!
!! version string, author info.
!!
     character(len=11), public, parameter :: V_AUTH = 'by li huang'

!!
!! @var V_INST
!!
!! version string, affiliation info.
!!
     character(len=36), public, parameter :: V_INST = 'China Academy of Engineering Physics'

!!
!! @var V_MAIL
!!
!! version string, email info.
!!
     character(len=22), public, parameter :: V_MAIL = 'lihuang.dmft@gmail.com'

!!
!! @var V_GPL3
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: V_GPL3 = 'GNU General Public License version 3'

  end module version
