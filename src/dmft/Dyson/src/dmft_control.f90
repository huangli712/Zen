!!!-----------------------------------------------------------------------
!!! project : dyson @ jacaranda
!!! program : control module
!!!           version module
!!! source  : dmft_control.f90
!!! type    : modules
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 02/23/2021 by li huang (created)
!!!           03/26/2025 by li huang (last modified)
!!! purpose : define the global control variables.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module control                                                   <<<
!!========================================================================

!!
!! @mod control
!!
!! define the control parameters and dimensional parameters.
!!
  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! code name of the current dynamical mean-field theory engine.
!!
     character(len = 09), public, save :: cname = 'JACARANDA'

!!
!! @var model
!!
!! name of strongly correlated materials or models.
!!
     character(len = 99), public, save :: model = 'SrVO3'

!!========================================================================
!!>>> integer variables: from dmft.in                                  <<<
!!========================================================================

!!
!! @var task
!!
!! control flag, determine the running mode of the code.
!!
!! if task == 1:
!!     the code will do the following jobs:
!!     (1) search the fermi level (depends on `lfermi`),
!!     (2) calculate the impurity level,
!!     (3) calculate the local green's function,
!!     (4) calculate the hybridization function,
!!     (5) calculate the local weiss's function,
!!     (6) write the above calculated results.
!!
!! if task == 2:
!!     the code will do the following jobs:
!!     (1) search the fermi level (depends on `lfermi`),
!!     (2) calculate correction for density matrix,
!!     (3) write the above calculated results.
!!
!! if task == 3:
!!     search the fermi level only (depends on `lfermi`).
!!
!! if task == 4:
!!     calculate the impurity level only.
!!
!! if task == 5:
!!     calculate all complex frequency dependent eigenvalues only.
!!     we can use this feature to calculate the fermi surface.
!!
!! if task == 6:
!!     calculate lattice green's functions.
!!     we can use this feature to calculate the spectral functions.
!!
!! if task == 7:
!!     calculate density of states.
!!
!! if task == 8:
!!     calculate unknown physical properties.
!!
     integer, public, save :: task   = 1

!!
!! @var axis
!!
!! control flag, determine the working axis for brillouin zone integration
!! and fermi level search.
!!
!! if axis == 1:
!!     imaginary axis.
!!
!! if axis == 2:
!!     real axis.
!!
     integer, public, save :: axis   = 1

!!========================================================================
!!>>> logical variables: from dmft.in                                  <<<
!!========================================================================

!!
!! @var lfermi
!!
!! control flag, determine whether the fermi level should be updated.
!!
!! if lfermi == .true.
!!     search the fermi level by using the bisection algorithm.
!!
!! if lfermi == .false.
!!     fix the fermi level. in other words, the dft fermi level is used.
!!
     logical, public, save :: lfermi = .true.

!!
!! @var ltetra
!!
!! control flag, determine whether the analytical tetrahedron method is
!! used to perform the brillouin zone integration. this feature has not
!! been implemented so far.
!!
!! if ltetra == .true.
!!     use the analytical tetrahedron algorithm.
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
!! inverse temperature, \beta = 1 / T.
!!
     real(dp), public, save :: beta  = 8.00_dp

!!
!! @var mc
!!
!! convergence criterion for fermi level search.
!!
     real(dp), public, save :: mc    = 0.0001_dp

!!========================================================================
!!>>> integer variables: from params.ir                                <<<
!!========================================================================

!!
!! @var nsort
!!
!! number of atomic sorts in the model.
!!
     integer, public, save :: nsort  = 3

!!
!! @var natom
!!
!! number of atoms in the model.
!!
     integer, public, save :: natom  = 5

!!
!! @var nband
!!
!! number of dft bands.
!!
     integer, public, save :: nband  = 30

!!
!! @var nkpt
!!
!! number of k-mesh points.
!!
     integer, public, save :: nkpt   = 729

!!
!! @var nspin
!!
!! number of spin orientations. for non-magnetic or paramagnetic systems
!! nspin = 1, while for magnetic systems, nspin = 2. here, spin-orbit
!! coupling has not been supported.
!!
     integer, public, save :: nspin  = 1

!!
!! @var ntet
!!
!! number of tetrahedra. note that ntet = 1 means that the tetrahedron
!! data are absent and `ltetra` must be .false.
!!
     integer, public, save :: ntet   = 1

!!
!! @var ngrp
!!
!! number of groups of projectors, which are used to create the Hilbert
!! subspace for correlated or non-correlated orbitals. note that `ngrp`
!! is always larger or equal to `nsite`. in other words, multiple groups
!! of projectors are permitted. but some of them might be non-correlated.
!!
     integer, public, save :: ngrp   = 1

!!
!! @var nwnd
!!
!! number of energy windows or band windows, which are used to restrict
!! how many Kohn-Sham states are included in the calculations. please
!! notice that nwnd is always equal to ngrp. in other words, there is
!! one-to-one match between each group of projectors and each window of
!! Kohn-Sham states.
!!
     integer, public, save :: nwnd   = 1

!!
!! @var nsite
!!
!! number of correlated electron problems, i.e, number of impurity sites
!! in which the correlated effect is considered. `nsite` should be smaller
!! or equal to `ngrp`.
!!
     integer, public, save :: nsite  = 1

!!
!! @var nmesh
!!
!! number of frequency points.
!!
     integer, public, save :: nmesh  = 8193

!!========================================================================
!!>>> real variables: from params.ir                                   <<<
!!========================================================================

!!
!! @var scale
!!
!! an universal scaling factor for the lattice constants.
!!
     real(dp), public, save :: scale = 4.00_dp

!!
!! @var fermi
!!
!! default fermi level, which is usually taken from the dft calculations.
!! when task = 1, 2, or 3, `fermi` might be updated.
!!
     real(dp), public, save :: fermi = 0.00_dp

!!
!! @var volt
!!
!! volume of a tetrahedron, which is used to renormalize the integration
!! by using the analytical tetrahedron algorithm.
!!
     real(dp), public, save :: volt  = 1.00_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

!!
!! @var nprocs
!!
!! number of processors: default value 1.
!!
     integer, public, save :: nprocs = 1

!!
!! @var myid
!!
!! the id of current process: default value 0.
!!
     integer, public, save :: myid   = 0

!!
!! @var master
!!
!! denote as the controller process: default value 0.
!!
     integer, public, save :: master = 0

!!
!! @var cid
!!
!! the id of current process in cartesian topology (cid == myid).
!!
     integer, public, save :: cid    = 0

!!
!! @var cx
!!
!! the x coordinates of current process in cartesian topology.
!!
     integer, public, save :: cx     = 0

!!
!! @var cy
!!
!! the y coordinates of current process in cartesian topology.
!!
     integer, public, save :: cy     = 0

  end module control

!!========================================================================
!!>>> module version                                                   <<<
!!========================================================================

!!
!! @mod version
!!
!! define the semantic version string.
!!
  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL = 'v0.7.5 @ 2025.03.26D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.7.5'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2025.03.26'

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
     character(len=22), public, parameter :: V_MAIL = 'huangli@caep.cn'

!!
!! @var V_GPL3
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: V_GPL3 = 'GNU General Public License version 3'

  end module version
