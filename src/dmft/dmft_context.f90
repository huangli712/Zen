!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_map      module
!!!           dmft_group    module
!!!           dmft_window   module
!!!           dmft_lattice  module
!!!           dmft_kmesh    module
!!!           dmft_tetra    module
!!!           dmft_eigen    module
!!!           dmft_projs    module
!!!           dmft_fmesh    module
!!!           dmft_eimps    module
!!!           dmft_sigma    module
!!!           dmft_green    module
!!!           dmft_weiss    module
!!!           dmft_delta    module
!!!           dmft_gamma    module
!!!           context       module
!!! source  : dmft_context.f90
!!! type    : modules
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           07/30/2021 by li huang (last modified)
!!! purpose : try to define the global modules and arrays, and implement
!!!           memory managment.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module dmft_map                                                  <<<
!!========================================================================

!!
!! @mod dmft_map
!!
!! define connections / mappings between the quantum impurity problems and
!! the groups of projectors (and band windows).
!!
  module dmft_map
     implicit none

!!
!! @var i_grp
!!
!! given a given quantum impurity problem, return the corresponding group
!! of projectors, impurity -> group.
!!
     integer, public, save, allocatable :: i_grp(:)

!!
!! @var i_wnd
!!
!! given a given quantum impurity problem, return the corresponding dft
!! band window, impurity -> window.
!!
     integer, public, save, allocatable :: i_wnd(:)

!!
!! @var g_imp
!!
!! given a given group of projectors, return the corresponding quantum
!! impurity problem, group -> impurity.
!!
!! 0 value means that this group of projectors is just for non-correlated
!! orbitals.
!!
     integer, public, save, allocatable :: g_imp(:)

!!
!! @var w_imp
!!
!! given a given dft band window, return the corresponding quantum impurity
!! problem, window -> impurity.
!!
!! 0 value means that this dft band window is for non-correlated orbitals.
!!
     integer, public, save, allocatable :: w_imp(:)

  end module dmft_map

!!========================================================================
!!>>> module dmft_group                                                <<<
!!========================================================================

!!
!! @mod dmft_group
!!
!! specify the traits of groups of projectors.
!!
  module dmft_group
     implicit none

!!
!! @var qdim
!!
!! maximum number of correlated orbitals for all groups. actually, it
!! should be equal to maxval(ndim).
!!
     integer, public, save :: qdim = -1

!!
!! @var shell
!!
!! specification of orbital shell.
!!
     character(len=5), public, save, allocatable :: shell(:)

!!
!! @var corr
!!
!! tell us this group is correlated or not.
!!
     logical, public, save, allocatable :: corr(:)

!!
!! @var site
!!
!! the corresponding atomic site of this group.
!!
     integer, public, save, allocatable :: site(:)

!!
!! @var l
!!
!! the corresponding angular momentum quantum number of this group.
!!
     integer, public, save, allocatable :: l(:)

!!
!! @var ndim
!!
!! number of projectors (orbitals) included in this group.
!!
     integer, public, save, allocatable :: ndim(:)

  end module dmft_group

!!========================================================================
!!>>> module dmft_window                                               <<<
!!========================================================================

!!
!! @mod dmft_window
!!
!! specify the dft band windows for projectors.
!!
  module dmft_window
     implicit none

!!
!! @var qbnd
!!
!! maximum number of dft bands for all the band windows. actually, it
!! should be equal to maxval(nbnd).
!!
     integer, public, save :: qbnd = -1

!!
!! @var bmin
!!
!! lower boundaries for the band windows.
!!
     integer, public, save, allocatable :: bmin(:)

!!
!! @var bmax
!!
!! upper boundaries for the band windows.
!!
     integer, public, save, allocatable :: bmax(:)

!!
!! @var nbnd
!!
!! number of bands for the band windows.
!!
     integer, public, save, allocatable :: nbnd(:)

!!
!! @var kwin
!!
!! momentum- and spin-dependent band windows.
!!
     integer, public, save, allocatable :: kwin(:,:,:,:)

  end module dmft_window

!!========================================================================
!!>>> module dmft_lattice                                              <<<
!!========================================================================

!!
!! @mod dmft_lattice
!!
!! contain crystallography information (i.e. crystal structures) for the
!! strongly correlated materials.
!!
  module dmft_lattice
     use constants, only : dp

     implicit none

!!
!! @var sorts
!!
!! sorts of atoms (chemical symbol).
!!
     character(len=2), public, save, allocatable :: sorts(:)

!!
!! @var atoms
!!
!! list of atoms (chemical symbol).
!!
     character(len=2), public, save, allocatable :: atoms(:)

!!
!! @var sortn
!!
!! number of atoms for each sort.
!!
     integer,  public, save, allocatable :: sortn(:)

!!
!! @var lvect
!!
!! three lattice vectors.
!!
     real(dp), public, save, allocatable :: lvect(:,:)

!!
!! @var coord
!!
!! atomic positions.
!!
     real(dp), public, save, allocatable :: coord(:,:)

  end module dmft_lattice

!!========================================================================
!!>>> module dmft_kmesh                                                <<<
!!========================================================================

!!
!! @mod dmft_kmesh
!!
!! contain the k-mesh and the corresponding integration weights.
!!
  module dmft_kmesh
     use constants, only : dp

     implicit none

!!
!! @var kmesh
!!
!! k-mesh in the brillouin zone.
!!
     real(dp), public, save, allocatable :: kmesh(:,:)

!!
!! @var weight
!!
!! integration weights for k-points. please pay attention to that they
!! have not been renormalized. we have to ensure sum(weight) = nkpt.
!!
     real(dp), public, save, allocatable :: weight(:)

  end module dmft_kmesh

!!========================================================================
!!>>> module dmft_tetra                                                <<<
!!========================================================================

!!
!! @mod dmft_tetra
!!
!! contain the tetrahedron information, which is used to carry out the
!! brillouin zone integration. this feature has not been implemented.
!!
  module dmft_tetra
     implicit none

!!
!! @var tetra
!!
!! contain tetrahedron information. tetra(itet,1:4) point to the four
!! vertices for given tetrahedron, while tetra(itet,5) denotes the
!! corresponding weight for this tetrahedron.
!!
     integer, public, save, allocatable :: tetra(:,:)

  end module dmft_tetra

!!========================================================================
!!>>> module dmft_eigen                                                <<<
!!========================================================================

!!
!! @mod dmft_eigen
!!
!! contain the Kohn-Sham eigenvalues and related occupations.
!!
  module dmft_eigen
     use constants, only : dp

     implicit none

!!
!! @var enk
!!
!! eigenvalues in the Kohn-Sham basis.
!!
     real(dp), public, save, allocatable :: enk(:,:,:)

!!
!! @var occupy
!!
!! occupations in the Kohn-Sham basis.
!!
     real(dp), public, save, allocatable :: occupy(:,:,:)

  end module dmft_eigen

!!========================================================================
!!>>> module dmft_projs                                                <<<
!!========================================================================

!!
!! @mod dmft_projs
!!
!! contain the local orbital projectors.
!!
  module dmft_projs
     use constants, only : dp

     implicit none

!!
!! @var chipsi
!!
!! overlap matrix between the local orbitals and the Kohn-Sham basis. its
!! definition is \langle \chi^{I}_{\alpha,k} | \psi_{b,k} \rangle, where
!! `I` means the index for correlated sites, \alpha means the index for
!! correlated orbitals. `b`, `k`, `s` are indices for dft bands, k-points,
!! and spins, respectively. of course, `b` is restricted by band windows,
!! and depends on k-points as well. so, chipsi is a rank-5 array.
!!
     complex(dp), public, save, allocatable :: chipsi(:,:,:,:,:)

!!
!! @var psichi
!!
!! overlap matrix between the Kohn-Sham basis and the local orbitals. its
!! definition is \langle \psi_{b,k} | \chi^{I}_{\alpha,k} \rangle.
!! actually, psichi can be obtained by chipsi via conjugate transpose.
!!
     complex(dp), public, save, allocatable :: psichi(:,:,:,:,:)

  end module dmft_projs

!!========================================================================
!!>>> module dmft_fmesh                                                <<<
!!========================================================================

!!
!! @mod dmft_fmesh
!!
!! contain the linear frequency mesh.
!!
  module dmft_fmesh
     use constants, only : dp

     implicit none

!!
!! @var fmesh
!!
!! linear frequency mesh. it can be defined on imaginary axis or real axis.
!!
     real(dp), public, save, allocatable :: fmesh(:)

  end module dmft_fmesh

!!========================================================================
!!>>> module dmft_eimps                                                <<<
!!========================================================================

!!
!! @mod dmft_eimps
!!
!! contain the local impurity levels.
!!
  module dmft_eimps
     use constants, only : dp

     implicit none

!!
!! @var eimps
!!
!! local impurity levels. eimps = \sum_k enk - mu.
!!
     complex(dp), public, save, allocatable :: eimps(:,:,:,:)

!!
!! @var eimpx
!!
!! local impurity levels shifted by double counting terms.
!! in other words, eimpx = eimps - sigdc.
!!
     complex(dp), public, save, allocatable :: eimpx(:,:,:,:)

  end module dmft_eimps

!!========================================================================
!!>>> module dmft_sigma                                                <<<
!!========================================================================

!!
!! @mod dmft_sigma
!!
!! contain the impurity self-energy functions.
!!
  module dmft_sigma
     use constants, only : dp

     implicit none

!!
!! @var sigdc
!!
!! dobule counting term for self-energy functions, which are determined
!! by the Zen framework. this code will read it from file sigma.dc.
!!
     complex(dp), public, save, allocatable :: sigdc(:,:,:,:)

!!
!! @var sigoo
!!
!! asymptotic values for bare self-energy functions (when \omega goes
!! to \infty). note that the double counting terms should be subtracted
!! from them.
!!
     complex(dp), public, save, allocatable :: sigoo(:,:,:,:)

!!
!! @var sigma
!!
!! impurity self-energy functions. they are usually taken from the output
!! of various quantum impurity solver. this code will read them from file
!! sigma.bare. note that the double counting terms should be subtracted
!! from them. in addition, sigma is frequency-dependent, but sigoo and
!! sigdc are not.
!!
     complex(dp), public, save, allocatable :: sigma(:,:,:,:,:)

  end module dmft_sigma

!!========================================================================
!!>>> module dmft_green                                                <<<
!!========================================================================

!!
!! @mod dmft_green
!!
!! contain the local green's functions.
!!
  module dmft_green
     use constants, only : dp

     implicit none

!!
!! @var green
!!
!! local green's functions. note that within the dynamical mean-field
!! theory, local green's functions should be equal to impurity green's
!! functions during self-consistent iterations.
!!
     complex(dp), public, save, allocatable :: green(:,:,:,:,:)

  end module dmft_green

!!========================================================================
!!>>> module dmft_weiss                                                <<<
!!========================================================================

!!
!! @mod dmft_weiss
!!
!! contain local weiss's functions.
!!
  module dmft_weiss
     use constants, only : dp

     implicit none

!!
!! @var weiss
!!
!! local weiss's functions.
!!
     complex(dp), public, save, allocatable :: weiss(:,:,:,:,:)

  end module dmft_weiss

!!========================================================================
!!>>> module dmft_delta                                                <<<
!!========================================================================

!!
!! @mod dmft_weiss
!!
!! contain local hybridization functions.
!!
  module dmft_delta
     use constants, only : dp

     implicit none

!!
!! @var delta
!!
!! local hybridization functions. the ct-qmc quantum impurity solver will
!! need them as input.
!!
     complex(dp), public, save, allocatable :: delta(:,:,:,:,:)

  end module dmft_delta

!!========================================================================
!!>>> module dmft_gamma                                                <<<
!!========================================================================

!!
!! @mod dmft_gamma
!!
!! contain dft + dmft correction for density matrix.
!!
  module dmft_gamma
     use constants, only : dp

     implicit none

!!
!! @var gamma
!!
!! dft + dmft correction for density matrix, which is used to carry out
!! fully charge self-consistent dft + dmft calculations.
!!
     complex(dp), public, save, allocatable :: gamma(:,:,:,:)

  end module dmft_gamma

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

!!
!! @mod context
!!
!! containing memory management subroutines, which initialize or destroy
!! all of the global variables and arrays.
!!
  module context
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : nsort, natom
     use control, only : nband, nkpt, nspin
     use control, only : ntet
     use control, only : ngrp, nwnd
     use control, only : nsite
     use control, only : nmesh

     use dmft_map
     use dmft_group
     use dmft_window
     use dmft_lattice
     use dmft_kmesh
     use dmft_tetra
     use dmft_eigen
     use dmft_projs
     use dmft_fmesh
     use dmft_eimps
     use dmft_sigma
     use dmft_green
     use dmft_weiss
     use dmft_delta
     use dmft_gamma

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
     public :: cat_alloc_map
     public :: cat_alloc_group
     public :: cat_alloc_window
     public :: cat_alloc_lattice
     public :: cat_alloc_kmesh
     public :: cat_alloc_tetra
     public :: cat_alloc_eigen
     public :: cat_alloc_projs
     public :: cat_alloc_fmesh
     public :: cat_alloc_eimps
     public :: cat_alloc_sigma
     public :: cat_alloc_green
     public :: cat_alloc_weiss
     public :: cat_alloc_delta
     public :: cat_alloc_gamma

     ! declaration of module procedures: deallocate memory
     public :: cat_free_map
     public :: cat_free_group
     public :: cat_free_window
     public :: cat_free_lattice
     public :: cat_free_kmesh
     public :: cat_free_tetra
     public :: cat_free_eigen
     public :: cat_free_projs
     public :: cat_free_fmesh
     public :: cat_free_eimps
     public :: cat_free_sigma
     public :: cat_free_green
     public :: cat_free_weiss
     public :: cat_free_delta
     public :: cat_free_gamma

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_map
!!
!! allocate memory for map-related variables.
!!
  subroutine cat_alloc_map()
     implicit none

!! [body

     ! allocate memory
     allocate(i_grp(nsite), stat = istat)
     allocate(i_wnd(nsite), stat = istat)
     allocate(g_imp(ngrp),  stat = istat)
     allocate(w_imp(nwnd),  stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_map','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     i_grp = 0
     i_wnd = 0
     g_imp = 0
     w_imp = 0

!! body]

     return
  end subroutine cat_alloc_map

!!
!! @sub cat_alloc_group
!!
!! allocate memory for group-related variables.
!!
  subroutine cat_alloc_group()
     implicit none

!! [body

     ! allocate memory
     allocate(shell(ngrp), stat = istat)
     allocate(corr(ngrp),  stat = istat)
     allocate(site(ngrp),  stat = istat)
     allocate(l(ngrp),     stat = istat)
     allocate(ndim(ngrp),  stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_group','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     shell = 's'
     corr  = .false.
     site  = 0
     l     = 0
     ndim  = 0

     ! special treatment for `qdim`.
     ! qdim should be initialized in dmft_setup_param().
     if ( qdim < 0 ) then
         call s_print_error('cat_alloc_group','qdim is less than 0')
     endif ! back if ( qdim < 0 ) block

!! body]

     return
  end subroutine cat_alloc_group

!!
!! @sub cat_alloc_window
!!
!! allocate memory for window-related variables.
!!
  subroutine cat_alloc_window()
     implicit none

!! [body

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

     ! special treatment for `qbnd`.
     ! qbnd should be initialized in dmft_setup_param().
     if ( qbnd < 0 ) then
         call s_print_error('cat_alloc_window','qbnd is less than 0')
     endif ! back if ( qbnd < 0 ) block

!! body]

     return
  end subroutine cat_alloc_window

!!
!! @sub cat_alloc_lattice
!!
!! allocate memory for lattice-related variables.
!!
  subroutine cat_alloc_lattice()
     implicit none

!! [body

     ! allocate memory
     allocate(sorts(nsort),   stat = istat)
     allocate(atoms(natom),   stat = istat)
     allocate(sortn(nsort),   stat = istat)
     allocate(lvect(3,3),     stat = istat)
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

!! body]

     return
  end subroutine cat_alloc_lattice

!!
!! @sub cat_alloc_kmesh
!!
!! allocate memory for kmesh-related variables.
!!
  subroutine cat_alloc_kmesh()
     implicit none

!! [body

     ! allocate memory
     allocate(kmesh(nkpt,3), stat = istat)
     allocate(weight(nkpt),  stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_kmesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     kmesh  = zero
     weight = zero

!! body]

     return
  end subroutine cat_alloc_kmesh

!!
!! @sub cat_alloc_tetra
!!
!! allocate memory for tetra-related variables.
!!
  subroutine cat_alloc_tetra()
     implicit none

!! [body

     ! allocate memory
     allocate(tetra(ntet,5), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_tetra','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     tetra = 0

!! body]

     return
  end subroutine cat_alloc_tetra

!!
!! @sub cat_alloc_eigen
!!
!! allocate memory for eigen-related variables.
!!
  subroutine cat_alloc_eigen()
     implicit none

!! [body

     ! allocate memory
     allocate(enk(nband,nkpt,nspin),    stat = istat)
     allocate(occupy(nband,nkpt,nspin), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_eigen','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     enk    = zero
     occupy = zero

!! body]

     return
  end subroutine cat_alloc_eigen

!!
!! @sub cat_alloc_projs
!!
!! allocate memory for projs-related variables.
!!
  subroutine cat_alloc_projs()
     implicit none

!! [body

     ! allocate memory
     allocate(chipsi(qdim,qbnd,nkpt,nspin,ngrp), stat = istat)
     allocate(psichi(qbnd,qdim,nkpt,nspin,ngrp), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_projs','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     chipsi = czero
     psichi = czero

!! body]

     return
  end subroutine cat_alloc_projs

!!
!! @sub cat_alloc_fmesh
!!
!! allocate memory for fmesh-related variables.
!!
  subroutine cat_alloc_fmesh()
     implicit none

!! [body

     ! allocate memory
     allocate(fmesh(nmesh), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fmesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     fmesh = zero

!! body]

     return
  end subroutine cat_alloc_fmesh

!!
!! @sub cat_alloc_eimps
!!
!! allocate memory for eimps-related variables.
!!
  subroutine cat_alloc_eimps()
     implicit none

!! [body

     ! allocate memory
     allocate(eimps(qdim,qdim,nspin,nsite), stat = istat)
     allocate(eimpx(qdim,qdim,nspin,nsite), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     eimps = czero
     eimpx = czero

!! body]

     return
  end subroutine cat_alloc_eimps

!!
!! @sub cat_alloc_sigma
!!
!! allocate memory for sigma-related variables.
!!
  subroutine cat_alloc_sigma()
     implicit none

!! [body

     ! allocate memory
     allocate(sigdc(qdim,qdim,nspin,nsite),       stat = istat)
     allocate(sigoo(qdim,qdim,nspin,nsite),       stat = istat)
     allocate(sigma(qdim,qdim,nmesh,nspin,nsite), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_sigma','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     sigdc = czero
     sigoo = czero
     sigma = czero

!! body]

     return
  end subroutine cat_alloc_sigma

!!
!! @sub cat_alloc_green
!!
!! allocate memory for green-related variables
!!
  subroutine cat_alloc_green()
     implicit none

!! [body

     ! allocate memory
     allocate(green(qdim,qdim,nmesh,nspin,nsite), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     green = czero

!! body]

     return
  end subroutine cat_alloc_green

!!
!! @sub cat_alloc_weiss
!!
!! allocate memory for weiss-related variables.
!!
  subroutine cat_alloc_weiss()
     implicit none

!! [body

     ! allocate memory
     allocate(weiss(qdim,qdim,nmesh,nspin,nsite), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_weiss','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     weiss = czero

!! body]

     return
  end subroutine cat_alloc_weiss

!!
!! @sub cat_alloc_delta
!!
!! allocate memory for delta-related variables.
!!
  subroutine cat_alloc_delta()
     implicit none

!! [body

     ! allocate memory
     allocate(delta(qdim,qdim,nmesh,nspin,nsite), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_delta','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     delta = czero

!! body]

     return
  end subroutine cat_alloc_delta

!!
!! @sub cat_alloc_gamma
!!
!! allocate memory for gamma-related variables.
!!
  subroutine cat_alloc_gamma()
     implicit none

!! [body

     ! allocate memory
     allocate(gamma(qbnd,qbnd,nkpt,nspin), stat = istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_gamma','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     gamma = czero

!! body]

     return
  end subroutine cat_alloc_gamma

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_map
!!
!! deallocate memory for map-related variables.
!!
  subroutine cat_free_map()
     implicit none

!! [body

     if ( allocated(i_grp) ) deallocate(i_grp)
     if ( allocated(i_wnd) ) deallocate(i_wnd)
     if ( allocated(g_imp) ) deallocate(g_imp)
     if ( allocated(w_imp) ) deallocate(w_imp)

!! body]

     return
  end subroutine cat_free_map

!!
!! @sub cat_free_group
!!
!! deallocate memory for group-related variables.
!!
  subroutine cat_free_group()
     implicit none

!! [body

     if ( allocated(shell) ) deallocate(shell)
     if ( allocated(corr)  ) deallocate(corr )
     if ( allocated(site)  ) deallocate(site )
     if ( allocated(l)     ) deallocate(l    )
     if ( allocated(ndim)  ) deallocate(ndim )

!! body]

     return
  end subroutine cat_free_group

!!
!! @sub cat_free_window
!!
!! deallocate memory for window-related variables.
!!
  subroutine cat_free_window()
     implicit none

!! [body

     if ( allocated(bmin) ) deallocate(bmin)
     if ( allocated(bmax) ) deallocate(bmax)
     if ( allocated(nbnd) ) deallocate(nbnd)
     if ( allocated(kwin) ) deallocate(kwin)

!! body]

     return
  end subroutine cat_free_window

!!
!! @sub cat_free_lattice
!!
!! deallocate memory for lattice-related variables.
!!
  subroutine cat_free_lattice()
     implicit none

!! [body

     if ( allocated(sorts) ) deallocate(sorts)
     if ( allocated(atoms) ) deallocate(atoms)
     if ( allocated(sortn) ) deallocate(sortn)
     if ( allocated(lvect) ) deallocate(lvect)
     if ( allocated(coord) ) deallocate(coord)

!! body]

     return
  end subroutine cat_free_lattice

!!
!! @sub cat_free_kmesh
!!
!! deallocate memory for kmesh-related variables.
!!
  subroutine cat_free_kmesh()
     implicit none

!! [body

     if ( allocated(kmesh)  ) deallocate(kmesh )
     if ( allocated(weight) ) deallocate(weight)

!! body]

     return
  end subroutine cat_free_kmesh

!!
!! @sub cat_free_tetra
!!
!! deallocate memory for tetra-related variables.
!!
  subroutine cat_free_tetra()
     implicit none

!! [body

     if ( allocated(tetra) ) deallocate(tetra)

!! body]

     return
  end subroutine cat_free_tetra

!!
!! @sub cat_free_eigen
!!
!! deallocate memory for eigen-related variables.
!!
  subroutine cat_free_eigen()
     implicit none

!! [body

     if ( allocated(enk)    ) deallocate(enk   )
     if ( allocated(occupy) ) deallocate(occupy)

!! body]

     return
  end subroutine cat_free_eigen

!!
!! @sub cat_free_projs
!!
!! deallocate memory for projs-related variables.
!!
  subroutine cat_free_projs()
     implicit none

!! [body

     if ( allocated(chipsi) ) deallocate(chipsi)
     if ( allocated(psichi) ) deallocate(psichi)

!! body]

     return
  end subroutine cat_free_projs

!!
!! @sub cat_free_fmesh
!!
!! deallocate memory for fmesh-related variables.
!!
  subroutine cat_free_fmesh()
     implicit none

!! [body

     if ( allocated(fmesh) ) deallocate(fmesh)

!! body]

     return
  end subroutine cat_free_fmesh

!!
!! @sub cat_free_eimps
!!
!! deallocate memory for eimps-related variables.
!!
  subroutine cat_free_eimps()
     implicit none

!! [body

     if ( allocated(eimps) ) deallocate(eimps)
     if ( allocated(eimpx) ) deallocate(eimpx)

!! body]

     return
  end subroutine cat_free_eimps

!!
!! @sub cat_free_sigma
!!
!! deallocate memory for sigma-related variables.
!!
  subroutine cat_free_sigma()
     implicit none

!! [body

     if ( allocated(sigdc) ) deallocate(sigdc)
     if ( allocated(sigoo) ) deallocate(sigoo)
     if ( allocated(sigma) ) deallocate(sigma)

!! body]

     return
  end subroutine cat_free_sigma

!!
!! @sub cat_free_green
!!
!! deallocate memory for green-related variables.
!!
  subroutine cat_free_green()
     implicit none

!! [body

     if ( allocated(green) ) deallocate(green)

!! body]

     return
  end subroutine cat_free_green

!!
!! @sub cat_free_weiss
!!
!! deallocate memory for weiss-related variables.
!!
  subroutine cat_free_weiss()
     implicit none

!! [body

     if ( allocated(weiss) ) deallocate(weiss)

!! body]

     return
  end subroutine cat_free_weiss

!!
!! @sub cat_free_delta
!!
!! deallocate memory for delta-related variables.
!!
  subroutine cat_free_delta()
     implicit none

!! [body

     if ( allocated(delta) ) deallocate(delta)

!! body]

     return
  end subroutine cat_free_delta

!!
!! @sub cat_free_gamma
!!
!! deallocate memory for gamma-related variables.
!!
  subroutine cat_free_gamma()
     implicit none

!! [body

     if ( allocated(gamma) ) deallocate(gamma)

!! body]

     return
  end subroutine cat_free_gamma

  end module context
