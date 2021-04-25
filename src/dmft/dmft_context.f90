!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_map
!!!           dmft_group
!!!           dmft_window
!!!           dmft_lattice
!!!           dmft_kmesh
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
!!!           04/25/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module dmft_map                                                  <<<
!!========================================================================

!!
!! @mod dmft_map
!!
!! define connections /mappings between the quantum impurity problems and
!! the groups of projectors
!!
  module dmft_map
     implicit none

!!
!! @var i_grp
!!
!! from a given quantum impurity problem, return the corresponding group
!! of projectors
!!
     integer, public, save, allocatable :: i_grp(:)

!!
!! @var g_imp
!!
!! from a given group of projectors, return the corresponding quantum
!! impurity problem
!!
     integer, public, save, allocatable :: g_imp(:)

  end module dmft_map

!!========================================================================
!!>>> module dmft_group                                                <<<
!!========================================================================

!!
!! @mod dmft_group
!!
!! specify the traits of groups of projectors
!!
  module dmft_group
     implicit none

!!
!! @var qdim
!!
!! maximum number of correlated orbitals in all groups
!!
     integer, public, save :: qdim = -1

!!
!! @var shell
!!
!! specification of orbital shell
!!
     character(len=4), public, save, allocatable :: shell(:)

!!
!! @var corr
!!
!! test which group is correlated
!!
     logical, public, save, allocatable :: corr(:)

!!
!! @var site
!!
!! the corresponding atomic site of group
!!
     integer, public, save, allocatable :: site(:)

!!
!! @var l
!!
!! the corresponding angular momentum quantum number of group
!!
     integer, public, save, allocatable :: l(:)

!!
!! @var ndim
!!
!! number of projectors (orbitals) of group
!!
     integer, public, save, allocatable :: ndim(:)

  end module dmft_group

!!========================================================================
!!>>> module dmft_window                                               <<<
!!========================================================================

!!
!! @mod dmft_window
!!
!! specify the band windows of groups of projectors
!!
  module dmft_window
     implicit none

!!
!! @var qbnd
!!
!! maximum number of bands for the band windows
!!
     integer, public, save :: qbnd = -1

!!
!! @var bmin
!!
!! lower boundaries for the band windows
!!
     integer, public, save, allocatable :: bmin(:)

!!
!! @var bmax
!!
!! upper boundaries for the band windows
!!
     integer, public, save, allocatable :: bmax(:)

!!
!! @var nbnd
!!
!! number of bands for the band windows
!!
     integer, public, save, allocatable :: nbnd(:)

!!
!! @var kwin
!!
!! momentum- and spin-dependent band windows
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
!! strongly correlated materials
!!
  module dmft_lattice
     use constants, only : dp

     implicit none

!!
!! @var sorts
!!
!! sorts of atoms
!!
     character(len=2), public, save, allocatable :: sorts(:)

!!
!! @var atoms
!!
!! list of atoms
!!
     character(len=2), public, save, allocatable :: atoms(:)

!!
!! @var sortn
!!
!! number of atoms of sorts
!!
     integer,  public, save, allocatable :: sortn(:)

!!
!! @var lvect
!!
!! three lattice vectors
!!
     real(dp), public, save, allocatable :: lvect(:,:)

!!
!! @var coord
!!
!! atomic positions
!!
     real(dp), public, save, allocatable :: coord(:,:)

  end module dmft_lattice

!!========================================================================
!!>>> module dmft_kmesh                                                <<<
!!========================================================================

!!
!! @mod dmft_kmesh
!!
!! contain the k-mesh and the corresponding integration weights  
!!
  module dmft_kmesh
     use constants, only : dp

     implicit none

!!
!! @var kmesh
!!
!! k-mesh in the brillouin zone
!!
     real(dp), public, save, allocatable :: kmesh(:,:)

!!
!! @var weight
!!
!! integration weights for k-points
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
!! brillouin zone integration
!!
  module dmft_tetra
     implicit none

!!
!! @var tetra
!!
!! contain tetrahedron information
!!
     integer, public, save, allocatable :: tetra(:,:)

  end module dmft_tetra

!!========================================================================
!!>>> module dmft_eigen                                                <<<
!!========================================================================

!!
!! @mod dmft_eigen
!!
!! contain the Kohn-Sham eigenvalues and occupations
!!
  module dmft_eigen
     use constants, only : dp

     implicit none

!!
!! @var enk
!!
!! eigenvalues in the Kohn-Sham basis
!!
     real(dp), public, save, allocatable :: enk(:,:,:)

!!
!! @var occupy
!!
!! occupations in the Kohn-Sham basis
!!
     real(dp), public, save, allocatable :: occupy(:,:,:)

  end module dmft_eigen

!!========================================================================
!!>>> module dmft_projs                                                <<<
!!========================================================================

!!
!! @mod dmft_projs
!!
!! contain the local orbital projections
!!
  module dmft_projs
     use constants, only : dp

     implicit none

!!
!! @var psichi
!!
!! overlap matrix between the Kohn-Sham basis and the local orbitals
!!
     complex(dp), public, save, allocatable :: psichi(:,:,:,:,:)

  end module dmft_projs

!!========================================================================
!!>>> module dmft_fmesh                                                <<<
!!========================================================================

!!
!! @mod dmft_fmesh
!!
!! contain the frequency mesh
!!
  module dmft_fmesh
     use constants, only : dp

     implicit none

!!
!! @var fmesh
!!
!! frequency mesh. it can be defined on imaginary axis or real axis
!!
     real(dp), public, save, allocatable :: fmesh(:)

  end module dmft_fmesh

!!========================================================================
!!>>> module dmft_sigma                                                <<<
!!========================================================================

!!
!! @mod dmft_sigma
!!
!! contain the self-energy functions
!!
  module dmft_sigma
     use constants, only : dp

     implicit none

!!
!! @var sigdc
!!
!! dobule counting term for self-energy functions
!!
     complex(dp), public, save, allocatable :: sigdc(:,:,:,:)

!!
!! @var sig_l
!!
!! local self-energy functions. they are usually taken from the output of
!! various quantum impurity solver
!!
     complex(dp), public, save, allocatable :: sig_l(:,:,:,:,:)

!!
!! @var sig_k
!!
!! self-energy functions embedded in k-space
!!
     complex(dp), public, save, allocatable :: sig_k(:,:,:,:,:)

  end module dmft_sigma

!!========================================================================
!!>>> module dmft_green                                                <<<
!!========================================================================

!!
!! @mod dmft_green
!!
!! contain the green's functions
!!
  module dmft_green
     use constants, only : dp

     implicit none

!!
!! @var grn_l
!!
!! local green's functions
!!
     complex(dp), public, save, allocatable :: grn_l(:,:,:,:,:)

!!
!! @var grn_k
!!
!! lattice green's functions
!!
     complex(dp), public, save, allocatable :: grn_k(:,:,:,:,:)

  end module dmft_green

!!========================================================================
!!>>> module dmft_weiss                                                <<<
!!========================================================================

!!
!! @mod dmft_weiss
!!
!! contain bath weiss functions and hybridization functions
!!
  module dmft_weiss
     use constants, only : dp

     implicit none

!!
!! @var wss_l
!!
!! local bath weiss functions
!!
     complex(dp), public, save, allocatable :: wss_l(:,:,:,:,:)

!!
!! @var wss_k
!!
!! lattice bath weiss functions
!!
     complex(dp), public, save, allocatable :: wss_k(:,:,:,:,:)

!!
!! @var hyb_l
!!
!! local hybridization functions
!!
     complex(dp), public, save, allocatable :: hyb_l(:,:,:,:,:)

!!
!! @var hyb_k
!!
!! lattice hybridization functions
!!
     complex(dp), public, save, allocatable :: hyb_k(:,:,:,:,:)

  end module dmft_weiss

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

!!
!! @mod context
!!
!! containing memory management subroutines, which initialize all of the
!! global variables and arrays
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
     use dmft_sigma
     use dmft_green
     use dmft_weiss

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
     public :: cat_alloc_sigma
     public :: cat_alloc_green
     public :: cat_alloc_weiss

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
     public :: cat_free_sigma
     public :: cat_free_green
     public :: cat_free_weiss

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_map
!!
!! allocate memory for map-related variables
!!
  subroutine cat_alloc_map()
     implicit none

! allocate memory
     allocate(i_grp(nsite), stat = istat)
     allocate(g_imp(ngrp),  stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_map','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     i_grp = 0
     g_imp = 0

     return
  end subroutine cat_alloc_map

!!
!! @sub cat_alloc_group
!!
!! allocate memory for group-related variables
!!
  subroutine cat_alloc_group()
     implicit none

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

! special treatment for qdim
! qdim should be initialized in dmft_setup_param() 
     if ( qdim < 0 ) then
         call s_print_error('cat_alloc_group','qdim is less than 0')
     endif ! back if ( qdim < 0 ) block

     return
  end subroutine cat_alloc_group

!!
!! @sub cat_alloc_window
!!
!! allocate memory for window-related variables
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

! special treatment for max_nbnd
! max_nbnd should be initialized in dmft_setup_param() 
     if ( max_nbnd < 0 ) then
         call s_print_error('cat_alloc_window','max_nbnd is less than 0')
     endif ! back if ( max_nbnd < 0 ) block

     return
  end subroutine cat_alloc_window

!!
!! @sub cat_alloc_lattice
!!
!! allocate memory for lattice-related variables
!!
  subroutine cat_alloc_lattice()
     implicit none

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

     return
  end subroutine cat_alloc_lattice

!!
!! @sub cat_alloc_kmesh
!!
!! allocate memory for kmesh-related variables
!!
  subroutine cat_alloc_kmesh()
     implicit none

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

     return
  end subroutine cat_alloc_kmesh

!!
!! @sub cat_alloc_tetra
!!
!! allocate memory for tetra-related variables
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
!! allocate memory for eigen-related variables
!!
  subroutine cat_alloc_eigen()
     implicit none

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

     return
  end subroutine cat_alloc_eigen

!!
!! @sub cat_alloc_projs
!!
!! allocate memory for projs-related variables
!!
  subroutine cat_alloc_projs()
     implicit none

! allocate memory
     allocate(psichi(qdim,max_nbnd,nkpt,nspin,ngrp), stat = istat)

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
!! allocate memory for fmesh-related variables
!!
  subroutine cat_alloc_fmesh()
     implicit none

! allocate memory
     allocate(fmesh(nmesh), stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fmesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     fmesh = zero

     return
  end subroutine cat_alloc_fmesh

!!
!! @sub cat_alloc_sigma
!!
!! allocate memory for sigma-related variables
!!
  subroutine cat_alloc_sigma()
     implicit none

! allocate memory
     allocate(sigdc(qdim,qdim,nspin,nsite),       stat = istat)
     allocate(sig_l(nmesh,qdim,qdim,nspin,nsite), stat = istat)
     allocate(sig_k(nmesh,max_nbnd,max_nbnd,nkpt,nspin),  stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_sigma','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sigdc = czero
     sig_l = czero
     sig_k = czero

     return
  end subroutine cat_alloc_sigma

!!
!! @sub cat_alloc_green
!!
!! allocate memory for green-related variables
!!
  subroutine cat_alloc_green()
     implicit none

! allocate memory
     allocate(grn_l(nmesh,qdim,qdim,nspin,nsite), stat = istat)
     allocate(grn_k(nmesh,max_nbnd,max_nbnd,nkpt,nspin),  stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     grn_l = czero
     grn_k = czero

     return
  end subroutine cat_alloc_green

!!
!! @sub cat_alloc_weiss
!!
!! allocate memory for weiss-related variables
!!
  subroutine cat_alloc_weiss()
     implicit none

! allocate memory
     allocate(wss_l(nmesh,qdim,qdim,nspin,nsite), stat = istat)
     allocate(wss_k(nmesh,max_nbnd,max_nbnd,nkpt,nspin),  stat = istat)
     allocate(hyb_l(nmesh,qdim,qdim,nspin,nsite), stat = istat)
     allocate(hyb_k(nmesh,max_nbnd,max_nbnd,nkpt,nspin),  stat = istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_weiss','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     wss_l = czero
     wss_k = czero
     hyb_l = czero
     hyb_k = czero

     return
  end subroutine cat_alloc_weiss

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_map
!!
!! deallocate memory for map-related variables
!!
  subroutine cat_free_map()
     implicit none

     if ( allocated(i_grp) ) deallocate(i_grp)
     if ( allocated(g_imp) ) deallocate(g_imp)

     return
  end subroutine cat_free_map

!!
!! @sub cat_free_group
!!
!! deallocate memory for group-related variables
!!
  subroutine cat_free_group()
     implicit none

     if ( allocated(shell) ) deallocate(shell)
     if ( allocated(corr)  ) deallocate(corr )
     if ( allocated(site)  ) deallocate(site )
     if ( allocated(l)     ) deallocate(l    )
     if ( allocated(ndim)  ) deallocate(ndim )

     return
  end subroutine cat_free_group

!!
!! @sub cat_free_window
!!
!! deallocate memory for window-related variables
!!
  subroutine cat_free_window()
     implicit none

     if ( allocated(bmin) ) deallocate(bmin)
     if ( allocated(bmax) ) deallocate(bmax)
     if ( allocated(nbnd) ) deallocate(nbnd)
     if ( allocated(kwin) ) deallocate(kwin)

     return
  end subroutine cat_free_window

!!
!! @sub cat_free_lattice
!!
!! deallocate memory for lattice-related variables
!!
  subroutine cat_free_lattice()
     implicit none

     if ( allocated(sorts) ) deallocate(sorts)
     if ( allocated(atoms) ) deallocate(atoms)
     if ( allocated(sortn) ) deallocate(sortn)
     if ( allocated(lvect) ) deallocate(lvect)
     if ( allocated(coord) ) deallocate(coord)

     return
  end subroutine cat_free_lattice

!!
!! @sub cat_free_kmesh
!!
!! deallocate memory for kmesh-related variables
!!
  subroutine cat_free_kmesh()
     implicit none

     if ( allocated(kmesh)  ) deallocate(kmesh )
     if ( allocated(weight) ) deallocate(weight)

     return
  end subroutine cat_free_kmesh

!!
!! @sub cat_free_tetra
!!
!! deallocate memory for tetra-related variables
!!
  subroutine cat_free_tetra()
     implicit none

     if ( allocated(tetra) ) deallocate(tetra)

     return
  end subroutine cat_free_tetra

!!
!! @sub cat_free_eigen
!!
!! deallocate memory for eigen-related variables
!!
  subroutine cat_free_eigen()
     implicit none

     if ( allocated(enk)    ) deallocate(enk   )
     if ( allocated(occupy) ) deallocate(occupy)

     return
  end subroutine cat_free_eigen

!!
!! @sub cat_free_projs
!!
!! deallocate memory for projs-related variables
!!
  subroutine cat_free_projs()
     implicit none

     if ( allocated(psichi) ) deallocate(psichi)

     return
  end subroutine cat_free_projs

!!
!! @sub cat_free_fmesh
!!
!! deallocate memory for fmesh-related variables
!!
  subroutine cat_free_fmesh()
     implicit none

     if ( allocated(fmesh) ) deallocate(fmesh)

     return
  end subroutine cat_free_fmesh

!!
!! @sub cat_free_sigma
!!
!! deallocate memory for sigma-related variables
!!
  subroutine cat_free_sigma()
     implicit none

     if ( allocated(sigdc) ) deallocate(sigdc)
     if ( allocated(sig_l) ) deallocate(sig_l)
     if ( allocated(sig_k) ) deallocate(sig_k)

     return
  end subroutine cat_free_sigma

!!
!! @sub cat_free_green
!!
!! deallocate memory for green-related variables
!!
  subroutine cat_free_green()
     implicit none

     if ( allocated(grn_l) ) deallocate(grn_l)
     if ( allocated(grn_k) ) deallocate(grn_k)

     return
  end subroutine cat_free_green

!!
!! @sub cat_free_weiss
!!
!! deallocate memory for weiss-related variables
!!
  subroutine cat_free_weiss()
     implicit none

     if ( allocated(wss_l) ) deallocate(wss_l)
     if ( allocated(wss_k) ) deallocate(wss_k)
     if ( allocated(hyb_l) ) deallocate(hyb_l)
     if ( allocated(hyb_k) ) deallocate(hyb_k)

     return
  end subroutine cat_free_weiss

  end module context
