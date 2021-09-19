!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_setup_tasks
!!!           dmft_setup_param
!!!           dmft_setup_system
!!!           dmft_input_map
!!!           dmft_input_group
!!!           dmft_input_window
!!!           dmft_input_lattice
!!!           dmft_input_kmesh
!!!           dmft_input_tetra
!!!           dmft_input_eigen
!!!           dmft_input_projs
!!!           dmft_input_sigdc
!!!           dmft_input_sigma
!!!           dmft_alloc_array
!!!           dmft_final_array
!!! source  : dmft_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           09/19/2021 by li huang (last modified)
!!! purpose : setup the configuration parameters and read the input data.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> config dynamical mean-field theory engine                        <<<
!!========================================================================

!!
!! @sub dmft_setup_tasks
!!
!! setup control parameters for the dynamical mean-field theory engine.
!! note that these parameters are extracted from the `dmft.in` file.
!! this code can run even without `dmft.in` file.
!!
  subroutine dmft_setup_tasks()
     use constants, only : dp

     use parser, only : p_create
     use parser, only : p_parse
     use parser, only : p_get
     use parser, only : p_destroy

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : task
     use control, only : axis
     use control, only : lfermi, ltetra
     use control, only : beta, mc
     use control, only : myid, master

     implicit none

!! local variables
     ! used to check whether the input file (dmft.in) exists
     logical :: exists

!! [body

     ! setup default parameters
     !--------------------------------------------------------------------
     task   = 1         ! computational task
     axis   = 1         ! type of frequency mesh
     !--------------------------------------------------------------------
     beta   = 8.00_dp   ! inverse temperature
     mc     = 0.0001_dp ! convergence criterion for fermi level search
     !--------------------------------------------------------------------
     lfermi = .true.    ! fermi level search
     ltetra = .true.    ! analytical tetrahedron algorithm
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

         ! inquire file status: dmft.in
         inquire (file = 'dmft.in', exist = exists)

         ! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             ! create the file parser
             call p_create()

             ! parse the config file
             call p_parse('dmft.in')

             ! extract parameters
             call p_get('task'  , task  )
             call p_get('axis'  , axis  )

             call p_get('beta'  , beta  )
             call p_get('mc'    , mc    )

             call p_get('lfermi', lfermi)
             call p_get('ltetra', ltetra)

             ! destroy the parser
             call p_destroy()
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes.
# if defined (MPI)

     call mp_bcast( task  , master )
     call mp_bcast( axis  , master )
     call mp_barrier()

     call mp_bcast( beta  , master )
     call mp_bcast( mc    , master )
     call mp_barrier()

     call mp_bcast( lfermi, master )
     call mp_bcast( ltetra, master )
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_setup_tasks

!!
!! @sub dmft_setup_param
!!
!! setup dimensional parameters and some real constants for the dynamical
!! mean-field theory engine. note that these parameters are extracted
!! from the `params.ir` file. but this code can run even without this
!! configuration file.
!!
  subroutine dmft_setup_param()
     use constants, only : dp
     use constants, only : zero
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : model
     use control, only : nsort, natom
     use control, only : nband, nkpt, nspin
     use control, only : ntet
     use control, only : ngrp, nwnd
     use control, only : nsite, nmesh
     use control, only : scale, fermi, volt
     use control, only : myid, master

     use context, only : qdim, qbnd

     implicit none

!! local variables
     ! used to check whether the input file (params.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! setup default parameters
     !--------------------------------------------------------------------
     model  = 'SrVO3'   ! name of system
     !--------------------------------------------------------------------
     nsort  = 3         ! number of atomic sorts
     natom  = 5         ! number of atoms
     nband  = 30        ! number of bands
     nkpt   = 729       ! number of k-points
     nspin  = 1         ! number of spins
     ntet   = 4374      ! number of tetrahedra
     !--------------------------------------------------------------------
     ngrp   = 1         ! number of groups for projectors
     qdim   = 3         ! maximum number of projectors in groups
     nwnd   = 1         ! number of windows for projectors
     qbnd   = 5         ! maximum number of bands in windows
     !--------------------------------------------------------------------
     nsite  = 1         ! number of impurity sites
     nmesh  = 8193      ! number of frequency points
     !--------------------------------------------------------------------
     scale  = 4.00_dp   ! scale factor for lattice constants
     fermi  = 0.00_dp   ! default fermi level
     volt   = 1.00_dp   ! volume of a tetrahedron
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

         ! inquire file status: params.ir
         inquire (file = 'params.ir', exist = exists)

         ! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then

             open(mytmp, file='params.ir', form='formatted', status='unknown')

             read(mytmp,*) ! skip header
             read(mytmp,*)
             read(mytmp,*)

             read(mytmp,*) ! for lattice block
             read(mytmp,*) chr1, chr2, model
             read(mytmp,*) chr1, chr2, scale
             read(mytmp,*) chr1, chr2, nsort
             read(mytmp,*) chr1, chr2, natom
             read(mytmp,*)

             read(mytmp,*) ! for eigen block
             read(mytmp,*) chr1, chr2, nband
             read(mytmp,*) chr1, chr2, nkpt
             read(mytmp,*) chr1, chr2, nspin
             read(mytmp,*) chr1, chr2, fermi
             read(mytmp,*)

             read(mytmp,*) ! for tetra block
             read(mytmp,*) chr1, chr2, ntet
             read(mytmp,*) chr1, chr2, volt
             read(mytmp,*)

             read(mytmp,*) ! for group block
             read(mytmp,*) chr1, chr2, ngrp
             read(mytmp,*) chr1, chr2, qdim
             read(mytmp,*)

             read(mytmp,*) ! for window block
             read(mytmp,*) chr1, chr2, nwnd
             read(mytmp,*) chr1, chr2, qbnd
             read(mytmp,*)

             read(mytmp,*) ! for sigma block
             read(mytmp,*) chr1, chr2, nsite
             read(mytmp,*) chr1, chr2, nmesh
             read(mytmp,*)

             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes.
# if defined (MPI)

     ! for lattice block
     call mp_bcast( model , master )
     call mp_bcast( scale , master )
     call mp_bcast( nsort , master )
     call mp_bcast( natom , master )
     call mp_barrier()

     ! for eigen block
     call mp_bcast( nband , master )
     call mp_bcast( nkpt  , master )
     call mp_bcast( nspin , master )
     call mp_bcast( fermi , master )
     call mp_barrier()

     ! for tetra block
     call mp_bcast( ntet  , master )
     call mp_bcast( volt  , master )
     call mp_barrier()

     ! for group block
     call mp_bcast( ngrp  , master )
     call mp_bcast( qdim  , master )
     call mp_barrier()

     ! for window block
     call mp_bcast( nwnd  , master )
     call mp_bcast( qbnd  , master )
     call mp_barrier()

     ! for sigma block
     call mp_bcast( nsite , master )
     call mp_bcast( nmesh , master )
     call mp_barrier()

# endif  /* MPI */

     ! reset fermi to zero, because it was calibrated by the adaptor.
     ! see ZenCore/src/plo.jl/plo_fermi() function for more details.
     fermi = zero

!! body]

     return
  end subroutine dmft_setup_param

!!
!! @sub dmft_setup_system
!!
!! setup correlated electron problem, Kohn-Sham dataset, and self-energy
!! functions for the dynamical mean-field theory engine. this subroutine
!! is only a dispatcher for the other individual subroutines.
!!
  subroutine dmft_setup_system()
     use control, only : ltetra

     implicit none

!! [body

     ! setup correlated electron problem
     call dmft_input_map()
     call dmft_input_group()
     call dmft_input_window()

     ! setup Kohn-Sham dataset
     call dmft_input_lattice()
     call dmft_input_kmesh()
     call dmft_input_eigen()
     call dmft_input_projs()
     !
     if ( ltetra .eqv. .true. ) then
         call dmft_input_tetra()
     endif ! back if ( ltetra .eqv. .true. ) block

     ! setup impurity self-energy functions and double counting terms
     call dmft_input_sigdc()
     call dmft_input_sigma()

!! body]

     return
  end subroutine dmft_setup_system

!!========================================================================
!!>>> setup correlated electron problem                                <<<
!!========================================================================

!!
!! @sub dmft_input_map
!!
!! read in connections / mappings between quantum impurity problems and
!! groups of projectors (see module dmft_map). the data can be used to
!! upfold or downfold the self-energy functions or green's functions.
!!
  subroutine dmft_input_map()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : ngrp, nwnd
     use control, only : nsite
     use control, only : myid, master

     use context, only : i_grp, i_wnd
     use context, only : g_imp, w_imp

     implicit none

!! local variables
     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (maps.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in mappings or connections between quantum impurity problems
     ! and groups of projectors.
     !
     ! this code can not run without the file `maps.ir`
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'maps.ir', exist = exists)

         ! file maps.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_map','file maps.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file maps.ir for reading
         open(mytmp, file='maps.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*)

         ! check nsite, ngrp, and nwnd
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nsite, 'nsite is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == ngrp, 'ngrp is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nwnd, 'nwnd is wrong')

         ! read data: i_grp
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) i_grp

         ! read data: i_wnd
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) i_wnd

         ! read data: g_imp
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) g_imp

         ! read data: w_imp
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) w_imp

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( i_grp, master )
     call mp_bcast( i_wnd, master )
     call mp_bcast( g_imp, master )
     call mp_bcast( w_imp, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     ! additional check for the data
     ! all of the impurity problems should share the same band window!
     call s_assert2(nsite <= ngrp, 'nsite must be smaller or equal to ngrp')
     !
     call s_assert2(nwnd == ngrp, 'nwnd must be equal to ngrp')

!! body]

     return
  end subroutine dmft_input_map

!!
!! @sub dmft_input_group
!!
!! read in groups of projectors (see module dmft_group). the data can
!! be used to upfold or downfold the self-energy functions and green's
!! functions.
!!
  subroutine dmft_input_group()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : ngrp
     use control, only : myid, master

     use context, only : qdim
     use context, only : shell, corr, site, l, ndim

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (groups.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in groups of projectors.
     ! apparently, this code can not run without the file `groups.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'groups.ir', exist = exists)

         ! file groups.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_group','file groups.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file groups.ir for reading
         open(mytmp, file='groups.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*)

         ! check ngrp
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == ngrp, 'ngrp is wrong')
         read(mytmp,*)

         ! read data
         do i=1,ngrp
             read(mytmp,*)
             read(mytmp,*) chr1, chr2, site(i)
             read(mytmp,*) chr1, chr2, l(i)
             read(mytmp,*) chr1, chr2, corr(i)
             read(mytmp,*) chr1, chr2, shell(i)
             read(mytmp,*) chr1, chr2, ndim(i)
             read(mytmp,*)
         enddo ! over i={1,ngrp} loop

         ! evaluate and check qdim
         itmp = maxval(ndim)
         call s_assert2(itmp == qdim, 'ndim or qdim is wrong')

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( shell, master )
     call mp_bcast( corr , master )
     call mp_bcast( site , master )
     call mp_bcast( l    , master )
     call mp_bcast( ndim , master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_group

!!
!! @sub dmft_input_window
!!
!! read in band windows of projectors (see module dmft_window). the data
!! are used to upfold or downfold the self-energy functions and green's
!! functions.
!!
  subroutine dmft_input_window()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nkpt, nspin
     use control, only : nwnd
     use control, only : myid, master

     use context, only : qbnd
     use context, only : xbnd
     use context, only : bmin, bmax, nbnd, qwin, kwin

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: s
     integer :: k

     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (windows.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in band windows of projectors.
     ! apparently, this code can not run without the file `windows.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'windows.ir', exist = exists)

         ! file windows.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_window','file windows.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file windows.ir for reading
         open(mytmp, file='windows.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*)

         ! check nwnd
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nwnd, 'nwnd is wrong')
         read(mytmp,*)

         ! read data
         do i=1,nwnd
             read(mytmp,*)
             read(mytmp,*) chr1, chr2, bmin(i)
             read(mytmp,*) chr1, chr2, bmax(i)
             read(mytmp,*) chr1, chr2, nbnd(i)
             read(mytmp,*) ! for kwin
             do s=1,nspin
                 do k=1,nkpt
                     read(mytmp,*) itmp, itmp, kwin(k,s,1,i), kwin(k,s,2,i)
                 enddo ! over k={1,nkpt} loop
             enddo ! over s={1,nspin} loop
             read(mytmp,*)
         enddo ! over i={1,nwnd} loop

         ! evaluate and check qbnd
         itmp = maxval(nbnd)
         call s_assert2(itmp == qbnd, 'nbnd or qbnd is wrong')

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( bmin , master )
     call mp_bcast( bmax , master )
     call mp_bcast( nbnd , master )
     call mp_bcast( kwin , master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     ! well, next we will try to build qwin and xbnd
     do s=1,nspin
         do k=1,nkpt
             qwin(k,s,1) = minval(kwin(k,s,1,:))
             qwin(k,s,2) = maxval(kwin(k,s,2,:))
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop
     !
     xbnd = maxval(bmax) - minval(bmin) + 1

!! body]

     return
  end subroutine dmft_input_window

!!========================================================================
!!>>> setup Kohn-Sham dataset                                          <<<
!!========================================================================

!!
!! @sub dmft_input_lattice
!!
!! read in key crystallography information (see module dmft_lattice).
!!
  subroutine dmft_input_lattice()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nsort, natom
     use control, only : myid, master

     use context, only : sorts, atoms, sortn, lvect, coord

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (lattice.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in crystallography information.
     ! this code can not run without the file `lattice.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'lattice.ir', exist = exists)

         ! file lattice.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_lattice','file lattice.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file lattice.ir for reading
         open(mytmp, file='lattice.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)

         ! check nsort and natom
         read(mytmp,*) ! empty line
         read(mytmp,*) ! skip _case
         read(mytmp,*) ! skip scale
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nsort, 'nsort is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == natom, 'natom is wrong')
         !
         read(mytmp,*) ! empty line

         ! read sorts
         read(mytmp,*) ! header
         read(mytmp,*) sorts
         read(mytmp,*) sortn
         read(mytmp,*) ! empty line

         ! read atoms
         read(mytmp,*) ! header
         read(mytmp,*) atoms
         read(mytmp,*) ! empty line

         ! read lvect
         read(mytmp,*) ! header
         do i=1,3
             read(mytmp,*) lvect(i,1:3)
         enddo ! over i={1,3} loop
         read(mytmp,*) ! empty line

         ! read coord
         read(mytmp,*) ! header
         do i=1,natom
             read(mytmp,*) coord(i,1:3)
         enddo ! over i={1,natom} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( sorts, master )
     call mp_bcast( atoms, master )
     call mp_bcast( sortn, master )
     call mp_bcast( lvect, master )
     call mp_bcast( coord, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_lattice

!!
!! @sub dmft_input_kmesh
!!
!! read in k-mesh and integration weights (see module dmft_kmesh).
!!
  subroutine dmft_input_kmesh()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nkpt
     use control, only : myid, master

     use context, only : kmesh, weight

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (kmesh.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in brillouin zone information.
     ! the code can not run without the file `kmesh.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'kmesh.ir', exist = exists)

         ! file kmesh.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_kmesh','file kmesh.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file kmesh.ir for reading
         open(mytmp, file='kmesh.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) ! empty line

         ! check nkpt and ndir
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nkpt, 'nkpt is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == 3, 'ndir is wrong')
         !
         read(mytmp,*) ! empty line

         ! read k-points and the corresponding weights
         do i=1,nkpt
             read(mytmp,*) kmesh(i,:), weight(i)
         enddo ! over i={1,nkpt} loop

         ! until now, `weight' has not been renormalized.
         ! their summations should be equal to nkpt.
         call s_assert2(sum(weight) == float(nkpt), 'weight is wrong')

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( kmesh, master )
     call mp_bcast(weight, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_kmesh

!!
!! @sub dmft_input_tetra
!!
!! read in data for tetrahedron integration (see module dmft_tetra).
!! note that this subroutine is only called when `ltetra` is true.
!!
  subroutine dmft_input_tetra()
     use constants, only : dp, mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : ntet, volt
     use control, only : myid, master

     use context, only : tetra

     implicit none

!! local variables
     ! loop index
     integer  :: i

     ! dummy integer variables
     integer  :: itmp

     ! used to check whether the input file (tetra.ir) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: rtmp

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in tetrahedron information if available
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'tetra.ir', exist = exists)

         ! file tetra.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_tetra','file tetra.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file tetra.ir for reading
         open(mytmp, file='tetra.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) ! empty line

         ! check ntet and volt
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == ntet, 'ntet is wrong')
         !
         read(mytmp,*) chr1, chr2, rtmp
         call s_assert2(rtmp == volt, 'volt is wrong')
         !
         read(mytmp,*) ! empty line

         ! read tetrahedron data
         do i=1,ntet
             read(mytmp,*) tetra(i,5), tetra(i,1:4)
         enddo ! over i={1,ntet} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( tetra, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_tetra

!!
!! @sub dmft_input_eigen
!!
!! read in band eigenvalues and band occupations (see module dmft_eigen).
!!
  subroutine dmft_input_eigen()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nband, nkpt, nspin
     use control, only : myid, master

     use context, only : enk, occupy

     implicit none

!! local variables
     ! loop index
     integer :: b
     integer :: k
     integer :: s

     ! dummy integer variables
     integer :: itmp

     ! used to check whether the input file (eigen.ir) exists
     logical :: exists

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in Kohn-Sham band structure information.
     ! this code can not run without the file `eigen.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'eigen.ir', exist = exists)

         ! file eigen.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_eigen','file eigen.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file eigen.ir for reading
         open(mytmp, file='eigen.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) ! empty line

         ! check nband, nkpt, and nspin
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nband, 'nband is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nkpt, 'nkpt is wrong')
         !
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nspin, 'nspin is wrong')
         !
         read(mytmp,*) ! empty line

         ! read eigenvalues and occupations data
         do s=1,nspin
             do k=1,nkpt
                 do b=1,nband
                     read(mytmp,*) enk(b,k,s), occupy(b,k,s)
                 enddo ! over b={1,nband} loop
             enddo ! over k={1,nkpt} loop
         enddo ! over s={1,nspin} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( enk  , master )
     call mp_bcast(occupy, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_eigen

!!
!! @sub dmft_input_projs
!!
!! read in overlap matrix between Kohn-Sham wavefunctions and local
!! orbitals, i.e., the local orbital projectors (see module dmft_projs).
!!
  subroutine dmft_input_projs()
     use constants, only : dp, mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : ngrp
     use control, only : nkpt, nspin
     use control, only : myid, master

     use context, only : ndim, nbnd
     use context, only : chipsi, psichi

     implicit none

!! local variables
     ! loop index
     integer  :: g
     integer  :: b
     integer  :: k
     integer  :: s
     integer  :: d

     ! dummy integer variables
     integer  :: itmp

     ! used to check whether the input file (projs.ir) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: re, im

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in local orbital projectors.
     ! this code can not run without the file `projs.ir`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'projs.ir', exist = exists)

         ! file projs.ir must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_projs','file projs.ir is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file projs.ir for reading
         open(mytmp, file='projs.ir', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)

         ! go through each group of projectors and read in the data
         do g=1,ngrp

             ! check group
             read(mytmp,*) ! empty line
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == g, 'group is wrong')

             ! check nproj
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == ndim(g), 'nproj is wrong')

             ! check nband
             !
             ! here, `nband` is not the number of bands used in the dft
             ! calculations. it is just the number of bands that contained
             ! in the band window for constructing the local orbital
             ! projection. for different groups of projectors, the related
             ! `nband` could be and might be different.
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == nbnd(g), 'nband is wrong')

             ! check nkpt
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == nkpt, 'nkpt is wrong')

             ! check nspin
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == nspin, 'nspin is wrong')
             read(mytmp,*) ! empty line

             ! parse the data
             do s=1,nspin
                 do k=1,nkpt
                     do b=1,nbnd(g)
                         do d=1,ndim(g)
                             read(mytmp,*) re, im
                             chipsi(d,b,k,s,g) = dcmplx(re,+im)
                             psichi(b,d,k,s,g) = dcmplx(re,-im)
                         enddo ! over d={1,ndim(g)} loop
                     enddo ! over b={1,nbnd(g)} loop
                 enddo ! over k={1,nkpt} loop
             enddo ! over s={1,nspin} loop

         enddo ! over g={1,ngrp} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast(chipsi, master )
     call mp_bcast(psichi, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_projs

!!========================================================================
!!>>> setup self-energy functions and double counting terms            <<<
!!========================================================================

!!
!! @sub dmft_input_sigdc
!!
!! read in double counting terms for the local self-energy functions (see
!! module dmft_sigma).
!!
  subroutine dmft_input_sigdc()
     use constants, only : dp, mytmp
     use constants, only : czero

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nspin
     use control, only : nsite
     use control, only : myid, master

     use context, only : i_grp
     use context, only : qdim, ndim
     use context, only : sigdc

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: s
     integer  :: m
     integer  :: n

     ! dummy integer variables
     integer  :: itmp

     ! used to check whether the input file (sigma.dc) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: re, im

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in double counting terms.
     ! the code can not run without the file `sigma.dc`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'sigma.dc', exist = exists)

         ! file sigma.dc must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_sigdc','file sigma.dc is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file sigma.dc for reading
         open(mytmp, file='sigma.dc', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) ! empty line

         ! check nsite
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nsite, 'nsite is wrong')

         ! check nspin
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nspin, 'nspin is wrong')

         ! check ndim
         do i=1,nsite
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == ndim(i_grp(i)), 'ndim is wrong')
             call s_assert2(itmp <= qdim, 'ndim is wrong')
         enddo ! over i={1,nsite} loop
         read(mytmp,*) ! empty line

         ! parse the data
         do i=1,nsite
             do s=1,nspin
                 !
                 read(mytmp,*) ! empty line
                 !
                 do m=1,ndim(i_grp(i))
                     do n=1,ndim(i_grp(i))
                         read(mytmp,*) re, im
                         sigdc(n,m,s,i_grp(i)) = dcmplx(re, im)
                     enddo ! over n={1,ndim(i_grp(i))} loop
                 enddo ! over m={1,ndim(i_grp(i))} loop
                 !
                 read(mytmp,*) ! empty line
                 !
             enddo ! over s={1,nspin} loop
         enddo ! over i={1,nsite} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( sigdc, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_sigdc

!!
!! @sub dmft_input_sigma
!!
!! read in bare self-energy functions from various quantum impurity
!! solvers (see module dmft_sigma).
!!
  subroutine dmft_input_sigma()
     use constants, only : dp, mytmp
     use constants, only : zero

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : beta
     use control, only : myid, master

     use context, only : i_grp
     use context, only : qdim, ndim
     use context, only : fmesh, sigma

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: s
     integer  :: m

     ! dummy integer variables
     integer  :: itmp

     ! used to check whether the input file (sigma.bare) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: rtmp
     real(dp) :: re, im

     ! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

!! [body

     ! read in local self-energy functions.
     ! this code can not run without the file `sigma.bare`.
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'sigma.bare', exist = exists)

         ! file sigma.bare must be present
         if ( exists .eqv. .false. ) then
             call s_print_error('dmft_input_sigma','file sigma.bare is absent')
         endif ! back if ( exists .eqv. .false. ) block

         ! open file sigma.bare for reading
         open(mytmp, file='sigma.bare', form='formatted', status='unknown')

         ! skip header
         read(mytmp,*)
         read(mytmp,*)
         read(mytmp,*) ! empty line

         ! check axis
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == axis, 'axis is wrong')

         ! check beta
         read(mytmp,*) chr1, chr2, rtmp
         call s_assert2(rtmp == beta, 'beta is wrong')

         ! check nsite
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nsite, 'nsite is wrong')

         ! check nmesh
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nmesh, 'nmesh is wrong')

         ! check nspin
         read(mytmp,*) chr1, chr2, itmp
         call s_assert2(itmp == nspin, 'nspin is wrong')

         ! check ndim
         do i=1,nsite
             read(mytmp,*) chr1, chr2, itmp
             call s_assert2(itmp == ndim(i_grp(i)), 'ndim is wrong')
             call s_assert2(itmp <= qdim, 'ndim is wrong')
         enddo ! over i={1,nsite} loop
         read(mytmp,*) ! empty line

         ! parse the data
         do i=1,nsite
             do s=1,nspin
                 !
                 read(mytmp,*) ! empty line
                 !
                 do m=1,nmesh
                     read(mytmp,*) chr2, itmp, fmesh(m)
                     do j=1,ndim(i_grp(i))
                         do k=1,ndim(i_grp(i))
                             read(mytmp,*) re, im
                             sigma(k,j,m,s,i_grp(i)) = dcmplx(re, im)
                         enddo ! over k={1,ndim(i_grp(i))} loop
                     enddo ! over j={1,ndim(i_grp(i))} loop
                 enddo ! over m={1,nmesh} loop
                 !
                 read(mytmp,*) ! empty line
                 !
             enddo ! over s={1,nspin} loop
         enddo ! over i={1,nsite} loop

         ! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

     ! block until all processes have reached here
     call mp_barrier()

     ! broadcast data
     call mp_bcast( fmesh, master )
     call mp_bcast( sigma, master )

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine dmft_input_sigma

!!========================================================================
!!>>> manage memory for dynamical mean-field theory engine             <<<
!!========================================================================

!!
!! @sub dmft_alloc_array
!!
!! allocate memory for global variables and then initialize them.
!!
  subroutine dmft_alloc_array()
     use context ! ALL

     implicit none

!! [body

     ! allocate memory for context module
     call cat_alloc_map()
     call cat_alloc_group()
     call cat_alloc_window()

     call cat_alloc_lattice()
     call cat_alloc_kmesh()
     call cat_alloc_tetra()
     call cat_alloc_eigen()
     call cat_alloc_projs()

     call cat_alloc_fmesh()
     call cat_alloc_eimps()
     call cat_alloc_sigma()
     call cat_alloc_green()
     call cat_alloc_weiss()
     call cat_alloc_delta()

     call cat_alloc_gcorr()

!! body]

     return
  end subroutine dmft_alloc_array

!!
!! @sub dmft_final_array
!!
!! garbage collection for this code, please refer to dmft_alloc_array.
!!
  subroutine dmft_final_array()
     use context ! ALL

     implicit none

!! [body

     ! deallocate memory for context module
     call cat_free_map()
     call cat_free_group()
     call cat_free_window()

     call cat_free_lattice()
     call cat_free_kmesh()
     call cat_free_tetra()
     call cat_free_eigen()
     call cat_free_projs()

     call cat_free_fmesh()
     call cat_free_eimps()
     call cat_free_sigma()
     call cat_free_green()
     call cat_free_weiss()
     call cat_free_delta()

     call cat_free_gcorr()

!! body]

     return
  end subroutine dmft_final_array
