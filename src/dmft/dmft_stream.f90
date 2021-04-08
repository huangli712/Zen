!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_setup_tasks
!!!           dmft_setup_param
!!!           dmft_setup_system
!!!           dmft_alloc_array
!!!           dmft_final_array
!!! source  : dmft_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           04/08/2021 by li huang (last modified)
!!! purpose :
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
!! note that these parameters are extracted from the dmft.in file
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

! local variables
! used to check whether the input file (dmft.in) exists
     logical :: exists

! setup default parameters
!-------------------------------------------------------------------------
     task   = 1         ! computational task
     axis   = 1         ! frequency mesh
!-------------------------------------------------------------------------
     lfermi = .true.    ! fermi level search
     ltetra = .true.    ! tetrahedron algorithm
!-------------------------------------------------------------------------
     beta   = 8.00_dp   ! inverse temperature
     mc     = 0.0001_dp ! convergence criterion for fermi level search
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

             call p_get('lfermi', lfermi)
             call p_get('ltetra', ltetra)

             call p_get('beta'  , beta  )
             call p_get('mc'    , mc    )

! destroy the parser
             call p_destroy()
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( task  , master )
     call mp_bcast( axis  , master )
     call mp_barrier()

     call mp_bcast( lfermi, master )
     call mp_bcast( ltetra, master )
     call mp_barrier()

     call mp_bcast( beta  , master )
     call mp_bcast( mc    , master )
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine dmft_setup_tasks

!!
!! @sub dmft_setup_param
!!
!! setup dimensional parameters and some real constants for the dynamical
!! mean-field theory engine. note that these parameters are extracted
!! from the dmft.in file
!!
  subroutine dmft_setup_param()
     use constants, only : dp
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

     implicit none

! local variables
! used to check whether the input file (params.ir) exists
     logical :: exists

! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

! setup default parameters
!-------------------------------------------------------------------------
     model  = 'SrVO3'   ! name of system
!-------------------------------------------------------------------------
     nsort  = 3         ! number of atomic sorts
     natom  = 5         ! number of atoms
     nband  = 30        ! number of bands
     nkpt   = 729       ! number of k-points
     nspin  = 1         ! number of spins
     ntet   = 4374      ! number of tetrahedra
!-------------------------------------------------------------------------
     ngrp   = 1         ! number of groups for projectors
     nwnd   = 1         ! number of windows for projectors
!-------------------------------------------------------------------------
     nsite  = 1         ! number of impurity sites
     nmesh  = 8193      ! number of frequency points
!-------------------------------------------------------------------------
     scale  = 4.00_dp   ! scale factor for lattice constants
     fermi  = 0.00_dp   ! default fermi level
     volt   = 1.00_dp   ! volume of a tetrahedron
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: params.in
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
             read(mytmp,*)

             read(mytmp,*) ! for window block
             read(mytmp,*) chr1, chr2, nwnd
             read(mytmp,*)

             read(mytmp,*) ! for sigma block
             read(mytmp,*) chr1, chr2, nsite
             read(mytmp,*) chr1, chr2, nmesh
             read(mytmp,*)

             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes
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
     call mp_barrier()

! for window block
     call mp_bcast( nwnd  , master )
     call mp_barrier()

! for sigma block
     call mp_bcast( nsite , master )
     call mp_bcast( nmesh , master )
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine dmft_setup_param

!!========================================================================
!!>>> parse Kohn-Sham data for dynamical mean-field theory engine      <<<
!!========================================================================

!!
!! @sub dmft_setup_system
!!
!! parse Kohn-Sham data for dynamical mean-field theory engine. this is
!! an entry for the other individual subroutines
!!
  subroutine dmft_setup_system()
     implicit none

! get descriptions of correlated subspace
     call dmft_input_group()
     call dmft_input_window()

! get Kohn-Sham data
     call dmft_input_lattice()
     call dmft_input_bzone()
     call dmft_input_tetra()
     call dmft_input_eigen()
     call dmft_input_projs()

! get impurity self-energy functions
     call dmft_input_sigdc()
     call dmft_input_sig_l()

     return
  end subroutine dmft_setup_system

!!
!! @sub dmft_input_group
!!
!! read in groups of projectors (see module dmft_group). the data are
!! used to embed or downfold the self-energy functions
!!
  subroutine dmft_input_group()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : ngrp
     use control, only : myid, master

     use context, only : max_ndim
     use context, only : shell, corr, site, l, ndim

     implicit none

! local variables
! loop index
     integer :: i

! dummy integer variables
     integer :: itmp

! used to check whether the input file (groups.ir) exists
     logical :: exists

! dummy character variables
     character(len = 5) :: chr1
     character(len = 2) :: chr2

! read in groups of projectors if available
!-------------------------------------------------------------------------
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

! check ngrp
         read(mytmp,*)
         read(mytmp,*) chr1, chr2, itmp
         read(mytmp,*)
         call s_assert2(itmp == ngrp, "ngrp is wrong")

! read data
         do i=1,ngrp
             read(mytmp,*)
             read(mytmp,*) chr1, chr2, site(i)
             read(mytmp,*) chr1, chr2, l(i)
             read(mytmp,*) chr1, chr2, corr(i)
             read(mytmp,*) chr1, chr2, shell(i)
             read(mytmp,*) chr1, chr2, ndim(i)
         enddo ! over i={1,ngrp} loop

! evaluate max_ndim
         max_ndim = maxval(ndim)

! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast( max_ndim, master )

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

     return
  end subroutine dmft_input_group

!!
!! @sub dmft_input_window
!!
!! read in windows of projectors (see module dmft_window). the data are
!! used to embed or downfold the self-energy functions
!!
  subroutine dmft_input_window()
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : nkpt, nspin
     use control, only : nwnd
     use control, only : myid, master

     use context, only : max_nbnd
     use context, only : bmin, bmax, nbnd, kwin

     implicit none

! local variables
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

! read in windows of projectors if available
!-------------------------------------------------------------------------
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

! check nwnd
         read(mytmp,*)
         read(mytmp,*) chr1, chr2, itmp
         read(mytmp,*)
         call s_assert2(itmp == nwnd, "nwnd is wrong")

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
         enddo ! over i={1,ngrp} loop

! evaluate max_nbnd
         max_nbnd = maxval(nbnd)

! close file handler
         close(mytmp)

     endif ! back if ( myid == master ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast data from master node to all children nodes
# if defined (MPI)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast( max_nbnd, master )

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

     return
  end subroutine dmft_input_window

!!
!! @sub dmft_input_lattice
!!
  subroutine dmft_input_lattice()
     implicit none

     return
  end subroutine dmft_input_lattice

!!
!! @sub dmft_input_bzone
!!
  subroutine dmft_input_bzone()
     implicit none

     return
  end subroutine dmft_input_bzone

!!
!! @sub dmft_input_tetra
!!
  subroutine dmft_input_tetra()
     implicit none

     return
  end subroutine dmft_input_tetra

!!
!! @sub dmft_input_eigen
!!
  subroutine dmft_input_eigen()
     implicit none

     return
  end subroutine dmft_input_eigen

!!
!! @sub dmft_input_projs
!!
  subroutine dmft_input_projs()
     implicit none

     return
  end subroutine dmft_input_projs

!!
!! @sub dmft_input_sigdc
!!
  subroutine dmft_input_sigdc()
     implicit none

     return
  end subroutine dmft_input_sigdc

!!
!! @sub dmft_input_sig_l
!!
  subroutine dmft_input_sig_l()
     implicit none

     return
  end subroutine dmft_input_sig_l

!!========================================================================
!!>>> manage memory for dynamical mean-field theory engine             <<<
!!========================================================================

!!
!! @sub dmft_alloc_array
!!
!! allocate memory for global variables and then initialize them
!! 
  subroutine dmft_alloc_array()
     use context ! ALL

     implicit none

! allocate memory for context module
     call cat_alloc_group()
     call cat_alloc_window()

     call cat_alloc_lattice()
     call cat_alloc_bzone()
     call cat_alloc_tetra()
     call cat_alloc_eigen()
     call cat_alloc_projs()

     call cat_alloc_fmesh()

     call cat_alloc_sigma()
     call cat_alloc_green()
     call cat_alloc_weiss()

     return
  end subroutine dmft_alloc_array

!!
!! @sub dmft_final_array
!!
!! garbage collection for this code, please refer to dmft_alloc_array
!!
  subroutine dmft_final_array()
     use context ! ALL

     implicit none

! deallocate memory for context module
     call cat_free_group()
     call cat_free_window()

     call cat_free_lattice()
     call cat_free_bzone()
     call cat_free_tetra()
     call cat_free_eigen()
     call cat_free_projs()

     call cat_free_fmesh()

     call cat_free_sigma()
     call cat_free_green()
     call cat_free_weiss()

     return
  end subroutine dmft_final_array
