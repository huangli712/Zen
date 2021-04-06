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
!!!           04/06/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

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

             print *, model, scale, nband, fermi, ntet
             

             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

     print *, 'hehe'

     return
  end subroutine dmft_setup_param

!!
!! @sub dmft_setup_system
!!
  subroutine dmft_setup_system()
     implicit none

     return
  end subroutine dmft_setup_system

!!
!! @sub dmft_input_group
!!
  subroutine dmft_input_group()
     implicit none

     return
  end subroutine dmft_input_group

!!
!! @sub dmft_input_window
!!
  subroutine dmft_input_window()
     implicit none

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
!! @sub dft_input_bzone
!!
  subroutine dft_input_bzone()
     implicit none

     return
  end subroutine dft_input_bzone

!!
!! @sub dft_input_tetra
!!
  subroutine dft_input_tetra()
     implicit none

     return
  end subroutine dft_input_tetra

!!
!! @sub dft_input_eigen
!!
  subroutine dft_input_eigen()
     implicit none

     return
  end subroutine dft_input_eigen

!!
!! @sub dft_input_projs
!!
  subroutine dft_input_projs()
     implicit none

     return
  end subroutine dft_input_projs

!!
!! @sub dft_input_sigdc
!!
  subroutine dft_input_sigdc()
     implicit none

     return
  end subroutine dft_input_sigdc

!!
!! @sub dft_input_sig_l
!!
  subroutine dft_input_sig_l()
     implicit none

     return
  end subroutine dft_input_sig_l

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
