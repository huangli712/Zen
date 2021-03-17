!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : 
!!! source  : dmft_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           03/17/2021 by li huang (last modified)
!!! purpose : 
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  subroutine dmft_setup_param()
     use parser, only : p_create
     use parser, only : p_parse
     use parser, only : p_get
     use parser, only : p_destroy

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control

     implicit none

! local variables
! used to check whether the input file (dmft.in) exists
     logical :: exists

     task = 1
     axis = 1
     beta = 8.00_dp

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
             call p_get('task', task)
             call p_get('axis', axis)
             call p_get('beta', beta)

! destroy the parser
             call p_destroy()
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( task , master )
     call mp_bcast( axis , master )
     call mp_barrier()

     call mp_bcast( beta , master )
     call mp_barrier()

# endif  /* MPI */

     nsort  = 3
     natom  = 5
     nband  = 30
     nkpt   = 729
     nspin  = 1
     ntet   = 4374
     ngrp   = 1
     nwnd   = 1

     scal  = 4.00_dp
     volt  = 1.00_dp
     fermi = 0.00_dp

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: params.in
         inquire (file = 'params.ir', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then

             open(mytmp, file='params.ir', form='formatted', status='unknown')
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

     return 
  end subroutine dmft_setup_param

  subroutine dmft_alloc_array()
     implicit none

     return
  end subroutine dmft_alloc_array

  subroutine dmft_final_array()
     implicit none

     return
  end subroutine dmft_final_array
