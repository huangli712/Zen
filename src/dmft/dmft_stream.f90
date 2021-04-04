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
!!!           04/04/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  subroutine dmft_setup_tasks()
     use constants, only : dp

     use control, only : task
     use control, only : axis
     use control, only : beta
     use control, only : myid, master

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

     return
  end subroutine dmft_setup_tasks

  subroutine dmft_setup_param()
     use constants, only : dp

     use control, only : nsort, natom
     use control, only : nband, nkpt, nspin
     use control, only : ntet
     use control, only : ngrp, nwnd
     use control, only : scale, fermi, volt
     use control, only : myid, master

     implicit none

     nsort  = 3
     natom  = 5
     nband  = 30
     nkpt   = 729
     nspin  = 1
     ntet   = 4374
     ngrp   = 1
     nwnd   = 1

     scale  = 4.00_dp
     fermi = 0.00_dp
     volt  = 1.00_dp

     return
  end subroutine dmft_setup_param

  subroutine dmft_setup_param1()
     !use parser, only : p_create
     !use parser, only : p_parse
     !use parser, only : p_get
     !use parser, only : p_destroy

     !use mmpi, only : mp_bcast
     !use mmpi, only : mp_barrier

     !use constants, only : mytmp
     !use control

     implicit none


! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: params.in
         inquire (file = 'params.ir', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then

             open(mytmp, file='params.ir', form='formatted', status='unknown')
             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*)
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_setup_param1

  subroutine dmft_setup_system()
     implicit none

     return
  end subroutine dmft_setup_system

  subroutine dmft_input_group()
  end

  subroutine dmft_input_window()
  end

  subroutine dmft_input_lattice()
  end

  subroutine dft_input_bzone()
  end
 
  subroutine dmft_alloc_array()
     implicit none

     return
  end subroutine dmft_alloc_array

  subroutine dmft_final_array()
     implicit none

     return
  end subroutine dmft_final_array
