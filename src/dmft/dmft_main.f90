!!!=========+=========+=========+=========+=========+=========+=========+!
!!! ZEN @ JACARANDA                                                      !
!!!                                                                      !
!!! A highly optimized dynamical mean field theory self-consistent calc- !
!!! ulation engine, which will be used together with density functional  !
!!! theory codes to explore the realistic correlated electron materials. !
!!!                                                                      !
!!! author  : Li Huang (China Academy of Engineering Physics)            !
!!! status  : (WARNING) IN TESTING STAGE, USE IT IN YOUR RISK            !
!!! comment : any question, please contact with lihuang.dmft@gmail.com   !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!========================================================================
  PROGRAM DMFT_MAIN !                                                  <<<
!!========================================================================

     use mmpi, only : mp_init      ! init mpi environment
     use mmpi, only : mp_finalize  ! finalize mpi environment
     use mmpi, only : mp_barrier   ! barrier to synchronize the data
     use mmpi, only : mp_comm_rank ! get index of current process
     use mmpi, only : mp_comm_size ! get number of processes

     use control, only : nprocs    ! number of processes
     use control, only : myid      ! index of current process
     use control, only : master    ! index of master process

     implicit none

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

     DMFT_START: BLOCK

! print the welcome messages
         if ( myid == master ) then ! only master node can do it
             call dmft_print_header()
         endif ! back if ( myid == master ) block

! setup the parameters
         call dmft_setup_tasks()
         call dmft_setup_param()
         STOP

! allocate memory spaces
         call dmft_alloc_array()

! print the runtime parameters
         if ( myid == master ) then ! only master node can do it
             call dmft_print_summary()
         endif ! back if ( myid == master ) block

     END BLOCK DMFT_START

     DMFT_SLEEP: BLOCK

! deallocate memory spaces
         call dmft_final_array()

! print the ending messages
         if ( myid == master ) then ! only master node can do it
             call dmft_print_footer()
         endif ! back if ( myid == master ) block

     END BLOCK DMFT_SLEEP

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

!!========================================================================
  END PROGRAM DMFT_MAIN !                                              <<<
!!========================================================================
