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

     if ( myid == master ) then
         print *, myid, nprocs
     endif

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
