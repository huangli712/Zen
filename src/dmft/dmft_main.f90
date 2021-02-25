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

     implicit none

     print *, "Hello world!"

!!========================================================================
  END PROGRAM DMFT_MAIN !                                              <<<
!!========================================================================
