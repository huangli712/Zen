!!!=========+=========+=========+=========+=========+=========+=========+!
!!! ZEN @ JACARANDA                                                      !
!!!                                                                      !
!!! A highly optimized hybridization expansion version continuous time   !
!!! quantum Monte Carlo (CTQMC) quantum impurity solver plus a classic   !
!!! dynamical mean field theory (DMFT) self-consistent engine            !
!!!                                                                      !
!!! author  : Li Huang (China Academy of Engineering Physics)            !
!!! status  : (WARNING) IN TESTING STAGE, USE IT IN YOUR RISK            !
!!! comment : this impurity solver is based on segment picture formalism !
!!!           any question, please contact with lihuang.dmft@gmail.com   !
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
