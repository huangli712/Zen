!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_dump_grn_l
!!!           dmft_dump_wss_l
!!!           dmft_dump_hyb_l
!!! source  : dmft_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           05/01/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_dump_grn_l
!!
  subroutine dmft_dump_grn_l(grn_l)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite, nmesh

     use context, only : qdim
     use context, only : ndim

     implicit none

! external arguments
! local green's function
     complex(dp), intent(in) :: grn_l(qdim,qdim,nmesh,nspin,nsite)

! local variables
! loop index
     integer :: t
     integer :: s
     integer :: m
     integer :: p, q

! open data file: dmft_grn_l.dat
     open(mytmp, file='dmft_grn_l.dat', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin
             write(mytmp,'(a7,i4,2X,a5,i4)') "# site:", t, "spin:", s
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_grn_l

!!
!! @sub dmft_dump_wss_l
!!
  subroutine dmft_dump_wss_l()
     implicit none

     return
  end subroutine dmft_dump_wss_l

!!
!! @sub dmft_dump_hyb_l
!!
  subroutine dmft_dump_hyb_l()
     implicit none

     return
  end subroutine dmft_dump_hyb_l
