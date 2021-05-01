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
!! write out local green's function in matsubara frequency space
!!
  subroutine dmft_dump_grn_l(grn_l)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite, nmesh

     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

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

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(a7,i4,2X,a5,i4)') "# site:", t, "spin:", s
             do m=1,nmesh
                 write(mytmp,'(a2,f16.8)') 'w:', fmesh(m)
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, grn_l(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

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
