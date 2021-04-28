!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_driver
!!!           map_chi_psi
!!!           map_psi_chi
!!! source  : dmft_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           04/28/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_driver
!!
  subroutine dmft_driver()
     implicit none

     call cal_grn_k()

     return
  end subroutine dmft_driver

!!
!! @sub cal_sig_k
!!
  subroutine cal_sig_k()
     implicit none

     return
  end subroutine cal_sig_k

!!
!! @sub cal_grn_k
!!
  subroutine cal_grn_k(w)
     use constants, only : dp
     use constants, only : czero, czi

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : fermi

     use context, only : qbnd
     use context, only : kwin
     use context, only : enk
     use context, only : fmesh
     use context, only : grn_k

     implicit none

! external arguments
     integer, intent(in) :: w

! local variables
! loop index
     integer :: s
     integer :: k
     integer :: m
     integer :: p
     integer :: q

     integer :: cbnd
     integer :: bs, be

! dummy arrays
     complex(dp) :: T(qbnd,qbnd)
     complex(dp) :: H(qbnd)

     grn_k = czero

     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,w)
             be = kwin(k,s,2,w)
             cbnd = be - bs + 1
             print *, k - 1, s, cbnd

             do m=1,nmesh
                 T = czero
                 H = czero

                 H(1:cbnd) = czi * fmesh(m) + fermi - enk(bs:be,k,s)
                 call s_diag_z(cbnd, H(1:cbnd), T(1:cbnd,1:cbnd))

                 if (m == 1025) then
                     print *, "before:" 
                     do p=1,qbnd
                         print *, p, T(p,p)
                     enddo
                 endif
                 call s_inv_z(cbnd, T(1:cbnd,1:cbnd))
                 if (m == 1025) then
                     print *, "after:"
                     do p=1,qbnd
                         print *, p, T(p,p)
                     enddo
                 endif

                 grn_k(1:cbnd,1:cbnd,m,k,s) = T(1:cbnd,1:cbnd)
             enddo ! over m={1,nmesh} loop
             STOP
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop
 
     return
  end subroutine cal_grn_k

  subroutine cal_grn_l()
     implicit none

     return
  end subroutine cal_grn_l

  subroutine cal_hyb_l()
     implicit none

     return
  end subroutine cal_hyb_l

  subroutine cal_wss_l()
     implicit none

     return
  end subroutine cal_wss_l

  subroutine map_chi_psi()
     implicit none

     return
  end subroutine map_chi_psi

  subroutine map_psi_chi()
     implicit none

     return
  end subroutine map_psi_chi
