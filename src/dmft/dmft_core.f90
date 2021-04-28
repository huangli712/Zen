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

     call cal_grn_k(1)
     call cal_grn_l(1)

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
             do m=1,nmesh
                 T = czero
                 H = czero

                 H(1:cbnd) = czi * fmesh(m) + fermi - enk(bs:be,k,s)
                 call s_diag_z(cbnd, H(1:cbnd), T(1:cbnd,1:cbnd))

! add self-energy function here

                 call s_inv_z(cbnd, T(1:cbnd,1:cbnd))

                 grn_k(1:cbnd,1:cbnd,m,k,s) = T(1:cbnd,1:cbnd)
             enddo ! over m={1,nmesh} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop
 
     return
  end subroutine cal_grn_k

  subroutine cal_grn_l(t)
     use constants, only : dp
     use constants, only : czero

     use control, only : nkpt, nspin, nmesh

     use context, only : kwin
     use context, only : grn_l, grn_k
     use context, only : qbnd, qdim, ndim
     use context, only : psichi

     implicit none

! external arguments
     integer, intent(in) :: t

! loop index
     integer :: s
     integer :: k
     integer :: m

     integer :: cbnd, cdim
     integer :: bs, be

     complex(dp) :: G(qbnd,qbnd)
     complex(dp) :: P(qdim,qbnd)

     grn_l(:,:,:,:,t) = czero

     do s=1,nspin
         do k=1,nkpt
             !!print *, s, k
             bs = kwin(k,s,1,1)
             be = kwin(k,s,2,1)
             cbnd = be - bs + 1
             cdim = ndim(t)

             P = czero
             P(1:cdim,1:cbnd) = psichi(1:cdim,1:cbnd,k,s,t)
             do m=1,nmesh
                 G = czero
                 G(1:cbnd,1:cbnd) = grn_k(1:cbnd,1:cbnd,m,k,s)
                 grn_l(1:cdim,1:cdim,m,s,t) = grn_l(1:cdim,1:cdim,m,s,t) + &
                 matmul(matmul(P(1:cdim,1:cbnd), G(1:cbnd,1:cbnd)), transpose(dconjg(P(1:cdim,1:cbnd))))
             enddo ! over m={1,nmesh} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     grn_l = grn_l / float(nkpt)
     print *, grn_l(1:ndim(t),1:ndim(t), 1, 1, t) 

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
