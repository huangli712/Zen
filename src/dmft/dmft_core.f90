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

  subroutine dmft_driver()
     use constants, only : dp, czero

     use control, only : nspin
     use control, only : nkpt, nband
     use control, only : nmesh

     use context, only : kwin, enk, qbnd
     use context, only : grn_k

     implicit none

! loop index
     integer :: s
     integer :: k
     integer :: m
     integer :: p
     integer :: q

     integer :: cbnd, bs, be

! dummy array
     complex(dp) :: T(qbnd,qbnd)
     complex(dp) :: H(qbnd,qbnd)
     complex(dp) :: hopping(nband)

     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,1)
             be = kwin(k,s,2,1)
             cbnd = be - bs + 1
             print *, k - 1, s, cbnd
             do m=1,nmesh
                 T = czero
                 H = czero
                 hopping = czero
                 hopping(bs:be) = enk(bs:be,k,s)
                 call s_diag_z(cbnd, hopping(bs:be), H(1:cbnd,1:cbnd))
                 T(1:cbnd, 1:cbnd) = -H(1:cbnd,1:cbnd)
                 if (m == 1) then
                     print *, "before:" 
                     do p=1,qbnd
                         print *, p, T(p,p)
                     enddo
                 endif
                 call s_inv_z(cbnd, T(1:cbnd,1:cbnd))
                 if (m == 1) then
                     print *, "after:"
                     do p=1,qbnd
                         print *, p, T(p,p)
                     enddo
                 endif
                 grn_k(1:cbnd,1:cbnd,m,k,s) = T(1:cbnd,1:cbnd)
             enddo ! over m={1,nmesh} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop
 
     return
  end subroutine dmft_driver

  subroutine map_chi_psi()
     implicit none

     return
  end subroutine map_chi_psi

  subroutine map_psi_chi()
     implicit none

     return
  end subroutine map_psi_chi
