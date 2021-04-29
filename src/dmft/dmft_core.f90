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

     !call cal_grn_k(1)
     call cal_grn_l(1)

     return
  end subroutine dmft_driver

  subroutine cal_grn_l(t)
     use constants, only : dp
     use constants, only : czero, czi

     use control, only : nkpt, nspin, nmesh, fermi

     use context, only : kwin
     use context, only : grn_l
     use context, only : qbnd, qdim, ndim
     use context, only : enk, fmesh

     implicit none

! external arguments
     integer, intent(in) :: t

! loop index
     integer :: s
     integer :: k
     integer :: m

     integer :: cbnd, cdim
     integer :: bs, be

     !complex(dp) :: Gm(qdim,qdim)
     !complex(dp) :: Tm(qbnd,qbnd)
     !complex(dp) :: Hm(qbnd)

     !Gm = czero
     !Tm = czero
     !Hm = czero

     complex(dp), allocatable :: Gm(:,:)
     complex(dp), allocatable :: Tm(:,:)
     complex(dp), allocatable :: Hm(:)

     grn_l(:,:,:,:,t) = czero
     cbnd = 0
     cdim = ndim(t)
     allocate(Gm(cdim,cdim))

     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,t)
             be = kwin(k,s,2,t)
             cbnd = be - bs + 1

             allocate(Hm(cbnd))
             allocate(Tm(cbnd,cbnd))
             do m=1,nmesh

                 Hm = czi * fmesh(m) + fermi - enk(bs:be,k,s)
                 call s_diag_z(cbnd, Hm, Tm)
                 !Hm(1:cbnd) = czi * fmesh(m) + fermi - enk(bs:be,k,s)
                 !call s_diag_z(cbnd, Hm(1:cbnd), Tm(1:cbnd,1:cbnd))

! add self-energy function here

                 call s_inv_z(cbnd, Tm)

                 call map_psi_chi(cbnd, cdim, k, s, t, Tm, Gm)
                 grn_l(1:cdim,1:cdim,m,s,t) = grn_l(1:cdim,1:cdim,m,s,t) + Gm
             enddo ! over m={1,nmesh} loop
             deallocate(Hm)
             deallocate(Tm)
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     grn_l = grn_l / float(nkpt)
     print *, grn_l(1:ndim(t),1:ndim(t), 1, 1, t) 
     deallocate(Gm)

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

!!
!! @sub map_psi_chi
!!
  subroutine map_psi_chi(cbnd, cdim, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : psichi, chipsi

     implicit none

! external arguments
     integer, intent(in) :: cbnd
     integer, intent(in) :: cdim
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(in) :: Mp(cbnd,cbnd)
     complex(dp), intent(out) :: Mc(cdim,cdim)

     Mc = matmul( matmul(psichi(1:cdim,1:cbnd,k,s,t), Mp), &
                         chipsi(1:cbnd,1:cdim,k,s,t) )

     return
  end subroutine map_psi_chi
