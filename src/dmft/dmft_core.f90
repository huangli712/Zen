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

     call cal_grn_l(1)

     return
  end subroutine dmft_driver

!!
!! @sub cal_grn_l
!!
!! try to calculate local green's function for given impurity site
!!
  subroutine cal_grn_l(t)
     use constants, only : dp
     use constants, only : czero, czi

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : fermi

     use context, only : i_grp
     use context, only : ndim
     use context, only : kwin
     use context, only : enk
     use context, only : fmesh
     use context, only : sig_l, sigdc
     use context, only : grn_l

     implicit none

! external arguments
! index for impurity sites
     integer, intent(in) :: t

! local variables
! loop index for spin
     integer :: s

! loop index for k-points
     integer :: k

! loop index for frequency mesh
     integer :: m

     integer :: cdim
     integer :: cbnd
     integer :: bs, be

     complex(dp), allocatable :: Gm(:,:)
     complex(dp), allocatable :: Tm(:,:)
     complex(dp), allocatable :: Hm(:)
     complex(dp), allocatable :: Sm(:,:)

     grn_l(:,:,:,:,t) = czero
     cbnd = 0
     cdim = ndim(t)

     sigdc = czero

     allocate(Gm(cdim,cdim))

     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,i_grp(t))
             be = kwin(k,s,2,i_grp(t))
             cbnd = be - bs + 1
             print *, s, k, t, cbnd, cdim

             allocate(Hm(cbnd))
             allocate(Tm(cbnd,cbnd))
             allocate(Sm(cbnd,cbnd))

             do m=1,nmesh

                 Hm = czi * fmesh(m) + fermi - enk(bs:be,k,s)
                 call s_diag_z(cbnd, Hm, Tm)

! add self-energy function here
                 Gm = sig_l(1:cdim,1:cdim,m,s,t) - sigdc(1:cdim,1:cdim,s,t)
                 if ( m == 1 ) then
                     Gm(1,1) = dcmplx(1.200, 5.0)
                     Gm(2,2) = dcmplx(1.000, -0.14) 
                     Gm(3,3) = dcmplx(2.200, 3.0)
                     Gm(4,4) = dcmplx(0.800, -0.1)
                     Gm(5,5) = dcmplx(1.255, 2.0_dp)

                     Gm(2,4) = dcmplx(-1.0, 0.34)
                     Gm(1,4) = dcmplx(-1.0, 0.34)
                     Gm(3,1) = dcmplx(-1.0, 0.34)
                     Gm(1,5) = dcmplx(-1.0, 0.34)
                 endif

                 call map_chi_psi(cdim, cbnd, k, s, t, Gm, Sm)

                 Tm = Tm - Sm

                 call s_inv_z(cbnd, Tm)

                 call map_psi_chi(cbnd, cdim, k, s, t, Tm, Gm)

                 grn_l(1:cdim,1:cdim,m,s,t) = grn_l(1:cdim,1:cdim,m,s,t) + Gm
             enddo ! over m={1,nmesh} loop

             deallocate(Hm)
             deallocate(Tm)
             deallocate(Sm)
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     grn_l = grn_l / float(nkpt)
     do s=1,ndim(t)
         print *, s,grn_l(s,s, 1, 1, t)
     enddo

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

!!
!! @sub map_chi_psi
!!
  subroutine map_chi_psi(cdim, cbnd, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : psichi
     use context, only : chipsi

     implicit none

! external arguments
     integer, intent(in) :: cdim
     integer, intent(in) :: cbnd
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(in)  :: Mc(cdim,cdim)
     complex(dp), intent(out) :: Mp(cbnd,cbnd)

     Mp = matmul( matmul(chipsi(1:cbnd,1:cdim,k,s,i_grp(t)), Mc), &
                         psichi(1:cdim,1:cbnd,k,s,i_grp(t)) )

     return
  end subroutine map_chi_psi

!!
!! @sub map_psi_chi
!!
  subroutine map_psi_chi(cbnd, cdim, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : psichi
     use context, only : chipsi

     implicit none

! external arguments
     integer, intent(in) :: cbnd
     integer, intent(in) :: cdim
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(in)  :: Mp(cbnd,cbnd)
     complex(dp), intent(out) :: Mc(cdim,cdim)

     Mc = matmul( matmul(psichi(1:cdim,1:cbnd,k,s,i_grp(t)), Mp), &
                         chipsi(1:cbnd,1:cdim,k,s,i_grp(t)) )

     return
  end subroutine map_psi_chi
