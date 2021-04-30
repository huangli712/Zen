!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_driver
!!!           map_chi_psi
!!!           map_psi_chi
!!! source  : dmft_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           04/30/2021 by li huang (last modified)
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
     use constants, only : mystd

     use control, only : nkpt, nspin

     use context, only : i_grp
     use context, only : ndim
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

! number of correlated orbitals for given impurity site
     integer :: cdim

! dummy array: for band dispersion (diagonal matrix)
     !!complex(dp), allocatable :: Tm(:,:)

! dummy array: for self-energy function (projected to Kohn-Sham basis)
     !!complex(dp), allocatable :: Sm(:,:)

! init cbnd and cdim
! cbnd will be k-dependent. it will be updated later
     cbnd = 0
     cdim = ndim(t)

! reset grn_l
     grn_l(:,:,:,:,t) = czero

! allocate memory for Gm
     !!allocate(Gm(cdim,cdim), stat = istat)
     !!if ( istat /= 0 ) then
     !!    call s_print_error('cal_grn_l','can not allocate enough memory')
     !!endif ! back if ( istat /= 0 ) block



     write(mystd,'(2X,a,i4)') 'calculate grn_l for site:', t
     write(mystd,'(2X,a)')  'add contributions from ...'

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt

! evaluate band window for the current k-point and spin
             bs = kwin(k,s,1,i_grp(t))
             be = kwin(k,s,2,i_grp(t))
             cbnd = be - bs + 1

! provide some information
             write(mystd,'(4X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)') 'window: ', bs, be, cbnd

! downfolding: Tm (Kohn-Sham basis) -> Gm (local basis)
                 !!call map_psi_chi(cbnd, cdim, k, s, t, Tm, Gm)

! save the results
                 !!grn_l(1:cdim,1:cdim,m,s,t) = grn_l(1:cdim,1:cdim,m,s,t) + Gm


         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! renormalize local green's function
     grn_l = grn_l / float(nkpt)

     !!do s=1,cdim
     !!    print *, s, grn_l(s,s,1,1,1)
     !!enddo

! deallocate memory
     !!deallocate(Gm)

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
!! @sub cal_sk_gk
!!
  subroutine cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)
     use constants, only : dp
     use constants, only : czi

     use control, only : axis
     use control, only : nmesh
     use control, only : fermi

     use context, only : enk
     use context, only : fmesh

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! band window: start index and end index for bands
     integer, intent(in) :: bs, be

     integer, intent(in) :: k
     integer, intent(in) :: s

     complex(dp), intent(in)  :: Sk(cbnd,cbnd,nmesh)
     complex(dp), intent(out) :: Gk(cbnd,cbnd,nmesh)

! local variables
! loop index for frequency mesh
     integer :: m

! status flag
     integer :: istat

! dummy array: for band dispersion (vector)
     complex(dp), allocatable :: Em(:), Hm(:)

! dummy array: for local green's function 
     complex(dp), allocatable :: Gm(:,:)

! allocate memory for Em, Hm, and Gm
     allocate(Em(cbnd),      stat = istat)
     allocate(Hm(cbnd),      stat = istat)
     allocate(Gm(cbnd,cbnd), stat = istat)

! evaluate Em, which is k-dependent, but frequency-independent
     Em = fermi - enk(bs:be,k,s)

     FREQ_LOOP: do m=1,nmesh

! consider imaginary axis or real axis
         if ( axis == 1 ) then
             Hm = czi * fmesh(m) + Em
         else
             Hm = fmesh(m) + Em
         endif

! convert Hm (vector) to Gm (diagonal matrix)
         call s_diag_z(cbnd, Hm, Gm)

! substract self-energy function from the hamiltonian
         Gk(:,:,m) = Gm - Sk(:,:,m)

! calculate lattice green's function
         call s_inv_z(cbnd, Gk(:,:,m))

     enddo FREQ_LOOP ! over m={1,nmesh} loop

! deallocate memory
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)
     if ( allocated(Gm) ) deallocate(Gm)

     return
  end subroutine cal_sk_gk

!!
!! @sub cal_sl_sk
!!
  subroutine cal_sl_sk(cdim, cbnd, k, s, t, Sk)
     use constants, only : dp
     use constants, only : czero

     use control, only : nmesh

     use context, only : sigdc, sig_l

     implicit none

! external arguments
     integer, intent(in) :: cdim
     integer, intent(in) :: cbnd
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(out) :: Sk(cbnd,cbnd,nmesh)

! local variables
     integer :: m
     integer :: istat

     complex(dp), allocatable :: Sl(:,:,:)

     allocate(Sl(cdim,cdim,nmesh), stat = istat)

! here we use Sl to save sig_l - sigdc
     do m=1,nmesh
         Sl(:,:,m) = sig_l(1:cdim,1:cdim,m,s,t) - sigdc(1:cdim,1:cdim,s,t)
     enddo

! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call map_chi_psi(cdim, cbnd, nmesh, k, s, t, Sl, Sk)

     deallocate(Sl)

     return
  end subroutine cal_sl_sk

!!
!! @sub map_chi_psi
!!
  subroutine map_chi_psi(cdim, cbnd, cmsh, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : psichi
     use context, only : chipsi

     implicit none

! external arguments
     integer, intent(in) :: cdim
     integer, intent(in) :: cbnd
     integer, intent(in) :: cmsh
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(in)  :: Mc(cdim,cdim,cmsh)
     complex(dp), intent(out) :: Mp(cbnd,cbnd,cmsh)

! local variables
     integer :: f

     complex(dp) :: Pc(cdim,cbnd)
     complex(dp) :: Cp(cbnd,cdim)

     Pc = psichi(1:cdim,1:cbnd,k,s,i_grp(t))
     Cp = chipsi(1:cbnd,1:cdim,k,s,i_grp(t))

     do f=1,cmsh
         Mp(:,:,f) = matmul( matmul( Cp, Mc(:,:,f) ), Pc )
     enddo ! over f={1,cmsh} loop

     return
  end subroutine map_chi_psi

!!
!! @sub map_psi_chi
!!
  subroutine map_psi_chi(cbnd, cdim, cmsh, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : psichi
     use context, only : chipsi

     implicit none

! external arguments
     integer, intent(in) :: cbnd
     integer, intent(in) :: cdim
     integer, intent(in) :: cmsh
     integer, intent(in) :: k
     integer, intent(in) :: s
     integer, intent(in) :: t

     complex(dp), intent(in)  :: Mp(cbnd,cbnd,cmsh)
     complex(dp), intent(out) :: Mc(cdim,cdim,cmsh)

! local variables
     integer :: f

     complex(dp) :: Pc(cdim,cbnd)
     complex(dp) :: Cp(cbnd,cdim)

     Pc = psichi(1:cdim,1:cbnd,k,s,i_grp(t))
     Cp = chipsi(1:cbnd,1:cdim,k,s,i_grp(t))

     do f=1,cmsh
         Mc(:,:,f) = matmul( matmul( Pc, Mp(:,:,f) ), Cp )
     enddo ! over f={1,cmsh} loop

     return
  end subroutine map_psi_chi
