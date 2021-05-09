!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_driver
!!!           dmft_try0
!!!           dmft_try1
!!!           dmft_try2
!!!           cal_fermi
!!!           cal_eimps
!!!           cal_grn_l
!!!           cal_wss_l
!!!           cal_hyb_l
!!!           cal_sl_sk
!!!           cal_sk_hk
!!!           cal_hk_ek
!!!           cal_sl_so
!!!           cal_so_ho
!!!           cal_ho_eo
!!!           cal_sk_so
!!!           cal_sk_gk
!!!           cal_gk_gl
!!!           dichotomy
!!!           cal_nelect
!!!           cal_occupy
!!!           cal_eigsys
!!!           fermi_dirac
!!!           map_chi_psi
!!!           map_psi_chi
!!!           one_chi_psi
!!!           one_psi_chi
!!! source  : dmft_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           05/08/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> driver subroutines: layer 1                                      <<<
!!========================================================================

!!
!! @sub dmft_driver
!!
!! core subroutine, dispatch the computational task
!!
  subroutine dmft_driver()
     use control, only : task

     implicit none

     DISPATCHER: select case ( task )

! task = 0, search the fermi level
         case (0)
             call dmft_try0()

! task = 1, calculate the local green's function
         case (1)
             call dmft_try1()

! task = 2, calculate density correction
         case (2)
             call dmft_try2()

         case default
             call s_print_error('dmft_driver','this feature is not supported')

     end select DISPATCHER

     return
  end subroutine dmft_driver

!!========================================================================
!!>>> driver subroutines: layer 2                                      <<<
!!========================================================================

!!
!! @sub dmft_try0
!!
!! to determine the fermi level, the global variable `fermi` may be
!! updated in this subroutine
!!
  subroutine dmft_try0()
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : myid, master

     implicit none

! check lfermi
     call s_assert2(lfermi .eqv. .true., 'lfermi must be true')

! call the computational subroutine to do this job
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     call cal_fermi()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try0

!!
!! @sub dmft_try1
!!
!! to calculate the local green's function, generate key inputs for the
!! quantum impurity solvers. the fermi level may be updated, depending
!! on the configuration parameter  
!!
  subroutine dmft_try1()
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : nsite
     use control, only : myid, master

     use context, only : grn_l

     implicit none

! local variables
! loop index for impurity sites
     integer :: t

! call the computational subroutine to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     if ( lfermi .eqv. .true. ) then
         call cal_fermi()
     endif ! back if ( lfermi .eqv. .true. ) block
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! call the computational subroutine to compute the local impurity levels
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Level'
     endif ! back if ( myid == master ) block
     !
     do t=1,nsite
         call cal_eimps(t)
     enddo ! over t={1,nsite} loop
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! call the computational subroutine to compute the local green's function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Green'
     endif ! back if ( myid == master ) block
     !
     do t=1,nsite
         call cal_grn_l(t)
     enddo ! over t={1,nsite} loop
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! call the computational subroutine to compute the local hybridization function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Hybri'
     endif ! back if ( myid == master ) block
     !
     do t=1,nsite
         call cal_hyb_l(t)
     enddo ! over t={1,nsite} loop
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save grn_l...'
         call dmft_dump_grn_l(grn_l)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try1

!!
!! @sub dmft_try2
!!
  subroutine dmft_try2()
     implicit none

     return
  end subroutine dmft_try2

!!========================================================================
!!>>> driver subroutines: layer 3                                      <<<
!!========================================================================

!!
!! @sub cal_fermi
!!
!! try to determine the fermi level
!!
  subroutine cal_fermi()
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : myid, master

     use context, only : qbnd

     implicit none

! local variables
! status flag
     integer  :: istat

! desired charge density
     real(dp) :: ndens

! dummy arrays, used to save the eigenvalues of H + \Sigma
     complex(dp), allocatable :: eigs(:,:,:,:)
     complex(dp), allocatable :: einf(:,:,:)

! allocate memory
     allocate(eigs(qbnd,nmesh,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_fermi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(einf(qbnd,nkpt,nspin),       stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_fermi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! calculate the nominal charge density according to the dft eigenvalues
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculating desired charge density'
     endif ! back if ( myid == master ) block
     !
     call cal_nelect(ndens)

! construct H + \Sigma, diagonalize it to obtain the dft + dmft eigenvalues
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculating dft + dmft eigenvalues'
     endif ! back if ( myid == master ) block
     !
     call cal_eigsys(eigs, einf)

! search the fermi level using bisection algorithm
! the global variable `fermi` will be updated within `dichotomy()`
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'searching fermi level'
     endif ! back if ( myid == master ) block
     !
     call dichotomy(ndens, eigs, einf)

! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

     return
  end subroutine cal_fermi

!!
!! @sub cal_eimps
!!
!! try to calculate local energy levels for given impurity site
!!
  subroutine cal_eimps(t)
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nkpt, nspin
     use control, only : myid, master

     use context, only : i_wnd
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
     use context, only : enk
     use context, only : eimps

     implicit none

! external arguments
! index for impurity sites
     integer, intent(in) :: t

! local variables
! loop index for spin
     integer :: s

! loop index for k-points
     integer :: k

! number of dft bands for given k-point and spin
     integer :: cbnd

! number of correlated orbitals for given impurity site
     integer :: cdim

! band window: start index and end index for bands
     integer :: bs, be

! status flag
     integer :: istat

! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

     complex(dp), allocatable :: Eimp(:,:)

! init cbnd and cdim
! cbnd will be k-dependent. it will be updated later
     cbnd = 0
     cdim = ndim(t)

! allocate memory
     allocate(Eimp(cdim,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! reset eimps
     eimps(:,:,:,t) = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,i4)') 'calculate eimps for site:', t
         write(mystd,'(4X,a)')  'add contributions from ...'
     endif ! back if ( myid == master ) block

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

! determine cbnd
             cbnd = be - bs + 1

! provide some useful information
             if ( myid == master ) then
                 write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
                 write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
                 write(mystd,'(2X,a,3i3)') 'window: ', bs, be, cbnd
             endif ! back if ( myid == master ) block

! allocate memory
             allocate(Em(cbnd),      stat = istat)
             allocate(Hm(cbnd,cbnd), stat = istat)
             if ( istat /= 0 ) then
                 call s_print_error('cal_eimps','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

! evaluate Em, which is just some dft eigenvalues 
             Em = enk(bs:be,k,s)

! convert `Em` to diagonal matrix `Hm`
             call s_diag_z(cbnd, Em, Hm)

! project hamiltonian to local basis
             call one_psi_chi(cbnd, cdim, k, s, t, Hm, Eimp)

! save the final results
             eimps(1:cdim,1:cdim,s,t) = eimps(1:cdim,1:cdim,s,t) + Eimp * weight(k)

! deallocate memory
             if ( allocated(Em) ) deallocate(Em)
             if ( allocated(Hm) ) deallocate(Hm)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! renormalize impurity levels
     eimps(:,:,:,t) = eimps(:,:,:,t) / float(nkpt)

! deallocate memory
     if ( allocated(Eimp) ) deallocate(Eimp)

     return
  end subroutine cal_eimps

!!
!! @sub cal_grn_l
!!
!! try to calculate local green's function for given impurity site
!!
  subroutine cal_grn_l(t)
     use constants, only : dp, mystd
     use constants, only : czero, czi

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : myid, master

     use context, only : i_wnd
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
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

! number of dft bands for given k-point and spin
     integer :: cbnd

! number of correlated orbitals for given impurity site
     integer :: cdim

! band window: start index and end index for bands
     integer :: bs, be

! status flag
     integer :: istat

! dummy array: for self-energy function (projected to Kohn-Sham basis)
     complex(dp), allocatable :: Sk(:,:,:)

! dummy array: for lattice green's function
     complex(dp), allocatable :: Gk(:,:,:)

! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:,:)

! init cbnd and cdim
! cbnd will be k-dependent. it will be updated later
     cbnd = 0
     cdim = ndim(t)

! allocate memory for Gl
     allocate(Gl(cdim,cdim,nmesh), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_grn_l','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! reset grn_l
     grn_l(:,:,:,:,t) = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,i4)') 'calculate grn_l for site:', t
         write(mystd,'(4X,a)')  'add contributions from ...'
     endif ! back if ( myid == master ) block

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

! determine cbnd
             cbnd = be - bs + 1

! provide some useful information
             if ( myid == master ) then
                 write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
                 write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
                 write(mystd,'(2X,a,3i3)') 'window: ', bs, be, cbnd
             endif ! back if ( myid == master ) block

! allocate memories for Sk and Gk. their sizes are k-dependent
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Gk(cbnd,cbnd,nmesh), stat = istat)
             if ( istat /= 0 ) then
                 call s_print_error('cal_grn_l','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

! build self-energy function, and then embed it into Kohn-Sham basis
             call cal_sl_sk(cdim, cbnd, k, s, t, Sk)

! calculate lattice green's function
             call cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)

! project lattice green's function to obtain local green's function
             call cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl)

! save the final results
             grn_l(1:cdim,1:cdim,:,s,t) = grn_l(1:cdim,1:cdim,:,s,t) + Gl * weight(k)

! deallocate memories
             if ( allocated(Sk) ) deallocate(Sk)
             if ( allocated(Gk) ) deallocate(Gk)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! renormalize local green's function
     grn_l(:,:,:,:,t) = grn_l(:,:,:,:,t) / float(nkpt)

! deallocate memory
     if ( allocated(Gl) ) deallocate(Gl)

!! DEBUG CODE
!<     do s=1,cdim
!<         print *, s, grn_l(s,s,1,1,1)
!<     enddo ! over s={1,cdim} loop
!! DEBUG CODE

     return
  end subroutine cal_grn_l

!!
!! @sub cal_wss_l
!!
!! try to calculate local weiss's function for given impurity site
!!
  subroutine cal_wss_l()
     implicit none

     return
  end subroutine cal_wss_l

!!
!! @sub cal_hyb_l
!!
!! try to calculate local hybridization function for given impurity site
!!
  subroutine cal_hyb_l(t)
     use constants, only : dp, mystd
     use constants, only : czi, czero

     use control, only : nspin
     use control, only : nmesh
     use control, only : fermi
     use control, only : myid, master

     use context, only : ndim
     use context, only : fmesh
     use context, only : sigdc, sig_l
     use context, only : grn_l
     use context, only : hyb_l
 
     implicit none

! external arguments
! index for impurity sites
     integer, intent(in) :: t

! local variables
! loop index for spin
     integer :: s

! loop index for frequency mesh
     integer :: m

! number of correlated orbitals for given impurity site
     integer :: cdim
     integer :: istat

     complex(dp) :: caux
     complex(dp), allocatable :: Tm(:,:)

     cdim = ndim(t)
     allocate(Tm(cdim,cdim), stat = istat)

     hyb_l(:,:,:,:,t) = czero

     SPIN_LOOP: do s=1,nspin
         MESH_LOOP: do m=1,nmesh
             caux = czi * fmesh(m) + fermi
             Tm = grn_l(1:cdim,1:cdim,m,s,t)
             call s_inv_z(cdim, Tm)
             Tm = caux - eimps(1:cdim,1:cdim,s,t) - sig_l(1:cdim,1:cdim,m,s,t) - Tm
             hyb_l(1:cdim,1:cdim,m,s,t) = Tm
         enddo MESH_LOOP ! over m={1,nmesh} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

     deallocate(Tm)

     return
  end subroutine cal_hyb_l

!!========================================================================
!!>>> service subroutines: set 1                                       <<<
!!========================================================================

!!
!! @sub cal_sl_sk
!!
!! try to substract the double-counting term from the local self-energy
!! function, and then map it from local basis to Kohn-Sham basis
!!
  subroutine cal_sl_sk(cdim, cbnd, k, s, t, Sk)
     use constants, only : dp
     use constants, only : czero

     use control, only : nmesh

     use context, only : sigdc, sig_l

     implicit none

! external arguments
! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! self-energy function in Kohn-Sham basis
     complex(dp), intent(out) :: Sk(cbnd,cbnd,nmesh)

! local variables
! loop index for frequency mesh
     integer :: m

! status flag
     integer :: istat

! dummy array: for local self-energy function
     complex(dp), allocatable :: Sl(:,:,:)

! allocate memory
     allocate(Sl(cdim,cdim,nmesh), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_sk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! here we use Sl to save sig_l - sigdc
     do m=1,nmesh
         Sl(:,:,m) = sig_l(1:cdim,1:cdim,m,s,t) - sigdc(1:cdim,1:cdim,s,t)
     enddo ! over m={1,nmesh} loop

! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call map_chi_psi(cdim, cbnd, nmesh, k, s, t, Sl, Sk)

! deallocate memory
     if ( allocated(Sl) ) deallocate(Sl)

     return
  end subroutine cal_sl_sk

!!
!! @sub cal_sk_hk
!!
!! try to build H(k) + \Sigma(i\omega_n)
!!
  subroutine cal_sk_hk(cbnd, bs, be, k, s, Sk, Hk)
     use constants, only : dp

     use control, only : nmesh

     use context, only : enk

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! band window: start index and end index for bands
     integer, intent(in) :: bs, be

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! self-energy function at Kohn-Sham basis
     complex(dp), intent(in)  :: Sk(cbnd,cbnd,nmesh)

! lattice green's function at given k-point and spin
     complex(dp), intent(out) :: Hk(cbnd,cbnd,nmesh)

! local variables
! loop index for frequency mesh
     integer :: m

! status flag
     integer :: istat

! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

! allocate memory
     allocate(Em(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_hk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Hm(cbnd,cbnd), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_hk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! evaluate Em, which is just some dft eigenvalues 
     Em = enk(bs:be,k,s)

! convert `Em` to diagonal matrix `Hm`
     call s_diag_z(cbnd, Em, Hm)

! combine `Hm` and `Sk` to build the effective hamiltonian 
     FREQ_LOOP: do m=1,nmesh
         Hk(:,:,m) = Hm + Sk(:,:,m)
     enddo FREQ_LOOP ! over m={1,nmesh} loop

! deallocate memory
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

     return
  end subroutine cal_sk_hk

!!
!! @sub cal_hk_ek
!!
!! try to diagonalize H(k) + \Sigma(i\omega_n), get all the eigenvalues
!!
  subroutine cal_hk_ek(cbnd, Hk, Ek)
     use constants, only : dp

     use control, only : nmesh

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! effective hamiltonian: H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in)  :: Hk(cbnd,cbnd,nmesh)

! resulting eigenvalues
     complex(dp), intent(out) :: Ek(cbnd,nmesh)

! local variables
! loop index for frequency mesh
     integer :: m

     FREQ_LOOP: do m=1,nmesh
         call s_eigvals_zg(cbnd, cbnd, Hk(:,:,m), Ek(:,m))
     enddo FREQ_LOOP ! over m={1,nmesh} loop

     return
  end subroutine cal_hk_ek

!!========================================================================
!!>>> service subroutines: set 2                                       <<<
!!========================================================================

!!
!! @sub cal_sl_so
!!
!! try to substract the double-counting term from the local self-energy
!! function. and then evaluate its asymptotic values at \omega = \infty.
!! finally, map it from local basis to Kohn-Sham basis
!!
  subroutine cal_sl_so(cdim, cbnd, k, s, t, So)
     use constants, only : dp
     use constants, only : czero

     use control, only : nmesh

     use context, only : sigdc, sig_l

     implicit none

! external arguments
! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! self-energy function at \omega = \infty in Kohn-Sham basis
     complex(dp), intent(out) :: So(cbnd,cbnd)

! local parameters
! how many frequency points are included to calculate the asymptotic
! values of self-energy function
     integer, parameter :: mcut = 16

! local variables
! loop index for frequency mesh
     integer :: m

! status flag
     integer :: istat

! dummy array: for local self-energy function
     complex(dp), allocatable :: Sl(:,:,:)

! dummy array: for lattice self-energy function
     complex(dp), allocatable :: Sk(:,:,:)

! allocate memory
! the last elements of Sl and Sk are used to store the averaged values
     allocate(Sl(cdim,cdim,mcut+1), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_so','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Sk(cbnd,cbnd,mcut+1), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_so','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! here we use Sl to save sig_l - sigdc
     do m=1,mcut
         Sl(:,:,m) = sig_l(1:cdim,1:cdim,nmesh+1-m,s,t) - sigdc(1:cdim,1:cdim,s,t)
     enddo ! over m={1,mcut} loop

! then the averaged values are stored at Sl(:,:,mcut+1)
     Sl(:,:,mcut+1) = czero
     do m=1,mcut
         Sl(:,:,mcut+1) = Sl(:,:,mcut+1) + Sl(:,:,m) / real(mcut)
     enddo ! over m={1,mcut} loop

! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call map_chi_psi(cdim, cbnd, mcut + 1, k, s, t, Sl, Sk)

! Sk(:,:,mcut + 1) is what we want, copy it to `So`
     So = Sk(:,:,mcut+1)

! deallocate memory
     if ( allocated(Sl) ) deallocate(Sl)
     if ( allocated(Sk) ) deallocate(Sk)

     return
  end subroutine cal_sl_so

!!
!! @sub cal_so_ho
!!
!! try to build H(k) + \Sigma(\infty)
!!
  subroutine cal_so_ho(cbnd, bs, be, k, s, So, Ho)
     use constants, only : dp

     use context, only : enk

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! band window: start index and end index for bands
     integer, intent(in) :: bs, be

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! self-energy function at Kohn-Sham basis, \omega = \infty
     complex(dp), intent(in)  :: So(cbnd,cbnd)

! lattice green's function at given k-point and spin
     complex(dp), intent(out) :: Ho(cbnd,cbnd)

! local variables
! status flag
     integer :: istat

! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

! allocate memory
     allocate(Em(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_so_ho','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Hm(cbnd,cbnd), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_so_ho','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! evaluate Em, which is just some dft eigenvalues 
     Em = enk(bs:be,k,s)

! convert `Em` to diagonal matrix `Hm`
     call s_diag_z(cbnd, Em, Hm)

! combine `Hm` and `So` to build the effective hamiltonian 
     Ho = Hm + So

! deallocate memory
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

     return
  end subroutine cal_so_ho

!!
!! @sub cal_ho_eo
!!
!! try to diagonalize H(k) + \Sigma(\infty), get all the eigenvalues
!!
  subroutine cal_ho_eo(cbnd, Ho, Eo)
     use constants, only : dp

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! effective hamiltonian: H(k) + \Sigma(\infty)
     complex(dp), intent(in)  :: Ho(cbnd,cbnd)

! resulting eigenvalues
     complex(dp), intent(out) :: Eo(cbnd)

     call s_eigvals_zg(cbnd, cbnd, Ho, Eo)

     return
  end subroutine cal_ho_eo

!!
!! @sub cal_sk_so
!!
!! try to evaluate \Sigma(i\omega_n \to \infty). it's function is similar
!! to `cal_sl_so()`
!!
  subroutine cal_sk_so(cbnd, Sk, So)
     use constants, only : dp
     use constants, only : czero

     use control, only : nmesh

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! self-energy function at Kohn-Sham basis
     complex(dp), intent(in)  :: Sk(cbnd,cbnd,nmesh)

! asymptotic values of self-energy function at Kohn-Sham basis 
     complex(dp), intent(out) :: So(cbnd,cbnd)

! local parameters
! how many frequency points are included to calculate the asymptotic
! values of self-energy function
     integer, parameter :: mcut = 16

! local variables
! loop index for frequency mesh
     integer :: m

! count the final `mcut` frequency points, and calculate the averaged value
     So = czero
     !
     do m=1,mcut
         So = So + Sk(:,:,nmesh + 1 - m)
     enddo ! over m={1,mcut} loop
     !
     So = So / real(mcut)

     return
  end subroutine cal_sk_so

!!========================================================================
!!>>> service subroutines: set 3                                       <<<
!!========================================================================

!!
!! @sub cal_sk_gk
!!
!! try to calculate lattice green's function at given k-point and spin.
!! note that this lattice green's function is not the nominal one. it is
!! connected with the impurity site. in other words, it is the lattice
!! green's function for the given impurity site. this subroutine needs
!! the self-energy function at Kohn-Sham basis (Sk), that is the reason
!! why it is called `cal_sk_gk`
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

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! self-energy function at Kohn-Sham basis
     complex(dp), intent(in)  :: Sk(cbnd,cbnd,nmesh)

! lattice green's function at given k-point and spin
     complex(dp), intent(out) :: Gk(cbnd,cbnd,nmesh)

! local variables
! loop index for frequency mesh
     integer :: m

! status flag
     integer :: istat

! dummy array: for band dispersion (vector)
     complex(dp), allocatable :: Em(:)

! dummy array: for effective hamiltonian (diagonal matrix)
     complex(dp), allocatable :: Hm(:,:)

! allocate memory
     allocate(Em(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_gk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Hm(cbnd,cbnd), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_gk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! evaluate Em, which is k-dependent, but frequency-independent
! if you want to consider magnetic field, you can add your codes here
     Em = fermi - enk(bs:be,k,s)

     FREQ_LOOP: do m=1,nmesh

! consider imaginary axis or real axis
         if ( axis == 1 ) then
             Em = czi * fmesh(m) + Em
         else
             Em = fmesh(m) + Em
         endif ! back if ( axis == 1 ) block

! convert Em (vector) to Hm (diagonal matrix)
         call s_diag_z(cbnd, Em, Hm)

! substract self-energy function from the hamiltonian
         Gk(:,:,m) = Hm - Sk(:,:,m)

! calculate lattice green's function by direct inversion
         call s_inv_z(cbnd, Gk(:,:,m))

     enddo FREQ_LOOP ! over m={1,nmesh} loop

! deallocate memory
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

     return
  end subroutine cal_sk_gk

!!
!! @sub cal_gk_gl
!!
!! try to calculate local green's function from lattice green's function
!! via downfolding method
!!
  subroutine cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl)
     use constants, only : dp

     use control, only : nmesh

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! lattice green's function at given k-point and spin
     complex(dp), intent(in)  :: Gk(cbnd,cbnd,nmesh)

! local green's function (contributions from the given k-point and spin)
     complex(dp), intent(out) :: Gl(cdim,cdim,nmesh)

     call map_psi_chi(cbnd, cdim, nmesh, k, s, t, Gk, Gl)

     return
  end subroutine cal_gk_gl

!!========================================================================
!!>>> service subroutines: set 4, fermi level search                   <<<
!!========================================================================

!!
!! @sub dichotomy
!!
!! try to locate the fermi level with the bisection method
!!
  subroutine dichotomy(desired, eigs, einf)
     use constants, only : dp, mystd

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : fermi
     use control, only : mc
     use control, only : myid, master

     use context, only : qbnd

     implicit none

! external arguments
! desired charge density
     real(dp), intent(in)    :: desired

! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(qbnd,nmesh,nkpt,nspin)

! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(in) :: einf(qbnd,nkpt,nspin)

! local parameters
! maximum number of loops for the bisection algorithm
     integer, parameter  :: max_loops = 100

! step for locating the energy boundary
     real(dp), parameter :: delta = 0.5_dp

! local variables
! loop index for the bisection algorithm 
     integer  :: loop

! left boundary, right boundary, and the final result for the fermi level
     real(dp) :: mu1, mu2, mu3

! the corresponding charge density
     real(dp) :: occ1, occ2, occ3

! sign
     real(dp) :: sign

! initialization, determine mu1, mu2, occ1, occ2, and sign
! if sign < 0, it means occ1 < desired, we should push mu2 to higher
! energy. if sign > 0, it means occ1 > desired. then mu1 will be
! the right boundary, and we should push mu2 to lower energy
     mu1 = fermi
     call cal_occupy(mu1, occ1, eigs, einf)
     !
     mu2 = mu1
     occ2 = occ1
     !
     sign = abs( occ1 - desired ) / ( occ1 - desired )

! first, determin the left and right boundaries
     loop = 1
     !
     if ( myid == master ) then
         write(mystd,'(6X,a,i2)',advance = 'no') 'iter: ', loop
         write(mystd,'(2X,a,f12.8)',advance = 'no') 'EF: ', mu1
         write(mystd,'(2X,a,f12.8)',advance = 'no') 'density: ', occ1
         write(mystd,'(2X,a,f12.8)',advance = 'no') 'desired: ', desired
         write(mystd,'(2X,a,f12.8)') 'diff: ', abs(occ1 - desired)
     endif ! back if ( myid == master ) block
     !
     do while ( loop <= max_loops .and. ( occ2 - desired ) * sign > 0 .and. abs( occ2 - desired ) > mc )
         loop = loop + 1
         mu2 = mu2 - sign * delta
         call cal_occupy(mu2, occ2, eigs, einf)
         if ( myid == master ) then
             write(mystd,'(6X,a,i2)',advance = 'no') 'iter: ', loop
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'EF: ', mu2
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'density: ', occ2
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'desired: ', desired
             write(mystd,'(2X,a,f12.8)') 'diff: ', abs(occ2 - desired)
         endif ! back if ( myid == master ) block
     enddo ! over do while loop

! exchange the left and right boundaries
     if ( mu1 > mu2 ) then
         mu3 = mu1; mu1 = mu2; mu2 = mu3
         occ3 = occ1; occ1 = occ2; occ2 = occ3
     endif ! back if ( mu1 > mu2 ) block

! now the fermi level should lie in the regine [mu1, mu2]
! refine the boundary to locate the fermi level
     if ( abs(occ1 - desired) < abs(occ2 - desired) ) then
         mu3 = mu1
         occ3 = occ1
     else
         mu3 = mu2
         occ3 = occ2
     endif
     !
     do while ( loop <= max_loops .and. abs( occ3 - desired ) > mc )
         loop = loop + 1
         mu3 = mu1 + ( mu2 - mu1 ) * ( desired - occ1 ) / ( occ2 - occ1 )
         call cal_occupy(mu3, occ3, eigs, einf)
         if ( ( occ1 - desired ) * ( occ3 - desired ) > 0 ) then
             mu1 = mu3
             occ1 = occ3
         else
             mu2 = mu3
             occ2 = occ3
         endif ! back if block
         if ( myid == master ) then
             write(mystd,'(6X,a,i2)',advance = 'no') 'iter: ', loop
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'EF: ', mu3
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'density: ', occ3
             write(mystd,'(2X,a,f12.8)',advance = 'no') 'desired: ', desired
             write(mystd,'(2X,a,f12.8)') 'diff: ', abs(occ3 - desired)
         endif ! back if ( myid == master ) block
     enddo ! over do while loop
     !
     if ( abs(occ3 - desired) < mc ) then
         write(mystd,'(6X,a)',advance = 'no') 'final results ->'
         write(mystd,'(2(2X,a,f12.8))') 'EF: ', mu3, 'density: ', occ3
     else
         call s_print_error('dichotomy', 'fail to locate the fermi level')
     endif ! back if ( abs(occ3 - desired) < mc ) block

! well, finally, we have to update the global variable `fermi`
     fermi = mu3

     return
  end subroutine dichotomy

!!
!! @sub cal_nelect
!!
!! try to calculate the number of valence electrons by dft occupations.
!! actually, what we obtain is the occupation numbers in the selected
!! band window
!!
  subroutine cal_nelect(nelect)
     use constants, only : dp
     use constants, only : zero, two

     use control, only : nkpt, nspin

     use context, only : kwin
     use context, only : weight
     use context, only : occupy

     implicit none

! external arguments
! number of relevant electrons
     real(dp), intent(out) :: nelect

! local variables
! index for k-points
     integer :: k

! index for spin
     integer :: s

! band window: start index and end index for bands
     integer :: bs, be

! basically, now we only support single band window. so the last index
! for kwin is always 1
     nelect = zero
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt
             bs = kwin(k,s,1,1)
             be = kwin(k,s,2,1)
             nelect = nelect + sum( occupy(bs:be,k,s) ) * weight(k)
         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! normalize `nelect`
     nelect = nelect / float(nkpt)

! consider the spins
     if ( nspin == 1 ) then
         nelect = nelect * two
     endif ! back if ( nspin == 1 ) block

     return
  end subroutine cal_nelect

!!
!! @sub cal_occupy
!!
!! for given fermi level, try to calculate the corresponding occupations
!!
  subroutine cal_occupy(fermi, val, eigs, einf)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czi, czero

     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : beta

     use context, only : qbnd
     use context, only : kwin
     use context, only : fmesh

     implicit none

! external arguments
! assumed fermi level
     real(dp), intent(in)  :: fermi

! occupation number
     real(dp), intent(out) :: val

! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(qbnd,nmesh,nkpt,nspin)

! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(in) :: einf(qbnd,nkpt,nspin)

! local variables
! loop index for bands
     integer  :: b

! loop index for k-points
     integer  :: k

! loop index for spins
     integer  :: s

! loop index for frequency mesh
     integer  :: m

! band window: start index and end index for bands
     integer  :: bs, be

! number of dft bands for given k-point and spin
     integer  :: cbnd

! status flag
     integer  :: istat

! complex(dp) dummy variable
     complex(dp) :: caux

! density matrix
     complex(dp), allocatable :: zocc(:,:)

! local green's function
     complex(dp), allocatable :: gloc(:,:,:)

! external functions
! used to calculate fermi-dirac function
     real(dp), external :: fermi_dirac

! allocate memory
     allocate(zocc(qbnd,nspin),       stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_occupy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(gloc(qbnd,nmesh,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_occupy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! calculate local green's function
! here, the asymptotic part is substracted
     gloc = czero
     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,1)
             be = kwin(k,s,2,1)
             cbnd = be - bs + 1
             do m=1,nmesh
                 caux = czi * fmesh(m) + fermi
                 do b=1,cbnd
                     gloc(b,m,s) = gloc(b,m,s) + one / ( caux - eigs(b,m,k,s) )
                     gloc(b,m,s) = gloc(b,m,s) - one / ( caux - einf(b,k,s) )
                 enddo ! over b={1,cbnd} loop
             enddo ! over m={1,nmesh} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

! calculate summation of the local green's function
     do s=1,nspin
         do b=1,qbnd
             zocc(b,s) = sum( gloc(b,:,s) ) / real(nkpt) * ( two / beta )
         enddo ! over b={1,cbnd} loop
     enddo ! over s={1,nspin} loop

! consider the contribution from asymptotic part 
     do s=1,nspin
         do k=1,nkpt
             bs = kwin(k,s,1,1)
             be = kwin(k,s,2,1)
             cbnd = be - bs + 1
             do b=1,cbnd
                 caux = einf(b,k,s) - fermi
                 zocc(b,s) = zocc(b,s) + fermi_dirac(real(caux)) / real(nkpt)
             enddo ! over b={1,cbnd} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

! sum up the density matrix
     val = real( sum(zocc) ) 

! consider the spins
     if ( nspin == 1 ) then
         val = val * two
     endif ! back if ( nspin == 1 ) block

! deallocate memory
     if ( allocated(zocc) ) deallocate(zocc)
     if ( allocated(gloc) ) deallocate(gloc)

     return
  end subroutine cal_occupy

!!
!! @sub cal_eigsys
!!
!! try to diagonalize H(k) + \Sigma(i\omega_n) and H(k) + \Sigma(\infty)
!! to obtain the corresponding eigenvalues 
!!
  subroutine cal_eigsys(eigs, einf)
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : i_wnd
     use context, only : ndim
     use context, only : qbnd
     use context, only : kwin

     implicit none

! external arguments
! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(out) :: eigs(qbnd,nmesh,nkpt,nspin)

! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(out) :: einf(qbnd,nkpt,nspin)

! local variables
! loop index for k-points
     integer :: k

! loop index for spins
     integer :: s

! loop index for impurity sites
     integer :: t

! number of dft bands for given k-point and spin
     integer :: cbnd

! number of correlated orbitals for given impurity site
     integer :: cdim

! band window: start index and end index for bands
     integer :: bs, be

! status flag
     integer :: istat

! self-energy functions, \Sigma(i\omega_n)
     complex(dp), allocatable :: Sk(:,:,:)

! H(k) + \Sigma(i\omega_n) 
     complex(dp), allocatable :: Hk(:,:,:)

! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), allocatable :: Ek(:,:)

! self-energy functions, \Sigma(\infty)
     complex(dp), allocatable :: So(:,:)

! H(k) + \Sigma(\infty) 
     complex(dp), allocatable :: Ho(:,:)

! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), allocatable :: Eo(:)

! well, here the number of impurity sites is restricted to be one 
! later we will remove this bug
     call s_assert2(nsite == 1, 'nsite should be 1')
     t = 1
     cdim = ndim(t)

! initialization
     eigs = czero
     einf = czero

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

! determine cbnd
             cbnd = be - bs + 1

! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)') 'window: ', bs, be, cbnd

! allocate memory
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Hk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Ek(cbnd,nmesh),      stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_eigsys','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block
             !
             allocate(So(cbnd,cbnd),       stat = istat)
             allocate(Ho(cbnd,cbnd),       stat = istat)
             allocate(Eo(cbnd),            stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_occupy','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

! construct H(k) + \Sigma(i\omega_n) and diagonalize it
             call cal_sl_sk(cdim, cbnd, k, s, t, Sk)
             !
             call cal_sk_hk(cbnd, bs, be, k, s, Sk, Hk)
             !
             call cal_hk_ek(cbnd, Hk, Ek)
             !
             eigs(1:cbnd,:,k,s) = Ek

! construct H(k) + \Sigma(\infty) and diagonalize it
!<             call cal_sk_so(cbnd, Sk, So)
             call cal_sl_so(cdim, cbnd, k, s, t, So)
             !
             call cal_so_ho(cbnd, bs, be, k, s, So, Ho)
             !
             call cal_ho_eo(cbnd, Ho, Eo)
             !
             einf(1:cbnd,k,s) = Eo

! deallocate memory
             if ( allocated(So) ) deallocate(So)
             if ( allocated(Ho) ) deallocate(Ho)
             if ( allocated(Eo) ) deallocate(Eo)

             if ( allocated(Sk) ) deallocate(Sk)
             if ( allocated(Hk) ) deallocate(Hk)
             if ( allocated(Ek) ) deallocate(Ek)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

     return
  end subroutine cal_eigsys

!!========================================================================
!!>>> service subroutines: set 5, fermi-dirac function                 <<<
!!========================================================================

!!
!! @fun fermi_dirac
!!
!! try to calculate the fermi-dirac function
!!
  function fermi_dirac(omega) result(value)
     use constants, only : dp
     use constants, only : zero, one

     use control, only : beta

     implicit none

! external arguments
! frequency point, \omega
     real(dp), intent(in) :: omega

! result value, return this
     real(dp) :: value

! check the range of omega to avoid numerical instability
     if      ( beta * omega >=  600.0_dp ) then
         value = zero
     else if ( beta * omega <= -600.0_dp ) then
         value = one
     else
         value = one / ( one + exp( beta * omega ) )
     endif ! back if ( beta * omega >=  600.0_dp ) block

     return
  end function fermi_dirac

!!========================================================================
!!>>> service subroutines: set 6, upfolding and downfolding            <<<
!!========================================================================

!!
!! @sub map_chi_psi
!!
!! service subroutine. map a function from local basis to Kohn-Sham
!! basis. you can call this procedure `embedding` or `upfold`
!!
  subroutine map_chi_psi(cdim, cbnd, nfrq, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

! external arguments
! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! number of frequency points
     integer, intent(in) :: nfrq

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! input array defined at local orbital (\chi) basis
     complex(dp), intent(in)  :: Mc(cdim,cdim,nfrq)

! output array defined at Kohn-Sham (\psi) basis
     complex(dp), intent(out) :: Mp(cbnd,cbnd,nfrq)

! local variables
! loop index for frequency mesh
     integer :: f

! status flag
     integer :: istat

! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

! allocate memory
     allocate(Cp(cdim,cbnd), stat = istat)
     allocate(Pc(cbnd,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('map_chi_psi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! copy data
     Cp = chipsi(1:cdim,1:cbnd,k,s,i_grp(t))
     Pc = psichi(1:cbnd,1:cdim,k,s,i_grp(t))

! upfolding or embedding
     do f=1,nfrq
         Mp(:,:,f) = matmul( matmul( Pc, Mc(:,:,f) ), Cp )
     enddo ! over f={1,nfrq} loop

! deallocate memory
     if ( allocated(Cp) ) deallocate(Cp)
     if ( allocated(Pc) ) deallocate(Pc)

     return
  end subroutine map_chi_psi

!!
!! @sub map_psi_chi
!!
!! service subroutine. map a function from Kohn-Sham basis to local
!! basis. you can call this procedure `projection` or `downfold`
!!
  subroutine map_psi_chi(cbnd, cdim, nfrq, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! number of frequency points
     integer, intent(in) :: nfrq

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! input array defined at Kohn-Sham (\psi) basis
     complex(dp), intent(in)  :: Mp(cbnd,cbnd,nfrq)

! output array defined at local orbital (\chi) basis
     complex(dp), intent(out) :: Mc(cdim,cdim,nfrq)

! local variables
! loop index for frequency mesh
     integer :: f

! status flag
     integer :: istat

! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

! allocate memory
     allocate(Cp(cdim,cbnd), stat = istat)
     allocate(Pc(cbnd,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('map_psi_chi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! copy data
     Cp = chipsi(1:cdim,1:cbnd,k,s,i_grp(t))
     Pc = psichi(1:cbnd,1:cdim,k,s,i_grp(t))

! downfolding or projection
     do f=1,nfrq
         Mc(:,:,f) = matmul( matmul( Cp, Mp(:,:,f) ), Pc )
     enddo ! over f={1,nfrq} loop

! deallocate memory
     if ( allocated(Cp) ) deallocate(Cp)
     if ( allocated(Pc) ) deallocate(Pc)

     return
  end subroutine map_psi_chi

!!
!! @sub one_chi_psi
!!
!! service subroutine. map a matrix from local basis to Kohn-Sham
!! basis. you can call this procedure `embedding` or `upfold`
!!
  subroutine one_chi_psi(cdim, cbnd, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

! external arguments
! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! input matrix defined at local orbital (\chi) basis
     complex(dp), intent(in)  :: Mc(cdim,cdim)

! output matrix defined at Kohn-Sham (\psi) basis
     complex(dp), intent(out) :: Mp(cbnd,cbnd)

! local variables
! status flag
     integer :: istat

! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

! allocate memory
     allocate(Cp(cdim,cbnd), stat = istat)
     allocate(Pc(cbnd,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('one_chi_psi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! copy data
     Cp = chipsi(1:cdim,1:cbnd,k,s,i_grp(t))
     Pc = psichi(1:cbnd,1:cdim,k,s,i_grp(t))

! upfolding or embedding
     Mp = matmul( matmul( Pc, Mc ), Cp )

! deallocate memory
     if ( allocated(Cp) ) deallocate(Cp)
     if ( allocated(Pc) ) deallocate(Pc)

     return
  end subroutine one_chi_psi

!!
!! @sub one_psi_chi
!!
!! service subroutine. map a matrix from Kohn-Sham basis to local
!! basis. you can call this procedure `projection` or `downfold`
!!
  subroutine one_psi_chi(cbnd, cdim, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

! external arguments
! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

! number of correlated orbitals for given impurity site
     integer, intent(in) :: cdim

! index for k-points
     integer, intent(in) :: k

! index for spin
     integer, intent(in) :: s

! index for impurity sites
     integer, intent(in) :: t

! input matrix defined at Kohn-Sham (\psi) basis
     complex(dp), intent(in)  :: Mp(cbnd,cbnd)

! output matrix defined at local orbital (\chi) basis
     complex(dp), intent(out) :: Mc(cdim,cdim)

! local variables
! status flag
     integer :: istat

! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

! allocate memory
     allocate(Cp(cdim,cbnd), stat = istat)
     allocate(Pc(cbnd,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('one_psi_chi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! copy data
     Cp = chipsi(1:cdim,1:cbnd,k,s,i_grp(t))
     Pc = psichi(1:cbnd,1:cdim,k,s,i_grp(t))

! downfolding or projection
     Mc = matmul( matmul( Cp, Mp ), Pc )

! deallocate memory
     if ( allocated(Cp) ) deallocate(Cp)
     if ( allocated(Pc) ) deallocate(Pc)

     return
  end subroutine one_psi_chi
