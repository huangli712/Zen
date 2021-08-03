!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : cal_sl_sk
!!!           cal_sk_hk
!!!           cal_hk_ek
!!!           cal_sl_so
!!!           cal_so_ho
!!!           cal_ho_eo
!!!           cal_sk_gk
!!!           cal_gk_gl
!!!           dichotomy
!!!           correction
!!!           cal_nelect
!!!           cal_occupy
!!!           cal_denmat
!!!           cal_eigsys
!!! source  : dmft_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           08/01/2021 by li huang (last modified)
!!! purpose : provide the core service subroutines for the work flow of
!!!           the dft + dmft calculations.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> service subroutines: set 1                                       <<<
!!========================================================================

!!
!! @sub cal_sl_sk
!!
!! try to upfold the self-energy function from local basis to Kohn-Sham
!! basis. here, we don't care whether the double counting terms have been
!! substracted from the self-energy functions.
!!
  subroutine cal_sl_sk(cdim, cbnd, k, s, t, Sk)
     use constants, only : dp
     use constants, only : czero

     use control, only : nmesh

     use context, only : sigma

     implicit none

!! external arguments
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

!! local variables
     ! loop index for frequency mesh
     integer :: m

     ! status flag
     integer :: istat

     ! dummy array: for local self-energy function
     complex(dp), allocatable :: Sl(:,:,:)

!! [body

     ! allocate memory
     allocate(Sl(cdim,cdim,nmesh), stat = istat)
     !
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_sk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! we use Sl to store parts of sigma.
     ! note that actually sigdc has been substracted from sigma before
     ! hand. please see cal_sigma() for more details.
     do m=1,nmesh
         Sl(:,:,m) = sigma(1:cdim,1:cdim,m,s,t)
     enddo ! over m={1,nmesh} loop

     ! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call map_chi_psi(cdim, cbnd, nmesh, k, s, t, Sl, Sk)

     ! deallocate memory
     if ( allocated(Sl) ) deallocate(Sl)

!! body]

     return
  end subroutine cal_sl_sk

!!
!! @sub cal_sk_hk
!!
!! try to build H(k) + \Sigma(i\omega_n). here, \Sigma should contain
!! contributions from all impurity sites. so, only when nsite = 1, we
!! can use the output of cal_sl_sk() as the input of this subroutine.
!!
  subroutine cal_sk_hk(cbnd, bs, be, k, s, Sk, Hk)
     use constants, only : dp

     use control, only : nmesh

     use context, only : enk

     implicit none

!! external arguments
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

     ! frequency-dependent effective hamiltonian at given k-point and spin
     complex(dp), intent(out) :: Hk(cbnd,cbnd,nmesh)

!! local variables
     ! loop index for frequency mesh
     integer :: m

     ! status flag
     integer :: istat

     ! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

!! [body

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

!! body]

     return
  end subroutine cal_sk_hk

!!
!! @sub cal_hk_ek
!!
!! try to diagonalize the effective hamiltonian: H(k) + \Sigma(i\omega_n),
!! get all the complex eigenvalues. here, we just assumed the effective
!! hamiltonian is a general complex matrix.
!!
  subroutine cal_hk_ek(cbnd, Hk, Ek)
     use constants, only : dp

     use control, only : nmesh

     implicit none

!! external arguments
     ! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

     ! effective hamiltonian: H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in)  :: Hk(cbnd,cbnd,nmesh)

     ! resulting eigenvalues
     complex(dp), intent(out) :: Ek(cbnd,nmesh)

!! local variables
     ! loop index for frequency mesh
     integer :: m

!! [body

     FREQ_LOOP: do m=1,nmesh
         call s_eigvals_zg(cbnd, cbnd, Hk(:,:,m), Ek(:,m))
     enddo FREQ_LOOP ! over m={1,nmesh} loop

!! body]

     return
  end subroutine cal_hk_ek

!!========================================================================
!!>>> service subroutines: set 2                                       <<<
!!========================================================================

!!
!! @sub cal_sl_so
!!
!! try to upfold the asymptotic values of self-energy functions at high
!! frequency (i.e `sigoo`) from local basis to Kohn-Sham basis. note that
!! the double-counting terms have been substracted from sigoo before hand.
!! please see cal_sigoo() for more details.
!!
  subroutine cal_sl_so(cdim, cbnd, k, s, t, So)
     use constants, only : dp
     use constants, only : czero

     use context, only : sigoo

     implicit none

!! external arguments
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

!! local variables
     ! status flag
     integer :: istat

     ! dummy array: for asymptotic self-energy function
     complex(dp), allocatable :: Sl(:,:)

!! [body

     ! allocate memory
     allocate(Sl(cdim,cdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sl_so','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! here we use Sl to save parts of sigoo
     Sl = sigoo(1:cdim,1:cdim,s,t)

     ! upfolding: Sl (local basis) -> Sk (Kohn-Sham basis)
     call one_chi_psi(cdim, cbnd, k, s, t, Sl, So)

     ! deallocate memory
     if ( allocated(Sl) ) deallocate(Sl)

!! body]

     return
  end subroutine cal_sl_so

!!
!! @sub cal_so_ho
!!
!! try to build H(k) + \Sigma(\infty). here, \Sigma should contain full
!! contributions from all impurity sites. so, only when nsite = 1, we
!! can use the output of cal_sl_so() as the input of this subroutine.
!!
  subroutine cal_so_ho(cbnd, bs, be, k, s, So, Ho)
     use constants, only : dp

     use context, only : enk

     implicit none

!! external arguments
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

!! local variables
     ! status flag
     integer :: istat

     ! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

!! [body

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

!! body]

     return
  end subroutine cal_so_ho

!!
!! @sub cal_ho_eo
!!
!! try to diagonalize the effective hamiltonian: H(k) + \Sigma(\infty),
!! get all the complex eigenvalues. here, we just assumed the effective
!! hamiltonian is a general complex matrix.
!!
  subroutine cal_ho_eo(cbnd, Ho, Eo)
     use constants, only : dp

     implicit none

!! external arguments
     ! number of dft bands for given k-point and spin
     integer, intent(in) :: cbnd

     ! effective hamiltonian: H(k) + \Sigma(\infty)
     complex(dp), intent(in)  :: Ho(cbnd,cbnd)

     ! resulting eigenvalues
     complex(dp), intent(out) :: Eo(cbnd)

!! [body

     call s_eigvals_zg(cbnd, cbnd, Ho, Eo)

!! body]

     return
  end subroutine cal_ho_eo

!!========================================================================
!!>>> service subroutines: set 3                                       <<<
!!========================================================================

!!
!! @sub cal_sk_gk
!!
!! try to calculate lattice green's function at given k-point and spin.
!! this subroutine needs the self-energy function at Kohn-Sham basis (i.e
!! `Sk`), that is the reason why it is named as `cal_sk_gk`. note that Sk
!! have to contain the full contributions from all impurity sites.
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

!! external arguments
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

!! local variables
     ! loop index for frequency mesh
     integer :: m

     ! status flag
     integer :: istat

     ! dummy array: for band dispersion (vector)
     complex(dp), allocatable :: Fm(:)
     complex(dp), allocatable :: Em(:)

     ! dummy array: for effective hamiltonian (diagonal matrix)
     complex(dp), allocatable :: Hm(:,:)

!! [body

     ! allocate memory
     allocate(Fm(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_gk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Em(cbnd),      stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_gk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Hm(cbnd,cbnd), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sk_gk','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! evaluate Fm, which is k-dependent, but frequency-independent.
     !
     ! if you want to consider magnetic or the other external fields,
     ! you can insert your codes here.
     Fm = fermi - enk(bs:be,k,s)

     ! loop over frequency points
     FREQ_LOOP: do m=1,nmesh

         ! consider imaginary axis or real axis
         if ( axis == 1 ) then
             Em = czi * fmesh(m) + Fm
         else
             Em = fmesh(m) + Fm
         endif ! back if ( axis == 1 ) block

         ! convert Em (vector) to Hm (diagonal matrix)
         call s_diag_z(cbnd, Em, Hm)

         ! substract self-energy function from the hamiltonian
         Gk(:,:,m) = Hm - Sk(:,:,m)

         ! calculate lattice green's function by direct inversion
         call s_inv_z(cbnd, Gk(:,:,m))

     enddo FREQ_LOOP ! over m={1,nmesh} loop

     ! deallocate memory
     if ( allocated(Fm) ) deallocate(Fm)
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

!! body]

     return
  end subroutine cal_sk_gk

!!
!! @sub cal_gk_gl
!!
!! try to calculate local green's function from lattice green's function
!! via downfolding procedure.
!!
  subroutine cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl)
     use constants, only : dp

     use control, only : nmesh

     implicit none

!! external arguments
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

!! [body

     call map_psi_chi(cbnd, cdim, nmesh, k, s, t, Gk, Gl)

!! body]

     return
  end subroutine cal_gk_gl

!!========================================================================
!!>>> service subroutines: set 4, fermi level and density matrix       <<<
!!========================================================================

!!
!! @sub dichotomy
!!
!! try to locate the fermi level with the bisection method.
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

!! external arguments
     ! desired charge density
     real(dp), intent(in)    :: desired

     ! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(qbnd,nmesh,nkpt,nspin)

     ! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(in) :: einf(qbnd,nkpt,nspin)

!! local parameters
     ! maximum number of loops for the bisection algorithm
     integer, parameter  :: max_loops = 100

     ! step for locating the energy boundary
     real(dp), parameter :: delta = 0.5_dp

!! local variables
     ! loop index for the bisection algorithm
     integer  :: loop

     ! left boundary, right boundary,
     ! and the final result for the fermi level.
     real(dp) :: mu1, mu2, mu3

     ! the corresponding charge density
     real(dp) :: occ1, occ2, occ3

     ! sign
     real(dp) :: sign

!! [body

     ! print message in the terminal
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'searching fermi level'
     endif ! back if ( myid == master ) block

     ! initialization, determine mu1, mu2, occ1, occ2, and sign.
     !
     ! (1) if sign < 0, it means occ1 < desired, we should push
     !     mu2 to higher energy.
     !
     ! (2) if sign > 0, it means occ1 > desired. then mu1 will be
     !     the right boundary, and we should push mu2 to lower energy.
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
     do while ( loop <= max_loops .and. &
                ( occ2 - desired ) * sign > 0 .and. &
                abs( occ2 - desired ) > mc )
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

     ! now the fermi level should lie in the regine [mu1, mu2].
     ! we should refine the boundary to locate the fermi level.
     if ( abs(occ1 - desired) < abs(occ2 - desired) ) then
         mu3 = mu1
         occ3 = occ1
     else
         mu3 = mu2
         occ3 = occ2
     endif ! back if block
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
         if ( myid == master ) then
             write(mystd,'(6X,a)',advance = 'no') 'final results ->'
             write(mystd,'(2(2X,a,f12.8))') 'EF: ', mu3, 'density: ', occ3
         endif ! back if ( myid == master ) block
     else
         call s_print_error('dichotomy', 'fail to locate the fermi level')
     endif ! back if ( abs(occ3 - desired) < mc ) block

     ! well, finally, we have to update the global variable `fermi`.
     fermi = mu3

     return
  end subroutine dichotomy

!!
!! @sub correction
!!
!! try to evaluate the difference between the dft density matrix and the
!! dft + dmft density matrix.
!!
  subroutine correction(kocc, gamma, ecorr)
     use constants, only : dp, mystd
     use constants, only : zero, two, czero

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qbnd
     use context, only : kwin
     use context, only : weight
     use context, only : enk, occupy

     implicit none

!! external arguments
     ! density matrix from dft + dmft calculations
     complex(dp), intent(in)  :: kocc(qbnd,qbnd,nkpt,nspin)

     ! correction for density matrix
     complex(dp), intent(out) :: gamma(qbnd,qbnd,nkpt,nspin)

     ! correction for band energy
     real(dp), intent(out)    :: ecorr

!! local variables
     ! index for spins
     integer :: s

     ! index for k-points
     integer :: k

     ! index for orbitals
     integer :: p, q

     ! index for impurity sites
     integer :: t

     ! band window: start index and end index for bands
     integer :: bs, be

     ! number of dft bands for given k-point and spin
     integer :: cbnd

     ! status flag
     integer :: istat

     ! dummy variable, used to perform mpi reduce operation for ecorr
     real(dp)    :: ecorr_mpi
     complex(dp) :: tr

     ! dummy arrays, used to build effective hamiltonian
     complex(dp), allocatable :: Em(:)
     complex(dp), allocatable :: Hm(:,:)

     ! dummy array, used to perform mpi reduce operation for gamma
     complex(dp), allocatable :: gamma_mpi(:,:,:,:)

!! [body

     ! allocate memory
     allocate(gamma_mpi(qbnd,qbnd,nkpt,nspin), stat = istat)
     !
     if ( istat /= 0 ) then
         call s_print_error('correction','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! reset gamma
     gamma = czero
     gamma_mpi = czero

     ! reset ecorr
     ecorr = zero
     ecorr_mpi = zero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate correction for density matrix'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     ! loop over spins and k-points
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect() for more details.
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

             ! determine cbnd
             cbnd = be - bs + 1

             ! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)',advance='no') 'window: ', bs, be, cbnd
             write(mystd,'(2X,a,i2)') 'proc: ', myid

             ! allocate memory
             allocate(Em(cbnd),      stat = istat)
             allocate(Hm(cbnd,cbnd), stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('correction','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! calculate the difference between dft + dmft density matrix
             ! `kocc` and the dft density matrix `occupy`. the results are
             ! saved at `gamma.
             do p = 1,cbnd
                 do q = 1,cbnd
                     if ( p /= q ) then
                         gamma(q,p,k,s) = kocc(q,p,k,s)
                     else
                         gamma(q,p,k,s) = kocc(q,p,k,s) - occupy(bs + q - 1,k,s)
                     endif
                 enddo ! over q={1,cbnd} loop
             enddo ! over p={1,cbnd} loop

             ! Now gamma(:,:,k,s) is ready, we would like to use it to
             ! calculate its contribution to band energy.

             ! evaluate Em, which is the eigenvalues.
             Em = enk(bs:be,k,s)

             ! convert `Em` to diagonal matrix `Hm`
             call s_diag_z(cbnd, Em, Hm)

             ! evaluate correction to band energy
             Hm = matmul(gamma(:,:,k,s), Hm)
             call s_trace_z(cbnd, Hm, tr)
             ecorr = ecorr + real(tr) * weight(k)

             ! deallocate memory
             if ( allocated(Em) ) deallocate(Em)
             if ( allocated(Hm) ) deallocate(Hm)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! collect data from all mpi processes
# if defined (MPI)
     !
     call mp_barrier()
     !
     call mp_allreduce(gamma, gamma_mpi)
     call mp_allreduce(ecorr, ecorr_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     gamma_mpi = gamma
     ecorr_mpi = ecorr

# endif /* MPI */

     ! get the final correction for density matrix
     gamma = gamma_mpi
     ecorr = ecorr_mpi / float(nkpt)
     if ( nspin == 1 ) then
         ecorr = ecorr * two
     endif ! back if ( nspin == 1 ) block

     ! deallocate memory
     deallocate(gamma_mpi)

!! body]

     return
  end subroutine correction

!!
!! @sub cal_nelect
!!
!! try to calculate the number of valence electrons by dft occupations
!! and weights. actually, what we obtain is the nominal occupation
!! numbers in the selected band window.
!!
  subroutine cal_nelect(nelect)
     use constants, only : dp, mystd
     use constants, only : zero, two

     use control, only : nkpt, nspin
     use control, only : myid, master

     use context, only : i_wnd
     use context, only : kwin
     use context, only : weight
     use context, only : occupy

     implicit none

!! external arguments
     ! number of relevant electrons
     real(dp), intent(out) :: nelect

!! local variables
     ! index for k-points
     integer :: k

     ! index for spin
     integer :: s

     ! band window: start index and end index for bands
     integer :: bs, be

!! [body

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculating desired charge density'
     endif ! back if ( myid == master ) block

!
! important remarks:
!
! there is an important question. which band window shall we used?
! multiple band windows are always possible since we don't have enough
! reason to restrict `nwnd` to be 1.
!
! the answer is the band window that is connected to the quantum
! impurity problems. basically, now we require that all the quantum
! impurity problems share the same band window. so in this subroutine,
! we only need to consider ONE band window, which is defined by
! i_wnd(1) or i_wnd(2). whatever, they should point to the same band
! window according to our assumption.
!

     ! reset nelect
     nelect = zero

     ! loop over spins and k-points to perform summation
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt
             bs = kwin(k,s,1,i_wnd(1))
             be = kwin(k,s,2,i_wnd(1))
             nelect = nelect + sum( occupy(bs:be,k,s) ) * weight(k)
         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

     ! don't forget to normalize `nelect`
     nelect = nelect / float(nkpt)

     ! consider the spins
     if ( nspin == 1 ) then
         nelect = nelect * two
     endif ! back if ( nspin == 1 ) block

!! body]

     return
  end subroutine cal_nelect

!!
!! @sub cal_occupy
!!
!! for given fermi level, try to calculate the corresponding occupations.
!! note that this subroutine only works in imaginary axis.
!!
  subroutine cal_occupy(fermi, val, eigs, einf)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czi, czero

     use control, only : axis
     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : beta

     use context, only : i_wnd
     use context, only : qbnd
     use context, only : kwin
     use context, only : weight
     use context, only : fmesh

     implicit none

!! external arguments
     ! assumed fermi level
     real(dp), intent(in)  :: fermi

     ! occupation number
     real(dp), intent(out) :: val

     ! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(qbnd,nmesh,nkpt,nspin)

     ! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(in) :: einf(qbnd,nkpt,nspin)

!! local variables
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
     complex(dp) :: caux, ctmp

     ! density matrix
     complex(dp), allocatable :: zocc(:,:)

     ! local green's function
     complex(dp), allocatable :: gloc(:,:,:)

!! external functions
     ! used to calculate fermi-dirac function
     real(dp), external :: fermi_dirac

!! [body

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

     ! check axis
     call s_assert2(axis == 1, 'axis is wrong')

     ! reset zocc and gloc
     zocc = czero
     gloc = czero

     ! loop over spins and k-points to perform k-summation to determine
     ! the local green's function
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=1,nkpt

             ! determine the band window
             ! see remarks in cal_nelect()
             bs = kwin(k,s,1,i_wnd(1))
             be = kwin(k,s,2,i_wnd(1))
             cbnd = be - bs + 1

             ! here, the asymptotic part is substracted
             do m=1,nmesh
                 caux = czi * fmesh(m) + fermi
                 do b=1,cbnd
                     ctmp = one / ( caux - eigs(b,m,k,s) ) - one / ( caux - einf(b,k,s) )
                     gloc(b,m,s) = gloc(b,m,s) + ctmp * weight(k)
                 enddo ! over b={1,cbnd} loop
             enddo ! over m={1,nmesh} loop

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

     ! calculate summation of the local green's function
     do s=1,nspin
         do b=1,qbnd
             zocc(b,s) = sum( gloc(b,:,s) ) / real(nkpt) * ( two / beta )
         enddo ! over b={1,cbnd} loop
     enddo ! over s={1,nspin} loop

     ! consider the contribution from asymptotic part
     do s=1,nspin
         do k=1,nkpt

             ! determine the band window
             ! see remarks in cal_nelect()
             bs = kwin(k,s,1,i_wnd(1))
             be = kwin(k,s,2,i_wnd(1))
             cbnd = be - bs + 1

             do b=1,cbnd
                 caux = einf(b,k,s) - fermi
                 zocc(b,s) = zocc(b,s) + fermi_dirac( real(caux) ) * weight(k) / real(nkpt)
             enddo ! over b={1,cbnd} loop

         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

     !
     ! actually, we should consider the correction due to finite
     ! frequency point here. later we will implement it.
     !
     ! TO_BE_DONE
     !

     ! sum up the density matrix
     val = real( sum(zocc) )

     ! consider the spins
     if ( nspin == 1 ) then
         val = val * two
     endif ! back if ( nspin == 1 ) block

     ! deallocate memory
     if ( allocated(zocc) ) deallocate(zocc)
     if ( allocated(gloc) ) deallocate(gloc)

!! body]

     return
  end subroutine cal_occupy

!!
!! @sub cal_denmat
!!
!! try to calculate the dft + dmft density matrix for given fermi level.
!!
  subroutine cal_denmat(kocc)
     use constants, only : dp, mystd
     use constants, only : one, czero

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : beta
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qbnd
     use context, only : ndim
     use context, only : kwin
     use context, only : fmesh

     implicit none

!! external arguments
     ! dft + dmft density matrix
     complex(dp), intent(out) :: kocc(qbnd,qbnd,nkpt,nspin)

!! local variables
     ! loop index for spin
     integer :: s

     ! loop index for k-points
     integer :: k

     ! loop index for orbitals
     integer :: p, q

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

     ! orbital density for given k and spin
     complex(dp) :: density

     ! dummy array: for self-energy function (upfolded to Kohn-Sham basis)
     complex(dp), allocatable :: Sk(:,:,:)
     complex(dp), allocatable :: Xk(:,:,:)

     ! dummy array: for lattice green's function
     complex(dp), allocatable :: Gk(:,:,:)

     ! dummy array: used to perform mpi reduce operation for kocc
     complex(dp), allocatable :: kocc_mpi(:,:,:,:)

!! [body

     ! allocate memory
     allocate(kocc_mpi(qbnd,qbnd,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_denmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! reset cbnd and cdim. they will be updated later.
     ! cbnd should be k-dependent and cdim should be impurity-dependent.
     cbnd = 0
     cdim = 0

     ! reset kocc and kocc_mpi
     kocc = czero
     kocc_mpi = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculate dft + dmft density matrix'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     ! loop over spins and k-points
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect() for more details.
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

             ! determine cbnd
             cbnd = be - bs + 1

             ! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)',advance='no') 'window: ', bs, be, cbnd
             write(mystd,'(2X,a,i2)') 'proc: ', myid

             ! allocate memories Sk, Xk, and Gk.
             ! their sizes are k-dependent.
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Xk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Gk(cbnd,cbnd,nmesh), stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_denmat','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! build self-energy function, and then upfold it into
             ! Kohn-Sham basis. Sk should contain contributions from
             ! all impurity sites
             Sk = czero
             do t=1,nsite
                 Xk = czero ! reset Xk
                 cdim = ndim(t)
                 call cal_sl_sk(cdim, cbnd, k, s, t, Xk)
                 Sk = Sk + Xk
             enddo ! over t={1,nsite} loop

             ! calculate lattice green's function
             call cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)

             ! try to calculate momentum-dependent and spin-dependent
             ! density matrix.
             do p=1,cbnd
                 do q=1,cbnd
                     call s_fft_density(nmesh, fmesh, Gk(q,p,:), density, beta)
                     if ( p == q ) then
                         kocc(q,p,k,s) = one + density
                     else
                         kocc(q,p,k,s) = density
                     endif
                 enddo ! over q={1,cbnd} loop
             enddo ! over p={1,cbnd} loop

             ! deallocate memories
             if ( allocated(Sk) ) deallocate(Sk)
             if ( allocated(Xk) ) deallocate(Xk)
             if ( allocated(Gk) ) deallocate(Gk)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! collect data from all mpi processes
# if defined (MPI)
     !
     call mp_barrier()
     !
     call mp_allreduce(kocc, kocc_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     kocc_mpi = kocc

# endif /* MPI */

     ! renormalize density matrix
     kocc = kocc_mpi

     ! deallocate memory
     if ( allocated(kocc_mpi) ) deallocate(kocc_mpi)

!! body]

     return
  end subroutine cal_denmat

!!
!! @sub cal_eigsys
!!
!! try to diagonalize H(k) + \Sigma(i\omega_n) and H(k) + \Sigma(\infty)
!! to obtain the corresponding eigenvalues.
!!
  subroutine cal_eigsys(eigs, einf)
     use constants, only : dp, mystd
     use constants, only : czero

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : ndim
     use context, only : qbnd
     use context, only : kwin

     implicit none

!! external arguments
     ! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(out) :: eigs(qbnd,nmesh,nkpt,nspin)

     ! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), intent(out) :: einf(qbnd,nkpt,nspin)

!! local variables
     ! loop index for spins
     integer :: s

     ! loop index for k-points
     integer :: k

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
     complex(dp), allocatable :: Xk(:,:,:)

     ! H(k) + \Sigma(i\omega_n)
     complex(dp), allocatable :: Hk(:,:,:)

     ! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), allocatable :: Ek(:,:)

     ! self-energy functions, \Sigma(\infty)
     complex(dp), allocatable :: So(:,:)
     complex(dp), allocatable :: Xo(:,:)

     ! H(k) + \Sigma(\infty)
     complex(dp), allocatable :: Ho(:,:)

     ! eigenvalues for H(k) + \Sigma(\infty)
     complex(dp), allocatable :: Eo(:)

     ! dummy array, used to perform mpi reduce operation for eigs
     complex(dp), allocatable :: eigs_mpi(:,:,:,:)

     ! dummy array, used to perform mpi reduce operation for einf
     complex(dp), allocatable :: einf_mpi(:,:,:)

!! [body

     ! allocate memory
     allocate(eigs_mpi(qbnd,nmesh,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eigsys','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(einf_mpi(qbnd,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eigsys','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialization
     eigs = czero
     einf = czero

     eigs_mpi = czero
     einf_mpi = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculating dft + dmft eigenvalues'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     ! loop over spins and k-points
     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect() subroutine.
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

             ! determine cbnd
             cbnd = be - bs + 1

             ! provide some useful information
             write(mystd,'(6X,a,i2)',advance='no') 'spin: ', s
             write(mystd,'(2X,a,i5)',advance='no') 'kpnt: ', k
             write(mystd,'(2X,a,3i3)',advance='no') 'window: ', bs, be, cbnd
             write(mystd,'(2X,a,i2)') 'proc: ', myid

             ! allocate memory
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Xk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Hk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Ek(cbnd,nmesh),      stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_eigsys','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block
             !
             allocate(So(cbnd,cbnd),       stat = istat)
             allocate(Xo(cbnd,cbnd),       stat = istat)
             allocate(Ho(cbnd,cbnd),       stat = istat)
             allocate(Eo(cbnd),            stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_eigsys','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! construct H(k) + \Sigma(i\omega_n) and diagonalize it
             Sk = czero
             !
             do t=1,nsite ! add contributions from all impurity sites
                 Xk = czero
                 cdim = ndim(t)
                 call cal_sl_sk(cdim, cbnd, k, s, t, Xk)
                 Sk = Sk + Xk
             enddo ! over t={1,nsite} loop
             !
             call cal_sk_hk(cbnd, bs, be, k, s, Sk, Hk)
             !
             call cal_hk_ek(cbnd, Hk, Ek)
             !
             eigs(1:cbnd,:,k,s) = Ek

             ! construct H(k) + \Sigma(\infty) and diagonalize it
             So = czero
             !
             do t=1,nsite ! add contributions from all impurity sites
                 Xo = czero
                 cdim = ndim(t)
                 call cal_sl_so(cdim, cbnd, k, s, t, Xo)
                 So = So + Xo
             enddo ! over t={1,nsite} loop
             !
             call cal_so_ho(cbnd, bs, be, k, s, So, Ho)
             !
             call cal_ho_eo(cbnd, Ho, Eo)
             !
             einf(1:cbnd,k,s) = Eo

             ! deallocate memory
             if ( allocated(So) ) deallocate(So)
             if ( allocated(Xo) ) deallocate(Xo)
             if ( allocated(Ho) ) deallocate(Ho)
             if ( allocated(Eo) ) deallocate(Eo)
             !
             if ( allocated(Sk) ) deallocate(Sk)
             if ( allocated(Xk) ) deallocate(Xk)
             if ( allocated(Hk) ) deallocate(Hk)
             if ( allocated(Ek) ) deallocate(Ek)

         enddo KPNT_LOOP ! over k={1,nkpt} loop
     enddo SPIN_LOOP ! over s={1,nspin} loop

! collect data from all mpi processes
# if defined (MPI)
     !
     call mp_barrier()
     !
     call mp_allreduce(eigs, eigs_mpi)
     call mp_allreduce(einf, einf_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     eigs_mpi = eigs
     einf_mpi = einf

# endif /* MPI */

     ! save the final results
     eigs = eigs_mpi
     einf = einf_mpi

     ! deallocate memory
     deallocate(eigs_mpi)
     deallocate(einf_mpi)

!! body]

     return
  end subroutine cal_eigsys
