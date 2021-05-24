!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_driver
!!!           dmft_try1
!!!           dmft_try2
!!!           dmft_try3
!!!           dmft_try4
!!!           dmft_try5
!!!           dmft_try999
!!!           cal_sigoo
!!!           cal_sig_l
!!!           cal_fermi
!!!           cal_eimps
!!!           cal_eimpx
!!!           cal_grn_l
!!!           cal_wss_l
!!!           cal_hyb_l
!!!           cal_sl_sk
!!!           cal_sk_hk
!!!           cal_hk_ek
!!!           cal_sl_so
!!!           cal_so_ho
!!!           cal_ho_eo
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
!!!           05/24/2021 by li huang (last modified)
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
!! core subroutine, dispatch all the computational tasks
!!
  subroutine dmft_driver()
     use control, only : task

     implicit none

! we have to preprocess the self-energy functions at first
! calculate sig_l -> sigoo -> sigoo - sigdc
     call cal_sigoo()

! calculate sig_l -> sig_l - sigdc (sig_l is updated)
     call cal_sig_l()

     DISPATCHER: select case ( task )

! task = 1, calculate the hybridization function, for one-shot calculation
         case (1)
             call dmft_try1()

! task = 2, calculate density correction, for self-consistent calculation
         case (2)
             call dmft_try2()

! task = 3, search the fermi level
         case (3)
             call dmft_try3()

! task = 4, calculate the impurity levels
         case (4)
             call dmft_try4()

! task = 5, calculate complex dft + dmft eigenvalues
         case (5)
             call dmft_try5()

! task = 999, only for test
         case (999)
             call dmft_try999()

         case default
             call s_print_error('dmft_driver','this feature is not supported')

     end select DISPATCHER

     return
  end subroutine dmft_driver

!!========================================================================
!!>>> driver subroutines: layer 2                                      <<<
!!========================================================================

!!
!! @sub dmft_try1
!!
!! to calculate the local green's function, generate key inputs for the
!! quantum impurity solvers. the fermi level may be updated, depending
!! on the configuration parameter. this subroutine is suitable for the
!! one-shot dft + dmft calculations.
!!
  subroutine dmft_try1()
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : fermi
     use control, only : myid, master

     use context, only : eimps, eimpx
     use context, only : grn_l
     use context, only : wss_l, hyb_l

     implicit none

! try to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     if ( lfermi .eqv. .true. ) then
         call cal_fermi()
     else
         if ( myid == master ) then
             write(mystd,'(4X,a)') 'SKIP'
         endif ! back if ( myid == master ) block
     endif ! back if ( lfermi .eqv. .true. ) block
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local impurity levels
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Level'
     endif ! back if ( myid == master ) block
     !
     call cal_eimps()
     !
     call cal_eimpx()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local green's function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Green'
     endif ! back if ( myid == master ) block
     !
     call cal_grn_l()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the hybridization function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Hybri'
     endif ! back if ( myid == master ) block
     !
     call cal_hyb_l()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! try to compute the local weiss's function
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Weiss'
     endif ! back if ( myid == master ) block
     !
     call cal_wss_l()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save fermi...'
         call dmft_dump_fermi(fermi)
         !
         write(mystd,'(4X,a)') 'save eimps...'
         call dmft_dump_eimps(eimps)
         !
         write(mystd,'(4X,a)') 'save eimpx...'
         call dmft_dump_eimpx(eimpx)
         !
         write(mystd,'(4X,a)') 'save grn_l...'
         call dmft_dump_grn_l(grn_l)
         !
         write(mystd,'(4X,a)') 'save hyb_l...'
         call dmft_dump_hyb_l(hyb_l)
         !
         write(mystd,'(4X,a)') 'save wss_l...'
         call dmft_dump_wss_l(wss_l)
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

!!
!! @sub dmft_try3
!!
!! try to search the fermi level, the global variable `fermi` may be
!! updated in this subroutine. it is just for testing purpose.
!!
  subroutine dmft_try3()
     use constants, only : mystd

     use control, only : cname
     use control, only : lfermi
     use control, only : fermi
     use control, only : myid, master

     implicit none

! check lfermi at first
     call s_assert2(lfermi .eqv. .true., 'lfermi must be true')

! try to search the fermi level
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Fermi'
     endif ! back if ( myid == master ) block
     !
     call cal_fermi()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save fermi...'
         call dmft_dump_fermi(fermi)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try3

!!
!! @sub dmft_try4
!!
!! try to determine the local impurity levels. this subroutine is just
!! for testing purpose.
!!
  subroutine dmft_try4()
     use constants, only : mystd

     use control, only : cname
     use control, only : myid, master

     use context, only : eimps, eimpx

     implicit none

! try to calculate the local impurity levels
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Level'
     endif ! back if ( myid == master ) block
     !
     call cal_eimps()
     !
     call cal_eimpx()
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save eimps...'
         call dmft_dump_eimps(eimps)
         !
         write(mystd,'(4X,a)') 'save eimpx...'
         call dmft_dump_eimpx(eimpx)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine dmft_try4

!!
!! @sub dmft_try5
!!
!! try to calculate all the complex dft + dmft eigenvalues. the subroutine
!! can be used in the postprocessing procedure.
!!
  subroutine dmft_try5()
     use constants, only : dp, mystd

     use control, only : cname
     use control, only : nkpt, nspin
     use control, only : nmesh
     use control, only : myid, master

     use context, only : qbnd

     implicit none

! local variables
! status flag
     integer  :: istat

! dummy array, used to save the eigenvalues of H + \Sigma(i\omega_n)
     complex(dp), allocatable :: eigs(:,:,:,:)

! dummy array, used to save the eigenvalues of H + \Sigma(ioo)
     complex(dp), allocatable :: einf(:,:,:)

! allocate memory
     allocate(eigs(qbnd,nmesh,nkpt,nspin), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('dmft_try5','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(einf(qbnd,nkpt,nspin),       stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('dmft_try5','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! try to diagonalize the effective hamiltonian
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Eigen'
     endif ! back if ( myid == master ) block
     !
     call cal_eigsys(eigs, einf)
     !
     if ( myid == master ) then
         write(mystd,*)
     endif ! back if ( myid == master ) block

! write the calculated results, only the master node can do it
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname // ' >>> Task : Write'
         !
         write(mystd,'(4X,a)') 'save eigen...'
         call dmft_dump_eigen(eigs)
         !
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

     return
  end subroutine dmft_try5

!!
!! @sub dmft_try999
!!
  subroutine dmft_try999()
     implicit none

     return
  end subroutine dmft_try999

!!========================================================================
!!>>> driver subroutines: layer 3                                      <<<
!!========================================================================

!!
!! @sub cal_sigoo
!!
!! try to calculate the asymptotic values for self-energy functions. the
!! double-counting terms will be removed as well. this function works for
!! Matsubara self-energy functions (bare) only.
!!
  subroutine cal_sigoo()
     use constants, only : dp
     use constants, only : czero

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : qdim
     use context, only : sigdc, sigoo, sig_l

     implicit none

! local parameters
! how many frequency points are included to calculate the asymptotic
! values of Matsubara self-energy function
     integer, parameter :: mcut = 16

! local variables
! loop index for frequency mesh
     integer :: m

! loop index for spins
     integer :: s

! loop index for impurity sites
     integer :: t

! status flag
     integer :: istat

! dummy array for the Matsubara self-energy functions
     complex(dp), allocatable :: Sm(:,:)

! check working axis
     call s_assert2(axis == 1, 'axis is wrong')

! allocate memory
     allocate(Sm(qdim,qdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sigoo','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! reset sigoo
     sigoo = czero

! loop over nsite and nspin
!
! we count the last `mcut` frequency points, then we try to calculate
! the averaged values. up to now, the double counting terms have not
! been substracted from sig_l. in other words, sig_l is still bare.
     do t=1,nsite
         do s=1,nspin
             Sm = czero
             !
             do m=1,mcut
                 Sm = Sm + sig_l(:,:,nmesh + 1 - m,s,t)
             enddo ! over m={1,mcut} loop
             !
             sigoo(:,:,s,t) = Sm / float(mcut)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! we substract the double counting terms from sigoo
     sigoo = sigoo - sigdc

! deallocate memory
     if ( allocated(Sm) ) deallocate(Sm)

     return
  end subroutine cal_sigoo

!!
!! @sub cal_sig_l
!!
!! try to substract the double counting terms from the bare Matsubara
!! self-energy functions. this function works for Matsubara self-energy
!! functions (bare) only.
!!
  subroutine cal_sig_l()
     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : sigdc, sig_l

     implicit none

! local variables
! loop index for frequency mesh
     integer :: m

! loop index for spins
     integer :: s

! loop index for impurity sites
     integer :: t

! check working axis
     call s_assert2(axis == 1, 'axis is wrong')

! substract the double counting terms
     do t=1,nsite
         do s=1,nspin
             do m=1,nmesh
                 sig_l(:,:,m,s,t) = sig_l(:,:,m,s,t) - sigdc(:,:,s,t)
             enddo ! over m={1,nmesh} loop
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

     return
  end subroutine cal_sig_l

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

     use context, only : qbnd

     implicit none

! local variables
! status flag
     integer  :: istat

! desired charge density
     real(dp) :: ndens

! dummy array, used to save the eigenvalues of H + \Sigma(i\omega_n)
     complex(dp), allocatable :: eigs(:,:,:,:)

! dummy array, used to save the eigenvalues of H + \Sigma(oo)
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
     call cal_nelect(ndens)

! construct H + \Sigma, diagonalize it to obtain the dft + dmft eigenvalues
     call cal_eigsys(eigs, einf)

! search the fermi level using bisection algorithm
! the global variable `fermi` will be updated within `dichotomy()`
     call dichotomy(ndens, eigs, einf)

! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

     return
  end subroutine cal_fermi

!!
!! @sub cal_eimps
!!
!! try to calculate local energy levels for all impurity sites. here,
!! eimps is defined as \sum_k \epsilon_{n,k} - \mu.
!!
  subroutine cal_eimps()
     use constants, only : dp, mystd
     use constants, only : czero

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : fermi
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qdim
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
     use context, only : enk
     use context, only : eimps

     implicit none

! local variables
! loop index for spins
     integer :: s

! loop index for k-points
     integer :: k

! index for impurity sites
     integer :: t

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

! dummy array, used to build site-dependent impurity level
     complex(dp), allocatable :: Xe(:,:)

! dummy array, used to perform mpi reduce operation for eimps
     complex(dp), allocatable :: eimps_mpi(:,:,:,:)

! allocate memory
     allocate(Xe(qdim,qdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(eimps_mpi(qdim,qdim,nspin,nsite), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_eimps','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! init cbnd and cdim
! cbnd will be k-dependent and cdim will be impurity-dependent. they will
! be updated later
     cbnd = 0
     cdim = 0

! reset eimps
     eimps = czero
     eimps_mpi = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate eimps for', nsite, 'sites'
         write(mystd,'(4X,a,2X,i4,2X,a)') 'add contributions from', nkpt, 'kpoints'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
! see remarks in cal_nelect()
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
                 call s_print_error('cal_eimps','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

! evaluate Em, which is just some dft eigenvalues
             Em = enk(bs:be,k,s) - fermi

! convert `Em` to diagonal matrix `Hm`
             call s_diag_z(cbnd, Em, Hm)

! project effective hamiltonian to local basis
             do t=1,nsite
                 Xe = czero
                 cdim = ndim(t)
                 call one_psi_chi(cbnd, cdim, k, s, t, Hm, Xe(1:cdim,1:cdim))
                 eimps(:,:,s,t) = eimps(:,:,s,t) + Xe * weight(k)
             enddo ! over t={1,nsite} loop

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
     call mp_allreduce(eimps, eimps_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     eimps_mpi = eimps

# endif /* MPI */

! renormalize the impurity levels
     eimps = eimps_mpi / float(nkpt)

! deallocate memory
     if ( allocated(Xe) ) deallocate(Xe)
     if ( allocated(eimps_mpi) ) deallocate(eimps_mpi)

     return
  end subroutine cal_eimps

!!
!! @sub cal_eimpx
!!
!! try to calculate local energy levels for all impurity sites. here,
!! eimpx is equal to eimps - sigdc.
!!
  subroutine cal_eimpx()
     use control, only : nspin
     use control, only : nsite

     use context, only : ndim
     use context, only : eimps, eimpx
     use context, only : sigdc

     implicit none

! local variables
! loop index for spins
     integer :: s

! index for impurity sites
     integer :: t

! number of correlated orbitals for given impurity site
     integer :: cdim

! substract the double counting terms from eimps to build eimpx
     do t=1,nsite
         do s=1,nspin
             cdim = ndim(t)
             eimpx(1:cdim,1:cdim,s,t) = eimps(1:cdim,1:cdim,s,t) - sigdc(1:cdim,1:cdim,s,t)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

     return
  end subroutine cal_eimpx

!!
!! @sub cal_grn_l
!!
!! try to calculate local green's function for all the impurity sites
!!
  subroutine cal_grn_l()
     use constants, only : dp, mystd
     use constants, only : czero, czi

     use mmpi, only : mp_barrier
     use mmpi, only : mp_allreduce

     use control, only : nkpt, nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : myid, master, nprocs

     use context, only : i_wnd
     use context, only : qdim
     use context, only : ndim
     use context, only : kwin
     use context, only : weight
     use context, only : grn_l

     implicit none

! local variables
! loop index for spin
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

! dummy array: for self-energy function (upfolded to Kohn-Sham basis)
     complex(dp), allocatable :: Sk(:,:,:)
     complex(dp), allocatable :: Xk(:,:,:)

! dummy array: for lattice green's function
     complex(dp), allocatable :: Gk(:,:,:)

! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:,:)

! dummy array: used to perform mpi reduce operation for grn_l
     complex(dp), allocatable :: grn_l_mpi(:,:,:,:,:)

! allocate memory for Gl
     allocate(Gl(qdim,qdim,nmesh), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_grn_l','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(grn_l_mpi(qdim,qdim,nmesh,nspin,nsite), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_grn_l','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! init cbnd and cdim
! cbnd will be k-dependent and cdim will be impurity-dependent. we will
! update them later.
     cbnd = 0
     cdim = 0

! reset grn_l
     grn_l = czero
     grn_l_mpi = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate grn_l for', nsite, 'sites'
         write(mystd,'(4X,a,2X,i4,2X,a)') 'add contributions from', nkpt, 'kpoints'
     endif ! back if ( myid == master ) block

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
! see remarks in cal_nelect() for more details
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

! allocate memories Sk, Xk, and Gk. their sizes are k-dependent
             allocate(Sk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Xk(cbnd,cbnd,nmesh), stat = istat)
             allocate(Gk(cbnd,cbnd,nmesh), stat = istat)
             !
             if ( istat /= 0 ) then
                 call s_print_error('cal_grn_l','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

! build self-energy function, and then upfold it into Kohn-Sham basis
! Sk should contain contributions from all impurity sites
             Sk = czero
             do t=1,nsite
                 Xk = czero
                 cdim = ndim(t)
                 call cal_sl_sk(cdim, cbnd, k, s, t, Xk)
                 Sk = Sk + Xk
             enddo ! over t={1,nsite} loop

! calculate lattice green's function
             call cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)
             print *, Gk(:,:,3)
             STOP

! downfold the lattice green's function to obtain local green's function,
! then we have to save the final results
             do t=1,nsite
                 Gl = czero
                 cdim = ndim(t)
                 call cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl(1:cdim,1:cdim,:))
                 grn_l(:,:,:,s,t) = grn_l(:,:,:,s,t) + Gl * weight(k)
             enddo ! over t={1,nsite} loop

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
     call mp_allreduce(grn_l, grn_l_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     grn_l_mpi = grn_l

# endif /* MPI */

! renormalize local green's function
     grn_l = grn_l_mpi / float(nkpt)

! deallocate memory
     if ( allocated(Gl) ) deallocate(Gl)
     if ( allocated(grn_l_mpi) ) deallocate(grn_l_mpi)

!! DEBUG CODE
     do s=1,qdim
         print *, s, grn_l(s,s,1,1,1)
     enddo ! over s={1,qdim} loop
!! DEBUG CODE

     return
  end subroutine cal_grn_l

!!
!! @sub cal_wss_l
!!
!! try to calculate local weiss's function for all impurity sites
!!
  subroutine cal_wss_l()
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nspin
     use control, only : nsite, nmesh
     use control, only : myid, master

     use context, only : ndim
     use context, only : sig_l
     use context, only : grn_l
     use context, only : wss_l

     implicit none

! local variables
! loop index for frequency mesh
     integer :: m

! loop index for spins
     integer :: s

! loop index for impurity sites
     integer :: t

! number of correlated orbitals for given impurity site
     integer :: cdim

! status flag
     integer :: istat

! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:)

! reset wss_l
     wss_l = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate wss_l for', nsite, 'sites'
     endif ! back if ( myid == master ) block

!
! try to calculate bath weiss's function using the following equation:
!     G^{-1}_{0} = G^{-1} + \Sigma
! please be aware that the double counting terms have been substracted
! from the self-energy function. see subroutine cal_sig_l().
!
     SITE_LOOP: do t=1,nsite
! get size of orbital space
         cdim = ndim(t)

! allocate memory
         allocate(Gl(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_wss_l','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block
         !
         Gl = czero

! loop over spins and frequency mesh
         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

! back local green's function to Gl
                 Gl = grn_l(1:cdim,1:cdim,m,s,t)

! inverse local green's function. now Gl is G^{-1}
                 call s_inv_z(cdim, Gl)

! plus the self-energy function. now Gl is G^{-1} + \Sigma
                 Gl = Gl + sig_l(1:cdim,1:cdim,m,s,t)

! inverse it again to obtain bath weiss's function. now Gl is G_0
                 call s_inv_z(cdim, Gl)

! save the final resuls to wss_l
                 wss_l(1:cdim,1:cdim,m,s,t) = Gl

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

! deallocate memory
         if ( allocated(Gl) ) deallocate(Gl)

     enddo SITE_LOOP ! over t={1,nsite} loop

     return
  end subroutine cal_wss_l

!!
!! @sub cal_hyb_l
!!
!! try to calculate hybridization function for all impurity sites
!!
  subroutine cal_hyb_l()
     use constants, only : dp, mystd
     use constants, only : czi, czero

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh
     use control, only : myid, master

     use context, only : ndim
     use context, only : fmesh
     use context, only : eimps
     use context, only : sig_l
     use context, only : grn_l
     use context, only : hyb_l

     implicit none

! local variables
! loop index for spin
     integer :: s

! index for impurity sites
     integer :: t

! loop index for frequency mesh
     integer :: m

! number of correlated orbitals for given impurity site
     integer :: cdim

! status flag
     integer :: istat

! dummy variables
     complex(dp) :: caux

! dummy arrays
     complex(dp), allocatable :: Im(:,:)
     complex(dp), allocatable :: Tm(:,:)
     complex(dp), allocatable :: Em(:,:)
     complex(dp), allocatable :: Sm(:,:)

! reset hyb_l
     hyb_l = czero

! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate hyb_l for', nsite, 'sites'
     endif ! back if ( myid == master ) block

     SITE_LOOP: do t=1,nsite

! determine dimensional parameter
         cdim = ndim(t)

! allocate memory
         allocate(Im(cdim,cdim), stat = istat)
         allocate(Tm(cdim,cdim), stat = istat)
         allocate(Em(cdim,cdim), stat = istat)
         allocate(Sm(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_hyb_l','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block

         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

! build identify
                 call s_identity_z(cdim, Im)

! get frequency point. note that the fermi level (chemical potential) is
! already included in the impurity levels `eimps`. so here we just ignore
! the fermi level.
                 if ( axis == 1 ) then
                     caux = czi * fmesh(m)
                 else
                     caux = fmesh(m)
                 endif ! back if ( axis == 1 ) block

! calculate G^{-1}
                 Tm = grn_l(1:cdim,1:cdim,m,s,t)
                 call s_inv_z(cdim, Tm)

! get self-energy function
! be aware that the double counting terms have been removed from the
! self-energy functions. see cal_sig_l() subroutine for more details.
                 Sm = sig_l(1:cdim,1:cdim,m,s,t)

! get local impurity levels. the local impurity levels are actually equal
! to \sum e_{nk} - \mu. see cal_eimps() subroutine for more details.
                 Em = eimps(1:cdim,1:cdim,s,t)

! assemble the hybridization function. actually, Sm + Tm is G^{-1}_0.
! please see cal_wss_l() subroutine for more details.
                 hyb_l(1:cdim,1:cdim,m,s,t) = caux * Im - Em - Sm - Tm

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

! deallocate memory
         if ( allocated(Im) ) deallocate(Im)
         if ( allocated(Tm) ) deallocate(Tm)
         if ( allocated(Em) ) deallocate(Em)
         if ( allocated(Sm) ) deallocate(Sm)

     enddo SITE_LOOP ! over t={1,nsite} loop

     do m=1,nmesh
         write(100,'(i5, 5f16.8)') m, fmesh(m), grn_l(1,1,m,1,1), hyb_l(1,1,m,1,1)
     enddo

     return
  end subroutine cal_hyb_l

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

     use context, only : sig_l

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

! we use Sl to store parts of sig_l. note that actually sigdc has been
! substracted from sig_l beforehand, see cal_sig_l() for more details.
     do m=1,nmesh
         Sl(:,:,m) = sig_l(1:cdim,1:cdim,m,s,t)
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
!! try to build H(k) + \Sigma(i\omega_n). here, \Sigma should contain
!! contributions from all impurity sites. so, only when nsite = 1, we
!! can use the output of cal_sl_sk() as the input of this subroutine.
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

! frequency-dependent effective hamiltonian at given k-point and spin
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
!! try to diagonalize the effective hamiltonian: H(k) + \Sigma(i\omega_n),
!! get all the complex eigenvalues. here, we just assumed the effective
!! hamiltonian is a general complex matrix.
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
!! try to upfold the asymptotic values of self-energy functions at high
!! frequency (i.e `sigoo`) from local basis to Kohn-Sham basis. note that
!! the double-counting terms have been substracted from sigoo beforehand.
!! please see cal_sigoo() for more details.
!!
  subroutine cal_sl_so(cdim, cbnd, k, s, t, So)
     use constants, only : dp
     use constants, only : czero

     use context, only : sigoo

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

! local variables
! status flag
     integer :: istat

! dummy array: for asymptotic self-energy function
     complex(dp), allocatable :: Sl(:,:)

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
!! try to diagonalize the effective hamiltonian: H(k) + \Sigma(\infty),
!! get all the complex eigenvalues. here, we just assumed the effective
!! hamiltonian is a general complex matrix.
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

!!========================================================================
!!>>> service subroutines: set 3                                       <<<
!!========================================================================

!!
!! @sub cal_sk_gk
!!
!! try to calculate lattice green's function at given k-point and spin.
!! this subroutine needs the self-energy function at Kohn-Sham basis (i.e
!! `Sk`), that is the reason why it is called `cal_sk_gk`. note that Sk
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
     complex(dp), allocatable :: Fm(:)
     complex(dp), allocatable :: Em(:)

! dummy array: for effective hamiltonian (diagonal matrix)
     complex(dp), allocatable :: Hm(:,:)

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

! evaluate Fm, which is k-dependent, but frequency-independent
! if you want to consider magnetic field, you can add your codes here
     Fm = fermi - enk(bs:be,k,s)

     FREQ_LOOP: do m=1,nmesh

! consider imaginary axis or real axis
         if ( axis == 1 ) then
             Em = czi * fmesh(m) + Fm
         else
             Em = fmesh(m) + Fm
         endif ! back if ( axis == 1 ) block

! convert Em (vector) to Hm (diagonal matrix)
         call s_diag_z(cbnd, Em, Hm)
         if ( m == 3) then
             print *, Hm
             STOP
         endif

! substract self-energy function from the hamiltonian
         Gk(:,:,m) = Hm - Sk(:,:,m)

! calculate lattice green's function by direct inversion
         call s_inv_z(cbnd, Gk(:,:,m))

     enddo FREQ_LOOP ! over m={1,nmesh} loop

! deallocate memory
     if ( allocated(Fm) ) deallocate(Fm)
     if ( allocated(Em) ) deallocate(Em)
     if ( allocated(Hm) ) deallocate(Hm)

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

     if ( myid == master ) then
         write(mystd,'(4X,a)') 'searching fermi level'
     endif ! back if ( myid == master ) block

! initialization, determine mu1, mu2, occ1, occ2, and sign
!
! (1) if sign < 0, it means occ1 < desired, we should push mu2 to higher
! energy. (2) if sign > 0, it means occ1 > desired. then mu1 will be the
! right boundary, and we should push mu2 to lower energy
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

! well, finally, we have to update the global variable `fermi`
     fermi = mu3

     return
  end subroutine dichotomy

!!
!! @sub cal_nelect
!!
!! try to calculate the number of valence electrons by dft occupations
!! and weights. actually, what we obtain is the occupation numbers in
!! the selected band window.
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

! perform summation
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

! check axis
     call s_assert2(axis == 1, 'axis is wrong')

! calculate local green's function
     gloc = czero
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
                     gloc(b,m,s) = gloc(b,m,s) + one / ( caux - eigs(b,m,k,s) )
                     gloc(b,m,s) = gloc(b,m,s) - one / ( caux - einf(b,k,s) )
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
! see remarks in cal_nelect()
             bs = kwin(k,s,1,i_wnd(1))
             be = kwin(k,s,2,i_wnd(1))
             cbnd = be - bs + 1
             do b=1,cbnd
                 caux = einf(b,k,s) - fermi
                 zocc(b,s) = zocc(b,s) + fermi_dirac( real(caux) ) / real(nkpt)
             enddo ! over b={1,cbnd} loop
         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

! actually, we should consider the correction due to finite frequency
! point here. later we will implement it.
!
! TO_BE_DONE

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

     if ( myid == master ) then
         write(mystd,'(4X,a)') 'calculating dft + dmft eigenvalues'
     endif ! back if ( myid == master ) block

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

! mpi barrier. waiting all processes reach here.
# if defined (MPI)
     !
     call mp_barrier()
     !
# endif /* MPI */

     SPIN_LOOP: do s=1,nspin
         KPNT_LOOP: do k=myid+1,nkpt,nprocs

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
!
! see remarks in cal_nelect() subroutine
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
!! basis. you can call this procedure `embedding` or `upfolding`.
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
     if ( istat /= 0 ) then
         call s_print_error('map_chi_psi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
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
     if ( istat /= 0 ) then
         call s_print_error('map_psi_chi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
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
     if ( istat /= 0 ) then
         call s_print_error('one_chi_psi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
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
     if ( istat /= 0 ) then
         call s_print_error('one_psi_chi','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
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
