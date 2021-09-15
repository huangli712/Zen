!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : cal_sigoo
!!!           cal_sigma
!!!           cal_fermi
!!!           cal_eimps
!!!           cal_eimpx
!!!           cal_green
!!!           cal_weiss
!!!           cal_delta
!!!           cal_gamma
!!! source  : dmft_flow.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/29/2021 by li huang (created)
!!!           09/15/2021 by li huang (last modified)
!!! purpose : implement the main work flow of dft + dmft calculation.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub cal_sigoo
!!
!! try to calculate the asymptotic values for self-energy functions. then
!! the double-counting terms will be removed as well. this function works
!! for Matsubara self-energy functions (bare) only.
!!
  subroutine cal_sigoo()
     use constants, only : dp
     use constants, only : czero

     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : qdim
     use context, only : sigdc, sigoo, sigma

     implicit none

!! local parameters
     ! how many frequency points are included to calculate the asymptotic
     ! values of Matsubara self-energy function.
     integer, parameter :: mcut = 16

!! local variables
     ! loop index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! status flag
     integer :: istat

     ! dummy array for the indices of frequency points
     integer, allocatable :: ip(:)

     ! dummy array for the Matsubara self-energy functions
     complex(dp), allocatable :: Sm(:,:)

!! [body

     ! allocate memory
     allocate(ip(mcut), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sigoo','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(Sm(qdim,qdim), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_sigoo','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! check working axis.
     ! now only the Matsubara frequency axis is supported.
     call s_assert2(axis == 1, 'axis is wrong')

     ! reset sigoo
     sigoo = czero

     ! build integer array for indices of frequency points
     call s_linspace_i(nmesh + 1 - mcut, nmesh, mcut, ip)

     ! loop over quantum impurities, spins, and frequency points.
     !
     ! we count the last `mcut` frequency points, then we try to
     ! calculate the averaged values. up to now, the double counting
     ! terms have not been substracted from sigma. in other words,
     ! sigma is still bare.
     do t=1,ngrp
         do s=1,nspin
             Sm = czero
             !
             do m=1,mcut
                 Sm = Sm + sigma(:,:,ip(m),s,t)
             enddo ! over m={1,mcut} loop
             !
             sigoo(:,:,s,t) = Sm / float(mcut)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,ngrp} loop

     ! we substract the double counting terms from sigoo
     sigoo = sigoo - sigdc

     ! deallocate memory
     if ( allocated(ip) ) deallocate(ip)
     if ( allocated(Sm) ) deallocate(Sm)

!! body]

     return
  end subroutine cal_sigoo

!!
!! @sub cal_sigma
!!
!! try to substract the double counting terms from the bare Matsubara
!! self-energy functions. this function works for Matsubara self-energy
!! functions (bare) only.
!!
  subroutine cal_sigma()
     use control, only : axis
     use control, only : nspin
     use control, only : nsite
     use control, only : nmesh

     use context, only : sigdc, sigma

     implicit none

!! local variables
     ! loop index for frequency mesh
     integer :: m

     ! loop index for spins
     integer :: s

     ! loop index for impurity sites
     integer :: t

!! [body

     ! check working axis.
     ! now only the Matsubara frequency axis is supported.
     call s_assert2(axis == 1, 'axis is wrong')

     ! loop over quantum impurities, spins, and frequency points.
     !
     ! substract the double counting terms: new sigma = sigma - sigdc.
     do t=1,nsite
         do s=1,nspin
             do m=1,nmesh
                 sigma(:,:,m,s,t) = sigma(:,:,m,s,t) - sigdc(:,:,s,t)
             enddo ! over m={1,nmesh} loop
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_sigma

!!
!! @sub cal_fermi
!!
!! try to determine the fermi level. in addition, the lattice occupancy
!! will be calculated at the same time.
!!
  subroutine cal_fermi(occup)
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nkpt, nspin
     use control, only : nmesh

     use context, only : qbnd

     implicit none

!! external arguments
     ! lattice occupancy
     real(dp), intent(out) :: occup

!! local variables
     ! desired charge density
     real(dp) :: ndens

     ! status flag
     integer  :: istat

     ! dummy array, used to save the eigenvalues of H + \Sigma(i\omega_n)
     complex(dp), allocatable :: eigs(:,:,:,:)

     ! dummy array, used to save the eigenvalues of H + \Sigma(oo)
     complex(dp), allocatable :: einf(:,:,:)

!! [body

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

     ! calculate the nominal charge density according to the raw
     ! dft eigenvalues.
     call cal_nelect(ndens); occup = ndens

     ! construct H + \Sigma, then diagonalize it to obtain the
     ! dft + dmft eigenvalues.
     call cal_eigsys(eigs, einf)

     ! search the fermi level using bisection algorithm.
     ! the global variable `fermi` will be updated within `dichotomy()`.
     call dichotomy(ndens, eigs, einf)

     ! deallocate memory
     if ( allocated(eigs) ) deallocate(eigs)
     if ( allocated(einf) ) deallocate(einf)

!! body]

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

!! local variables
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

!! [body

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

     ! reset cbnd and cdim. they will be updated later.
     ! cbnd should be k-dependent and cdim should be impurity-dependent.
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

             ! evaluate band window for the current k-point and spin.
             !
             ! i_wnd(t) returns the corresponding band window for given
             ! impurity site t. see remarks in cal_nelect().
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

             ! evaluate `Em`, which is the eigenvalues substracted
             ! by the fermi level.
             Em = enk(bs:be,k,s) - fermi

             ! convert `Em` to diagonal matrix `Hm`
             call s_diag_z(cbnd, Em, Hm)

             ! project effective hamiltonian from the Kohn-Sham basis
             ! to the local basis, and then sum it up.
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

!! body]

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

!! local variables
     ! index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! number of correlated orbitals for given impurity site
     integer :: cdim

!! [body

     ! substract the double counting terms from eimps to build eimpx
     do t=1,nsite
         do s=1,nspin
             cdim = ndim(t)
             eimpx(1:cdim,1:cdim,s,t) = eimps(1:cdim,1:cdim,s,t) - sigdc(1:cdim,1:cdim,s,t)
         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_eimpx

!!
!! @sub cal_green
!!
!! try to calculate local green's function for all the impurity sites.
!!
  subroutine cal_green()
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
     use context, only : green

     implicit none

!! local variables
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

     ! dummy array: used to perform mpi reduce operation for green
     complex(dp), allocatable :: green_mpi(:,:,:,:,:)

!! [body

     ! allocate memory for Gl and green_mpi
     allocate(Gl(qdim,qdim,nmesh), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
     !
     allocate(green_mpi(qdim,qdim,nmesh,nspin,nsite), stat = istat)
     if ( istat /= 0 ) then
         call s_print_error('cal_green','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! reset cbnd and cdim. they will be updated later.
     ! cbnd should be k-dependent and cdim should be impurity-dependent.
     cbnd = 0
     cdim = 0

     ! reset green
     green = czero
     green_mpi = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate green for', nsite, 'sites'
         write(mystd,'(4X,a,2X,i4,2X,a)') 'add contributions from', nkpt, 'kpoints'
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
                 call s_print_error('cal_green','can not allocate enough memory')
             endif ! back if ( istat /= 0 ) block

             ! build self-energy function, and then upfold it into
             ! Kohn-Sham basis. Sk should contain contributions from
             ! all impurity sites.
             Sk = czero
             do t=1,nsite
                 Xk = czero ! reset Xk
                 cdim = ndim(t)
                 call cal_sl_sk(cdim, cbnd, k, s, t, Xk)
                 Sk = Sk + Xk
             enddo ! over t={1,nsite} loop

             ! calculate lattice green's function
             call cal_sk_gk(cbnd, bs, be, k, s, Sk, Gk)

             ! downfold the lattice green's function to obtain local
             ! green's function, then we have to perform k-summation.
             do t=1,nsite
                 Gl = czero
                 cdim = ndim(t)
                 call cal_gk_gl(cbnd, cdim, k, s, t, Gk, Gl(1:cdim,1:cdim,:))
                 green(:,:,:,s,t) = green(:,:,:,s,t) + Gl * weight(k)
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
     call mp_allreduce(green, green_mpi)
     !
     call mp_barrier()
     !
# else  /* MPI */

     green_mpi = green

# endif /* MPI */

     ! renormalize local green's function
     green = green_mpi / float(nkpt)

     ! deallocate memory
     if ( allocated(Gl) ) deallocate(Gl)
     if ( allocated(green_mpi) ) deallocate(green_mpi)

!! body]

     return
  end subroutine cal_green

!!
!! @sub cal_weiss
!!
!! try to calculate local weiss's function for all impurity sites.
!!
  subroutine cal_weiss()
     use constants, only : dp, mystd
     use constants, only : czero

     use control, only : nspin
     use control, only : nsite, nmesh
     use control, only : myid, master

     use context, only : ndim
     use context, only : sigma
     use context, only : green
     use context, only : weiss

     implicit none

!! local variables
     ! loop index for impurity sites
     integer :: t

     ! loop index for spins
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! status flag
     integer :: istat

     ! dummy array: for local green's function
     complex(dp), allocatable :: Gl(:,:)

!! [body

     ! reset weiss
     weiss = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate weiss for', nsite, 'sites'
     endif ! back if ( myid == master ) block

!
! remarks:
!
! try to calculate bath weiss's function using the following equation:
!
!     G^{-1}_{0} = G^{-1} + \Sigma
!
! please be aware that the double counting terms have been substracted
! from the self-energy function. see subroutine cal_sigma().
!

     ! loop over quantum impurities
     SITE_LOOP: do t=1,nsite
         ! get size of orbital space
         cdim = ndim(t)

         ! allocate memory
         allocate(Gl(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_weiss','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block
         !
         Gl = czero

         ! loop over spins and frequency mesh
         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

                 ! copy local green's function to Gl
                 Gl = green(1:cdim,1:cdim,m,s,t)

                 ! inverse local green's function.
                 ! now Gl is G^{-1}.
                 call s_inv_z(cdim, Gl)

                 ! plus the self-energy function.
                 ! now Gl is G^{-1} + \Sigma.
                 Gl = Gl + sigma(1:cdim,1:cdim,m,s,t)

                 ! inverse it again to obtain bath weiss's function.
                 ! now Gl is G_0.
                 call s_inv_z(cdim, Gl)

                 ! save the final resuls to weiss
                 weiss(1:cdim,1:cdim,m,s,t) = Gl

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

         ! deallocate memory
         if ( allocated(Gl) ) deallocate(Gl)

     enddo SITE_LOOP ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_weiss

!!
!! @sub cal_delta
!!
!! try to calculate hybridization function for all impurity sites.
!!
  subroutine cal_delta()
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
     use context, only : sigma
     use context, only : green
     use context, only : delta

     implicit none

!! local variables
     ! index for impurity sites
     integer :: t

     ! loop index for spin
     integer :: s

     ! loop index for frequency mesh
     integer :: m

     ! number of correlated orbitals for given impurity site
     integer :: cdim

     ! status flag
     integer :: istat

     ! dummy variables
     complex(dp) :: caux

     ! dummy arrays
     ! for identity matrix
     complex(dp), allocatable :: Im(:,:)

     ! for local impurity levels matrix
     complex(dp), allocatable :: Em(:,:)

     ! for green's function matrix
     complex(dp), allocatable :: Tm(:,:)

     ! for self-energy function matrix
     complex(dp), allocatable :: Sm(:,:)

!! [body

     ! reset delta
     delta = czero

     ! print some useful information
     if ( myid == master ) then
         write(mystd,'(4X,a,2X,i2,2X,a)') 'calculate delta for', nsite, 'sites'
     endif ! back if ( myid == master ) block

     ! loop over quantum impurities
     SITE_LOOP: do t=1,nsite

         ! determine dimensional parameter
         cdim = ndim(t)

         ! allocate memory
         allocate(Im(cdim,cdim), stat = istat)
         allocate(Em(cdim,cdim), stat = istat)
         allocate(Tm(cdim,cdim), stat = istat)
         allocate(Sm(cdim,cdim), stat = istat)
         !
         if ( istat /= 0 ) then
             call s_print_error('cal_delta','can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block

         ! loop over spins and frequency meshes
         SPIN_LOOP: do s=1,nspin
             MESH_LOOP: do m=1,nmesh

                 ! build identity matrix
                 call s_identity_z(cdim, Im)

                 ! get frequency point.
                 ! note that the fermi level (chemical potential) is
                 ! already included in the impurity levels `eimps`.
                 ! so here we just ignore the fermi level.
                 if ( axis == 1 ) then
                     caux = czi * fmesh(m)
                 else
                     caux = fmesh(m)
                 endif ! back if ( axis == 1 ) block

                 ! get local impurity levels.
                 ! the local impurity levels are actually equal to
                 !
                 !     \sum e_{nk} - \mu.
                 !
                 ! see cal_eimps() subroutine for more details.
                 Em = eimps(1:cdim,1:cdim,s,t)

                 ! calculate G^{-1}
                 Tm = green(1:cdim,1:cdim,m,s,t)
                 call s_inv_z(cdim, Tm)

                 ! get self-energy function.
                 ! be aware that the double counting terms have been
                 ! removed from the self-energy functions. see the
                 ! cal_sigma() subroutine for more details.
                 Sm = sigma(1:cdim,1:cdim,m,s,t)

                 ! assemble the hybridization function.
                 ! actually, Tm + Sm is G^{-1}_0.
                 ! please see cal_weiss() subroutine for more details.
                 delta(1:cdim,1:cdim,m,s,t) = caux * Im - Em - Tm - Sm

             enddo MESH_LOOP ! over m={1,nmesh} loop
         enddo SPIN_LOOP ! over s={1,nspin} loop

         ! deallocate memory
         if ( allocated(Im) ) deallocate(Im)
         if ( allocated(Em) ) deallocate(Em)
         if ( allocated(Tm) ) deallocate(Tm)
         if ( allocated(Sm) ) deallocate(Sm)

     enddo SITE_LOOP ! over t={1,nsite} loop

!! body]

     return
  end subroutine cal_delta

!!
!! @sub cal_gamma
!!
!! try to calculate correlation-induced correction for density matrix.
!!
  subroutine cal_gamma(ecorr)
     use constants, only : dp

     use control, only : nkpt, nspin

     use context, only : qbnd
     use context, only : gamma

     implicit none

!! external arguments
     ! correction to band energy
     real(dp), intent(out) :: ecorr

!! local variables
     ! status flag
     integer  :: istat

     ! dft + dmft density matrix
     complex(dp), allocatable :: kocc(:,:,:,:)

!! [body

     ! allocate memory
     allocate(kocc(qbnd,qbnd,nkpt,nspin), stat = istat)
     !
     if ( istat /= 0 ) then
         call s_print_error('cal_gamma','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! calculate new density matrix based
     ! on the dft + dmft eigenvalues.
     call cal_denmat(kocc)

     ! calculate the difference between dft
     ! and dft + dmft density matrices.
     call correction(kocc, gamma, ecorr)

     ! deallocate memory
     if ( allocated(kocc) ) deallocate(kocc)

!! body]

     return
  end subroutine cal_gamma
