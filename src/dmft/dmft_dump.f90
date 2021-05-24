!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : dmft_dump_fermi
!!!           dmft_dump_eimps
!!!           dmft_dump_eimpx
!!!           dmft_dump_eigen
!!!           dmft_dump_grn_l
!!!           dmft_dump_wss_l
!!!           dmft_dump_hyb_l
!!! source  : dmft_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 02/23/2021 by li huang (created)
!!!           05/23/2021 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dmft_dump_fermi
!!
!! write out fermi level
!!
  subroutine dmft_dump_fermi(fermi)
     use constants, only : dp
     use constants, only : mytmp

     implicit none

! external arguments
! fermi level
     real(dp), intent(in) :: fermi

! open data file: dmft.fermi
     open(mytmp, file='dmft.fermi', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,f16.8)') '# fermi: ', fermi

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_fermi

!!
!! @sub dmft_dump_eimps
!!
!! write out local impurity levels, eimps. note that the double counting
!! terms have NOT been substracted from them.
!!
  subroutine dmft_dump_eimps(eimps)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite

     use context, only : qdim
     use context, only : ndim

     implicit none

! external arguments
! local impurity levels
     complex(dp), intent(in) :: eimps(qdim,qdim,nspin,nsite)

! local variables
! loop index for impurity sites
     integer :: t

! loop index for spins
     integer :: s

! loop index for correlated orbitals
     integer :: p, q

! open data file: dmft.eimps
     open(mytmp, file='dmft.eimps', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', t, 'spin:', s, 'dims:', ndim(t)
             do q=1,ndim(t)
                 do p=1,ndim(t)
                     write(mytmp,'(2i4,2f16.8)') p, q, eimps(p,q,s,t)
                 enddo ! over p={1,ndim(t)} loop
             enddo ! over q={1,ndim(t)} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_eimps

!!
!! @sub dmft_dump_eimpx
!!
!! write out local impurity levels, eimps. note that the double counting
!! terms have been substracted from them.
!!
  subroutine dmft_dump_eimpx(eimpx)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite

     use context, only : qdim
     use context, only : ndim

     implicit none

! external arguments
! local impurity levels. eimpx = eimps - sigdc
     complex(dp), intent(in) :: eimpx(qdim,qdim,nspin,nsite)

! local variables
! loop index for impurity sites
     integer :: t

! loop index for spins
     integer :: s

! loop index for correlated orbitals
     integer :: p, q

! open data file: dmft.eimpx
     open(mytmp, file='dmft.eimpx', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', t, 'spin:', s, 'dims:', ndim(t)
             do q=1,ndim(t)
                 do p=1,ndim(t)
                     write(mytmp,'(2i4,2f16.8)') p, q, eimpx(p,q,s,t)
                 enddo ! over p={1,ndim(t)} loop
             enddo ! over q={1,ndim(t)} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_eimpx

!!
!! @sub dmft_dump_eigen
!!
!! write out complex dft + dmft eigenvalues
!!
  subroutine dmft_dump_eigen(eigs)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nkpt, nspin
     use control, only : nmesh

     use context, only : i_wnd
     use context, only : qbnd
     use context, only : kwin
     use context, only : fmesh

     implicit none

! external arguments
! eigenvalues for H(k) + \Sigma(i\omega_n)
     complex(dp), intent(in) :: eigs(qbnd,nmesh,nkpt,nspin)

! local variables
! loop index for k-points
     integer :: k

! loop index for spins
     integer :: s

! loop index for impurity sites
     integer :: t

! loop index for frequency grid
     integer :: m

! loop index for bands in band window
     integer :: q

! number of dft bands for given k-point and spin
     integer :: cbnd

! band window: start index and end index for bands
     integer :: bs, be

! open data file: dmft.eigen
     open(mytmp, file='dmft.eigen', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nkpt : ', nkpt
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qbnd : ', qbnd

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do s=1,nspin
         do k=1,nkpt

! evaluate band window for the current k-point and spin
! i_wnd(t) returns the corresponding band window for given impurity site t
!
! see remarks in cal_nelect() subroutine
             t = 1 ! t is fixed to 1
             bs = kwin(k,s,1,i_wnd(t))
             be = kwin(k,s,2,i_wnd(t))

! determine cbnd
             cbnd = be - bs + 1

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# kpt:', k, 'spin:', s, 'cbnd:', cbnd
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 do q=1,cbnd
                     write(mytmp,'(i4,2f16.8)') q, eigs(q,m,k,s)
                 enddo ! over q={1,cbnd} loop
             enddo ! over m={1,nmesh} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over k={1,nkpt} loop
     enddo ! over s={1,nspin} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_eigen

!!
!! @sub dmft_dump_grn_l
!!
!! write out local green's function in matsubara frequency space
!!
  subroutine dmft_dump_grn_l(grn_l)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite, nmesh

     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

! external arguments
! local green's function
     complex(dp), intent(in) :: grn_l(qdim,qdim,nmesh,nspin,nsite)

! local variables
! loop index for impurity sites
     integer :: t

! loop index for spins
     integer :: s

! loop index for frequency grid
     integer :: m

! loop index for correlated orbitals
     integer :: p, q

! open data file: dmft.grn_l
     open(mytmp, file='dmft.grn_l', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', t, 'spin:', s, 'dims:', ndim(t)
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, grn_l(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_grn_l

!!
!! @sub dmft_dump_wss_l
!!
!! write out local weiss's function in matsubara frequency space
!!
  subroutine dmft_dump_wss_l(wss_l)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite, nmesh

     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

! external arguments
! local weiss's function
     complex(dp), intent(in) :: wss_l(qdim,qdim,nmesh,nspin,nsite)

! local variables
! loop index for impurity sites
     integer :: t

! loop index for spins
     integer :: s

! loop index for frequency grid
     integer :: m

! loop index for correlated orbitals
     integer :: p, q

! open data file: dmft.wss_l
     open(mytmp, file='dmft.wss_l', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', t, 'spin:', s, 'dims:', ndim(t)
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, wss_l(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_wss_l

!!
!! @sub dmft_dump_hyb_l
!!
!! write out local hybridization function in matsubara frequency space
!!
  subroutine dmft_dump_hyb_l(hyb_l)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nspin
     use control, only : nsite, nmesh

     use context, only : qdim
     use context, only : ndim
     use context, only : fmesh

     implicit none

! external arguments
! local hybridization function
     complex(dp), intent(in) :: hyb_l(qdim,qdim,nmesh,nspin,nsite)

! local variables
! loop index for impurity sites
     integer :: t

! loop index for spins
     integer :: s

! loop index for frequency grid
     integer :: m

! loop index for correlated orbitals
     integer :: p, q

! open data file: dmft.hyb_l
     open(mytmp, file='dmft.hyb_l', form='formatted', status='unknown')

! write parameters
     write(mytmp,'(a9,i4)') '# nsite: ', nsite
     write(mytmp,'(a9,i4)') '# nspin: ', nspin
     write(mytmp,'(a9,i4)') '# nmesh: ', nmesh
     write(mytmp,'(a9,i4)') '# qdim : ', qdim

! write separators
     write(mytmp,*)
     write(mytmp,*)

! write body
     do t=1,nsite
         do s=1,nspin

! write data for given spin and site
             write(mytmp,'(3(a,i4,2X))') '# site:', t, 'spin:', s, 'dims:', ndim(t)
             do m=1,nmesh
                 write(mytmp,'(a2,i6,f16.8)') 'w:', m, fmesh(m)
                 do q=1,ndim(t)
                     do p=1,ndim(t)
                         write(mytmp,'(2i4,2f16.8)') p, q, hyb_l(p,q,m,s,t)
                     enddo ! over p={1,ndim(t)} loop
                 enddo ! over q={1,ndim(t)} loop
             enddo ! over m={1,nmesh} loop

! write separators
             write(mytmp,*)
             write(mytmp,*)

         enddo ! over s={1,nspin} loop
     enddo ! over t={1,nsite} loop

! close data file
     close(mytmp)

     return
  end subroutine dmft_dump_hyb_l
