!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : map_chi_psi
!!!           one_chi_psi
!!!           map_psi_chi
!!!           one_psi_chi
!!! source  : dmft_core.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/29/2021 by li huang (created)
!!!           07/30/2021 by li huang (last modified)
!!! purpose : service subroutines for upfolding and downfolding.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> service subroutines: upfolding                                   <<<
!!========================================================================

!!
!! @sub map_chi_psi
!!
!! service subroutine. map a matrix function from local basis to Kohn-Sham
!! basis. you can call this procedure `embedding` or `upfolding`.
!!
  subroutine map_chi_psi(cdim, cbnd, nfrq, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

!! external arguments
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

!! local variables
     ! loop index for frequency mesh
     integer :: f

     ! status flag
     integer :: istat

     ! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

!! [body

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

!! body]

     return
  end subroutine map_chi_psi

!!
!! @sub one_chi_psi
!!
!! service subroutine. map a single matrix from local basis to Kohn-Sham
!! basis. you can call this procedure `embedding` or `upfold`.
!!
  subroutine one_chi_psi(cdim, cbnd, k, s, t, Mc, Mp)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

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

     ! input matrix defined at local orbital (\chi) basis
     complex(dp), intent(in)  :: Mc(cdim,cdim)

     ! output matrix defined at Kohn-Sham (\psi) basis
     complex(dp), intent(out) :: Mp(cbnd,cbnd)

!! local variables
     ! status flag
     integer :: istat

     ! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

!! [body

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

!! body]

     return
  end subroutine one_chi_psi

!!========================================================================
!!>>> service subroutines: downfolding                                 <<<
!!========================================================================

!!
!! @sub map_psi_chi
!!
!! service subroutine. map a matrix function from Kohn-Sham basis to local
!! basis. you can call this procedure `projection` or `downfold`.
!!
  subroutine map_psi_chi(cbnd, cdim, nfrq, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

     implicit none

!! external arguments
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

!! local variables
     ! loop index for frequency mesh
     integer :: f

     ! status flag
     integer :: istat

     ! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

!! [body

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

!! body]

     return
  end subroutine map_psi_chi

!!
!! @sub one_psi_chi
!!
!! service subroutine. map a single matrix from Kohn-Sham basis to local
!! basis. you can call this procedure `projection` or `downfold`.
!!
  subroutine one_psi_chi(cbnd, cdim, k, s, t, Mp, Mc)
     use constants, only : dp

     use context, only : i_grp
     use context, only : chipsi
     use context, only : psichi

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

     ! input matrix defined at Kohn-Sham (\psi) basis
     complex(dp), intent(in)  :: Mp(cbnd,cbnd)

     ! output matrix defined at local orbital (\chi) basis
     complex(dp), intent(out) :: Mc(cdim,cdim)

!! local variables
     ! status flag
     integer :: istat

     ! overlap matrix between local orbitals and Kohn-Sham wave-functions
     complex(dp), allocatable :: Cp(:,:)
     complex(dp), allocatable :: Pc(:,:)

!! [body

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

!! body]

     return
  end subroutine one_psi_chi
