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
     use constants, only : dp
     use control, only : nspin

     implicit none

     integer :: s
     integer :: k
     integer :: p
     integer :: q
     integer :: f

     do s=1,nspin
         do k=1,nkpt

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
