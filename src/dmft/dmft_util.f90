!!!-----------------------------------------------------------------------
!!! project : jacaranda
!!! program : fermi_dirac
!!! source  : dmft_util.f90
!!! type    : functions
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/29/2021 by li huang (created)
!!!           07/30/2021 by li huang (last modified)
!!! purpose : define some useful utility functions and subroutines.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @fun fermi_dirac
!!
!! try to calculate the fermi-dirac function.
!!
  function fermi_dirac(omega) result(val)
     use constants, only : dp
     use constants, only : zero, one

     use control, only : beta

     implicit none

!! external arguments
     ! frequency point, \omega
     real(dp), intent(in) :: omega

     ! result value, return this
     real(dp) :: val

!! [body

     ! check the range of omega to avoid numerical instability
     if      ( beta * omega >=  600.0_dp ) then
         val = zero
     !
     else if ( beta * omega <= -600.0_dp ) then
         val = one
     !
     else
         val = one / ( one + exp( beta * omega ) )
     !
     endif ! back if ( beta * omega >=  600.0_dp ) block

!! body]

     return
  end function fermi_dirac
