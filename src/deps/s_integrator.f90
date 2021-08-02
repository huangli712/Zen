!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : s_int_trapezoid
!!!           s_int_simpson
!!! source  : s_integrator.f90
!!! type    : functions
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/20/2014 by li huang (created)
!!!           07/29/2022 by li huang (last modified)
!!! purpose : the purpose of these functions is to implement the composite
!!!           trapezoid or composite simpson integration algorithms.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @fun s_int_trapezoid
!!
!! numerical integration with trapezoid algorithm.
!!
  function s_int_trapezoid(f, a, b, n) result(val)
     use constants, only : dp
     use constants, only : zero, two

     implicit none

!! external arguments
     ! number of data points
     integer, intent(in)  :: n

     ! left and right boundries for numerical integration
     real(dp), intent(in) :: a
     real(dp), intent(in) :: b

!! external function
     ! f means the integrand
     real(dp), external   :: f

!! local variables
     ! loop index
     integer  :: i

     ! return value
     real(dp) :: val

     ! step for integration
     real(dp) :: h

     ! sum for trapezoid rule
     real(dp) :: trapSum

!! [body

     ! evaluate the step
     h = ( b - a ) / dble(n)

     ! calculate trapezoid sum
     trapSum = zero
     do i=1,n-1
         trapSum = trapSum + f(a+dble(i)*h)
     enddo ! over i={1,n-1} loop

     ! calculate the final value
     val = ( h / two ) * ( f(a) + f(b) + two*trapSum )

!! body]

     return
  end function s_int_trapezoid

!!
!! @fun s_int_simpson
!!
!! numerical integration with simpson algorithm.
!!
  function s_int_simpson(f, a, b, n) result(val)
     use constants, only : dp
     use constants, only : zero

     implicit none

!! external arguments
     ! number of data points
     integer, intent(in)  :: n

     ! boundries for numerical integration
     real(dp), intent(in) :: a
     real(dp), intent(in) :: b

!! external function
     ! f means the integrand
     real(dp), external   :: f

!! local variables
     ! loop index
     integer  :: i

     ! return value
     real(dp) :: val

     ! step for integration
     real(dp) :: h

     ! sum for trapezoid rule
     real(dp) :: oddSum
     real(dp) :: evenSum

!! [body

     ! evaluate the step
     h = ( b - a ) / dble(n)

     ! calculate simpson sum
     evenSum = zero
     oddSum = zero
     !
     do i=1,n-1
         if ( mod(i,2) == 0 ) then
             evenSum = evenSum + f(a+dble(i)*h)
         else
             oddSum = oddSum + f(a+dble(i)*h)
         endif ! back if ( mod(i,2) == 0 ) block
     enddo ! over i={1,n-1} loop

     ! calculate the final value
     val = ( h / 3.0_dp ) * ( f(a) + f(b) + 2.0_dp * evenSum + 4.0_dp * oddsum )

!! body]

     return
  end function s_int_simpson
