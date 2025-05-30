!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : s_fft_tails
!!!           s_fft_forward
!!!           s_fft_backward
!!!           s_fft_density
!!! source  : s_fourier.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 07/10/2014 by li huang (created)
!!!           12/30/2024 by li huang (last modified)
!!! purpose : these subroutines are used to do fast fourier transformation
!!!           for green's or hybridization functions.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub s_fft_tails
!!
!! calculate high frequency tails using K. Haule's trick.
!!
  subroutine s_fft_tails(rtail, mfreq, rmesh, green)
     use constants, only : dp
     use constants, only : zero, one

     implicit none

!! external arguments
     ! number of matsubara frequency points
     integer, intent(in)     :: mfreq

     ! the desired high frequency tail
     real(dp), intent(out)   :: rtail

     ! matsubara frequency grid
     real(dp), intent(in)    :: rmesh(mfreq)

     ! function on matsubara frequency space
     complex(dp), intent(in) :: green(mfreq)

!! local parameters
     ! number of matsubara frequency points.
     ! these points are used to calculate high frequency tail.
     integer, parameter :: ntail = 128

!! local variables
     ! loop index
     integer  :: j

     ! dummy variables
     real(dp) :: Sn, Sx, Sy
     real(dp) :: Sxx, Sxy

!! [body

     Sn = zero
     Sx = zero
     Sy = zero
     !
     Sxx = zero
     Sxy = zero
     !
     do j=mfreq - ntail, mfreq
         Sn = Sn + one
         Sx = Sx + one / rmesh(j)**2
         Sy = Sy + aimag(green(j)) * rmesh(j)
         Sxx = Sxx + one / rmesh(j)**4
         Sxy = Sxy + aimag(green(j)) * rmesh(j) / rmesh(j)**2
     enddo ! over j={mfreq - nfreq, mfreq} loop

     rtail = (Sx * Sxy - Sxx * Sy) / (Sn * Sxx - Sx * Sx)

!! body]

     return
  end subroutine s_fft_tails

!!
!! @sub s_fft_forward
!!
!! fourier from imaginary time space forward to matsubara frequency space
!! using linear fourier algorithm.
!!
  subroutine s_fft_forward(ntime, tmesh, ftau, mfreq, rmesh, fmat)
     use constants, only : dp
     use constants, only : zero

     implicit none

!! external arguments
     ! number of imaginary time points
     integer, intent(in)  :: ntime

     ! number of matsubara frequency points
     integer, intent(in)  :: mfreq

     ! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

     ! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

     ! function on imaginary time axis
     real(dp), intent(in) :: ftau(ntime)

     ! function on matsubara frequency axis
     complex(dp), intent(out) :: fmat(mfreq)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! dummy variables
     real(dp) :: sre, sim

     real(dp) :: c0, c1
     real(dp) :: s0, s1
     real(dp) :: g0, g1, dg

!! [body

     do i=1,mfreq
         sre = zero
         sim = zero

         do j=1,ntime-1
             c0 = cos( tmesh(j)   * rmesh(i) )
             c1 = cos( tmesh(j+1) * rmesh(i) )
             s0 = sin( tmesh(j)   * rmesh(i) )
             s1 = sin( tmesh(j+1) * rmesh(i) )
             !
             g0 = ftau(j)
             g1 = ftau(j+1)
             !
             dg = ( g1 - g0 ) / ( tmesh(j+1) - tmesh(j) )
             !
             sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / rmesh(i) ) / rmesh(i)
             sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / rmesh(i) ) / rmesh(i)
         enddo ! over j={1,ntime-1} loop

         fmat(i) = dcmplx(sre, sim)
     enddo ! over i={1,mfreq} loop

!! body]

     return
  end subroutine s_fft_forward

!!
!! @sub s_fft_backward
!!
!! fourier from matsubara frequency space backward to imaginary time space.
!!
  subroutine s_fft_backward(mfreq, rmesh, fmat, ntime, tmesh, ftau, beta)
     use constants, only : dp
     use constants, only : pi, zero, two, half

     implicit none

!! external arguments
     ! number of matsubara frequency points
     integer, intent(in)   :: mfreq

     ! number of imaginary time points
     integer, intent(in)   :: ntime

     ! inverse temperature
     real(dp), intent(in)  :: beta

     ! matsubara frequency mesh
     real(dp), intent(in)  :: rmesh(mfreq)

     ! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

     ! function on imaginary time axis
     real(dp), intent(out) :: ftau(ntime)

     ! function on matsubara frequency axis
     complex(dp), intent(in) :: fmat(mfreq)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! real(dp) dummy variables
     real(dp) :: raux
     real(dp) :: tail

!! [body

     ! calculate high frequency tails need to be subtracted
     call s_fft_tails(tail, mfreq, rmesh, fmat)

     ! perform infourier transformation
     do i=1,ntime
         raux = zero
         do j=1,mfreq
             raux = raux + cos( rmesh(j) * tmesh(i) ) *   real( fmat(j) )
             raux = raux + sin( rmesh(j) * tmesh(i) ) * (aimag( fmat(j) ) + tail / rmesh(j))
         enddo ! over j={1,mfreq} loop
         ftau(i) = two * raux / beta - half * tail
     enddo ! over i={1,ntime} loop

     ! corrections for the boundary point
     raux = real( fmat(mfreq) ) * rmesh(mfreq) / pi
     ftau(1) = ftau(1) + raux
     ftau(ntime) = ftau(ntime) - raux

     ! additional corrections, may be useful for lda + dmft calculations
     ftau(1) = 3.0_dp * ftau(2) - 3.0_dp * ftau(3) + ftau(4)
     ftau(ntime) = 3.0_dp * ftau(ntime-1) - 3.0_dp * ftau(ntime-2) + ftau(ntime-3)

!! body]

     return
  end subroutine s_fft_backward

!!
!! @sub s_fft_density
!!
!! try to evaluate the occupation number via fourier transformation.
!!
  subroutine s_fft_density(mfreq, rmesh, fmat, density, beta)
     use constants, only : dp
     use constants, only : pi, zero, two, half

     implicit none

!! external arguments
     ! number of matsubara frequency points
     integer, intent(in)   :: mfreq

     ! inverse temperature
     real(dp), intent(in)  :: beta

     ! matsubara frequency mesh
     real(dp), intent(in)  :: rmesh(mfreq)

     ! function on matsubara frequency axis
     complex(dp), intent(in)  :: fmat(mfreq)

     ! complex occupation number
     complex(dp), intent(out) :: density

!! local variables
     ! loop index
     integer  :: j

     ! real(dp) dummy variables
     real(dp) :: raux
     real(dp) :: tail

!! [body

     ! calculate high frequency tails need to be subtracted
     call s_fft_tails(tail, mfreq, rmesh, fmat)

     ! perform inverted fourier transformation
     raux = zero
     !
     do j=1,mfreq
         raux = raux + real( fmat(j) )
     enddo ! over j={1,mfreq} loop
     !
     density = two * raux / beta - half * tail

     ! corrections for the boundary point
     raux = real( fmat(mfreq) ) * rmesh(mfreq) / pi
     density = density + raux

!! body]

     return
  end subroutine s_fft_density
