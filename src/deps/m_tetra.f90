!!!-----------------------------------------------------------------------
!!! project : flink @ sakura
!!! program : mtetra
!!! source  : m_tetra.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 06/07/2006 by li huang (created)
!!!           07/29/2021 by li huang (last modified)
!!! purpose : this module is devoted to compute the analytical integration
!!!           weights for brillouin zone sampling.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module mtetra
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

!! module parameters
     ! kind for double precision
     integer, private, parameter  :: dp   = kind(1.0d0)

     ! well-known $\pi$
     real(dp), private, parameter :: pi   = 3.141592653589793238462643383279_dp

     ! 0.0 in double precision form
     real(dp), private, parameter :: zero = 0.00_dp

     ! $\epsilon$ in double precision form
     real(dp), private, parameter :: eps8 = 1.00E-8

     ! broadening parameter for smearing algorithm
     real(dp), private, parameter :: gamm = 0.15_dp

!!========================================================================
!!>>> declare common variables                                         <<<
!!========================================================================

!! module variables
     ! blochl corrections for dweight
     real(dp), save, private :: cweight

     ! density of states weights at the four corners of a given tetrahedron
     real(dp), save, private :: dweight(4)

     ! integration weights at the four corners of a given tetrahedron
     real(dp), save, private :: tweight(4)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! Blochl algorithm for (integrated) density of states
     public  :: tetra_blochl_weight1
     public  :: tetra_blochl_weight2

     ! Labmbin-Vigneron algorithm for crystal Green's function
     public  :: tetra_lambin_weight

     ! Gaussian broadening algorithm for (integrated) density of states
     public  :: smearing_gauss_weight1
     public  :: smearing_gauss_weight2
     public  :: smearing_gauss_weight3

     ! Fermi-Dirac broadening algorithm for (integrated) density of states
     public  :: smearing_fermi_weight1
     public  :: smearing_fermi_weight2
     public  :: smearing_fermi_weight3

     ! Marzari-Vanderbilt broadening algorithm for (integrated) density of states
     public  :: smearing_marzari_weight1
     public  :: smearing_marzari_weight2
     public  :: smearing_marzari_weight3

     ! private procedures invoked by tetra_blochl_weightX() subroutine
     private :: tetra_p_ek1
     private :: tetra_p_ek12
     private :: tetra_p_ek23
     private :: tetra_p_ek34
     private :: tetra_p_ek4

     ! private procedures invoked by tetra_lambin_weightX() subroutine
     private :: tetra_lv
     private :: tetra_lv_ekarb
     private :: tetra_lv_ek2idn
     private :: tetra_lv_ek22idn
     private :: tetra_lv_ek3idn
     private :: tetra_lv_ek4idn
     private :: tetra_lv_creorder
     private :: tetra_lv_prime
     private :: tetra_lv_setmap

  contains ! encapsulated functionality

!!========================================================================
!!>>> Blochl algorithm, driver layer                                   <<<
!!========================================================================

!!
!! @sub tetra_blochl_weight1
!!
!! Peter E. Blochl algorithm for (integrated) density of states and
!! relevant integration weights.
!!
!! note: s_qsorter() is implemented in s_util.f90.
!!
  subroutine tetra_blochl_weight1(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)    :: z

     ! density of states
     real(dp), intent(out)   :: dos

     ! integrated density of states
     real(dp), intent(out)   :: tos

     ! one-particle energies at the corners of the tetrahedron
     real(dp), intent(inout) :: e(4)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

     ! initialize tos and dos
     dos = zero
     tos = zero

     ! initialize dweight and tweight
     dweight = zero
     tweight = zero

     ! initialize cweight
     cweight = zero

     ! sort the corner energy according to increasing values
     call s_qsorter(4, e)

     ! remove degenerate
     do i=1,3
         if ( abs( e(i) - e(i+1) ) < eps8 ) then
             e(i) = e(i) + eps8 / real(i)
         endif
     enddo ! over i={1,3} loop
     !
     do i=1,4
         if ( abs( e(i) - z ) < eps8 / 10.0_dp ) then
             e(i) = e(i) + eps8 / 10.0_dp / real(i)
         endif
     enddo ! over i={1,4} loop

     ! find the case, to update dweight, tweight, and cweight.
     ! case 1, fully unoccupied tetrahedron
     if      ( z < e(1) ) then
         call tetra_p_ek1 (    )

     ! case 2, partially occupied tetrahedron
     else if ( z < e(2) .and. z > e(1) ) then
         call tetra_p_ek12(z, e)

     ! case 3, partially occupied tetrahedron
     else if ( z < e(3) .and. z > e(2) ) then
         call tetra_p_ek23(z, e)

     ! case 4, partially occupied tetrahedron
     else if ( z < e(4) .and. z > e(3) ) then
         call tetra_p_ek34(z, e)

     ! case 5, fully occupied tetrahedron
     else if ( z > e(4) ) then
         call tetra_p_ek4 (    )

     endif ! back if ( z < e(1) ) block

     ! add up Blochl corrections for density of states weights,
     ! apply equation (22).
     do i=1,4
         do j=1,4
             dweight(i) = dweight(i) + ( e(j) - e(i) ) * cweight * 0.025_dp
         enddo ! over j={1,4} loop
     enddo ! over i={1,4} loop

     ! compute density of states
     dos = sum( dweight )

     ! add up Blochl corrections for integration weights,
     ! apply equation (22).
     do i=1,4
         do j=1,4
             tweight(i) = tweight(i) + ( e(j) - e(i) ) * dos * 0.025_dp
         enddo ! over j={1,4} loop
     enddo ! over i={1,4} loop

     ! compute integrated density of states
     tos = sum( tweight )

!! body]

     return
  end subroutine tetra_blochl_weight1

!!
!! @sub tetra_blochl_weight2
!!
!! Peter E. Blochl algorithm for (integrated) density of states and
!! relevant integration weights.
!!
!! note: s_qsorter() is implemented in s_util.f90.
!!
  subroutine tetra_blochl_weight2(z, e, ddd, ttt)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)    :: z

     ! integration weight for density of states
     real(dp), intent(out)   :: ddd(4)

     ! integration weight for integrated density of states
     real(dp), intent(out)   :: ttt(4)

     ! one-particle energies at the corners of the tetrahedron
     real(dp), intent(inout) :: e(4)

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! density of states
     real(dp) :: dos

     ! integrated density of states
     real(dp) :: tos

!! [body

     ! initialize tos and dos
     dos = zero
     tos = zero

     ! initialize dweight and tweight
     dweight = zero
     tweight = zero

     ! initialize cweight
     cweight = zero

     ! sort the corner energy according to increasing values
     call s_qsorter(4, e)

     ! remove degenerate
     do i=1,3
         if ( abs( e(i) - e(i+1) ) < eps8 ) then
             e(i) = e(i) + eps8 / real(i)
         endif
     enddo ! over i={1,3} loop
     !
     do i=1,4
         if ( abs( e(i) - z ) < eps8 / 10.0_dp ) then
             e(i) = e(i) + eps8 / 10.0_dp / real(i)
         endif
     enddo ! over i={1,4} loop

     ! find the case, to update dweight, tweight, and cweight.
     ! case 1, fully unoccupied tetrahedron
     if      ( z < e(1) ) then
         call tetra_p_ek1 (    )

     ! case 2, partially occupied tetrahedron
     else if ( z < e(2) .and. z > e(1) ) then
         call tetra_p_ek12(z, e)

     ! case 3, partially occupied tetrahedron
     else if ( z < e(3) .and. z > e(2) ) then
         call tetra_p_ek23(z, e)

     ! case 4, partially occupied tetrahedron
     else if ( z < e(4) .and. z > e(3) ) then
         call tetra_p_ek34(z, e)

     ! case 5, fully occupied tetrahedron
     else if ( z > e(4) ) then
         call tetra_p_ek4 (    )

     endif ! back if ( z < e(1) ) block

     ! add up Blochl corrections for density of states weights,
     ! apply equation (22).
     do i=1,4
         do j=1,4
             dweight(i) = dweight(i) + ( e(j) - e(i) ) * cweight * 0.025_dp
         enddo ! over j={1,4} loop
     enddo ! over i={1,4} loop

     ! compute density of states
     dos = sum( dweight )
     ddd = dweight

     ! add up Blochl corrections for integration weights,
     ! apply equation (22).
     do i=1,4
         do j=1,4
             tweight(i) = tweight(i) + ( e(j) - e(i) ) * dos * 0.025_dp
         enddo ! over j={1,4} loop
     enddo ! over i={1,4} loop

     ! compute integrated density of states
     tos = sum( tweight )
     ttt = tweight

!! body]

     return
  end subroutine tetra_blochl_weight2

!!========================================================================
!!>>> Blochl algorithm, service layer                                  <<<
!!========================================================================

!!
!! @sub tetra_p_ek1
!!
!! Blochl algorithm, case 1, for fully unoccupied tetrahedron.
!!
  subroutine tetra_p_ek1()
     implicit none

!! [body

     ! integration weights, apply equation (B1)
     tweight = zero

     ! density of states weights
     dweight = zero

     ! corrections for dweight
     cweight = zero

     return
  end subroutine tetra_p_ek1

!!
!! @sub tetra_p_ek12
!!
!! Blochl algorithm, case 2, for partially occupied tetrahedron.
!!
  subroutine tetra_p_ek12(z, e)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in) :: z

     ! one-particle energies at the corners of the tetrahedron
     real(dp), intent(in) :: e(4)

!! local variables
     ! intermediated variables
     real(dp) :: c, dc

     ! ze_{i} = e - e_{i}
     real(dp) :: ze1

     ! e_{ij} = e_{i} - e_{j}
     real(dp) :: e21, e31, e41

!! [body

     ! setup common variables
     ze1 = z - e(1)
     e21 = e(2) - e(1)
     e31 = e(3) - e(1)
     e41 = e(4) - e(1)

     ! intermediate variable, apply equation (B6)
     c = ze1 * ze1 * ze1 / ( 4.0_dp * e21 * e31 * e41 )
     dc = 3.0_dp * ze1 * ze1 / ( 4.0_dp * e21 * e31 * e41 )

     ! integration weights
     ! apply equation (B2)
     tweight(1) = c * ( 4.0_dp - ze1 * ( 1.0_dp / e21 + 1.0_dp / e31 + 1.0_dp / e41 ) )

     ! apply equation (B3)
     tweight(2) = c * ze1 / e21

     ! apply equation (B4)
     tweight(3) = c * ze1 / e31

     ! apply equation (B5)
     tweight(4) = c * ze1 / e41

     ! density of states weights
     dweight(1) = 4.0_dp * dc - ( dc * ze1 + c ) * ( 1.0_dp / e21 + 1.0_dp / e31 + 1.0_dp / e41 )
     dweight(2) = dc * ze1 / e21 + c / e21
     dweight(3) = dc * ze1 / e31 + c / e31
     dweight(4) = dc * ze1 / e41 + c / e41

     ! corrections for dweight
     cweight = 6.0_dp * ze1 / ( e21 * e31 * e41 )

!! body]

     return
  end subroutine tetra_p_ek12

!!
!! @sub tetra_p_ek23
!!
!! Blochl algorithm, case 3, for partially occupied tetrahedron.
!!
  subroutine tetra_p_ek23(z, e)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in) :: z

     ! one-particle energies at the corners of the tetrahedron
     real(dp), intent(in) :: e(4)

!! local variables
     ! intermediated variables
     real(dp) ::  c1,  c2,  c3
     real(dp) :: dc1, dc2, dc3

     ! ze_{i} = e - e_{i}
     real(dp) :: ze1, ze2, ze3, ze4

     ! e_{ij} = e_{i} - e_{j}
     real(dp) :: e21, e31, e32, e41, e42

!! [body

     ! setup common variables
     ze1 = z - e(1)
     ze2 = z - e(2)
     ze3 = z - e(3)
     ze4 = z - e(4)
     !
     e21 = e(2) - e(1)
     e31 = e(3) - e(1)
     e32 = e(3) - e(2)
     e41 = e(4) - e(1)
     e42 = e(4) - e(2)

     ! intermediate variables
     ! apply equation (B11)
     c1  = ze1 * ze1 / ( 4.0_dp * e41 * e31 )
     dc1 = ze1 / ( 2.0_dp * e41 * e31 )

     ! apply equation (B12)
     c2  = - ze1 * ze2 * ze3 / ( 4.0_dp * e41 * e32 * e31 )
     dc2 = ( - ze2 * ze3 - ze1 * ze3 - ze1 * ze2 ) / ( 4.0_dp * e41 * e32 * e31 )

     ! apply equation (B13)
     c3  = - ze2 * ze2 * ze4 / ( 4.0_dp * e42 * e32 * e41 )
     dc3 = ( - 2.0_dp * ze2 * ze4 - ze2 * ze2 ) / ( 4.0_dp * e42 * e32 * e41 )

     ! integration weights
     ! apply equation (B7)
     tweight(1) = c1 - ( c1 + c2 ) * ze3 / e31 - ( c1 + c2 + c3 ) * ze4 / e41

     ! apply equation (B8)
     tweight(2) = c1 + c2 + c3 - ( c2 + c3 ) * ze3 / e32 - c3 * ze4 / e42

     ! apply equation (B9)
     tweight(3) = ( c1 + c2 ) * ze1 / e31 + ( c2 + c3 ) * ze2 / e32

     ! apply equation (B10)
     tweight(4) = ( c1 + c2 + c3 ) * ze1 / e41 + c3 * ze2 / e42

     ! density of states weights
     dweight(1) = dc1 - ( ( dc1 + dc2 ) * ze3 + c1 + c2 ) / e31 - ( ( dc1 + dc2 + dc3 ) * ze4 + c1 + c2 + c3 ) / e41
     dweight(2) = dc1 + dc2 + dc3 - ( ( dc2 + dc3 ) * ze3 + c2 + c3 ) / e32 - ( dc3 * ze4 + c3 ) / e42
     dweight(3) = ( ( dc1 + dc2 ) * ze1 + c1 + c2 ) / e31 + ( ( dc2 + dc3 ) * ze2 + c2 + c3 ) / e32
     dweight(4) = ( ( dc1 + dc2 + dc3 ) * ze1 + c1 + c2 + c3 ) / e41 + ( dc3 * ze2 + c3 ) / e42

     ! corrections for dweight
     cweight = 6.0_dp * ( 1.0_dp - ( e31 + e42 ) * ze2 / ( e32 * e42 ) ) / ( e31 * e41 )

!! body]

     return
  end subroutine tetra_p_ek23

!!
!! @sub tetra_p_ek34
!!
!! Blochl algorithm, case 4, for partially occupied tetrahedron.
!!
  subroutine tetra_p_ek34(z, e)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in) :: z

     ! one-particle energies at the corners of the tetrahedron
     real(dp), intent(in) :: e(4)

!! local variables
     ! intermediated variables
     real(dp) :: c, dc

     ! ze_{i} = e - e_{i}
     real(dp) :: ze4

     ! e_{ij} = e_{i} - e_{j}
     real(dp) :: e41, e42, e43

!! [body

     ! setup common variables
     ze4 = z - e(4)
     !
     e41 = e(4) - e(1)
     e42 = e(4) - e(2)
     e43 = e(4) - e(3)

     ! intermediate variables, apply equation (B18)
     c = - ze4 * ze4 * ze4 / ( 4.0_dp * e41 * e42 * e43 )
     dc = - 3.0_dp * ze4 * ze4 / ( 4.0_dp * e41 * e42 * e43 )

     ! integration weights
     ! apply equation (B14)
     tweight(1) = 0.25_dp + c * ze4 / e41

     ! apply equation (B15)
     tweight(2) = 0.25_dp + c * ze4 / e42

     ! apply equation (B16)
     tweight(3) = 0.25_dp + c * ze4 / e43

     ! apply equation (B17)
     tweight(4) = 0.25_dp - c * ( 4.0_dp + ( 1.0_dp / e41 + 1.0_dp / e42 + 1.0_dp / e43 ) * ze4 )

     ! density of states weights
     dweight(1) = ( dc * ze4 + c ) / e41
     dweight(2) = ( dc * ze4 + c ) / e42
     dweight(3) = ( dc * ze4 + c ) / e43
     dweight(4) = - 4.0_dp * dc - ( 1.0_dp / e41 + 1.0_dp / e42 + 1.0_dp / e43) * ( dc * ze4 + c )

     ! corrections for dweight
     cweight = 6.0_dp * ze4 / ( e41 * e42 * e43 )

!! body]

     return
  end subroutine tetra_p_ek34

!!
!! @sub tetra_p_ek4
!!
!! Blochl algorithm, case 5, for fully occupied tetrahedron.
!!
  subroutine tetra_p_ek4()
     implicit none

!! [body

     ! integration weights, apply equation (B19)
     tweight = 0.25_dp

     ! density of states weights
     dweight = zero

     ! corrections for dweight
     cweight = zero

!! body]

     return
  end subroutine tetra_p_ek4

!!========================================================================
!!>>> Lambin-Vigneron algorithm, driver layer                          <<<
!!========================================================================

!!
!! @sub tetra_lambin_weight
!!
!! Lambin-Vigneron algorithm for crystal Green's function. Here we have
!! implemented the equations presented in G. Palsson's Ph.D thesis.
!!
  subroutine tetra_lambin_weight(z, e, weights)
     implicit none

!! external arguments
     ! current complex frequency
     complex(dp), intent(in)  :: z

     ! energy eigenvalues for four corners at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integration weights
     complex(dp), intent(out) :: weights(4)

!! local variables
     ! corner degeneracy index
     integer :: flag

     ! corner reordering index
     integer :: isrt(4)

     ! local copy of energy values
     complex(dp) :: ecp(4)

     ! linear weight factors
     complex(dp) :: res(4)

!! [body

     ! make local copy of corner energies
     ecp = e

     ! reorder the corner energies
     call tetra_lv_creorder(ecp, isrt, flag)

     ! select the appropriate routine to compute the integral
     select case (flag)

         ! arbitrary eigenvalues
         case (1)
             call tetra_lv_ekarb  (z, ecp, res)

         ! two eigenvalues are identical, e1 = e2
         case (2)
             call tetra_lv_ek2idn (z, ecp, res)

         ! two pairs of energy are identical, e1 = e2 and e3 = e4
         case (3)
             call tetra_lv_ek22idn(z, ecp, res)

         ! three eigenvalues are identical, e1 = e2 = e3
         case (4)
             call tetra_lv_ek3idn (z, ecp, res)

         ! four (all) eigenvalues are identical
         case (5)
             call tetra_lv_ek4idn (z, ecp, res)

     end select

     ! put integration weights in correct order
     weights( isrt(1) ) = res(1)
     weights( isrt(2) ) = res(2)
     weights( isrt(3) ) = res(3)
     weights( isrt(4) ) = res(4)

!! body]

     return
  end subroutine tetra_lambin_weight

!!========================================================================
!!>>> Lambin-Vigneron algorithm, service layer                         <<<
!!========================================================================

!!
!! @sub tetra_lv
!!
!! We define the Lambin-Vigneron function as follows:
!!    lv(x,y) = z ( 1 - z * log[1 + 1/z] )  where z = x / y
!! or
!!    lv(x,y) = ( 1 - log[1+w] / w ) / w   where w = y / x.
!! for small z or w we use the series expansions:
!!    lv(x,y) = z ( 1 + z ( log(z) - z ( 1 - z ( 1/2 - z ( 1/3 - z ( 1/4 - z ...
!!    lv(x,y) = 1/2 - w ( 1/3 - w ( 1/4 - w ( 1/5 - ...
!!
  function tetra_lv(x, y) result(lv)
     implicit none

!! external arguments
     complex(dp), intent(in) :: x
     complex(dp), intent(in) :: y

!! local variables
     complex(dp) :: z
     complex(dp) :: w

     ! function type
     complex(dp) :: lv

!! local parameters
     complex(dp), parameter :: i2 = dcmplx(1.0_dp/2.0_dp)
     complex(dp), parameter :: i3 = dcmplx(1.0_dp/3.0_dp)
     complex(dp), parameter :: i4 = dcmplx(1.0_dp/4.0_dp)
     complex(dp), parameter :: i5 = dcmplx(1.0_dp/5.0_dp)
     complex(dp), parameter :: i6 = dcmplx(1.0_dp/6.0_dp)
     complex(dp), parameter :: i7 = dcmplx(1.0_dp/7.0_dp)
     complex(dp), parameter :: i8 = dcmplx(1.0_dp/8.0_dp)

!! [body

     if ( abs(y) > abs(x) ) then
         z = x / y
         lv = z * ( 1.0_dp - z * ( log(x+y) - log(x) ) )
     else
         w = y / x
         if ( abs(w) > 0.1_dp ) then
             lv = ( 1.0_dp - ( log(x+y) - log(x) ) / w ) / w
         else
             lv = i2 - w * ( i3 - w * ( i4 - w * ( i5 - w * ( i6 - w * ( i7 - w * i8 ) ) ) ) )
         endif ! back if ( abs(w) > 0.1_dp ) block
     endif ! back if ( abs(y) > abs(x) ) block

!! body]

     return
  end function tetra_lv

!!
!! @sub tetra_lv_ekarb
!!
!! evaluate integrated weights when eigenvalues at each corner
!! are explicit.
!!
  subroutine tetra_lv_ekarb(z, e, res)
     implicit none

!! external arguments
     ! current complex energy
     complex(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integrated weights
     complex(dp), intent(out) :: res(4)

!! local variables
     ! ze_i = z - e_i
     complex(dp) :: ze1, ze2, ze3, ze4

     ! e_ij = e_i - e_j
     complex(dp) :: e12, e13, e14
     complex(dp) :: e21, e23, e24
     complex(dp) :: e31, e32, e34
     complex(dp) :: e41, e42, e43

     ! p_ij = ze_i * tetra_lv(ze_i, e_ij)
     complex(dp) :: p12, p13, p14
     complex(dp) :: p21, p23, p24
     complex(dp) :: p31, p32, p34
     complex(dp) :: p41, p42, p43

!! [body

     ! build intermediate variables ze_{i}
     ze1 = z - e(1)
     ze2 = z - e(2)
     ze3 = z - e(3)
     ze4 = z - e(4)

     ! build intermediate variables e_{ij}
     e12 = e(1) - e(2)
     e13 = e(1) - e(3)
     e14 = e(1) - e(4)

     e21 = e(2) - e(1)
     e23 = e(2) - e(3)
     e24 = e(2) - e(4)

     e31 = e(3) - e(1)
     e32 = e(3) - e(2)
     e34 = e(3) - e(4)

     e41 = e(4) - e(1)
     e42 = e(4) - e(2)
     e43 = e(4) - e(3)

     ! build intermediate variables p_{ij}
     p12 = ze1 * tetra_lv( ze1, e12 )
     p13 = ze1 * tetra_lv( ze1, e13 )
     p14 = ze1 * tetra_lv( ze1, e14 )

     p21 = ze2 * tetra_lv( ze2, e21 )
     p23 = ze2 * tetra_lv( ze2, e23 )
     p24 = ze2 * tetra_lv( ze2, e24 )

     p31 = ze3 * tetra_lv( ze3, e31 )
     p32 = ze3 * tetra_lv( ze3, e32 )
     p34 = ze3 * tetra_lv( ze3, e34 )

     p41 = ze4 * tetra_lv( ze4, e41 )
     p42 = ze4 * tetra_lv( ze4, e42 )
     p43 = ze4 * tetra_lv( ze4, e43 )

     ! build res(1)
     res(1) = p21 / ( e32 * e42 ) + p31 / ( e23 * e43 ) + p41 / ( e24 * e34 )

     ! build res(2)
     res(2) = p12 / ( e31 * e41 ) + p32 / ( e13 * e43 ) + p42 / ( e14 * e34 )

     ! build res(3)
     res(3) = p13 / ( e21 * e41 ) + p23 / ( e12 * e42 ) + p43 / ( e14 * e24 )

     ! build res(4)
     res(4) = p14 / ( e21 * e31 ) + p24 / ( e12 * e32 ) + p34 / ( e13 * e23 )

!! body]

     return
  end subroutine tetra_lv_ekarb

!!
!! @sub tetra_lv_ek2idn
!!
!! evaluate integrated weights when two corner energies are identical.
!!
  subroutine tetra_lv_ek2idn(z, e, res)
     implicit none

!! external arguments
     ! current complex energy
     complex(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integrated weights
     complex(dp), intent(out) :: res(4)

!! local variables
     ! ze_i = z - e_i
     complex(dp) :: ze2, ze3, ze4

     ! e_ij = e_i - e_j
     complex(dp) :: e23, e24
     complex(dp) :: e32, e34
     complex(dp) :: e42, e43

     complex(dp) :: p2, p3, p4
     complex(dp) :: v2
     complex(dp) :: v31, v32, v41, v42

!! [body

     ! build intermediate variables ze_{i}
     ze2 = z - e(2)
     ze3 = z - e(3)
     ze4 = z - e(4)

     ! build intermediate variables e_{ij}
     e23 = e(2) - e(3)
     e24 = e(2) - e(4)

     e32 = e(3) - e(2)
     e34 = e(3) - e(4)

     e42 = e(4) - e(2)
     e43 = e(4) - e(3)

     ! build intermediate variables p_{i}
     p2 = ze2 / ( e32 * e42 )
     p3 = ze3 / ( e23 * e43 )
     p4 = ze4 / ( e24 * e34 )

     ! build intermediate variables v_{ij}
     v2  = ze2 / ( e23 * e24 )
     v31 = ze3 / ( e23 * e23 )
     v32 = ze3 / ( e23 * e24 )
     v41 = ze4 / ( e24 * e24 )
     v42 = ze4 / ( e23 * e24 )

     ! build res(2)
     res(2) = 0.5_dp * p2 + p3 * tetra_lv( ze3, e32 ) + p4 * tetra_lv( ze4, e42 )

     ! build res(3)
     res(3) = v2 - ( 2.0_dp * v32 + v41 ) * tetra_lv( ze2, e23 ) + v41 * tetra_lv( ze4, e43 )

     ! build res(4)
     res(4) = v2 - ( 2.0_dp * v42 + v31 ) * tetra_lv( ze2, e24 ) + v31 * tetra_lv( ze3, e34 )

     ! build res(1)
     res(1) = res(2)

!! body]

     return
  end subroutine tetra_lv_ek2idn

!!
!! @sub tetra_lv_ek22idn
!!
!! evaluate integrated weights when e1 = e2 and e3 = e4.
!!
  subroutine tetra_lv_ek22idn(z, e, res)
     implicit none

!! external arguments
     ! current complex energy
     complex(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integrated weights
     complex(dp), intent(out) :: res(4)

!! local variables
     ! ze_i = z - e_i
     complex(dp) :: ze2, ze3

     ! e_ij = e_i - e_j
     complex(dp) :: e23, e32

     complex(dp) :: p21, p22, p31, p32

     ! build intermediate variables ze_{i}
     ze2 = z - e(2)
     ze3 = z - e(3)

     ! build intermediate variables e_{ij}
     e23 = e(2) - e(3)
     e32 = e(3) - e(2)

     ! build intermediate variables p_{ij}
     p21 = ze2 / ( e23 * e23 )
     p22 = ze2 / ( e32 * e32 )
     p31 = ze3 / ( e23 * e23 )
     p32 = ze3 / ( e32 * e32 )

     ! build res(2)
     res(2) = p32 + 0.5_dp * p22 - 3.0_dp * p22 * tetra_lv( ze3, e32 )

     ! build res(3)
     res(3) = p21 + 0.5_dp * p31 - 3.0_dp * p31 * tetra_lv( ze2, e23 )

     ! build res(1)
     res(1) = res(2)

     ! build res(4)
     res(4) = res(3)

!! body]

     return
  end subroutine tetra_lv_ek22idn

!!
!! @sub tetra_lv_ek3idn
!!
!! evaluate integrated weights when e1 = e2 = e3.
!!
  subroutine tetra_lv_ek3idn(z, e, res)
     implicit none

!! external arguments
     ! current complex energy
     complex(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integrated weights
     complex(dp), intent(out) :: res(4)

!! local variables
     ! actually they are constants
     complex(dp) :: c1, c2

     ! ze_i = z - e_i
     complex(dp) :: ze3, ze4

     ! e_ij = e_i - e_j
     complex(dp) :: e34, e43

     complex(dp) :: p31, p32, p41, p42

!! [body

     ! build intermediate variables c1 and c2
     c1 = dcmplx(2.0_dp/6.0_dp)
     c2 = dcmplx(5.0_dp/6.0_dp)

     ! build intermediate variables ze_{i}
     ze3 = z - e(3)
     ze4 = z - e(4)

     ! build intermediate variables e_{ij}
     e34 = e(3) - e(4)
     e43 = e(4) - e(3)

     ! build intermediate variables p_{ij}
     p31 = ze3 / ( e34 * e34 )
     p32 = ze3 / ( e43 * e43 )
     p41 = ze4 / ( e34 * e34 )
     p42 = ze4 / ( e43 * e43 )

     ! build res(3)
     res(3) = c1 * p31 - c2 * p41 + p41 * tetra_lv( ze4, e43 )

     ! build res(4)
     res(4) = p42 + 0.5_dp * p32 - 3.0_dp * p32 * tetra_lv( ze4, e43 )

     ! build res(1)
     res(1) = res(3)

     ! build res(2)
     res(2) = res(3)

!! body]

     return
  end subroutine tetra_lv_ek3idn

!!
!! @sub tetra_lv_ek4idn
!!
!! evaluate integrated weights when eigenvalues at four corners are
!! the same.
!!
  subroutine tetra_lv_ek4idn(z,e,res)
     implicit none

!! external arguments
     ! current complex energy
     complex(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     complex(dp), intent(in)  :: e(4)

     ! integrated weights
     complex(dp), intent(out) :: res(4)

!! local variables
     ! loop index
     integer :: i

!! [body

     do i=1,4
         res(i) = 0.25_dp / ( z - e(i) )
     enddo ! over i={1,4} loop

!! body]

     return
  end subroutine tetra_lv_ek4idn

!!
!! @sub tetra_lv_creorder
!!
!! this subroutine reorders four complex energies, e(i), in such a way
!! that any identical energies are in the low index array slots. there are
!! a total of 15 different possibilities that can occur and to distinguish
!! between these we set up a unique array of prime numbers and then use the
!! product of the array elements to identify which case occurs.
!!
!! note that the way we check for identical elements could in principle
!! result in more than 15 possibilities if say element 1 and 2 are close
!! within tolerance and element 2 and 3 are close within tolerance it could
!! occur that element 1 and 3 would not be counted as equal if they were
!! not within the tolerance limit. this case is identified by a negative
!! value for the flag, flg, and should be treated (as error) in the calling
!! routine.
!!
  subroutine tetra_lv_creorder(e, isrt, flg)
     implicit none

!! external arguments
     ! flg == 1 if all energies different;
     ! flg == 2 if two energies are identical;
     ! flg == 3 if two pairs of energies are identical;
     ! flg == 4 if three energies are identical;
     ! flg == 5 if all energies are identical;
     ! flg ==-1 if no order is found.
     integer, intent(out) :: flg

     ! pointer array
     integer, intent(out) :: isrt(4)

     ! energy points (should be overwritten)
     complex(dp), intent(inout) :: e(4)

!! local variables
     ! loop index
     integer :: p
     integer :: q

     ! index pointer
     integer :: idx

     ! real and imaginary part of z
     real(dp) :: rez
     real(dp) :: imz

     ! z = e_p - e_q
     complex(dp) :: z

     ! dummy arrays for e
     complex(dp) :: enew(4)

!! [body

     ! set up prime number matrix
     idx = 1
     do p=1,3
         do q=p+1,4
             z = e(p) - e(q)
             rez = dabs( dreal(z) )
             imz = dabs( dimag(z) )
             !
             if ( rez <= eps8 .and. imz <= eps8 ) then
                 idx = idx * tetra_lv_prime( p, q )
             endif
         enddo ! over q={p+1,4} loop
     enddo ! over p={1,3} loop

     ! find flg
     flg = 1
     !
     if ( idx > 1    ) flg = flg + 1
     if ( idx > 13   ) flg = flg + 1
     if ( idx > 35   ) flg = flg + 1
     if ( idx > 1001 ) flg = flg + 1

     ! find new ordering
     if ( idx == 1 .or. idx == 2 .or. idx == 42 .or. idx == 30030 ) then
         call tetra_lv_setmap( isrt, 1, 2, 3, 4 )

     else if ( idx == 3 ) then
         call tetra_lv_setmap( isrt, 1, 3, 2, 4 )

     else if ( idx == 5 ) then
         call tetra_lv_setmap( isrt, 1, 4, 2, 3 )

     else if ( idx ==  7 .or. idx == 35 ) then
         call tetra_lv_setmap( isrt, 2, 3, 1, 4 )

     else if ( idx == 11 .or. idx == 33 ) then
         call tetra_lv_setmap( isrt, 2, 4, 1, 3 )

     else if ( idx == 13 .or. idx == 26 ) then
         call tetra_lv_setmap( isrt, 3, 4, 1, 2 )

     else if ( idx ==  195 ) then
         call tetra_lv_setmap( isrt, 1, 3, 4, 2 )

     else if ( idx ==  110 ) then
         call tetra_lv_setmap( isrt, 1, 2, 4, 3 )

     else if ( idx == 1001 ) then
         call tetra_lv_setmap( isrt, 2, 3, 4, 1 )

     else
         flg = -1
         call tetra_lv_setmap( isrt, 1, 2, 3, 4 )

     endif

     ! reorder energies
     do p=1,4
         enew(p) = e( isrt(p) )
     enddo ! over p={1,4} loop

     ! overwrite original energies
     do p=1,4
         e(p) = enew(p)
     enddo ! over p={1,4} loop

!! body]

     return
  end subroutine tetra_lv_creorder

!!
!! @sub tetra_lv_prime
!!
!! this function associates a prime number to pairs of small integers.
!! the association is the following:
!!    (1,2) = 2,
!!    (1,3) = 3,
!!    (1,4) = 5,
!!    (2,3) = 7,
!!    (2,4) = 11,
!!    (3,4) = 13.
!!
  function tetra_lv_prime(p, q) result(prime)
     implicit none

!! external arguments
     ! pair of integers
     integer, intent(in) :: p
     integer, intent(in) :: q

!! local variables
     ! function type, it means prime number
     integer :: prime

!! [body

     prime = p * q
     !
     if ( prime == 4 .or. prime == 6 .or. prime == 12 ) then
         prime = prime + 1
     else if ( prime == 8 ) then
         prime = prime + 3
     endif

!! body]

     return
  end function tetra_lv_prime

!!
!! @sub tetra_lv_setmap
!!
!! this subroutine takes in one vector, v, of integers and four numbers
!! a, b, c, d and assigns them the elements of the vector v, which is
!! of length four.
!!
  subroutine tetra_lv_setmap(v, a, b, c, d)
     implicit none

!! external arguments
     integer, intent(in) :: a
     integer, intent(in) :: b
     integer, intent(in) :: c
     integer, intent(in) :: d
     !
     integer, intent(out) :: v(4)

!! [body

     v(1) = a
     v(2) = b
     v(3) = c
     v(4) = d

!! body]

     return
  end subroutine tetra_lv_setmap

!!========================================================================
!!>>> Gaussian smearing algorithm                                      <<<
!!========================================================================

!!
!! @sub smearing_gauss_weight1
!!
!! standard Gaussian broadening algorithm for (integrated) density
!! of states.
!!
  subroutine smearing_gauss_weight1(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! energies at given k-point and band
     real(dp), intent(in)  :: e

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! dummy variables
     real(dp) :: dummy

!! [body

     dummy = ( z - e ) / gamm
     tos = 0.125_dp * ( 1.0_dp - erf(-dummy) )
     dos = 0.25_dp * exp(-dummy**2.0) / ( sqrt(pi) * gamm )

!! body]

     return
  end subroutine smearing_gauss_weight1

!!
!! @sub smearing_gauss_weight2
!!
!! standard Gaussian broadening algorithm for (integrated) density
!! of states.
!!
  subroutine smearing_gauss_weight2(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! loop index
     integer :: i

     ! dummy variables
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         tweight(i) = 0.125_dp * ( 1.0_dp - erf(-dummy) )
         dweight(i) = 0.25_dp * exp(-dummy**2.0) / ( sqrt(pi) * gamm )
     enddo ! over i={1,4} loop

     ! sum up weights to calculate the density of states
     dos = sum( dweight )

     ! sum up the weights to calculate the integrated density of states
     tos = sum( tweight )

!! body]

     return
  end subroutine smearing_gauss_weight2

!!
!! @sub smearing_gauss_weight3
!!
!! standard Gaussian broadening algorithm for (integrated) density
!! of states.
!!
  subroutine smearing_gauss_weight3(z, e, dd, tt)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! integration weight for density of states
     real(dp), intent(out) :: dd(4)

     ! integration weight for number of states
     real(dp), intent(out) :: tt(4)

!! local variables
     ! loop index
     integer  :: i

     ! dummy variables
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         tweight(i) = 0.125_dp * ( 1.0_dp - erf(-dummy) )
         dweight(i) = 0.25_dp * exp(-dummy**2.0) / ( sqrt(pi) * gamm )
     enddo ! over i={1,4} loop

     ! set up weights to calculate the density of states
     dd = dweight

     ! set up the weights to calculate the integrated density of states
     tt = tweight

!! body]

     return
  end subroutine smearing_gauss_weight3

!!========================================================================
!!>>> Fermi-Dirac smearing algorithm                                   <<<
!!========================================================================

!!
!! @sub smearing_fermi_weight1
!!
!! Fermi-Dirac broadening algorithm for (integrated) density of states.
!!
  subroutine smearing_fermi_weight1(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! energies at given k-point and band
     real(dp), intent(in)  :: e

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! dummy variables
     real(dp) :: dummy

!! [body

     dummy = ( z - e ) / gamm
     tos = 1.0_dp / (1.0_dp +  exp(-dummy) )
     dos = 0.5_dp / (1.0_dp + cosh( dummy) )

!! body]

     return
  end subroutine smearing_fermi_weight1

!!
!! @sub smearing_fermi_weight2
!!
!! Fermi-Dirac broadening algorithm for (integrated) density of states.
!!
  subroutine smearing_fermi_weight2(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! loop index
     integer :: i

     ! dummy variables
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         tweight(i) = 1.0_dp / ( 1.0_dp +  exp(-dummy) )
         dweight(i) = 0.5_dp / ( 1.0_dp + cosh( dummy) )
     enddo ! over i={1,4} loop

     ! sum up weights to calculate the density of states
     dos = sum( dweight )

     ! sum up the weights to calculate the integrated density of states
     tos = sum( tweight )

!! body]

     return
  end subroutine smearing_fermi_weight2

!!
!! @sub smearing_fermi_weight3
!!
!! Fermi-Dirac broadening algorithm for (integrated) density of states.
!!
  subroutine smearing_fermi_weight3(z, e, dd, tt)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! integration weight for density of states
     real(dp), intent(out) :: dd(4)

     ! integration weight for number of states
     real(dp), intent(out) :: tt(4)

!! local variables
     ! loop index
     integer :: i

     ! dummy variables
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         tweight(i) = 1.0_dp / (1.0_dp +  exp(-dummy) )
         dweight(i) = 0.5_dp / (1.0_dp + cosh( dummy) )
     enddo ! over i={1,4} loop

     ! set up weights to calculate the density of states
     dd = dweight

     ! set up the weights to calculate the integrated density of states
     tt = tweight

!! body]

     return
  end subroutine smearing_fermi_weight3

!!========================================================================
!!>>> cold smearing algorithm                                          <<<
!!========================================================================

!!
!! @sub smearing_marzari_weight1
!!
!! cold smearing algorithm, by Marzari and Vanderbilt, for (integrated)
!! density of states.
!!
  subroutine smearing_marzari_weight1(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! energies at given k-point and band
     real(dp), intent(in)  :: e

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! dummy variables
     real(dp) :: xp
     real(dp) :: dummy

!! [body

     dummy = ( z - e ) / gamm
     xp = dummy - 1.0_dp / sqrt(2.0_dp)
     tos = 0.5_dp * erf(xp) + 1.0_dp / sqrt(2.0_dp * pi) * exp(-xp**2.0) + 0.5_dp
     dos = 2.0_dp / sqrt(pi) * exp(-xp**2.0) * ( 2.0_dp - sqrt(2.0_dp) * dummy )

!! body]

     return
  end subroutine smearing_marzari_weight1

!!
!! @sub smearing_marzari_weight2
!!
!! cold smearing algorithm, by Marzari and Vanderbilt, for (integrated)
!! density of states.
!!
  subroutine smearing_marzari_weight2(z, e, dos, tos)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! density of states
     real(dp), intent(out) :: dos

     ! number of states
     real(dp), intent(out) :: tos

!! local variables
     ! loop index
     integer :: i

     ! dummy variables
     real(dp) :: xp
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         xp = dummy - 1.0_dp / sqrt(2.0_dp)
         tweight(i) = 0.5_dp * erf(xp) + 1.0_dp / sqrt(2.0_dp * pi) * exp(-xp**2.0) + 0.5_dp
         dweight(i) = 2.0_dp / sqrt(pi) * exp(-xp**2.0) * ( 2.0_dp - sqrt(2.0_dp) * dummy )
     enddo ! over i={1,4} loop

     ! sum up weights to calculate the density of states
     dos = sum( dweight )

     ! sum up the weights to calculate the integrated density of states
     tos = sum( tweight )

!! body]

     return
  end subroutine smearing_marzari_weight2

!!
!! @sub smearing_marzari_weight3
!!
!! cold smearing algorithm, by Marzari and Vanderbilt, for (integrated)
!! density of states.
!!
  subroutine smearing_marzari_weight3(z, e, dd, tt)
     implicit none

!! external arguments
     ! current energy
     real(dp), intent(in)  :: z

     ! corner energies at given tetrahedron
     real(dp), intent(in)  :: e(4)

     ! integration weight for density of states
     real(dp), intent(out) :: dd(4)

     ! integration weight for number of states
     real(dp), intent(out) :: tt(4)

!! local variables
     ! loop index
     integer :: i

     ! dummy variables
     real(dp) :: xp
     real(dp) :: dummy

!! [body

     do i=1,4
         dummy = ( z - e(i) ) / gamm
         xp = dummy - 1.0_dp / sqrt(2.0_dp)
         tweight(i) = 0.5_dp * erf(xp) + 1.0_dp / sqrt(2.0_dp * pi) * exp(-xp**2.0) + 0.5_dp
         dweight(i) = 2.0_dp / sqrt(pi) * exp(-xp**2.0) * ( 2.0_dp - sqrt(2.0_dp) * dummy )
     enddo ! over i={1,4} loop

     ! set up weights to calculate the density of states
     dd = dweight

     ! set up the weights to calculate the integrated density of states
     tt = tweight

!! body]

     return
  end subroutine smearing_marzari_weight3

  end module mtetra
