These subroutines are used to perform spline interpolation.

## Type

subroutines

## Source

`src/s_spline.f90`

## Usage

**(1)** Calculate 1-order derivates for a given function.

```fortran
subroutine s_spl_deriv1(ydim, xval, yval, d1y)
```

```fortran
     ! dimension of xval and yval
     integer, intent(in)   :: ydim

     ! old knots
     real(dp), intent(in)  :: xval(ydim)

     ! old function values to be interpolated
     real(dp), intent(in)  :: yval(ydim)

     ! 1-order derivates
     real(dp), intent(out) :: d1y(ydim)
```

**(2)** Calculate 2-order derivates for a given function.

```fortran
subroutine s_spl_deriv2(ydim, xval, yval, startu, startd, d2y)
```

```fortran
     ! dimension of xval and yval
     integer, intent(in)   :: ydim

     ! first-derivate at point 1
     real(dp), intent(in)  :: startu

     ! first-derivate at point ydim
     real(dp), intent(in)  :: startd

     ! old knots
     real(dp), intent(in)  :: xval(ydim)

     ! old function values to be interpolated
     real(dp), intent(in)  :: yval(ydim)

     ! 2-order derivates
     real(dp), intent(out) :: d2y(ydim)
```

**(3)** Evaluate function value at a given point.

```fortran
function s_spl_funct(xdim, xval, yval, d2y, x)
```

```fortran
     ! dimension of xval and yval
     integer, intent(in)  :: xdim

     ! new mesh point
     real(dp), intent(in) :: x

     ! old mesh
     real(dp), intent(in) :: xval(xdim)

     ! old function value
     real(dp), intent(in) :: yval(xdim)

     ! 2-order derviates of old function
     real(dp), intent(in) :: d2y(xdim)
```
