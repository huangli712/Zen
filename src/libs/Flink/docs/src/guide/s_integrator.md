It provides two subroutines to perform simple numerical integration.

## Type

subroutines

## Source

`src/s_integrator.f90`

## Usage

To use `s_int_trapezoid()` or `s_int_simpson()`, you have to define the integrand at first. For example:

```fortran
function f(x)
    use constants, only : dp

    implicit none

    real(dp) :: x
    real(dp) :: f

    f = x * x
end function f
```

Next, you have to determine the lower bound `a` and upper bound `b`, and the number of points `n`. Noted that now both the `s_int_trapezoid()` and `s_int_simpson()` functions only support the 1-D numerical integration.

**(1)** Use `s_int_trapezoid()`.

```fortran
procedure( real(dp) ) :: s_int_trapezoid
procedure( real(dp) ) :: f
real(dp) :: val

val = s_int_trapezoid(f, a, b, n)
```

**(2)** Use `s_int_simpson()`.

```fortran
procedure( real(dp) ) :: s_int_simpson
procedure( real(dp) ) :: f
real(dp) :: val

val = s_int_simpson(f, a, b, n)
```

## Theory

**Trapezoid rule**

```math
\int_{a}^{b} f(x) dx \approx
\frac{h}{2}
\left[
f(x_0) + 2f(x_1) + \cdots + 2f(x_{n-1}) + f(x_n)
\right],
```

where $h = (b-a)/n$.

**Simpson rule**

```math
\int_{a}^{b} f(x) dx \approx
\frac{h}{3}
\left[
f(x_0) + f(x_n) +
4 \sum_{i = 1}^{n/2} f(x_{2i-1}) +
2 \sum_{i = 1}^{n/2-1} f(x_{2i})
\right],
```

where $h = (b-a)/n$, and $n$ must be even number.
