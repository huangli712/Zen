It provides a lot of subroutines to manupulate the vectors.

## Type

subroutines

## Source

`src/s_vector.f90`

## Usage

**(1)** Mesh generation

```fortran
subroutine s_linspace_i(xmin, xmax, n, x)
subroutine s_linspace_d(xmin, xmax, n, x)
subroutine s_linspace_z(xmin, xmax, n, x)
```

To create a linear mesh `x` in interval [`xmin`, `xmax`], `n` is size of array `x`.

**(2)** Calculate sum of vector.

```fortran
subroutine s_cumsum_i(n, v, vsum)
subroutine s_cumsum_d(n, v, vsum)
subroutine s_cumsum_z(n, v, vsum)
```

Return the cumsum of an array `v`. `n` is size of array `v`, and `vsum` is cumsum of array `v`.

**(3)** Calculate product of vector.

```fortran
subroutine s_cumprod_i(n, v, vprod)
subroutine s_cumprod_d(n, v, vprod)
subroutine s_cumprod_z(n, v, vprod)
```

Return the cumproduct of an array `v`. `n` is size of array `v`, and `vprod` is cumprod of array `v`.

**(4)** Swap two vectors.

```fortran
subroutine s_swap_i(n, ix, iy)
subroutine s_swap_d(n, dx, dy)
subroutine s_swap_z(n, zx, zy)
```

Exchange two vectors, where `n` is dimension of the vectors.

**(5)** Linear mixing for vectors.

```fortran
subroutine s_mix_i(n, ix, iy, alpha)
subroutine s_mix_d(n, dx, dy, alpha)
subroutine s_mix_z(n, zx, zy, alpha)
```

Perform linear mixing for two vectors, where `n` is dimension of the vectors, and `alpha` is the mixing parameter $\alpha$.

**(6)** Convert diagonal elements of matrix to vector.

```fortran
subroutine s_vecadd_i(n, ix, iy, alpha)
subroutine s_vecadd_d(n, dx, dy, alpha)
subroutine s_vecadd_z(n, zx, zy, alpha)
```

Try to add diagonal elements of a matrix to a vector, where `n` is dimension of the vectors, and `alpha` is the prefactor. `*x` is vector, and `*y` is a `n`-by-`n` matrix.

!!! note

    `_i` means integer version, `_d` real(dp) version, and `_z` complex(dp) version.
