It provides a lot of subroutines to manipulate the matrix and implement the linear algebra.

## Type

subroutines

## Source

`src/s_matrix.f90`

## Usage

**(1)** Build constants (0) matrix.

```fortran
subroutine s_zeros_i(n, A)
subroutine s_zeros_d(n, A)
subroutine s_zeros_z(n, A)
```

To build a matrix with all elements are zero. Here `n` means size of matrix, `A` means input/output matrix.

**(2)** Build constants (1) matrix.

```fortran
subroutine s_ones_i(n, A)
subroutine s_ones_d(n, A)
subroutine s_ones_z(n, A)
```

To build a matrix with all elements are zero. Here `n` means size of matrix, `A` means input/output matrix.

**(3)** Build constants (any values) matrix.

```fortran
subroutine s_any_i(n, i, A)
subroutine s_any_d(n, d, A)
subroutine s_any_z(n, z, A)
```

To build a matrix with all elements are given by `i`, `d`, or `z`. Here `n` means size of matrix, `A` means input/output matrix. `i`, `d`, and `z` are initial values of matrix elements.

**(4)** Build diagonal matrix.

```fortran
subroutine s_eye_i(n, k, A)
subroutine s_eye_d(n, k, A)
subroutine s_eye_z(n, k, A)
```

To build a matrix with ones on the diagonal and zeros elsewhere. Here `n` means size of matrix, `A` means input/output matrix, `k` means  index of the diagonal.

**(5)** Build identity matrix.

```fortran
subroutine s_identity_i(n, A)
subroutine s_identity_d(n, A)
subroutine s_identity_z(n, A)
```

To build an identity matrix. Here `n` means size of matrix, `A` means input/output matrix.

**(6)** Build diagonal matrix from vector.

```fortran
subroutine s_diag_i(n, v, A)
subroutine s_diag_d(n, v, A)
subroutine s_diag_z(n, v, A)
```

To build a diagonal matrix from a vector. Here `n` means size of matrix, `A` means input/output matrix, `v` means input vector.

**(7)** Calculate trace for matrix.

```fortran
subroutine s_trace_d(n, A, tr)
subroutine s_trace_z(n, A, tr)
```

Return trace for a matrix. Here `n` means size of matrix, `A` means input/output matrix, `tr` means matrix's trace.

**(8)** Calculate determinant for matrix.

```fortran
subroutine s_det_d(ndim, dmat, ddet)
subroutine s_det_z(ndim, zmat, zdet)
```

where `ndim` means dimension of the matrix, `ddet` and `zdet` are the determinants, `dmat` and `zmat` are the input matrices.

**(9)** Calculate matrix inversion.

```fortran
subroutine s_inv_d(ndim, dmat)
subroutine s_inv_z(ndim, zmat)
```

where `ndim` means dimension of the matrix, `dmat` and `zmat` are the input/output matrices.

**(10)** General eigensystem problem.

```fortran
subroutine s_eig_dg(ldim, ndim, amat, eval, evec)
subroutine s_eig_zg(ldim, ndim, zmat, zeig, zvec)
subroutine s_eigvals_dg(ldim, ndim, amat, eval)
subroutine s_eigvals_zg(ldim, ndim, zmat, zeig)
```

Here `ldim` is leading dimension of matrix, `ndim` is the order of the matrix, `amat` and `zmat` are the input matrices, `eval` and `zeig` mean eigenvalues, `evec` and `zvec` mean eigenvectors.

**(11)** Symmetric eigensystem problem.

```fortran
subroutine s_eig_sy(ldim, ndim, amat, eval, evec)
subroutine s_eig_he(ldim, ndim, amat, eval, evec)
subroutine s_eigvals_sy(ldim, ndim, amat, eval)
subroutine s_eigvals_he(ldim, ndim, amat, eval)
```

Here `ldim` is leading dimension of matrix, `ndim` is the order of the matrix, `amat` is the input matrix, `eval` means eigenvalues, `evec` means eigenvectors.

**(12)** Linear equation solver.

```fortran
subroutine s_solve_dg(n, nrhs, A, B)
subroutine s_solve_zg(n, nrhs, A, B)
subroutine s_solve_sy(n, nrhs, A, B)
subroutine s_solve_he(n, nrhs, A, B)
```

Try to solve a linear system $AX = B$. Here, `n` is the number of linear equations, `nrhs` is the number of right-hand sides, `A` and `B` are input matrices. On exit, `B` is overwritten by the solution matrix $X$.

**(13)** General singular value decomposition.

```fortran
subroutine s_svd_dg(m, n, min_mn, amat, umat, svec, vmat)
subroutine s_svd_zg(m, n, min_mn, amat, umat, svec, vmat)
```

To perform the singular values decomposition for a general `m`-by-`n` matrix $A$, where $A = U \Sigma V'$, return its left vectors, right vectors, and singular values. Here, `m` is number of rows of $A$ matrix, `n` is number of columns of $A$ matrix, `min_mn` is minimal value of `m` and `n`, `amat` means $A$, `umat` means  left vectors $U$, `svec` means singular values $\Sigma$, `vmat` means right vectors $V'$.

!!! note

    `_i` means integer version, `_d` real(dp) version, and `_z` complex(dp) version. `_dg` means real(dp) general version, `_zg` complex(dp) general version, `_sy` real(dp) symmetric version, `_he` complex(dp) Hermitian version.
