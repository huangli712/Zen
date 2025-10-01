In this module, we implement some basic sparse matrix algebra. Now it supports both double precision real and complex numbers. The internal format of sparse matrix used in this module is CSR (compressed sparse row) format.

## Type

module

## Source

`src/m_sparse.f90`

## Usage

**(1)** Import sparse module support.

```fortran
use sparse
```

**(2)** Subroutines

* **csr_alloc**, allocate memory for sparse matrix.
* **csr_free**, free memory for sparse matrix.
* **csr_print**, print sparse matrix to terminal.
* **csr_dns**, convert sparse matrix to dense matrix.
* **dns_csr**, convert dense matrix to sparse matrix.
* **csr_csr**, copy sparse matrix to sparse matrix.
* **get_csr**, return the elements in sparse matrix.
* **csr_mv**, perform sparse matrix-vector multiplication.
* **csr_mm**, perform sparse matrix-matrix multiplication.
* **csr_md**, perform sparse matrix-diagonal matrix multiplication.
* **csr_dm**, perform diagonal matrix-sparse matrix multiplication.
* **csr_plus**, addition for two sparse matrices.

**(3)** Key parameters

* **nrows**, row dimension of dense matrix.
* **ncols**, column dimension of dense matrix.
* **nnz**, maximum number of nonzero elements allowed.

**(4)** Sparse matrix in CSR format

* Method 1

```fortran
! ia, ja, a, input matrix in compressed sparse row format
integer, intent(in) :: ia(nrows+1)
integer, intent(in) :: ja(nnz)
real(dp), intent(in) :: a(nnz)
```

or

```fortran
! ia, ja, a, input matrix in compressed sparse row format
integer, intent(in) :: ia(nrows+1)
integer, intent(in) :: ja(nnz)
complex(dp), intent(in) :: a(nnz)
```

* Method 2

```fortran
! sparse_t: abstract struct for sparse matrix
type, private, abstract :: sparse_t

    ! nrows: number of rows
    integer :: nrows = 0

    ! ncols: number of columns
    integer :: ncols = 0

    ! nnz: number of non-zero values
    integer :: nnz   = 0

end type sparse_t

! csr_t: sparse matrix in compressed sparse row format
type, private, extends(sparse_t) :: csr_t

    ! rowptr: matrix row pointer
    integer, allocatable :: rowptr(:)

    ! colptr: matrix column pointer
    integer, allocatable :: colptr(:)

end type csr_t

! csr_d: compressed sparse row format, real(dp) version
type, public, extends(csr_t) :: csr_d
    !
    real(dp), allocatable :: V(:)
    !
end type csr_d

! csr_z: compressed sparse row format, complex(dp) version
type, public, extends(csr_t) :: csr_z
    !
    complex(dp), allocatable :: V(:)
    !
end type csr_z
```

Note that only the `csr_d` and `csr_z` structs are public. The `csr_d` struct is for real sparse matrix, while `csr_z` struct is for complex sparse matrix.

**(5)** Sample code

```fortran
program test
    use constants, only : dp
    use sparse

    implicit none

    ! dimensional parameters
    integer, parameter :: nrows = 5
    integer, parameter :: ncols = 5
    integer, parameter :: nnz = 10

    ! dense matrix
    real(dp) :: AA(nrows,ncols)

    ! sparse matrix
    type (csr_d) :: TA

    ! initialize dense matrix
    AA = 0.0_dp
    AA(1,1) = 0.31_dp
    AA(1,4) = 4.07_dp
    AA(2,2) = 0.84_dp
    AA(2,5) = 2.60_dp
    AA(3,1) = 6.18_dp
    AA(3,2) = 7.12_dp
    AA(4,4) = 0.82_dp
    AA(5,3) = 1.12_dp

    ! allocate memory for sparse matrix
    call csr_alloc(nrows, ncols, nnz, TA)
    call csr_print(TA)

    ! convert dense matrix to sparse matrix
    call dns_csr(nrows, ncols, AA, TA)
    call csr_print(TA)

    ! deallocate memory for sparse matrix
    call csr_free(TA)
end program test
```
