It is a common module which defines some common used numerical or physical constants. They are as follows:

* Numerical precision: `sp`, `dp`
* File units: `mystd`, `myout`, `mytmp`
* Numerical constants (real number): `pi`, `zero`, `one`, `two`, `half`
* Numerical constants (small number): `eps6`, `eps8`, `epst`, `epss`
* Numerical constants (complex number): `czi`, `cone`, `czero`
* Physical constants: `ev2k`, `ry2e`, `ha2e`

We always need to import this module.

## Type

module

## Source

`src/m_constants.f90`

## Usage

**(1)** Import constants module completely.

```fortran
program test
    use constants

    implicit none

    real(dp) :: A
    A = one
end program test
```

**(2)** Import constants module partially.

```fortran
program test
    use constants, only : dp
    use constants, only : one

    implicit none

    real(dp) :: A
    A = one
end program test
```
