The following two random number generators (generates a random number between 0 and 1, real double precision) are supported by now.

* MT19937
```text
Mersenne Twister pseudorandom number generator
by Makoto Matsumoto and Takuji Nishimura
Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator
ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, January pp.3-30 (1998).
```

* SFMT
```text
SIMD-oriented Fast Mersenne Twister
by Mutsuo Saito and Makoto Matsumoto
SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number Generator
Monte Carlo and Quasi-Monte Carlo Methods 2006, Springer, 2008, pp. 607-622.
```

## Type

module

## Source

`src/m_spring.f90`

## Usage

**(1)** Use MT19937.

```fortran
use spring
call spring_mt_init(seed)
r = spring_mt_stream()
r = spring_mt_string()
```

**(2)** Use SFMT.

```fortran
use spring
call spring_sfmt_init(seed)
r = spring_sfmt_stream()
r = spring_sfmt_string()
```

The functions `spring_mt_stream()` and `sprint_sfmt_stream()` obtain a psuedo random real number in the range (0,1), i.e., a number
greater than 0 and less than 1. The functions `spring_mt_string()` and `sprint_sfmt_string()` obtain a psuedo random real number in the range [0,1], i.e., a number greater than or equal to 0 and less than or equal to 1.

!!! note

    Since SFMT has a better performance, it is preferable.
