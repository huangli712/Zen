It provides some subroutines to calculate the special functions.

## Type

subroutines

## Source

`src/s_function.f90`

## Usage

**(1)** Orthogonal polynomial basis.

```fortran
subroutine s_leg_basis(lemax, legrd, lmesh, rep_l)
subroutine s_che_basis(chmax, chgrd, cmesh, rep_c)
subroutine s_svd_basis(svmax, svgrd, smesh, rep_s, bose, beta)
```

`lemax` means maximum order for legendre orthogonal polynomial, `legrd` means number of mesh points for legendre orthogonal polynomial, `lmesh` means mesh for legendre orthogonal polynomial in [-1,1], `rep_l` saves legendre orthogonal polynomial defined on [-1,1].

`chmax` means maximum order for chebyshev orthogonal polynomial, `chgrd` means number of mesh points for chebyshev orthogonal polynomial, `cmesh` means mesh for chebyshev orthogonal polynomial in [-1,1], `rep_c` saves chebyshev orthogonal polynomial defined on [-1,1].

`svmax` means maximum order for svd orthogonal polynomial, `svgrd` means number of mesh points for svd orthogonal polynomial, `smesh` means mesh for svd orthogonal polynomial in [-1,1], `rep_s` saves svd orthogonal polynomial defined on [-1,1], `bose` determines whether the bosonic kernel is used, `beta` means inverse system temperature $\beta$.

**(2)** Spheric Bessel function.

```fortran
subroutine s_sph_jl(lmax, x, jl)
```

It computes the spherical Bessel functions of the first kind, $j_l(x)$, for argument $x$ and $l = 0, 1, \ldots, l_{max}$.

**(3)** Bernstein polynomial.

```fortran
subroutine s_bezier(n, x, bern)
```

`n` means the degree of the bernstein polynomials to be used. For any given $n$, there is a set of $n + 1$ bernstein polynomials, each of degree $n$, which form a basis for polynomials on [0,1]. `x` means the evaluation point. `bern` saves the values of the $n+1$ bernstein polynomials at $x$.

**(4)** Some helper functions for `s_svd_basis()`.

```fortran
function s_safe_exp(x)
```

It is a safe exp call to avoid data overflow.

```fortran
function s_f_kernel(tau, omega, beta)
function s_b_kernel(tau, omega, beta)
```

They are used to calculate fermionic or bosonic kernel functions. `tau` means $\tau$, `omega` means $\omega$, `beta` means inverse system temperature $\beta$.

## Theory

**Legendre orthogonal polynomial**

The Legendre orthogonal polynomials obey the three term recurrence relation, known as Bonnet’s recursion formula:

```math
P_0(x) = 1,
```

```math
P_1(x) = x,
```

```math
(n+1) P_{n+1}(x) = (2n+1) P_n(x) - n P_{n-1}(x).
```

**Chebyshev orthogonal polynomial**

The Chebyshev orthogonal polynomials of the second kind can be defined by the following recurrence relation:

```math
U_0(x) = 1,
```

```math
U_1(x) = 2x,
```

```math
U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x).
```

**SVD orthogonal polynomial**

In the imaginary-frequency domain, the Lehmann representation reads,

```math
G(i\nu) = \int_{-\infty}^\infty d\omega
\underbrace{\frac{1}{i\nu - \omega}}_{\equiv K(i\nu, \omega)}
A(\omega),
```

where $A(\omega)$ is a spectral function. $K(i\nu,\omega)$ is the so-called analytic continuation kernel. The Lehmann representation can be transformed to the imaginary-time domain as

```math
G(\tau) = -\int_{-\infty}^\infty
d\omega K(\tau,\omega) A(\omega),
```

where $0 < \tau < \beta$ and

```math
K(\tau,\omega) \equiv
-\frac{1}{\beta} \sum_{i\nu} e^{-i\nu \tau} K(i\nu,\omega)=
\begin{cases}
    \frac{e^{-\tau\omega}}{1+e^{-\beta\omega}} & (\mathrm{fermion}),\\
    \frac{e^{-\tau\omega}}{1-e^{-\beta\omega}} & (\mathrm{boson}).
\end{cases}
```

The minus sign originates from our convention $K(\tau, \omega) > 0$. To avoid the divergence of the bosonic kernel at $\omega=0$, we reformulate Equation of $G(\tau)$ as

```math
G(\tau)= -\int_{-\infty}^\infty d{\omega}
K^\mathrm{L}(\tau,\omega) \rho(\omega),
```

where $K^\mathrm{L}(\tau,\omega)$ is the *logistic kernel* defined as

```math
K^\mathrm{L}(\tau,\omega) =
\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}},
```

and $\rho(\omega)$ is the modified spectral function

```math
\rho(\omega) \equiv
\begin{cases}
    A(\omega) & (\mathrm{fermion}),\\
    \frac{A(\omega)}{\tanh(\beta \omega/2)} & (\mathrm{boson}).
\end{cases}
```

The singular value expnasion of the kernel reads

```math
K^\mathrm{L}(\tau, \omega) = \sum_{l=0}^\infty U_l(\tau) S_l V_l(\omega),
```

for $\omega \in [-\omega_{max}, \omega_{max}]$ with $\omega_{max}$ ($> 0$) being a cut-off frequency. $U_l(\tau)$ and $V_l(\omega)$ are left and right singular functions and $S_l$ is the singular values (with $S_0>S_1>S_2>...>0$). The two sets of singular functions $U$ and $V$ make up the basis functions of the so-called Intermediate Representation (IR), which depends on $\beta$ and the cutoff $\omega_{max}$. For the peculiar choice of the regularization for the bosonic kernel using $K^\mathrm{L}$, these basis functions do not depend on statistical properties. The basis functions $U_l(\tau)$ are transformed to the imaginary-frequency axis as

```math
U_l(i\nu) \equiv \int_0^\beta d \tau e^{i\nu\tau} U_l(\tau).
```

Some of the information regarding real-frequency properties of the system is often lost during transition into the imaginary-time domain, so that the imaginary-frequency Green's function does hold less information than the real-frequency Green's function. The reason for using IR lies within its compactness and ability to display that information in imaginary quantities.

The decay of the singular values depends on $\beta$ and $\omega_{max}$ only through the dimensionless parameter $\Lambda \equiv \beta\omega_{max}$.

**Spheric Bessel function**

The following recursion relation

```math
j_{l+1}(x)=\frac{2l+1}{x}j_l(x)-j_{l-1}(x),
```

is used either downwards for $x < l$ or upwards for $x \ge l$. For $x \ll 1$, the following asymtotic form is used:

```math
j_l(x) \approx \frac{x^l}{(2l+1)!!}.
```

This procedure is numerically stable and accurate to near this machine precision for $l \le 50$.

**Bernstein polynomial**

The bernstein polynomials are assumed to be based on [0,1]. The formula reads:

```math
B(N,I)(X) = \frac{N!}{I!(N-I)!} (1-X)^{(N-I)} X^I.
```

First values:

    B(0,0)(X) = 1
    B(1,0)(X) =      1-X
    B(1,1)(X) =                X
    B(2,0)(X) =     (1-X)**2
    B(2,1)(X) = 2 * (1-X)    * X
    B(2,2)(X) =                X**2
    B(3,0)(X) =     (1-X)**3
    B(3,1)(X) = 3 * (1-X)**2 * X
    B(3,2)(X) = 3 * (1-X)    * X**2
    B(3,3)(X) =                X**3
    B(4,0)(X) =     (1-X)**4
    B(4,1)(X) = 4 * (1-X)**3 * X
    B(4,2)(X) = 6 * (1-X)**2 * X**2
    B(4,3)(X) = 4 * (1-X)    * X**3
    B(4,4)(X) =                X**4

Special values:

* B(N,I)(X) has a unique maximum value at X = I/N.
* B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
* For X = 1/2,

```math
B(N,I)(1/2) = \frac{C(N,K)}{2^N}.
```

* For a fixed X and N, the polynomials add up to 1:

```math
\sum_{I=0}^{N} B(N,I)(X) = 1.
```
