In this module, two tetrahedron integration algorithms are implemented. They are:

* P. E. Blochl tetrahedron integration algorithm,
* Lambin-Vigneron tetrahedron integration algorithm.

Note that the former is used to calculate the DFT density of states, while the later is used to calculate the lattice Green's functions during the DFT + DMFT calculations. Besides, we also implement three smearing algorithm, which can be used to calculate the DFT density of states as well. They are:

* Gaussian smearing algorithm,
* Fermi-Dirac smearing algorithm,
* Marzari-Vanderbilt cold smearing algorithm.

This module depends on the `s_util.f90/s_qsorter()` subroutine.

## Type

module

## Source

`src/m_tetra.f90`

## Usage

**(1)** P. E. Blochl tetrahedron integration algorithm.

```fortran
call tetra_blochl_weight1(z, e, dos, tos)
```

or

```fortran
call tetra_blochl_weight2(z, e, ddd, ttt)
```

For the external energy `z` and the one-particle energies at the corners of the tetrahedron `e`, the subroutine `tetra_blochl_weight1()` returns total contributions of the given tetrahedron to the density of states `dos` and integrated density of states `tos`, the subroutine `tetra_blochl_weight2()` returns corner-dependent contributions of the given tetrahedron to the density of states `ddd` and integrated density of states `ttt`.

!!! warning

    The subroutines
    * `tetra_blochl_weight1()`
    * `tetra_blochl_weight2()`
    are private. So, please call the public interface `tetra_weight()`.

**(2)** Lambin-Vigneron tetrahedron integration algorithm.

```fortran
call tetra_lambin_weight(z, e, weights)
```

For the external energy `z` and the energy eigenvalues for four corners at given tetrahedron, the subroutine `tetra_lambin_weight()` returns corner-dependent integration weights `weights` for lattice Green's function. Note that now `z`, `e`, and `weights` are all complex.

!!! warning

    The subroutine `tetra_lambin_weight()` is private. So, please call the public interface `tetra_weight()`.

**(3)** Gaussian smearing algorithm.

```fortran
call smearing_gauss_weight1(z, e, dos, tos)
```

or

```fortran
call smearing_gauss_weight2(z, e, dos, tos)
```

or

```fortran
call smearing_gauss_weight3(z, e, dd, tt)
```

Here, `z` means external energy, `e` means energy at given k-point and band in `smearing_gauss_weight1()`, corner energies at given tetrahedron in `smearing_gauss_weight2()` and `smearing_gauss_weight3()`. In `smearing_gauss_weight1()` and `smearing_gauss_weight2()`, the total contributions of the given tetrahedron to the density of states `dos` and integrated density of states `tos` are returned. In `smearing_gauss_weight3()`, the corner-dependent contributions (integration weights) of the given tetrahedron to the density of states `dd` and integrated density of states `tt` are returned.

!!! warning

    The subroutines
    * `smearing_gauss_weight1()`
    * `smearing_gauss_weight2()`
    * `smearing_gauss_weight3()`
    are private. So, please call the public interface `smearing_gauss()`.

**(4)** Fermi-Dirac smearing algorithm.

```fortran
call smearing_fermi_weight1(z, e, dos, tos)
```

or

```fortran
call smearing_fermi_weight2(z, e, dos, tos)
```

or

```fortran
call smearing_fermi_weight3(z, e, dd, tt)
```

Here, `z` means external energy, `e` means energy at given k-point and band in `smearing_fermi_weight1()`, corner energies at given tetrahedron in `smearing_fermi_weight2()` and `smearing_fermi_weight3()`. In `smearing_fermi_weight1()` and `smearing_fermi_weight2()`, the total contributions of the given tetrahedron to the density of states `dos` and integrated density of states `tos` are returned. In `smearing_fermi_weight3()`, the corner-dependent contributions (integration weights) of the given tetrahedron to the density of states `dd` and integrated density of states `tt` are returned.

!!! warning

    The subroutines
    * `smearing_fermi_weight1()`
    * `smearing_fermi_weight2()`
    * `smearing_fermi_weight3()`
    are private. So, please call the public interface `smearing_fermi()`.

**(5)** Marzari-Vanderbilt cold smearing algorithm.

```fortran
call smearing_marzari_weight1(z, e, dos, tos)
```

or

```fortran
call smearing_marzari_weight2(z, e, dos, tos)
```

or

```fortran
call smearing_marzari_weight3(z, e, dd, tt)
```

Here, `z` means external energy, `e` means energy at given k-point and band in `smearing_marzari_weight1()`, corner energies at given tetrahedron in `smearing_marzari_weight2()` and `smearing_marzari_weight3()`. In `smearing_marzari_weight1()` and `smearing_marzari_weight2()`, the total contributions of the given tetrahedron to the density of states `dos` and integrated density of states `tos` are returned. In `smearing_marzari_weight3()`, the corner-dependent contributions (integration weights) of the given tetrahedron to the density of states `dd` and integrated density of states `tt` are returned.

!!! warning

    The subroutines
    * `smearing_marzari_weight1()`
    * `smearing_marzari_weight2()`
    * `smearing_marzari_weight3()`
    are private. So, please call the public interface `smearing_marzari()`.

## Theory

### Analytical Tetrahedron Integration (Blochl's Algorithm)

**Number of states**

Here we list the expressions for the integrated density of states or number of states function $n(\epsilon)$ from a given tetrahedron. The number of states function is used to determine the Fermi level. The expressions shown are similar to those given in previous papers. For the sake of simplicity we have omitted all band indices. The one-particle energies at the corners of the tetrahedron are $\epsilon_1$, $\epsilon_2$, $\epsilon_3$, $\epsilon_4$, which are ordered according to increasing values. $\epsilon_{ij}$ is a shorthand notation for $\epsilon_i - \epsilon_j$. $V_T$ is the reciprocal-space volume of the tetrahedron, and $V_G$ is the volume of the reciprocal unit cell.

```math
n(\epsilon) = 0
```

for $\epsilon < \epsilon_1$,

```math
n(\epsilon) = \frac{V_T}{V_G} \frac{(\epsilon - \epsilon_1)^3}{\epsilon_{21}\epsilon_{31}\epsilon_{41}}
```

for $\epsilon_1 < \epsilon < \epsilon_2$,

```math
n(\epsilon) = \frac{V_T}{V_G}
\frac{1}{\epsilon_{31}\epsilon_{41}}
\left[
\epsilon^2_{21} +
3\epsilon_{21}(\epsilon-\epsilon_2) +
3(\epsilon - \epsilon_2)^2 -
\frac{\epsilon_{31} + \epsilon_{42}}{\epsilon_{32}\epsilon_{42}}(\epsilon - \epsilon_2)^3
\right]
```

for $\epsilon_2 < \epsilon < \epsilon_3$,

```math
n(\epsilon) = \frac{V_T}{V_G}
\left(
1 - \frac{(\epsilon_4 - \epsilon)^3}{\epsilon_{41}\epsilon_{42}\epsilon_{43}}
\right)
```

for $\epsilon_3 < \epsilon < \epsilon_4$, and

```math
n(\epsilon) = \frac{V_T}{V_G}
```

for $\epsilon > \epsilon_4$.

**Integration weights**

Here we give the expressions that result in the integration weights $w_{nj}$. By $w_1$, $w_2$, $w_3$, and $w_4$ we denote the contribution to the integration weights at the four corners of a tetrahedron, which are again ordered according to increasing one-particle energies. Note that the band index $n$ is suppressed.

For a fully tetrahedron, i.e., $\epsilon_F < \epsilon_1$, the contributions vanish:

```math
w_1 = w_2 = w_3 = w_4 = 0.
```

For $\epsilon_1 < \epsilon_F < \epsilon_2$, we obtain

```math
w_1 = C
\left[
4 - (\epsilon_F - \epsilon_1)
\left(\frac{1}{\epsilon_{21}} + \frac{1}{\epsilon_{31}} + \frac{1}{\epsilon_{41}}\right)
\right],
```

```math
w_2 = C \frac{\epsilon_F - \epsilon_1}{\epsilon_{21}},
```

```math
w_3 = C \frac{\epsilon_F - \epsilon_1}{\epsilon_{31}},
```

```math
w_4 = C \frac{\epsilon_F - \epsilon_1}{\epsilon_{41}},
```

with

```math
C = \frac{V_T}{4V_G} \frac{(\epsilon_F - \epsilon_1)^3}{\epsilon_{21}\epsilon_{31}\epsilon_{41}}.
```

For $\epsilon_2 < \epsilon_F < \epsilon_3$, we obtain

```math
w_1 = C_1 +
(C_1 + C_2) \frac{\epsilon_3 - \epsilon_F}{\epsilon_{31}} +
(C_1 + C_2 + C_3)\frac{\epsilon_4 - \epsilon_F}{\epsilon_{41}},
```

```math
w_2 = C_1 + C_2 + C_3 +
(C_2 + C_3) \frac{\epsilon_3 - \epsilon_F}{\epsilon_{32}} +
C_3 \frac{\epsilon_4 - \epsilon_F}{\epsilon_{42}},
```

```math
w_3 = (C_1 + C_2) \frac{\epsilon_F - \epsilon_1}{\epsilon_{31}} +
(C_2 + C_3) \frac{\epsilon_F - \epsilon_2}{\epsilon_{32}},
```

```math
w_4 = (C_1 + C_2 + C_3)\frac{\epsilon_F - \epsilon_1}{\epsilon_{41}} +
C_3 \frac{\epsilon_F - \epsilon_2}{\epsilon_{42}},
```

with

```math
C_1 = \frac{V_T}{4V_G} \frac{(\epsilon_F - \epsilon_1)^2}{\epsilon_{41}\epsilon_{31}},
```

```math
C_2 = \frac{V_T}{4V_G}
\frac{(\epsilon_F - \epsilon_1)(\epsilon_F - \epsilon_2)(\epsilon_3 - \epsilon_F)}{\epsilon_{41}\epsilon_{32}\epsilon_{31}},
```

```math
C_3 = \frac{V_T}{4V_G}
\frac{(\epsilon_F - \epsilon_2)^2(\epsilon_4 - \epsilon_F)}{\epsilon_{42}\epsilon_{32}\epsilon_{41}}.
```

For $\epsilon_3 < \epsilon_F < \epsilon_4$, the weights are

```math
w_1 = \frac{V_T}{4V_G} - C \frac{\epsilon_4 - \epsilon_F}{\epsilon_{41}},
```

```math
w_2 = \frac{V_T}{4V_G} - C \frac{\epsilon_4 - \epsilon_F}{\epsilon_{42}},
```

```math
w_3 = \frac{V_T}{4V_G} - C \frac{\epsilon_4 - \epsilon_F}{\epsilon_{43}},
```

```math
w_4 = \frac{V_T}{4V_G} - C
\left[
4 - \left(\frac{1}{\epsilon_{41}} + \frac{1}{\epsilon_{42}} + \frac{1}{\epsilon_{43}}\right)
(\epsilon_4 - \epsilon_F)
\right],
```

with

```math
C = \frac{V_T}{4V_G} \frac{(\epsilon_4 - \epsilon_F)^3}{\epsilon_{41}\epsilon_{42}\epsilon_{43}}.
```

For a fully occupied tetrahedron the contribution for each corner is identical:

```math
w_1 = w_2 = w_3 = w_4 = \frac{V_T}{4V_G}.
```

The corresponding correction $dw_i$ to the integration weights has an extremely simple form:

```math
dw_i = \sum_T \frac{1}{40} D_T (\epsilon_F) \sum^4_{j=1} (\epsilon_j - \epsilon_i).
```

**Density of states**

The contribution of one tetrahedron to the density of states at energy $\epsilon$ is given by

```math
D_T(\epsilon) = 0
```

for $\epsilon < \epsilon_1$ and for $\epsilon_4 < \epsilon$,

```math
D_T(\epsilon) = \frac{V_T}{V_G}
\frac{3(\epsilon-\epsilon_1)^2}{\epsilon_{21}\epsilon_{31}\epsilon_{41}}
```

for $\epsilon_1 < \epsilon < \epsilon_2$,

```math
D_T(\epsilon) = \frac{V_T}{V_G}
\frac{1}{\epsilon_{31}\epsilon_{41}}
\left[
3\epsilon_{21} + 6(\epsilon - \epsilon_2)
- 3\frac{(\epsilon_{31} + \epsilon_{42})(\epsilon - \epsilon_2)^2}{\epsilon_{32}\epsilon_{42}}
\right]
```

for $\epsilon_2 < \epsilon < \epsilon_3$, and

```math
D_T(\epsilon) = \frac{V_T}{V_G}
\frac{3(\epsilon_4 - \epsilon)^2}{\epsilon_{41}\epsilon_{42}\epsilon_{43}}
```

for $\epsilon_3 < \epsilon < \epsilon_4$.

**References:**

[1] P. E. Blochl, O. Jepsen, and O. K. Anderson, *Phys. Rev. B* **49**, 16223 (1994)

### Analytical Tetrahedron Integration (Lambin-Vigneron's Algorithm)

In this part we will write up the equations of Lambin and Vigneron for the weight factors used in the analytical tetrahedron method. The purpose of this is to rewrite their formulas in a simpler, more readable and hopefully more
stable form.

Before we start writing out the equations let us define a few symbols. The external frequency is denoted by $z$ and the energy at corner number $i$ is denoted by $\epsilon_i$.  We define $z_i = z-\epsilon_i$ and also $\epsilon_{ij} = \epsilon_i-\epsilon_j$.

In the general case when all the tetrahedron corners have different energies the equation for the weight factor at corner $i$ is according to LV

```math
\begin{align*}
r_i
&\equiv
6\int_{0}^{1}dc\int_{0}^{1-c}db\int_{0}^{1-b-c}da
\frac{(1-a-b-c)\delta_{i1}+a\delta_{i2}+b\delta_{i3}+c\delta_{i4}}
{z_1-a\epsilon_{21}-b\epsilon_{31}-c\epsilon_{41}} \\
&=
\frac{z_i^2}{\prod_{k\neq i}\epsilon_{ki}}+\sum_{j\neq i}
\frac{z_j^3}{\prod_{k\neq j}\epsilon_{kj}}
\frac{\log\left(\frac{z_j}{z_i}\right)}{\epsilon_{ij}}.
\end{align*}
```

If we now use the identity:

```math
0 = \sum_{i=1}^{4}\frac{z_i^2}{\prod_{j\neq i}(z_i-z_j)} =
  \sum_{i=1}^{4}\frac{z_i^2}{\prod_{j\neq i}\epsilon_{ji}}.
```

the equation for the weight factor can be written as

```math
\begin{align*}
r_i
&=
\sum_{j \neq i}
\frac{z_j^2}{\prod_{k \neq j}\epsilon_{kj}}
\left[
\frac{z_j}{\epsilon_{ij}}\log\left(\frac{z_j}{z_i}\right)-1
\right] \\
&=
\sum_{j \neq i}
\frac{z_j}{\prod_{k \neq i,j}\epsilon_{kj}}
\Big[
\frac{z_j}{\epsilon_{ji}}-
\frac{z_j^2}{\epsilon_{ij}^2}\log\left(1+\frac{\epsilon_{ji}}{z_j}\right)
\Big] \\
&=
\sum_{j \neq i}
\frac{z_j}{\prod_{k \neq i,j}\epsilon_{kj}}lv\left(\frac{z_j}{\epsilon_{ji}}\right).
\end{align*}
```

Here we have defined the so called Lambin-Vigneron function as:

```math
lv(x) = x\left(1-x\log\left[1+\frac{1}{x}\right]\right).
```

This is a simple function which is easy to compute and also has the added benefit of having nice asymptotic limits, i.e.

```math
\lim_{x\rightarrow 0}lv(x) = 0,
```

and

```math
\lim_{x\rightarrow \infty}lv(x) = \frac{1}{2}.
```

For considering the cases where the some of the corner energies are identical it is useful to write out the equations above explicitly for each corner:

```math
\begin{align*}
r_1 &=
 \frac{z_2}{\epsilon_{32}\epsilon_{42}}lv\left(\frac{z_2}{\epsilon_{21}}\right)
+\frac{z_3}{\epsilon_{23}\epsilon_{43}}lv\left(\frac{z_3}{\epsilon_{31}}\right)
+\frac{z_4}{\epsilon_{24}\epsilon_{34}}lv\left(\frac{z_4}{\epsilon_{41}}\right),
 \\
r_2 &=
 \frac{z_1}{\epsilon_{31}\epsilon_{41}}lv\left(\frac{z_1}{\epsilon_{12}}\right)
+\frac{z_3}{\epsilon_{13}\epsilon_{43}}lv\left(\frac{z_3}{\epsilon_{32}}\right)
+\frac{z_4}{\epsilon_{14}\epsilon_{34}}lv\left(\frac{z_4}{\epsilon_{42}}\right),
 \\
r_3 &=
 \frac{z_1}{\epsilon_{21}\epsilon_{41}}lv\left(\frac{z_1}{\epsilon_{13}}\right)
+\frac{z_2}{\epsilon_{12}\epsilon_{42}}lv\left(\frac{z_2}{\epsilon_{23}}\right)
+\frac{z_4}{\epsilon_{14}\epsilon_{24}}lv\left(\frac{z_4}{\epsilon_{43}}\right),
 \\
r_4 &=
 \frac{z_1}{\epsilon_{21}\epsilon_{31}}lv\left(\frac{z_1}{\epsilon_{14}}\right)
+\frac{z_2}{\epsilon_{12}\epsilon_{32}}lv\left(\frac{z_2}{\epsilon_{24}}\right)
+\frac{z_3}{\epsilon_{13}\epsilon_{23}}lv\left(\frac{z_3}{\epsilon_{34}}\right).
\end{align*}
```

The next case we consider is when two corners have equal energies and we assume here that the corners have been ordered in such as way that $\epsilon_1 = \epsilon_2$.  Looking at the formulas above we see that the limit $\epsilon_1 \rightarrow \epsilon_2$ is trivial take for the weight factors for corners 1 and 2 but for corners 3 and 4 two terms must be taken together and using Mathematica we get:

```math
\begin{align*}
r_2 &=
 \frac{1}{2}\frac{z_2}{\epsilon_{32}\epsilon_{42}}
+\frac{z_3}{\epsilon_{23}\epsilon_{43}}lv\left(\frac{z_3}{\epsilon_{32}}\right)
+\frac{z_4}{\epsilon_{24}\epsilon_{34}}lv\left(\frac{z_4}{\epsilon_{42}}\right),
 \\
r_3 &=
 \frac{z_2}{\epsilon_{23}\epsilon_{24}}
 -\left[\frac{2z_3}{\epsilon_{23}\epsilon_{24}}+\frac{z_4}{\epsilon_{24}^2}\right]
 lv\left(\frac{z_2}{\epsilon_{23}}\right)
+\frac{z_4}{\epsilon_{24}^2}lv\left(\frac{z_4}{\epsilon_{43}}\right),
 \\
r_4 &=
 \frac{z_2}{\epsilon_{23}\epsilon_{24}}
 -\left[\frac{2z_4}{\epsilon_{23}\epsilon_{24}}+\frac{z_3}{\epsilon_{23}^2}\right]
 lv\left(\frac{z_2}{\epsilon_{24}}\right)
+\frac{z_3}{\epsilon_{23}^2}lv\left(\frac{z_3}{\epsilon_{34}}\right).
\end{align*}
```

The next case we consider is $\epsilon_1 = \epsilon_2$ and $\epsilon_4 = \epsilon_3$ and we see from the formulas in equation (see above) that the limit $\epsilon_4 \rightarrow \epsilon_3$ is trivial for corners 3 and 4 and the weight factor for corner 2 follows from symmetry. Hence we get:

```math
\begin{align*}
r_2 &=
 \frac{z_3}{\epsilon_{32}^2}+\frac{1}{2}\frac{z_2}{\epsilon_{32}^2}
 -3\frac{z_2}{\epsilon_{32}^2}lv\left(\frac{z_3}{\epsilon_{32}}\right),
 \\
r_3 &=
 \frac{z_2}{\epsilon_{23}^2}+\frac{1}{2}\frac{z_3}{\epsilon_{23}^2}
 -3\frac{z_3}{\epsilon_{23}^2}lv\left(\frac{z_2}{\epsilon_{23}}\right).
\end{align*}
```

Next case is $\epsilon_1 = \epsilon_2 = \epsilon_3$ and this case we can either obtain from the formulas in equation (44) or simply do the corresponding integral which is simple in this case since we have:

```math
\begin{align*}
  r_3 &= \int_{0}^{1}dc\int_{0}^{1-c}db\int_{0}^{1-b-c}da
  \frac{b}{z_3-c\epsilon_{43}},
  \\
  r_4 &= \int_{0}^{1}dc\int_{0}^{1-c}db\int_{0}^{1-b-c}da
  \frac{c}{z_3-c\epsilon_{43}}.
\end{align*}
```

Performing these integrals we obtain:

```math
\begin{align*}
r_3 &=
 \frac{2z_3-5z_4}{6\epsilon_{34}^2}+
 \frac{z_4}{\epsilon_{34}^2}lv\left(\frac{z_4}{\epsilon_{43}}\right), \\
r_4 &=
 \frac{z_4}{\epsilon_{43}^2}+\frac{1}{2}\frac{z_3}{\epsilon_{43}^2}
 -3\frac{z_3}{\epsilon_{43}^2}lv\left(\frac{z_4}{\epsilon_{43}}\right).
\end{align*}
```

The final case is of course the simplest one where all the corners have the same energy and in that case we obtain:

```math
r_4 = \frac{1}{4}\frac{1}{z_4}.
```

**References:**

[1] Gunnar Palsson, *PhD Thesis*, Rutgers University (2001)

[2] P. Lambin and J. P. Vigneron, *Phys. Rev. B* **29**, 3430 (1984)

### Smearning Algorithm

**Gaussian smearing**

```math
n(x) = \frac{1 - \text{erfc}{(-x)}} {8},
```

```math
D(x) = \frac{1}{4} \frac{\exp(-x^2/2)}{\sqrt{\pi}\gamma},
```

where $x = (z - e) / \gamma$. $\gamma$ is a adjustable parameter, $z$ is the external energy, $e$ is the energy at given $k$-point and band or corner energies at given tetrahedron.

**Fermi-Dirac smearing**

```math
n(x) = \frac{1}{1 + \exp{(-x)}},
```

```math
D(x) = \frac{1}{2} \frac{1}{1 + \cosh{(x)}},
```

where $x = (z - e) / \gamma$. $\gamma$ is a adjustable parameter, $z$ is the external energy, $e$ is the energy at given $k$-point and band or corner energies at given tetrahedron.

**Marzari-Vanderbilt cold smearing**

```math
n(x) = \frac{\text{erf}(\tilde{x})}{2} + \frac{1}{\sqrt{2\pi}}
\exp(-\tilde{x}^2) + \frac{1}{2},
```

```math
D(x) = \frac{2}{\sqrt{\pi}}
\exp{(-\tilde{x}^2)}( 2 - \sqrt{2}x ).
```

where $x = (z - e) / \gamma$ and $\tilde{x} = x - \frac{1}{\sqrt{2}}$. $\gamma$ is a adjustable parameter, $z$ is the external energy, $e$ is the energy at given $k$-point and band or corner energies at given tetrahedron.

**References:**

[1] N. Marzari, D. Vanderbilt, *et al.*, *Phys. Rev. Lett.* **82**, 3296 (1999)
