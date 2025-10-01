The major features of the ACTest toolkit are as follows:

* ACTest can randomly generate any number of ``A(\omega)``, as well as corresponding ``G(\tau)`` or ``G(i\omega_n)``.

* ACTest employs one or more parameterized Gaussian, Lorentzian, ``\delta``-like, rectangular, and Rise-And-Decay peaks to assemble ``A(\omega)``. ``A(\omega)`` can be either positive definite for fermionic Green's functions, or non-positive definite for bosonic Green's functions and matrix-valued Green's functions. The frequency grid, i.e., ``\omega``, can be either linear or non-linear (such as tangent, Lorentzian, and half-Lorentzian grids).

* ACTest supports fermionic, bosonic, and symmetric bosonic kernels to generate Green's functions on either imaginary time or Matsubara frequency axes.

* To mimic realistic quantum Monte Carlo simulation data, ACTest supports artificial noise. The noise could be imaginary time correlated or a simple normal distribution. The noise level is adjustable.

* ACTest includes a built-in testing dataset, ACT100, which can serve as a relatively fair standard for examining different analytic continuation methods and codes.

* ACTest is already interfaced with [ACFlow](https://github.com/huangli712/ACFlow), which is a full-fledged and open-source analytic continuation toolkit. ACTest can access various analytic continuation methods in the ACFlow toolkit, launch them to perform analytic continuation calculations, and provide benchmark reports on their accuracy and efficiency.

* ACTest is interfaced with the [MiniPole](https://github.com/Green-Phys/MiniPole) code, which implements the minimal pole method.

* ACTest also provides a plotting script, that can be used to visualize and compare the true and reconstructed spectral functions.

* ACTest is an open-source software developed in Julia language. It is quite easy to be extended to implement new features or support the other analytic continuation codes, such as [Nevanlinna.jl](https://github.com/SpM-lab/Nevanlinna.jl) and [SmoQyDEAC.jl](https://github.com/SmoQySuite/SmoQyDEAC.jl). Furthermore, ACTest offers extensive documentation and examples, making it user-friendly.

In the following text, we will elaborate on the technical details inside ACTest.
