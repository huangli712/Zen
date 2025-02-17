#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
#

"""
    ACFlow

ACFlow is a modern software toolkit for solving the many-body analytic
continuation problem. It is usually used to convert the single-particle
or two-particle correlators from imaginary axis to real axis. Now this
toolkit is under heavy development. **PLEASE USE IT AT YOUR OWN RISK**.

Nowadays the ACFlow toolkit supports the following algorithms:

* Maximum Entropy Method (`MaxEnt` solver, `recommended`)
* Barycentric Rational Function Method (`BarRat` solver, `recommended`)
* Nevanlinna Analytical Continuation (`NevanAC` solver, `experimental`)
* Stochastic Analytic Continuation (`StochAC` and `StochSK` solvers)
* Stochastic Optimization Method (`StochOM` solver)
* Stochastic Pole eXpansion (`StochPX` solver, `recommended`)

More algorithms will be implemented in the future.

Note that ACFlow toolkit has been designed to be integrated into the
`Zen` package. Actually, it is also compatible with the `iQIST` toolkit.
In the other words, the end user can use it to analytically continue the
imaginary time (or Matsubara frequency) data generated by the various
quantum impurity solvers in the `iQIST` toolkit. Of course, it is quite
easy to implement some kinds of interfaces for the other quantum impurity
solvers.

For more details about how to obtain, install and use the ACFlow toolkit,
please visit the following website:

* `https://huangli712.github.io/projects/acflow/index.html`

Any suggestions, comments, and feedbacks are welcome. Enjoy it!
"""
module ACFlow

#=
### *Using Standard Libraries*
=#

using Distributed
using LinearAlgebra
using Statistics
using Random
using Dates
using Printf
using DelimitedFiles
using TOML

#=
### *Using Third-Party Libraries*
=#

#=
The `Zygote` package is used to to calculate gradient via an automatic
differentiation approach. The `NevanAC` solver depends on it to perform
the Hardy basis function optimization. However, once you have trouble in
installing `Zygote`, we also provide an built-in gradient function,
which implements a finite difference approach (it is quite slow). Please
see the comments in `nac.jl` about how to activate the built-in gradient
function.
=#

#using Zygote

#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the ACFlow toolkit.

*Members* :

```text
I32, I64, API -> Numerical types (Integer).
F32, F64, APF -> Numerical types (Float).
C32, C64, APC -> Numerical types (Complex).
R32, R64, APR -> Numerical types (Union of Integer and Float).
N32, N64, APN -> Numerical types (Union of Integer, Float, and Complex).
#
__LIBNAME__   -> Name of this julia toolkit.
__VERSION__   -> Version of this julia toolkit.
__RELEASE__   -> Released date of this julia toolkit.
__AUTHORS__   -> Authors of this julia toolkit.
#
authors       -> Print the authors of ACFlow to screen.
```
=#

#
include("global.jl")
#
export I32, I64, API
export F32, F64, APF
export C32, C64, APC
export R32, R64, APR
export N32, N64, APN
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors

#=
### *Includes And Exports* : *types.jl*
=#

#=
*Summary* :

Define some dicts and structs, which are used to store the config
parameters or represent some essential data structures.

*Members* :

```text
DType           -> Customized type.
ADT             -> Customized type.
#
PBASE           -> Configuration dict for general setup.
PMaxEnt         -> Configuration dict for MaxEnt solver.
PBarRat         -> Configuration dict for BarRat solver.
PNevanAC        -> Configuration dict for NevanAC solver.
PStochAC        -> Configuration dict for StochAC solver.
PStochSK        -> Configuration dict for StochSK solver.
PStochOM        -> Configuration dict for StochOM solver.
PStochPX        -> Configuration dict for StochPX solver.
#
AbstractSolver  -> Abstract AC solver.
MaxEntSolver    -> It represents the MaxEnt solver.
BarRatSolver    -> It represents the BarRat solver.
NevanACSolver   -> It represents the NevanAC solver.
StochACSolver   -> It represents the StochAC solver.
StochSKSolver   -> It represents the StochSK solver.
StochOMSolver   -> It represents the StochOM solver.
StochPXSolver   -> It represents the StochPX solver.
#
AbstractData    -> Abstract input data in imaginary axis.
RawData         -> Raw input data.
GreenData       -> Preprocessed input data.
#
AbstractGrid    -> Abstract mesh for input data.
FermionicImaginaryTimeGrid -> Grid in fermionic imaginary time axis.
FermionicFragmentTimeGrid -> Grid in fermionic imaginary time axis (incomplete).
FermionicMatsubaraGrid -> Grid in fermionic Matsubara frequency axis.
FermionicFragmentMatsubaraGrid -> Grid in fermionic Matsubara frequency axis (incomplete).
BosonicImaginaryTimeGrid -> Grid in bosonic imaginary time axis.
BosonicFragmentTimeGrid -> Grid in bosonic imaginary time axis (incomplete).
BosonicMatsubaraGrid -> Grid in bosonic Matsubara frequency axis.
BosonicFragmentMatsubaraGrid -> Grid in bosonic Matsubara frequency axis (incomplete).
#
AbstractMesh    -> Abstract grid for calculated spectral function.
LinearMesh      -> Linear mesh.
TangentMesh     -> Tangent mesh.
LorentzMesh     -> Lorentzian mesh.
HalfLorentzMesh -> Lorentzian mesh at half-positive axis.
DynamicMesh     -> Dynamic (very fine) mesh for stochastic-like solvers.
#
AbstractMC      -> Abstract Monte Carlo engine.
StochACMC       -> Monte Carlo engine used in the StochAC solver.
StochSKMC       -> Monte Carlo engine used in the StochSK solver.
StochOMMC       -> Monte Carlo engine used in the StochOM solver.
StochPXMC       -> Monte Carlo engine used in the StochPX solver.
```
=#

#
include("types.jl")
#
export DType
export ADT
#
export PBASE
export PMaxEnt
export PBarRat
export PNevanAC
export PStochAC
export PStochSK
export PStochOM
export PStochPX
#
export AbstractSolver
export MaxEntSolver
export BarRatSolver
export NevanACSolver
export StochACSolver
export StochSKSolver
export StochOMSolver
export StochPXSolver
#
export AbstractData
export RawData
export GreenData
#
export AbstractGrid
export FermionicImaginaryTimeGrid
export FermionicFragmentTimeGrid
export FermionicMatsubaraGrid
export FermionicFragmentMatsubaraGrid
export BosonicImaginaryTimeGrid
export BosonicFragmentTimeGrid
export BosonicMatsubaraGrid
export BosonicFragmentMatsubaraGrid
#
export AbstractMesh
export LinearMesh
export TangentMesh
export LorentzMesh
export HalfLorentzMesh
export DynamicMesh
#
export AbstractMC
export StochACMC
export StochSKMC
export StochOMMC
export StochPXMC

#=
### *Includes And Exports* : *util.jl*
=#

#=
*Summary* :

To provide some useful utility macros and functions. They can be used
to colorize the output strings, query the environments, and parse the
input strings, etc.

*Members* :

```text
@cswitch      -> C-style switch.
@time_call    -> Evaluate a function call and print the elapsed time.
@pcs          -> Print colorful strings.
#
require       -> Check julia envirnoment.
setup_args    -> Setup ARGS manually.
query_args    -> Query program's arguments.
trace_error   -> Write exceptions or errors to terminal or external file.
catch_error   -> Catch the thrown exceptions or errors.
welcome       -> Print welcome message.
overview      -> Print runtime information of ACFlow.
goodbye       -> Say goodbye.
sorry         -> Say sorry.
prompt        -> Print some messages or logs to the output devices.
line_to_array -> Convert a line to a string array.
```
=#

#
include("util.jl")
#
export @cswitch
export @time_call
export @pcs
#
export require
export setup_args
export query_args
export trace_error
export catch_error
export welcome
export overview
export goodbye
export sorry
export prompt
export line_to_array

#=
### *Includes And Exports* : *math.jl*
=#

#=
*Summary* :

To provide some numerical algorithms, such as numerical integrations,
interpolations, Einstein summation notation, and optimization method.

*Members* :

```text
secant    -> Root finding secant algorithm.
newton    -> Root finding newton algorithm.
#
trapz     -> Numerical integration (composite trapezoidal rule).
simpson   -> Numerical integration (simpson rule).
#
second_derivative -> Calculate second derivative.
gradient_via_fd ->  Calculate gradient via finite difference algorithm.
#
AbstractInterpolation -> Abstract interpolation
LinearInterpolation -> Linear interpolation.
QuadraticInterpolation -> Quadratic interpolation.
CubicSplineInterpolation -> Cubic spline interpolation.
#
@einsum   -> Macro for Einstein summation notation.
#
curve_fit -> Try to fit the given (x,y) data to a predefined model.
#
optimize  -> Try to minimize f(x) via the BFGS algorithm.
```
=#

#
include("math.jl")
#
export secant
export newton
#
export trapz
export simpson
#
export second_derivative
export gradient_via_fd
#
export AbstractInterpolation
export LinearInterpolation
export QuadraticInterpolation
export CubicSplineInterpolation
#
export @einsum
#
export curve_fit
#
export optimize

#=
### *Includes And Exports* : *grid.jl*
=#

#=
*Summary* :

To implement various grids for the input data.

*Members* :

```text
rebuild! -> Rebuild the grid.
resize!  -> Change size of the grid.
reverse! -> Reverse the grid.
```
=#

#
include("grid.jl")
#
export rebuild!
export resize!
export reverse!

#=
### *Includes And Exports* : *mesh.jl*
=#

#=
*Summary* :

To implement various meshes for the calculated spectral functions.

*Members* :

```text
nearest -> Return index of the nearest point to a given number.
```
=#

#
include("mesh.jl")
#
export nearest

#=
### *Includes And Exports* : *config.jl*
=#

#=
*Summary* :

To extract, parse, verify, and print the configuration parameters.
They are stored in external files (case.toml) or dictionaries.

*Members* :

```text
inp_toml   -> Parse case.toml, return raw configuration information.
fil_dict   -> Fill dicts for configuration parameters.
see_dict   -> Display all the relevant configuration parameters.
rev_dict_b -> Update dict (PBASE) for configuration parameters.
rev_dict_m -> Update dict (PMaxEnt) for configuration parameters.
rev_dict_r -> Update dict (PBarRat) for configuration parameters.
rev_dict_n -> Update dict (PNevanAC) for configuration parameters.
rev_dict_a -> Update dict (PStochAC) for configuration parameters.
rev_dict_k -> Update dict (PStochSK) for configuration parameters.
rev_dict_s -> Update dict (PStochOM) for configuration parameters.
rev_dict_x -> Update dict (PStochPX) for configuration parameters.
chk_dict   -> Check dicts for configuration parameters.
_v         -> Verify dict's values.
get_b      -> Extract value from dict (PBASE dict), return raw value.
get_m      -> Extract value from dict (PMaxEnt dict), return raw value.
get_r      -> Extract value from dict (PBarRat dict), return raw value.
get_n      -> Extract value from dict (PNevanAC dict), return raw value.
get_a      -> Extract value from dict (PStochAC dict), return raw value.
get_k      -> Extract value from dict (PStochSK dict), return raw value.
get_s      -> Extract value from dict (PStochOM dict), return raw value.
get_x      -> Extract value from dict (PStochPX dict), return raw value.
```
=#

#
include("config.jl")
#
export inp_toml
export fil_dict
export see_dict
export rev_dict_b
export rev_dict_m
export rev_dict_r
export rev_dict_n
export rev_dict_a
export rev_dict_k
export rev_dict_s
export rev_dict_x
export chk_dict
export _v
export get_b
export get_m
export get_r
export get_n
export get_a
export get_k
export get_s
export get_x

#=
### *Includes And Exports* : *inout.jl*
=#

#=
*Summary* :

To read the input data or write the calculated results.

*Members* :

```text
read_real_data    -> Read data in imaginary time axis.
read_cmplx_data   -> Read data in Matsubara frequency axis.
#
write_spectrum    -> Write spectral functions.
write_backward    -> Write reproduced input data in imaginary axis.
write_complete    -> Write full data in real axis.
write_misfit      -> Write α-dependent χ².
write_goodness    -> Write Θ-dependent χ².
write_model       -> Write default model function.
write_prony       -> Write Prony approximation.
write_barycentric -> Write barycentric rational function approximation.
write_hamiltonian -> Write effective hamiltonian for StochAC solver.
write_passed      -> Write indices of selected solutions StochOM and StochPX solvers.
write_pole        -> Write details of generated poles for StochPX solver.
write_probability -> Write Bayesian a-posteriori probability.
write_statistics  -> Write statistics info. for StochAC/StochSK/StochOM/StochPX solver.
```
=#

#
include("inout.jl")
#
export read_real_data
export read_cmplx_data
#
export write_spectrum
export write_backward
export write_complete
export write_misfit
export write_goodness
export write_model
export write_prony
export write_barycentric
export write_hamiltonian
export write_passed
export write_pole
export write_probability
export write_statistics

#=
### *Includes And Exports* : *model.jl*
=#

#=
*Summary* :

To define some default model functions.

*Members* :

```text
build_flat_model         -> Construct a flat model.
build_gaussian_model     -> Construct a gaussian model.
build_1gaussian_model    -> Construct a shifted gaussian model.
build_2gaussians_model   -> Construct a two-gaussians model.
build_lorentzian_model   -> Construct a lorentzian model.
build_1lorentzian_model  -> Construct a shifted lorentzian model.
build_2lorentzians_model -> Construct a two-lorentzians model.
build_risedecay_model    -> Construct a rise-and-decay model.
build_file_model         -> Construct a model from external file.
build_func_model         -> Construct a model by user-defined function.
```
=#

#
include("model.jl")
#
export build_flat_model
export build_gaussian_model
export build_1gaussian_model
export build_2gaussians_model
export build_lorentzian_model
export build_1lorentzian_model
export build_2lorentzians_model
export build_risedecay_model
export build_file_model
export build_func_model

#=
### *Includes And Exports* : *kernel.jl*
=#

#=
*Summary* :

To define various kernel functions.

*Members* :

```text
build_kernel        -> Build kernel function.
build_kernel_symm   -> Build kernel function for symmetric case.
#
make_blur           -> Add preblur effect to the spectral functions.
make_singular_space -> Perform singular value decomposition for kernel.
make_gauss_peaks    -> Generate a series gaussian peaks in a linear mesh.
```
=#

#
include("kernel.jl")
#
export build_kernel
export build_kernel_symm
#
export make_blur
export make_singular_space
export make_gauss_peaks

#=
### *Includes And Exports* : *maxent.jl*
=#

#=
*Summary* :

To implement the MaxEnt solver for analytic continuation problem.

*Members* :

```text
MaxEntContext   -> Essential struct for the MaxEnt solver.
#
solve           -> Wrapper function for the MaxEnt solver.
init            -> Initialize maximum entropy simulation.
run             -> Perform maximum entropy simulation.
last            -> Postprocess the calculated results and write them.
#
historic        -> historic algorithm.
classic         -> Classic algorithm.
bryan           -> Bryan algorithm.
chi2kink        -> Chi2kink algorithm.
optimizer       -> Optimize the non-linear equation.
#
precompute      -> Precompute some key coefficients.
f_and_J         -> Define the function that need to be optimized.
f_and_J_od      -> Define the function that need to be optimized (offdiag version).
svd_to_real     -> From singular to real space.
svd_to_real_od  -> From singular to real space (offdiag version).
calc_entropy    -> Calculate entropy.
calc_entropy_od -> Calculate entropy (offdiag version).
calc_bayes      -> Calculate Bayesian probability.
calc_bayes_od   -> Calculate Bayesian probability (offdiag version).
calc_chi2       -> Calculate χ² function.
```
=#

#
include("maxent.jl")
#
export MaxEntContext
#
export solve
export init
export run
export last
#
export historic
export classic
export bryan
export chi2kink
export optimizer
#
export precompute
export f_and_J
export f_and_J_od
export svd_to_real
export svd_to_real_od
export calc_entropy
export calc_entropy_od
export calc_bayes
export calc_bayes_od
export calc_chi2

#=
### *Includes And Exports* : *rfa.jl*
=#

#=
*Summary* :

To implement the BarRat solver for analytic continuation problem.

*Members* :

```text
BarycentricFunction -> Barycentric representation of a rational function.
bc_nodes        -> Return nodes of the rational function.
bc_values       -> Return values of the rational function.
bc_weights      -> Return weights of the rational function.
bc_degree       -> Return degree of the rational function.
bc_poles        -> Return poles of the rational function.
aaa             -> Adaptive Antoulas-Anderson algorithm.
#
PronyApproximation -> Prony approximation of a complex-valued function.
prony_data      -> Preprocess the input Matsubara data.
prony_svd       -> Perform singular value decomposition.
prony_idx       -> Choose an optimal vector in the orthogonal matrix V.
prony_v         -> Extract a vector from the orthogonal matrix V.
prony_gamma     -> Evaluate Γₚ for Prony approximation.
prony_omega     -> Evaluate Ωₚ for Prony approximation.
#
BarRatContext   -> Essential struct for the BarRat solver.
#
solve           -> Wrapper function for the BarRat solver.
init            -> Initialize barycentric rational function simulation.
run             -> Perform barycentric rational function simulation.
last            -> Postprocess the calculated results and write them.
poles!          -> Get pole representation for the Matsubara Green's function.
```
=#

#
include("rfa.jl")
#
export BarycentricFunction
export bc_nodes
export bc_values
export bc_weights
export bc_degree
export bc_poles
export aaa
#
export PronyApproximation
export prony_data
export prony_svd
export prony_idx
export prony_v
export prony_gamma
export prony_omega
#
export BarRatContext
#
export solve
export init
export run
export last
export poles!

#=
### *Includes And Exports* : *nac.jl*
=#

#=
*Summary* :

To implement the NevanAC solver for analytic continuation problem.

*Members* :

```text
NevanACContext  -> Essential struct for the NevanAC solver.
#
solve           -> Wrapper function for the NevanAC solver.
init            -> Initialize Nevanlinna analytical continuation simulation.
run             -> Perform Nevanlinna analytical continuation simulation.
last            -> Postprocess the calculated results and write them.
#
precompute      -> Precompute some key arrays.
calc_mobius     -> Perform Mobius transformation.
calc_inv_mobius -> Perform inverse Mobius transformation.
calc_pick       -> Calculate the Pick matrix.
calc_phis       -> Try to calculate the Φ vector.
calc_abcd       -> Try to calculate the coefficients matrix abcd.
calc_hbasis     -> Try to calculate the Hardy basis.
calc_hmatrix    -> Try to calculate ``[f^k(z), f^k(z)^*]`` for 0 ≤ 𝑘 ≤ 𝐻-1.
calc_theta      -> Try to calculate the contractive function θ(z).
calc_green      -> Evaluate the Green's function via Nevanlinna interpolant.
calc_noptim     -> Evaluate the optimal value for the size of input data.
calc_hmin!      -> Evaluate the minimum value for the order of Hardy basis.
calc_hopt!      -> Evaluate the optimal value for the order of Hardy basis.
hardy_optimize! -> Try to perform Hardy basis optimization.
smooth_norm     -> Establish the smooth norm.
check_pick      -> Check whether the input data are valid.
check_causality -> Check causality of the Hardy coefficients `𝑎𝑏`.
```
=#

#
include("nac.jl")
#
export NevanACContext
#
export solve
export init
export run
export last
#
export precompute
export calc_mobius
export calc_inv_mobius
export calc_pick
export calc_phis
export calc_abcd
export calc_hbasis
export calc_hmatrix
export calc_theta
export calc_green
export calc_noptim
export calc_hmin!
export calc_hopt!
export hardy_optimize!
export smooth_norm
export check_pick
export check_causality

#=
### *Includes And Exports* : *sac.jl*
=#

#=
*Summary* :

To implement the StochAC solver for analytic continuation problem.

*Members* :

```text
StochACElement -> A struct that contains Monte Carlo field configurations.
StochACContext -> Essential struct for the StochAC solver.
#
solve          -> Wrapper function for the StochAC solver.
init           -> Initialize stochastic analytic continuation simulation.
run (prun)     -> Perform stochastic analytic continuation simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
warmup         -> Warmup Monte Carlo engine.
sample         -> Sample field configurations via metropolis algorithm.
measure        -> Measure spectral functions and internal energies.
#
init_iodata    -> Preprocess the input data.
init_mc        -> Create a StochACMC struct.
init_element   -> Create a StochACElement struct.
init_context   -> Create a StochACContext struct.
#
calc_fmesh     -> Build very dense mesh in [wmin,wmax].
calc_phi       -> Calculate ϕ function.
calc_delta     -> Precompute δ functions.
calc_hamil     -> Calculate α-resolved Hc.
calc_htau      -> Calculate α-resolved h(τ).
calc_alpha     -> Calculate α parameters.
constraints    -> Limit the position of δ functions.
#
try_move_s     -> Try to shift the position of single δ function.
try_move_p     -> Try to shift the positions of two δ functions.
try_move_a     -> Try to change the weights of two δ functions.
try_move_x     -> Try to exchange configurations between two adjacent layers.
```
=#

#
include("sac.jl")
#
export StochACElement
export StochACContext
#
export solve
export init
export run
export prun
export average
export last
#
export warmup
export sample
export measure
#
export init_iodata
export init_mc
export init_element
export init_context
#
export calc_fmesh
export calc_phi
export calc_delta
export calc_hamil
export calc_htau
export calc_alpha
export constraints
#
export try_move_s
export try_move_p
export try_move_a
export try_move_x

#=
### *Includes And Exports* : *san.jl*
=#

#=
*Summary* :

To implement the StochSK solver for analytic continuation problem.

*Members* :

```text
StochSKElement -> A struct that contains Monte Carlo field configurations.
StochSKContext -> Essential struct for the StochSK solver.
#
solve          -> Wrapper function for the StochSK solver.
init           -> Initialize stochastic analytic continuation simulation.
run (prun)     -> Perform stochastic analytic continuation simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
warmup         -> Warmup Monte Carlo engine.
sample         -> Sample field configurations via metropolis algorithm.
measure        -> Measure spectral functions.
shuffle        -> Shuffle field configurations.
#
init_iodata    -> Preprocess the input data.
init_mc        -> Create a StochSKMC struct.
init_element   -> Create a StochSKElement struct.
init_context   -> Create a StochSKContext struct.
#
calc_fmesh     -> Build very dense mesh in [wmin,wmax].
calc_correlator-> Calculate correlator function from field configuration.
calc_goodness  -> Calculate χ² function.
calc_theta     -> Search optimal Θ.
constraints    -> Limit the position of δ functions.
#
try_move_s     -> Try to shift the position of one δ function.
try_move_p     -> Try to shift the positions of two δ functions.
try_move_q     -> Try to shift the positions of four δ functions.
```
=#

#
include("san.jl")
#
export StochSKElement
export StochSKContext
#
export solve
export init
export run
export prun
export average
export last
#
export warmup
export sample
export measure
export shuffle
#
export init_iodata
export init_mc
export init_element
export init_context
#
export calc_fmesh
export calc_correlator
export calc_goodness
export calc_theta
export constraints
#
export try_move_s
export try_move_p
export try_move_q

#=
### *Includes And Exports* : *som.jl*
=#

#=
*Summary* :

To implement the StochOM solver for analytic continuation problem.

*Members* :

```text
Box            -> A struct for describing the field configuration.
StochOMElement -> A struct that contains Monte Carlo field configurations.
StochOMContext -> Essential struct for the StochOM solver.
#
solve          -> Wrapper function for the StochOM solver.
init           -> Initialize stochastic optimization simulation.
run (prun)     -> Perform stochastic optimization simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
update         -> Sample field configurations via metropolis algorithm.
#
init_iodata    -> Preprocess the input data.
init_mc        -> Create a StochOMMC struct.
init_element   -> Create a StochOMElement struct.
init_context   -> Prepare data for a StochOMContext struct.
#
eval_lambda    -> Build Λ function.
calc_error     -> Calculate χ² function.
calc_green     -> Reproduce Green's function via the field configurations.
calc_norm      -> Calculate norm of the field configurations.
constraints    -> Limit the position of δ functions.
#
try_insert     -> Try to insert a new box in the configuration.
try_remove     -> Try to remove a box.
try_shift      -> Try to shift a box.
try_width      -> Try to change width of a box.
try_height     -> Try to change height of a box.
try_split      -> Try to split a box.
try_merge      -> Try to merge two boxes.
#
Pdx            -> Try to calculate the probability density function.
```
=#

#
include("som.jl")
#
export Box
export StochOMElement
export StochOMContext
#
export solve
export init
export run
export prun
export average
export last
#
export update
#
export init_iodata
export init_mc
export init_element
export init_context
#
export eval_lambda
export calc_error
export calc_green
export calc_norm
export constraints
#
export try_insert
export try_remove
export try_shift
export try_width
export try_height
export try_split
export try_merge
#
export Pdx

#=
### *Includes And Exports* : *spx.jl*
=#

#=
*Summary* :

To implement the StochPX solver for analytic continuation problem.

*Members* :

```text
StochPXElement -> A struct that contains Monte Carlo field configurations.
StochPXContext -> Essential struct for the StochPX solver.
#
solve          -> Wrapper function for the StochPX solver.
init           -> Initialize stochastic pole expansion simulation.
run (prun)     -> Perform stochastic pole expansion simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
sample         -> Sample field configurations via simulated annealing algorithm.
measure        -> Record Monte Carlo field configurations.
#
init_iodata    -> Preprocess the input data.
init_mc        -> Create a StochPXMC struct.
init_element   -> Create a StochPXElement struct.
init_context   -> Prepare data for a StochPXContext struct.
#
reset_mc       -> Reset counters in StochPXMC struct.
reset_element  -> Reset Monte Carlo field configurations.
reset_context  -> Reset Green's function and goodness-of-fit function.
#
calc_fmesh     -> Build very dense mesh for poles.
calc_lambda    -> Precompute 1 / (iωₙ - ϵ).
calc_green     -> Reproduce Green's function via the field configurations.
calc_chi2      -> Calculate goodness-of-fit function.
constraints    -> Limit the position of poles.
#
try_move_s     -> Shift positions of single pole.
try_move_p     -> Shift positions of two poles.
try_move_a     -> Change amplitudes of two poles.
try_move_x     -> Exchange amplitudes of two poles.
```
=#

#
include("spx.jl")
#
export StochPXElement
export StochPXContext
#
export solve
export init
export run
export prun
export average
export last
#
export sample
export measure
#
export init_iodata
export init_mc
export init_element
export init_context
#
export reset_mc
export reset_element
export reset_context
#
export calc_fmesh
export calc_lambda
export calc_green
export calc_chi2
export constraints
#
export try_move_s
export try_move_p
export try_move_a
export try_move_x

#=
### *Includes And Exports* : *base.jl*
=#

#=
*Summary* :

To provide basic workflow for the users of the ACFlow toolkit.

*Members* :

```text
solve       -> Select solver to solve the analytic continuation problem.
#
reprod      -> Try to generate the input data via calculated spectrum.
kramers     -> Calculate real part of response function.
#
setup_param -> Setup parameters.
read_param  -> Read parameters from case.toml.
read_data   -> Read the input data.
#
make_data   -> Preprocess the input data.
make_grid   -> Generate grid for the input data.
make_mesh   -> Generate mesh for the calculated spectrum.
make_model  -> Generate default model function.
make_kernel -> Generate kernel function.
```
=#

#
include("base.jl")
#
export solve
#
export reprod
export kramers
#
export setup_param
export read_param
export read_data
#
export make_data
export make_grid
export make_mesh
export make_model
export make_kernel

#=
### *PreCompile*
=#

export _precompile

"""
    _precompile()

Here, we would like to precompile the whole `ACFlow` toolkit to reduce
the runtime latency and speed up the successive calculations.
"""
function _precompile()
    prompt("Loading...")

    # Get an array of the names exported by the `ACFlow` module
    nl = names(ACFlow)

    # Go through each name
    cf = 0 # Counter
    for i in eachindex(nl)
        # Please pay attention to that nl[i] is a Symbol, we need to
        # convert it into string and function, respectively.
        str = string(nl[i])
        fun = eval(nl[i])

        # For methods only (macros must be excluded)
        if fun isa Function && !startswith(str, "@")
            # Increase the counter
            cf = cf + 1

            # Extract the signature of the function
            # Actually, `types` is a Core.SimpleVector.
            types = nothing
            try
                types = typeof(fun).name.mt.defs.sig.types
            catch
                @printf("Function %15s (#%3i) is skipped.\r", str, cf)
                continue
            end

            # Convert `types` from SimpleVector into Tuple
            # If length(types) is 1, the method is without arguments.
            T = ()
            if length(types) > 1
                T = tuple(types[2:end]...)
            end

            # Precompile them one by one
            #println(i, " -> ", str, " -> ", length(types), " -> ", T)
            precompile(fun, T)
            @printf("Function %24s (#%4i) is compiled.\r", str, cf)
        end
    end

    prompt("Well, ACFlow is compiled and loaded ($cf functions).")
    prompt("We are ready to go!")
    println()
    flush(stdout)
end

"""
    __init__()

This function would be executed immediately after the module is loaded
at runtime for the first time. It works at the REPL mode only.
"""
__init__() = begin
    isinteractive() && _precompile()
end

end # END OF MODULE
