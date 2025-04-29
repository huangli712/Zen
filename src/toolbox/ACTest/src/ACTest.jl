#
# Project : Lily
# Source  : ACTest.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/28
#

"""
    ACTest

ACTest is an automatical generator for fermionic or bosonic correlation
functions. It is used to generate a lot of correlators and corresponding
spectral functions. They are then employed to test the newly developed
analytic continuation methods or tools.

Now ACTest is used together with the `ACFlow` and `MiniPole` toolkits.
But it is easily to be interfaced with the other analytic continuation
codes.

For more details about how to obtain, install and use the ACTest toolkit,
please visit the following website:

* `https://huangli712.github.io/projects/actest/index.html`

Any suggestions, comments, and feedbacks are welcome. Enjoy it!
"""
module ACTest

#=
### *Using Standard Libraries*
=#

using Distributed
using LinearAlgebra
using Statistics
using Random
using Dates
using Printf
using TOML

#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the ACTest toolkit.

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
authors       -> Print the authors of ACTest to screen.
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
PTEST           -> Configuration dict for general setup.
#
AbstractGrid    -> Abstract mesh for input data.
FermionicImaginaryTimeGrid -> Grid in fermionic imaginary time axis.
FermionicMatsubaraGrid -> Grid in fermionic Matsubara frequency axis.
BosonicImaginaryTimeGrid -> Grid in bosonic imaginary time axis.
BosonicMatsubaraGrid -> Grid in bosonic Matsubara frequency axis.
#
AbstractMesh    -> Abstract grid for calculated spectral function.
LinearMesh      -> Linear mesh.
TangentMesh     -> Tangent mesh.
LorentzMesh     -> Lorentzian mesh.
HalfLorentzMesh -> Lorentzian mesh at half-positive axis.
DynamicMesh     -> Dynamic mesh. It represents the mesh that read from file.
#
AbstractPeak    -> Abstract peak in spectral function.
GaussianPeak    -> Gaussian peak in spectral function.
LorentzianPeak  -> Lorentzian peak in spectral function.
DeltaPeak       -> δ-like peak in spectral function.
RectanglePeak   -> Rectangle peak in spectral function.
RiseDecayPeak   -> Rise-and-Decay peak in spectral function.
#
AbstractFunction -> Abstract function.
SpectralFunction -> Spectral function.
GreenFunction   -> Green's function.
```
=#

#
include("types.jl")
#
export DType
export ADT
#
export PTEST
#
export AbstractGrid
export FermionicImaginaryTimeGrid
export FermionicMatsubaraGrid
export BosonicImaginaryTimeGrid
export BosonicMatsubaraGrid
#
export AbstractMesh
export LinearMesh
export TangentMesh
export LorentzMesh
export HalfLorentzMesh
export DynamicMesh
#
export AbstractPeak
export GaussianPeak
export LorentzianPeak
export DeltaPeak
export RectanglePeak
export RiseDecayPeak
#
export AbstractFunction
export SpectralFunction
export GreenFunction

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
welcome       -> Print welcome message.
overview      -> Print runtime information of ACTest.
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

To provide some numerical algorithms, such as numerical integrations
and Einstein summation notation.

*Members* :

```text
trapz     -> Numerical integration (composite trapezoidal rule).
simpson   -> Numerical integration (simpson rule).
#
@einsum   -> Macro for Einstein summation notation.
```
=#

#
include("math.jl")
#
export trapz
export simpson
#
export @einsum

#=
### *Includes And Exports* : *grid.jl*
=#

#=
*Summary* :

To implement various grids for the correlation function data.

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
### *Includes And Exports* : *peak.jl*
=#

#=
*Summary* :

To generate various peaks, which should be used to construct the spectrum.

*Members* :

```text
𝑝(ω) -> Function call to peak generator.
```
=#

#
include("peak.jl")
#

#=
### *Includes And Exports* : *config.jl*
=#

#=
*Summary* :

To extract, parse, verify, and print the configuration parameters.
They are stored in external files (case.toml) or dictionaries.

*Members* :

```text
inp_toml -> Parse case.toml, return raw configuration information.
fil_dict -> Fill dicts for configuration parameters.
see_dict -> Display all the relevant configuration parameters.
rev_dict -> Update dict (PTEST) for configuration parameters.
chk_dict -> Check dicts for configuration parameters.
_v       -> Verify dict's values.
get_t    -> Extract value from dict (PTEST dict), return raw value.
```
=#

#
include("config.jl")
#
export inp_toml
export fil_dict
export see_dict
export rev_dict
export chk_dict
export _v
export get_t

#=
### *Includes And Exports* : *inout.jl*
=#

#=
*Summary* :

To read the input data or write the calculated results.

*Members* :

```text
write_spectrum -> Write spectral functions.
write_backward -> Write reproduced input data in imaginary axis.
```
=#

#
include("inout.jl")
#
export write_spectrum
export write_backward

#=
### *Includes And Exports* : *kernel.jl*
=#

#=
*Summary* :

To define various kernel functions.

*Members* :

```text
build_kernel      -> Build kernel function.
build_kernel_symm -> Build kernel function for symmetric case.
```
=#

#
include("kernel.jl")
#
export build_kernel
export build_kernel_symm

#=
### *Includes And Exports* : *dataset.jl*
=#

#=
*Summary* :

To define a standard dataset (namely ACT100), which contains 100 typical
spectral functions.

*Members* :

```text
STD_FG  -> Dicts for fermionic systems (continuum spectra).
STD_FD  -> Dicts for fermionic systems (discrete spectra).
STD_FRD -> Dicts for fermionic systems (non-positive definite spectra).
#
STD_BG  -> Dicts for bosonic systems (continuum spectra).
STD_BD  -> Dicts for bosonic systems (discrete spectra).
STD_BRD -> Dicts for bosonic systems (non-positive definite spectra).
```
=#

#
include("dataset.jl")
#
export STD_FG
export STD_FD
export STD_FRD
#
export STD_BG
export STD_BD
export STD_BRD

#=
### *Includes And Exports* : *base.jl*
=#

#=
*Summary* :

To provide basic workflow for the users of the ACTest toolkit.

*Members* :

```text
reprod      -> Try to generate the correlator via calculated spectrum.
#
setup_param -> Setup parameters.
read_param  -> Read parameters from case.toml.
#
make_data_std -> Generate standard dataset (namely ACT100).
make_data   -> Generate spectral functions and corresponding correlators.
make_peak   -> Generate various peaks.
make_spectrum -> Generate spectral function by peaks.
make_noise  -> Generate noise.
make_green  -> Generate Green's functions at imaginary axis.
make_grid   -> Generate grid for the Green's function data.
make_mesh   -> Generate mesh for the calculated spectrum.
make_kernel -> Generate kernel function.
```
=#

#
include("base.jl")
#
export reprod
#
export setup_param
export read_param
#
export make_data_std
export make_data
export make_peak
export make_spectrum
export make_noise
export make_green
export make_grid
export make_mesh
export make_kernel

#=
### *PreCompile*
=#

export _precompile

"""
    _precompile()

Here, we would like to precompile the whole `ACTest` toolkit to reduce
the runtime latency and speed up the successive calculations.
"""
function _precompile()
    prompt("Loading...")

    # Get an array of the names exported by the `ACTest` module
    nl = names(ACTest)

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

    prompt("Well, ACTest is compiled and loaded ($cf functions).")
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

end
