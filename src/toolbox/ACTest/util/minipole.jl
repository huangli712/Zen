#!/usr/bin/env julia

#
# This script is used to start analytic continuation simulations with the
# MiniPole toolkit. It will launch only 1 process.
#
# This script won't work when the green.data.* files are unavailable.
#
# Usage:
#
#     $ julia minipole.jl act.toml std=false inds=[]
#
# The arguments `std` and `inds` are optional.
#
# (1) Perform normal test.
#
#     $ julia minipole.jl act.toml
#
# (2) Perform standard test (using ACT100 dataset).
#
#     $ julia minipole.jl act.toml std=true
#
# (3) Perform normal test, only tests 11, 12, and 13 are treated.
#
#     $ julia minipole.jl act.toml std=false inds=[11,12,13]
#
# (4) Perform standard test (using ACT100 dataset), only tests 1~40 are used.
#
#     $ julia minipole.jl act.toml std=true inds=1:40
#
# In order to execute this script correctly, please read the following
# instructions carefully.
#
# Important notes:
#
# 1. Install MiniPole
#
# Please download the source codes of MiniPole from github:
#
#    https://github.com/Green-Phys/MiniPole
#
# And then execute the following commands:
#
#    $ python setup.py install
#
# 2. Install PyCall.jl
#
# The MiniPole toolkit is written by the Python language. In order to call
# it from Julia script, we need the PyCall.jl package. Note that PyCall.jl
# doesn't belong to the Julia's standard library, we have to install it
# manually:
#
#    $ julia> using Pkg
#    $ julia> Pkg.add("PyCall")
#
# 3. Setup PyCall.jl
#
# PyCall.jl is shipped with a built-in Python runtime environment (it is
# probably in your .julia/conda directory). By default, PyCall.jl will use
# its own Python interpreter. But, we should enfore PyCall.jl to adopt the
# system's Python interpreter. So, please execute the following commands
# in Julia's REPL:
#
#    julia> ENV["PYTHON"] = "... path of the python executable ..."
#    julia> using Pkg
#    julia> Pkg.build("PyCall")
#
# To force Julia to use its own Python distribution, via Conda, simply set
# ENV["PYTHON"] to the empty string "" and re-run Pkg.build("PyCall"). To
# check which Python distribution is being used by PyCall.jl, please input
#
#    julia> using PyCall
#    julia> PyCall.pyprogramname
#
# 4. Setup act.toml
#
# Usually there are two blocks in the `act.toml` file. In the [Test] block,
# we have to change the `solver` parameter to "MiniPole". In the [Solver]
# block, please set the parameters for the MiniPole solver. A typical
# [Solver] block is as follows:
#
#     [Solver]
#     n0 = "auto"
#     n0_shift = 0
#     err = 1e-4
#     compute_const = true
#
# Please see documentation of MiniPole for all the possible parameters. Be
# careful, all the parameters are optional.
#
# 5. Tests about this script
#
# See actest/test/B01 and actest/test/B02.
#

haskey(ENV,"ACTEST_HOME") && pushfirst!(LOAD_PATH, ENV["ACTEST_HOME"])

using ACTest

using PyCall
using Printf
using DelimitedFiles

"""
    get_dict()

Prepare configurations for the MiniPole toolkit.

### Returns
* B -> A dict containing basic parameters for the test.
* S -> A dict containing basic parameters for analytic continuation solver.
"""
function get_dict()
    # General setup
    B = Dict{String,Any}(
        "finput" => "green.data",
        "ktype"  => get_t("ktype"),
        "grid"   => get_t("grid"),
        "mesh"   => get_t("mesh"),
        "ngrid"  => get_t("ngrid"),
        "nmesh"  => get_t("nmesh"),
        "wmax"   => get_t("wmax"),
        "wmin"   => get_t("wmin"),
        "beta"   => get_t("beta"),
        "offdiag" => get_t("fnpd"),
        "pmesh" => get_t("pmesh"),
    )

    # For analytic continuation solver
    cfg = inp_toml(query_args(), true)
    S = cfg["Solver"]

    return B, S
end

"""
    fix_dict!(i::I64, B::Dict{String,Any})

Fix configuration dynamically. This is essential for standard test. Note
that for standard test, the correlation functions could be fermionic or
bosonic, diagonal or non-diagonal. We have to make sure the configurations
are consistent with the original setups.

It seems that the MiniPole solver doesn't care about this.

### Arguments
* i -> Index for the current test.
* B -> A dict containing basic parameters for the test.

### Returns
`B` will be modified.
"""
function fix_dict!(i::I64, B::Dict{String,Any})
    # Get dicts for the standard test (ACT100)
    ACT100 = union(STD_FG, STD_FD, STD_FRD, STD_BG, STD_BD, STD_BRD)

    # We have to make sure ntest == 100
    ntest = get_t("ntest")
    @assert ntest == length(ACT100) == 100

    # Fix ktype, grid, and mesh
    B["ktype"] = ACT100[i]["ktype"]
    B["grid"] = ACT100[i]["grid"]
    B["mesh"] = ACT100[i]["mesh"]
    B["offdiag"] = ACT100[i]["fnpd"]
end

"""
    get_error(
        i::I64,
        mesh::Vector{F64},
        Aout::Vector{F64},
        B::Dict{String,Any}
    )

Evaluate error for the current test. It just calculates the distance
between the true and calculated spectral function. Note that the exact
spectral function should be read from the image.data.i file.

### Arguments
* i -> Index for the current test.
* mesh -> Real frequency mesh, ω.
* Aout -> Calculated spectral function by the MiniPole toolkit, A(ω).
* B -> A dict containing basic parameters for the test.
"""
function get_error(
    i::I64,
    mesh::Vector{F64},
    Aout::Vector{F64},
    B::Dict{String,Any}
    )
    # Read true solution
    data = readdlm("image.data." * string(i))
    ω = data[:,1]
    Ainp = data[:,2]

    # If there is a bosonic system, Ainp is actually A(ω) / ω. We should
    # convert it to A(ω).
    if B["ktype"] != "fermi"
        @. Ainp = Ainp * ω
    end

    # Calculate the difference
    error = trapz(mesh, abs.(Ainp .- Aout)) / trapz(mesh, abs.(Ainp))

    # Sometimes Aout could be extremely high δ-like peaks. We have to
    # take care of these cases.
    if error > 1000.0
        error = Inf
    end

    return error
end

"""
    write_summary(
        inds::Vector{I64},
        error::Vector{F64},
        ctime::Vector{F64}
    )

Write summary for the tests to external file `summary.data`.

### Arguments
* inds -> Indices for the tests.
* error -> Errors for the tests.
* ctime -> Elapsed times for the tests.
"""
function write_summary(
    inds::Vector{I64},
    error::Vector{F64},
    ctime::Vector{F64}
    )
    # Get number of tests
    ntest = length(inds)

    open("summary.data", "a") do fout
        println(fout, "# index            error         time (s) passed")
        #
        for i in inds
            @printf(fout, "%7i %16.12f %16.12f", i, error[i], ctime[i])
            if error[i] == 0.0
                @printf(fout, "%7s\n", "false")
            else
                @printf(fout, "%7s\n", "true")
            end
        end
        #
        println(fout, "# Number of tests: ", ntest)
        println(fout, "# Failed tests: ", count(x -> iszero(x), error[inds]))
        println(fout, "# Abnormal tests: ", count(x -> isinf(x), error[inds]))
    end
end

"""
    make_test(std::Bool = false, inds::Vector{I64} = I64[])

Perform analytic continuation simulations using the MiniPole toolkit. if
`std` is true, then the ACT100 dataset is considered. if `inds` is not
empty, then only the selected tests are handled.

### Arguments
* std -> Is the ACT100 dataset is adopted?
* inds -> A collection of indices for the tests.
"""
function make_test(std::Bool = false, inds::Vector{I64} = I64[])
    # Get number of tests (ntest).
    # cinds is used to store the indices of tests.
    ntest = get_t("ntest")
    if isempty(inds)
        cinds = collect(1:ntest)
    else
        cinds = inds
    end

    # Prepare some counters and essential variables
    nfail = 0
    nsucc = 0
    error = zeros(F64, ntest)
    ctime = zeros(F64, ntest)

    # Prepare configurations
    B, S = get_dict()

    # Start the loop
    for i in cinds
        @printf("Test -> %4i / %4i\n", i, ntest)
        #
        # Setup configurations further for the current test
        B["finput"] = "green.data." * string(i)
        # If we want to perform standard test, we have to change `ktype`
        # and `grid` parameters dynamically.
        if std == true
            println("Note: the act100 dataset is being used!")
            fix_dict!(i, B)
        end
        #
        # Prepare real mesh for spectral function
        mesh = make_mesh(B["ktype"])
        #
        # Transfer parameters to the MiniPole toolkit
        py"setup_param"(B, S, mesh.mesh)
        #
        try
            # Solve the analytic continuation problem.
            # The elapsed time is recorded as well.
            start = time_ns()
            mesh, Aout = py"solve"()
            finish = time_ns()
            #
            # Backup the calculated results for further analytics
            cp("Aout.data", "Aout.data." * string(i), force = true)
            cp("Gout.data", "Gout.data." * string(i), force = true)
            cp("repr.data", "repr.data." * string(i), force = true)
            #
            # Calculate the accuracy / error
            error[i] = get_error(i, mesh, Aout, B)
            ctime[i] = (finish - start) / 1e9
            #
            # Increase the counter
            @printf("Accuracy: %16.12f\n", error[i])
            nsucc = nsucc + 1
        catch ex
            error[i] = 0.0
            ctime[i] = 0.0
            nfail = nfail + 1
            println("Something wrong for test case $i")
        end
        #
        println()
    end

    @assert nfail + nsucc == length(cinds)
    println("Only $nsucc / $(length(cinds)) tests can survive!")
    println("Please check summary.data for more details.")
    println()

    # Write summary for the test
    write_summary(cinds, error, ctime)
end

"""
    main()

Entry of this script. It will parse the command line arguments and call
the corresponding functions.
"""
function main()
    nargs = length(ARGS)

    # Besides the case.toml, no arguments.
    #
    # $ minipole.jl act.toml
    if nargs == 1
        make_test()
    end

    # Two arguments. Besides the case.toml, we can specify whether the
    # ACT100 dataset is used.
    #
    # $ minipole.jl act.toml std=true
    if nargs == 2
        std = parse(Bool, split(ARGS[2],"=")[2])
        make_test(std)
    end

    # Three arguments. We can specify whether the ACT100 dataset is used,
    # and the indices of selected tests.
    #
    # $ minipole.jl act.toml std=true inds=[11,12,13]
    # $ minipole.jl act.toml std=true inds=11:13
    if nargs == 3
        std = parse(Bool, split(ARGS[2],"=")[2])
        str = split(ARGS[3],"=")[2]
        if contains(str, ",")
            inds = parse.(Int, split(chop(str; head=1, tail=1), ','))
        else
            arr = parse.(Int, split(str, ':'))
            inds = collect(arr[1]:arr[2])
        end
        make_test(std, inds)
    end
end

"""
    python()

Define the interface to the MiniPole toolkit. Now it doesn't support the
DLR feature of the MiniPole toolkit.
"""
function python()
    py"""
    import numpy as np
    from mini_pole import MiniPole

    # Global variables
    #
    # For general configurations
    _B = None
    #
    # For real mesh
    _ω = None
    #
    # For MiniPole's parameters. See mini_pole.py for more details.
    _P = {
        'n0' : "auto",
        'n0_shift' : 0,
        'err' : None,
        'err_type' : "abs",
        'M' : None,
        'symmetry' : False,
        'G_symmetric' : False,
        'compute_const' : False,
        'plane' : None,
        'include_n0' : False,
        'k_max' : 999,
        'ratio_max' : 10
    }

    # Update the configurations
    def setup_param(B, S, ω):
        global _B, _ω, _P
        #
        _B = B
        _ω = ω
        #
        # We should scan all the keys in S, and check whether it is valid.
        for k in S.keys():
            if k not in _P:
                print("error: this parameter " + k + " is not supported")
                import sys
                sys.exit(-1)
            else:
                _P[k] = S[k]

    # Read Matsubara data from input file
    def read_data():
        iωₙ, Gᵣ, Gᵢ = np.loadtxt(
            _B["finput"],
            unpack = True,
            usecols = (0,1,2)
        )
        G = Gᵣ + Gᵢ * 1j
        return iωₙ, G

    # Write A(ω), G(ω), and G(iωₙ) to external files
    def write_data(iωₙ, Gout, Grep, Aout):
        with open("Gout.data", "w") as f:
            for i in range(_ω.size):
                print(_ω[i], Gout[i].real, Gout[i].imag, file = f)
        #
        with open("repr.data", "w") as f:
            for i in range(iωₙ.size):
                print(iωₙ[i], Grep[i].real, Grep[i].imag, file = f)
        #
        with open("Aout.data", "w") as f:
            for i in range(_ω.size):
                print(_ω[i], Aout[i], file = f)

    # Calculate Green's function by pole representation
    def calc_green(z, 𝔸, 𝕏):
        Gz = 0.0
        for i in range(𝕏.size):
            Gz += 𝔸[i] / (z - 𝕏[i])
        return Gz

    # Calculate spectral function
    def calc_spectrum(G):
        return -1.0 / np.pi * G.imag

    # Solve the analytic continuation problem by the MiniPole solver
    def solve():
        # Read Matsubara data
        iωₙ, G = read_data()
        #
        # Solve the problem
        p = MiniPole(
            G, iωₙ,
            n0 = _P["n0"],
            n0_shift = _P["n0_shift"],
            err = _P["err"],
            err_type = _P["err_type"],
            M = _P["M"],
            symmetry = _P["symmetry"],
            G_symmetric = _P["G_symmetric"],
            compute_const = _P["compute_const"],
            plane = _P["plane"],
            include_n0 = _P["include_n0"],
            k_max = _P["k_max"],
            ratio_max = _P["ratio_max"]
        )
        #
        # Get pole representation
        location = p.pole_location
        weight = p.pole_weight.reshape(-1)
        #
        # Calculate G(ω), G(iωₙ), and A(ω)
        Gout = calc_green(_ω, weight, location)
        Grep = calc_green(iωₙ * 1j, weight, location)
        Aout = calc_spectrum(Gout)
        #
        # Write analytic continuation results
        write_data(iωₙ, Gout, Grep, Aout)

        return _ω, Aout
    """
end

python()
welcome()
overview()
read_param()
main()
goodbye()
