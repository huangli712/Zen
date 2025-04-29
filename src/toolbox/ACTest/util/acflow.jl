#!/usr/bin/env julia

#
# This script is used to start analytic continuation simulations with the
# ACFlow toolkit. It will launch only 1 process.
#
# This script can be easily modified to support the other analytic
# configuration tools / methods, or support parallel calculations.
#
# Usage:
#
#     $ ./acflow.jl act.toml std=false inds=[]
#
# The arguments `std` and `inds` are optional.
#
# (1) Perform normal test.
#
#     $ ./acflow.jl act.toml
#
# (2) Perform standard test (using ACT100 dataset).
#
#     $ ./acflow.jl act.toml std=true
#
# (3) Perform normal test, only tests 11, 12, and 13 are treated.
#
#     $ ./acflow.jl act.toml std=false inds=[11,12,13]
#
# (4) Perform standard test (using ACT100 dataset), only tests 1~40 are used.
#
#     $ ./acflow.jl act.toml std=true inds=1:40
#

haskey(ENV,"ACTEST_HOME") && pushfirst!(LOAD_PATH, ENV["ACTEST_HOME"])
haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACTest
using ACFlow:setup_param
using ACFlow:read_data
using ACFlow:solve

using Printf
using DelimitedFiles

# Prepare configurations for the ACFlow toolkit
function get_dict()
    # General setup
    B = Dict{String,Any}(
        "finput" => "green.data",
        "solver" => get_t("solver"),
        "ktype"  => get_t("ktype"),
        "mtype"  => "flat",
        "grid"   => get_t("grid"),
        "mesh"   => get_t("mesh"),
        "ngrid"  => get_t("ngrid"),
        "nmesh"  => get_t("nmesh"),
        "wmax"   => get_t("wmax"),
        "wmin"   => get_t("wmin"),
        "beta"   => get_t("beta"),
        "offdiag" => get_t("offdiag"),
        "pmesh" => get_t("pmesh"),
    )

    # For analytic continuation solver
    cfg = inp_toml(query_args(), true)
    S = cfg["Solver"]

    return B, S
end

# Fix configuration dynamically. This is essential for standard test.
# Note that for standard test, the correlation functions could be
# fermionic or bosonic, diagonal or non-diagonal. We have to make sure
# the configurations are consistent with the original setups.
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
    B["offdiag"] = ACT100[i]["offdiag"]
end

# Evaluate error for the current test. It just calculates the distance
# between the true and calculated spectral function.
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
        # For the following solvers, their outputs (`Aout`) are A(ω) / ω
        # as well. We have to convert them to A(ω) too.
        if B["solver"] in ("MaxEnt", "StochAC", "StochSK", "StochOM")
            @. Aout = Aout * mesh
        end
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

# Write summary for the tests to external file `summary.data`
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

# Perform analytic continuation simulations using the ACFlow toolkit.
# if `std` is true, then the ACT100 dataset is considered.
# if `inds` is not empty, then only the selected tests are handled.
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
        # Transfer parameters to the ACFlow toolkit
        setup_param(B, S)
        #
        try
            # Solve the analytic continuation problem.
            # The elapsed time is recorded as well.
            start = time_ns()
            mesh, Aout, _ = solve(read_data())
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

# Entry of this script. It will parse the command line arguments and call
# the corresponding functions.
function main()
    nargs = length(ARGS)

    # Besides the case.toml, no arguments.
    #
    # $ acflow.jl act.toml
    if nargs == 1
        make_test()
    end

    # Two arguments. Besides the case.toml, we can specify whether the
    # ACT100 dataset is used.
    #
    # $ acflow.jl act.toml std=true
    if nargs == 2
        std = parse(Bool, split(ARGS[2],"=")[2])
        make_test(std)
    end

    # Three arguments. We can specify whether the ACT100 dataset is used,
    # and the indices of selected tests.
    #
    # $ acflow.jl act.toml std=true inds=[11,12,13]
    # $ acflow.jl act.toml std=true inds=11:13
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

welcome()
overview()
read_param()
main()
goodbye()
