#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/08
#

#
# Driver Functions
#

"""
    ready()

Examine whether all the conditions, including input files and working
directories, for DFT + DMFT calculations are ready.

See also: [`go`](@ref).
"""
function ready()
    # R1: Check the input files
    query_inps(get_d("engine"))

    # R2: Prepare the working directories
    make_trees()
end

"""
    go()

Dispatcher for DFT + DMFT calculations. Note that it can not call the
`cycle3()`-`cycle8()` functions.

See also: [`ready`](@ref).
"""
function go()
    # Get calculation mode
    mode = get_m("mode")

    # Choose suitable computational driver
    @cswitch mode begin
        # One-shot DFT + DMFT calculations
        @case 1
            cycle1()
            break

        # Fully self-consistent DFT + DMFT calculations
        @case 2
            cycle2()
            break

        # To be implemented
        @default
            sorry()
            break
    end
end

"""
    final()

Finalize the DFT + DMFT calculations.

See also: [`go`](@ref).
"""
function final()
    sorry()
end

#=
*Remarks 1*:

We would like to perform two successive DFT runs if `get_d("loptim")`` is
true. The purpose of the first DFT run is to evaluate the fermi level.
Then an energy window is determined. We will use this window to generate
optimal projectors in the second DFT run.

On the other hand, if `get_d("loptim")` is false, only the first DFT run
is enough.
=#

#=
*Remarks 2*:

We want better *optimal projectors*.

In the previous DFT run, `initial` fermi level = 0 -> `wrong` energy
window -> `wrong` optimial projectors. But at this point, the fermi
level is updated, so we have to generate the optimal projectors
again within this new window by doing addition DFT calculation.
=#

#=
*Remarks 3*:

The key Kohn-Sham data inclue lattice structures, 𝑘-mesh and its weights,
tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. At
first, the adaptor will read in these data from the output files of DFT
engine. And then it will process the raw projectors (such as parsing,
labeling, grouping, filtering, and rotatation). Finally, the adaptor will
write down the processed data to some specified files using the `IR`
format.
=#

#=
*Remarks 4*:

Now everything is ready. We are going to solve the DMFT self-consistent
equation iterately.
=#

"""
    cycle1()

Perform one-shot DFT + DMFT calculations. In other words, the charge
density won't be fed back to the DFT engine. The self-consistency is
only achieved at the DMFT level.

See also: [`cycle2`](@ref), [`go`](@ref).
"""
function cycle1()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

#
# Initialization (C01-C04)
#
    prompt("ZEN", "Initialization")

    # C01: Perform DFT calculation (for the first time)
    dft_run(it, lr)

    # C02: Perform DFT calculation (for the second time)
    if get_d("loptim")
        dft_run(it, lr)
    end

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    adaptor_run(it, lr, ai)

    # C04: Prepare default self-energy functions
    sigma_core(it, lr, ai, "reset")

#
# DFT + DMFT Iterations (C05-C12)
#
    prompt("ZEN", "Iterations")
    save_it(it, lr)

    for iter = 1:get_m("niter")
        # Print the log
        prompt("ZEN", "Cycle $iter")
        prompt(lr.log, "")
        prompt(lr.log, "< dft_dmft_cycle >")

        # Update IterInfo struct
        incr_it(it)

        # C05: Tackle with the double counting term
        sigma_core(it, lr, ai, "dcount")

        # C06: Perform DMFT calculation with `task` = 1
        dmft_run(it, lr, 1)

        # C07: Mix the hybridization functions
        mixer_core(it, lr, ai, "delta")

        # C08: Mix the local impurity levels
        mixer_core(it, lr, ai, "eimpx")

        # C09: Split and distribute the hybridization functions
        sigma_core(it, lr, ai, "split")

        # C10: Solve the quantum impurity problems
        solver_run(it, lr, ai)

        # C11: Gather and combine the impurity self-functions
        sigma_core(it, lr, ai, "gather")

        # C12: Mix the impurity self-energy functions
        mixer_core(it, lr, ai, "sigma")

        # Print the cycle info
        save_it(it, lr)
    end

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle2()

Perform fully self-consistent DFT + DMFT calculations. The self-consistency
is achieved at both DFT and DMFT levels.

See also: [`cycle1`](@ref), [`go`](@ref).
"""
function cycle2()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

#
# Initialization (C01-C04)
#
    prompt("ZEN", "Initialization")

    # C01: Perform DFT calculation (for the first time)
    dft_run(it, lr)

    # C02: Perform DFT calculation (for the second time)
    if get_d("loptim")
        dft_run(it, lr)
    end

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    adaptor_run(it, lr, ai)

    # C04: Prepare default self-energy functions
    sigma_core(it, lr, ai, "reset")

end

"""
    cycle3()

Perform DFT calculations only. If there are something wrong, then you
have chance to adjust the DFT input files manually (for example, you
can modify `vasp_incar()/vasp.jl` by yourself).

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle3()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C01: Perform DFT calculation (for the first time)
    dft_run(it, lr)

    # C02: Perform DFT calculation (for the second time)
    if get_d("loptim")
        dft_run(it, lr)
    end

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle4(task::I64)

Perform DMFT calculations only. The users can execute it in the REPL mode
to see whether the DMFT engine works properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle4(task::I64)
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C01: Execuate the DMFT engine
    dmft_run(it, lr, task)

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle5()

Perform calculations using quantum impurity solvers only. The users can
execute it in the REPL mode to see whether the quantum impurity solvers
work properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle5()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execuate the quantum impurity solvers
    solver_run(it, lr, ai)

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle6()

Perform calculations using Kohn-Sham adaptor only. The users can execute
it in the REPL mode to see whether the Kohn-Sham adaptor works properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle6()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execute the Kohn-Sham adaptor
    adaptor_run(it, lr, ai)

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle7(task::String = "reset")

Perform calculations using self-energy engine only. The users can execute
it in the REPL mode to see whether the self-energy engine works properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle7(task::String = "reset")
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execute the Kohn-Sham adaptor
    sigma_core(it, lr, ai, task)

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle8(task::String = "sigma")

Perform calculations using mixer engine only. The users can execute
it in the REPL mode to see whether the mixer engine works properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function cycle8(task::String = "sigma")
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Further setup the IterInfo struct
    it.I₃ = 1 # Please modify them to fit your requirements
    it.I₁ = 10

    # C02: Execute the Kohn-Sham adaptor
    mixer_core(it, lr, ai, task)

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

#
# Service Functions
#

#=
*Remarks*:

In order to terminate the `Zen` code, the following two conditions
should be fulfilled at the same time.

* The argument `force_exit` is true.

* The `case.stop` file exists (from `query_stop()``).

We usually use this functon to stop the whole DFT + DMFT iterations.
=#

"""
    monitor(force_exit::Bool = false)

Determine whether we need to terminate the Zen code.

See also: [`query_stop`](@ref).
"""
function monitor(force_exit::Bool = false)
    if force_exit && query_stop()
        exit(-1)
    end
end

"""
    incr_it(it::IterInfo)

Modify the internal counters in IterInfo struct.

See also: [`IterInfo`](@ref), [`Logger`](@ref).
"""
function incr_it(it::IterInfo)
    # Get calculation mode
    mode = get_m("mode")

    # For one-shot DFT + DMFT mode
    if mode == 1
        it.I₃ = 1
        it.I₁ = it.I₁ + 1
        it.I₄ = it.I₄ + 1
    # For fully charge self-consistent DFT + DMFT mode
    else
        sorry()
    end
end

"""
    save_it(it::IterInfo, lr::Logger)

Try to record the iteration information in the `case.cycle` file.

See also: [`IterInfo`](@ref), [`Logger`](@ref).
"""
function save_it(it::IterInfo, lr::Logger)
    # Extract parameter `nsite`
    nsite = get_i("nsite")
    @assert nsite == length(it.nf)

    # Write the header
    if it.I₄ == 0
        print(lr.cycle, "#   #   #   #   μ₀        μ₁        ")
        for t = 1:nsite
            print(lr.cycle, "dc$(subscript(t))       ")
        end
        for t = 1:nsite
            print(lr.cycle, "nf$(subscript(t))       ")
        end
        println(lr.cycle, "et")
        # Write separator
        println(lr.cycle, repeat('-', 46 + 20*nsite))
    # Write iteration information
    else
        @printf(lr.cycle, "%-4i", it.I₄)
        @printf(lr.cycle, "%-4i", it.I₃)
        @printf(lr.cycle, "%-4i", it.I₁)
        @printf(lr.cycle, "%-4i", it.I₂)
        if it.μ₀ < 0.0
            @printf(lr.cycle, "%-10.5f", it.μ₀)
        else
            @printf(lr.cycle, "+%-9.5f", it.μ₀)
        end
        if it.μ₁ < 0.0
            @printf(lr.cycle, "%-10.5f", it.μ₁)
        else
            @printf(lr.cycle, "+%-9.5f", it.μ₁)
        end
        for t = 1:nsite
            @printf(lr.cycle, "%-10.5f", it.dc[t])
        end
        for t = 1:nsite
            @printf(lr.cycle, "%-10.5f", it.nf[t])
        end
        @printf(lr.cycle, "%-10.5f", it.et)
        println(lr.cycle)
    end

    # Flush the IOStream
    flush(lr.cycle)
end

#=
*Remarks*:

The working directories include `dft`, `dmft1`, `dmft2`, and `impurity.i`.
If they exist already, it would be better to remove them at first.
=#

"""
    make_trees()

Prepare the working directories at advance.

See also: [`rm_trees`](@ref).
"""
function make_trees()
    # Build an array for folders
    dir_list = ["dft", "dmft1", "dmft2"]
    for i = 1:get_i("nsite")
        push!(dir_list, "impurity.$i")
    end

    # Go through these folders, create them one by one.
    for i in eachindex(dir_list)
        dir = dir_list[i]
        if isdir(dir)
            rm(dir, force = true, recursive = true)
        end
        mkdir(dir)
    end
end

"""
    rm_trees()

Remove the working directories finally.

See also: [`make_trees`](@ref).
"""
function rm_trees()
    # Build an array for folders
    dir_list = ["dft", "dmft1", "dmft2"]
    for i = 1:get_i("nsite")
        push!(dir_list, "impurity.$i")
    end

    # Go through these folders, remove them one by one.
    for i in eachindex(dir_list)
        dir = dir_list[i]
        if isdir(dir)
            rm(dir, force = true, recursive = true)
        end
    end
end

"""
    dft_run(it::IterInfo, lr::Logger)

Simple driver for DFT engine. It performs three tasks: (1) Examine
the runtime environment for the DFT engine. (2) Launch the DFT engine.
(3) Backup the output files by DFT engine for next iterations.

Now only the VASP engine is supported. If you want to support the other
DFT engine, this function must be adapted.

See also: [`adaptor_run`](@ref), [`dmft_run`](@ref), [`solver_run`](@ref).
"""
function dft_run(it::IterInfo, lr::Logger)
    # Determine the chosen engine
    engine = get_d("engine")

    # Print the log
    prompt("DFT")
    prompt(lr.log, engine)

    # Enter dft directory
    cd("dft")

    # Activate the chosen DFT engine
    @cswitch engine begin
        # For VASP
        @case "vasp"
            vasp_init(it)
            vasp_exec(it)
            vasp_save(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    dmft_run(it::IterInfo, lr::Logger, task::I64)

Simple driver for DMFT engine. It performs three tasks: (1) Examine
the runtime environment for the DMFT engine. (2) Launch the DMFT engine.
(3) Backup the output files by DMFT engine for next iterations.

The argument `task` is used to specify running mode of the DMFT code.

See also: [`adaptor_run`](@ref), [`dft_run`](@ref), [`solver_run`](@ref).
"""
function dmft_run(it::IterInfo, lr::Logger, task::I64)
    # Examine the argument `task`
    @assert task === 1 || task === 2

    # Print the log
    prompt("DMFT")
    prompt(lr.log, "dmft$task")

    # Enter dmft1 or dmft2 directory
    cd("dmft$task")

    # Activate the chosen DMFT engine
    @cswitch task begin
        # Solve the DMFT self-consistent equation
        @case 1
            dmft_init(it, 1)
            dmft_exec(it, 1)
            dmft_save(it, 1)
            break

        # Generate DMFT correction for DFT charge density
        @case 2
            dmft_init(it, 2)
            dmft_exec(it, 2)
            dmft_save(it, 2)
            break
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    solver_run(it::IterInfo, lr::Logger, ai::Array{Impurity,1})

Simple driver for quantum impurity solvers. It performs three tasks: (1)
Examine the runtime environment for quantum impurity solver. (2) Launch
the quantum impurity solver. (3) Backup output files by quantum impurity
solver for next iterations.

Now only the `ct_hyb1`, `ct_hyb2`, `hub1`, and `norg` quantum impurity
solvers are supported. If you want to support the other quantum impurity
solvers, this function must be adapted.

See also: [`adaptor_run`](@ref), [`dft_run`](@ref), [`dmft_run`](@ref).
"""
function solver_run(it::IterInfo, lr::Logger, ai::Array{Impurity,1})
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Determine the chosen solver
        engine = get_s("engine")

        # Extract the Impurity strcut
        imp = ai[i]

        # Print the log
        prompt("Solvers")
        prompt(lr.log, engine)

        # Enter impurity.i directory
        cd("impurity.$i")

        # Activate the chosen quantum impurity solver
        @cswitch engine begin
            @case "ct_hyb1"
                s_qmc1_init(it, imp)
                s_qmc1_exec(it)
                s_qmc1_save(it, imp)
                break

            @case "ct_hyb2"
                s_qmc2_init(it)
                s_qmc2_exec(it)
                s_qmc2_save(it)
                break

            @case "hub1"
                s_hub1_init(it)
                s_hub1_exec(it)
                s_hub1_save(it)
                break

            @case "norg"
                s_norg_init(it)
                s_norg_exec(it)
                s_norg_save(it)
                break

            @default
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end

    # Monitor the status
    monitor(true)
end

"""
    adaptor_run(it::IterInfo, lr::Logger, ai::Array{Impurity,1})

Simple driver for the adaptor. It performs three tasks: (1) Initialize
the adaptor, to check whether the essential files exist. (2) Parse the
Kohn-Sham data output by the DFT engine, try to preprocess them, and
then transform them into IR format. (3) Backup the files by adaptor.

For the first task, only the VASP adaptor is supported. While for the
second task, only the PLO adaptor is supported. If you want to support
more adaptors, please adapt this function.

See also: [`dft_run`](@ref), [`dmft_run`](@ref), [`solver_run`](@ref).
"""
function adaptor_run(it::IterInfo, lr::Logger, ai::Array{Impurity,1})
    # Enter dft directory
    cd("dft")

    #
    # A0: Create a dict named DFTData
    #
    # This dictionary is for storing the Kohn-Sham band structure and
    # related data. The key-value pairs would be inserted into this
    # dict dynamically.
    #
    DFTData = Dict{Symbol,Any}()

    #
    # A1: Parse the original Kohn-Sham data
    #
    # Choose suitable driver function according to DFT engine. The
    # Kohn-Sham data will be stored in the DFTData dict.
    #
    engine = get_d("engine")
    prompt("Adaptor")
    prompt(lr.log, "adaptor::$engine")
    @cswitch engine begin
        # For VASP
        @case "vasp"
            vasp_files()
            vasp_adaptor(DFTData)
            break

        @default
            sorry()
            break
    end

    #
    # A2: Process the original Kohn-Sham data
    #
    # Well, now we have the Kohn-Sham data. But they can not be used
    # directly. We have to check and process them carefully. Please
    # pay attention to that the DFTData dict will be modified in
    # the `plo_adaptor()` function.
    #
    # The plo_adaptor() function also has the ability to calculate
    # some selected physical quantities (such as overlap matrix and
    # density of states) to check the correctness of the Kohn-Sham
    # data. This feature will be activated automatically if you are
    # using the `src/tools/test.jl` tool to examine the DFT data.
    #
    projtype = get_d("projtype")
    prompt("Adaptor")
    prompt(lr.log, "adaptor::$projtype")
    @cswitch projtype begin
        # For projected local orbital scheme
        @case "plo"
            plo_adaptor(DFTData, ai)
            break

        # For maximally localized wannier function scheme
        @case "wannier"
            sorry()
            break

        @default
            sorry()
            break
    end

    #
    # A3: Output the processed Kohn-Sham data
    #
    # Ok, now the Kohn-Sham data are ready. We would like to write them
    # to some specified files with the IR format. Then these files will
    # be saved immediately.
    #
    prompt("Adaptor")
    prompt(lr.log, "adaptor::ir")
    ir_adaptor(DFTData)
    ir_save(it)

    #
    # A4: Clear the DFTData dict
    #
    empty!(DFTData)
    @assert isempty(DFTData)

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    sigma_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1}, task::String = "reset")

Simple driver for functions for processing the self-energy functions
and hybridization functions.

Now it supports four tasks: `reset`, `dcount`, `split`, `gather`. It
won't change the current directory.

See also: [`mixer_core`](@ref).
"""
function sigma_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1}, task::String = "reset")
    # Check the given task
    @assert task in ("reset", "dcount", "split", "gather")

    # Print the log
    prompt("Sigma")
    prompt(lr.log, "sigma::$task")

    # Launch suitable subroutine
    @cswitch task begin
        # Generate default self-energy functions and store them
        @case "reset"
            sigma_reset(ai)
            break

        # Calculate the double counting term and store it
        @case "dcount"
            sigma_dcount(it, ai)
            break

        # Split the hybridization functions and store them
        @case "split"
            sigma_split(ai)
            break

        # Collect impurity self-energy functions and combine them
        @case "gather"
            sigma_gather(it, ai)
            break

        @default
            sorry()
            break
    end

    # Monitor the status
    monitor(true)
end

"""
    mixer_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1}, task::String = "sigma")

Simple driver for the mixer. It will try to mix the self-energy functions
or hybridization functions and generate a new one.

Now it supports four tasks: `sigma`, `delta`, `eimpx`, `gamma`. It
won't change the current directory.

See also: [`sigma_core`](@ref).
"""
function mixer_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1}, task::String = "sigma")
    # Check the given task
    @assert task in ("sigma", "delta", "eimpx", "gamma")

    # Check iteration number
    if it.I₃ == 1
        if it.I₁ < 2
            return
        end
    end

    # Print the log
    prompt("Mixer")
    prompt(lr.log, "mixer::$task")

    # Launch suitable subroutine
    @cswitch task begin
        # Try to mix the self-energy functions
        @case "sigma"
            mixer_sigma(it, ai)
            break

        # Try to mix the hybridization functions
        @case "delta"
            mixer_delta(it, ai)
            break

        # Try to mix the local impurity levels
        @case "eimpx"
            mixer_eimpx(it, ai)
            break

        # Try to mix the charge density
        @case "gamma"
            mixer_gamma(it)
            break

        @default
            sorry()
            break
    end

    # Monitor the status
    monitor(true)
end
