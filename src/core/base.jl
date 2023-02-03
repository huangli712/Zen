#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/20
#

#=
### *Driver Functions* : *Layer 1*
=#

"""
    ready(mkdir::Bool = true)

Examine whether all the conditions (including input files and working
directories) for DFT + DMFT calculations are ready.

See also: [`go`](@ref).
"""
function ready(mkdir::Bool = true)
    # Print the header
    prompt("Ready")

    # R1: Setup _engine_, the DFT backend.
    engine = get_d("engine")
    global _engine_ = str_to_struct(engine, "Engine")
    println("DFT backend: $engine -> $(typeof(_engine_))")

    # R2: Setup _solver_, the quentum impurity solver.
    solver = get_s("engine")
    global _solver_ = str_to_struct(solver, "Solver")
    println("Quantum impurity solver: $solver -> $(typeof(_solver_))")

    # R3: Setup _adaptor_, the DFT-DMFT adaptor.
    adaptor = get_d("projtype")
    global _adaptor_ = str_to_struct(adaptor, "Adaptor")
    println("Projection scheme: $adaptor -> $(typeof(_adaptor_))")

    # R4: Check the input files
    query_inps(_engine_)
    println("Input files: ready")

    # R5: Prepare the working directories
    mkdir && begin
        build_trees()
        println("Build directories: ready")
    end

    # Print the footer
    println()
    flush(stdout)
end

"""
    go()

Dispatcher for DFT + DMFT calculations. Note that it can not call the
`try_dft()`-`try_mixer()` series functions. These functions are designed
only for testing purpose.

See also: [`ready`](@ref).
"""
function go()
    # Get calculation mode
    mode = get_m("mode")

    # Choose suitable computational driver
    # One-shot DFT + DMFT calculations
    mode == 1 && cycle1()

    # Fully self-consistent DFT + DMFT calculations
    mode == 2 && cycle2()
end

"""
    final()

Finalize the DFT + DMFT calculations.

See also: [`go`](@ref).
"""
function final()
    prompt("Thank you for using $(__LIBNAME__). See you later.")
end

#=
### *Driver Functions* : *Layer 2*
=#

#=
*Some Explanations for the DFT + DMFT Algorithm*

*Remarks 1* :

If the DFT backend is the `vasp` code, we would like to carry out two
successive DFT runs. The purpose of the first DFT run is to evaluate
the fermi level. Then an energy window is determined. We will use this
window to generate optimal projectors in the second DFT run.

If the DFT backend is the `quantum espresso` (`pwscf`) code, we also
carry out two successive DFT runs. The first runs is used to generate
converged charge density, while the second runs is to produce accurate
Kohn-Sham states at given 𝑘-mesh.

*Remarks 2* :

We want better *optimal projectors* for the `vasp` code.

In the first DFT run, `initial` fermi level = 0 -> `wrong` energy
window -> `wrong` optimial projectors. But at this point, the fermi
level is updated, so we have to generate the optimal projectors
again within this new window by doing addition DFT calculation.

*Remarks 3* :

The key Kohn-Sham data inclue lattice structures, 𝑘-mesh and its weights,
tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. At
first, the adaptor will read in these data from the output files of DFT
engine. And then it will process the raw projectors (such as parsing,
labeling, grouping, filtering, and rotatation). Finally, the adaptor will
write down the processed data to some specified files using the internal
`IR` format.

*Remarks 4* :

Once everything is ready, we are going to solve the DMFT self-consistent
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
    prompt("Initialization")
    it.sc = 0 # In preparation mode

    # C01: Perform DFT calculation
    @time_call dft_core(it, lr)

    # C02: Perform DFT calculation again (for the vasp code only)
    is_vasp() && @time_call dft_core(it, lr)

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    adaptor_core(it, lr, ai)

    # C04: Prepare default self-energy functions
    @time_call sigma_core(it, lr, ai, "reset")

#
# DFT + DMFT Iterations (C05-C12)
#
    prompt("Iterations")
    it.sc = 1 # In one-shot DFT + DMFT mode

    # Print the cycle info
    show_it(it, lr)

    for iter = 1:it.M₃
        # Print the log
        prompt("Cycle $iter")
        prompt(lr.log, "")
        prompt(lr.log, "< dft_dmft_cycle >")

        # Update configuration parameters
        refresh(ai)

        # Update IterInfo struct
        incr_it(it)

        # C05: Tackle with the double counting term
        @time_call sigma_core(it, lr, ai, "dcount")

        # C06: Perform DMFT calculation with `task` = 1
        @time_call dmft_core(it, lr, 1)

        # C07: Mix the hybridization functions
        @time_call mixer_core(it, lr, ai, "Δ")

        # C08: Mix the local impurity levels
        @time_call mixer_core(it, lr, ai, "E")

        # C09: Split and distribute the hybridization functions
        @time_call sigma_core(it, lr, ai, "split")

        # C10: Solve the quantum impurity problems
        @time_call solver_core(it, lr, ai)

        # C11: Gather and combine the impurity self-functions
        @time_call sigma_core(it, lr, ai, "gather")

        # C12: Mix the impurity self-energy functions
        @time_call mixer_core(it, lr, ai, "Σ")

        # C18: Check the convergence for total energy
        energy_core(it)

        # Print the cycle info
        show_it(it, lr)

        # If the convergence has been achieved, then break the cycle.
        conv_it(it) && break
    end

    # C97: Kill the DFT engine if it is still alive
    suicide(it)

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
is achieved at both DFT and DMFT levels. So far this function only supports
the `vasp` + `plo` mode.

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
    prompt("Initialization")
    it.sc = 0 # In preparation mode

    # C01: Perform DFT calculation
    @time_call dft_core(it, lr)

    # C02: Perform DFT calculation again (for the vasp code only)
    is_vasp() && @time_call dft_core(it, lr)

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    adaptor_core(it, lr, ai)

    # C04: Prepare default self-energy functions
    @time_call sigma_core(it, lr, ai, "reset")

#
# DFT + DMFT Iterations (C05-C12)
#
    prompt("Iterations")
    it.sc = 2 # In self-consistent DFT + DMFT mode

    # Print the cycle info
    show_it(it, lr)

    # C05: Launch the self-consistent engine
    dft_core(it, lr)

    # Wait the DFT engine to finish its job and sleep
    suspend(4) # Apply a larger time interval

    # Outer: DFT + DMFT LOOP
    for iter = 1:it.M₃

        # Print the log
        prompt("Cycle $iter")
        prompt(lr.log, "")
        prompt(lr.log, "< dft_dmft_cycle >")

        # Update configuration parameters
        refresh(ai)

        # Update IterInfo struct, fix it.I₃.
        incr_it(it, 3, iter)

        # Inner: Adaptor BLOCK
        # Try to extract Kohn-Sham dataset
        begin

            # C06: Apply the adaptor to extract new Kohn-Sham dataset
            adaptor_core(it, lr, ai)

        end

        # Inner: DMFT₁ BLOCK
        # Try to solve the quantum impurity problems
        for iter1 = 1:it.M₁

            # Update IterInfo struct, fix it.I₁.
            incr_it(it, 1, iter1)

            # C07: Tackle with the double counting term
            @time_call sigma_core(it, lr, ai, "dcount")

            # C08: Perform DMFT calculation with `task` = 1
            @time_call dmft_core(it, lr, 1)

            # C09: Mix the hybridization functions
            @time_call mixer_core(it, lr, ai, "Δ")

            # C10: Mix the local impurity levels
            @time_call mixer_core(it, lr, ai, "E")

            # C11: Split and distribute the hybridization functions
            @time_call sigma_core(it, lr, ai, "split")

            # C12: Solve the quantum impurity problems
            @time_call solver_core(it, lr, ai)

            # C13: Gather and combine the impurity self-functions
            @time_call sigma_core(it, lr, ai, "gather")

            # C14: Mix the impurity self-energy functions
            @time_call mixer_core(it, lr, ai, "Σ")

            # Print the cycle info
            show_it(it, lr)
            show_it("dmft1", iter1, it.M₁)

        end # END OF ITER1 LOOP

        # Inner: DMFT₂ BLOCK
        # Try to generate update for density matrix
        begin

            # Update IterInfo struct, fix it.I₂.
            incr_it(it, 2, 1)

            # C15: Perform DMFT calculation with `task` = 2
            @time_call dmft_core(it, lr, 2) # Generate correction for density matrix

            # C16: Mix the correction for density matrix
            @time_call mixer_core(it, lr, ai, "Γ")

            # Print the cycle info
            show_it("dmft2", 1, 1)

        end

        # Inner: DFT BLOCK
        # Try DFT engine with a fixed charge density update
        for iter2 = 1:it.M₂

            # Update IterInfo struct, fix it.I₂.
            incr_it(it, 2, iter2)

            # C17: Reactivate the DFT engine
            dft_core(it, lr, true)

            # Print the cycle info
            show_it(it, lr)
            show_it("dft", iter2, it.M₂)

        end # END OF ITER2 LOOP
        #
        # Reset the counter in IterInfo: I₁, I₂.
        zero_it(it)

        # C18: Check the convergence for total energy
        energy_core(it)

        # If the convergence has been achieved, then break the cycle.
        conv_it(it) && break
    end # END OF ITER LOOP

    # C97: Kill the DFT engine if it is still alive
    suicide(it)

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
    try_dft()

Perform DFT calculations only. If there are something wrong, then you
have chance to adjust the DFT input files manually (for example, you
can modify `vaspc_incar()` in `vasp.jl` by yourself).

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_dft()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C01: Perform DFT calculation
    @time_call dft_core(it, lr)

    # C02: Perform DFT calculation again (for the vasp code only)
    is_vasp() && @time_call dft_core(it, lr)

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
    try_dmft(task::I64)

Perform DMFT calculations only. The users can execute it in the REPL mode
to see whether the DMFT engine works properly. If `task = 1`, it means to
generate a new hybridization function for quantum impurity solver. On the
other hand, it `task = 2`, it means to generate a DMFT correction for the
density matrix, which will be used by the DFT engine.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_dmft(task::I64)
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C01: Execuate the DMFT engine
    @time_call dmft_core(it, lr, task)

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
    try_solver()

Perform calculations using quantum impurity solvers only. The users can
execute it in the REPL mode to see whether the quantum impurity solvers
work properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_solver()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execuate the quantum impurity solvers
    @time_call solver_core(it, lr, ai)

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
    try_adaptor()

Perform calculations using Kohn-Sham adaptor only. The users can execute
it in the REPL mode to see whether the Kohn-Sham adaptor works properly.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_adaptor()
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execute the Kohn-Sham adaptor
    adaptor_core(it, lr, ai)

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
    try_sigma(task::String = "reset")

Perform calculations using self-energy engine only. The users can execute
it in the REPL mode to see whether the self-energy engine works properly.
The argument `task` can be `reset`, `dcount`, `gather`, and `split`.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_sigma(task::String = "reset")
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Execute the Kohn-Sham adaptor
    @time_call sigma_core(it, lr, ai, task)

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
    try_mixer(task::String = "Σ")

Perform calculations using mixer engine only. The users can execute
it in the REPL mode to see whether the mixer engine works properly.
The argument `task` can be `Σ`, `Δ`, `E`, and `Γ`.

In order to run this function correctly, users should try to modify
the predefined parameters in step `C01`.

See also: [`cycle1`](@ref), [`cycle2`](@ref).
"""
function try_mixer(task::String = "Σ")
    # C-2: Create IterInfo struct
    it = IterInfo()

    # C-1: Create Logger struct
    lr = Logger(query_case())

    # C00: Initialize the quantum impurity problems
    ai = GetImpurity()

    # C01: Further setup the IterInfo struct
    #
    # Please modify the I₃ and I₁ parameters to fit your requirements
    it.I₃ = 2; it.I₁ = 0; it.I₂ = 1; it.sc = 2 # Test mixer_gcorr()

    # C02: Execute the Kohn-Sham adaptor
    @time_call mixer_core(it, lr, ai, task)

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

#=
### *Service Functions* : *Layer 1*
=#

"""
    refresh(ai::Array{Impurity,1})

Update the configuration parameters dynamically during the
DFT + DMFT iterations.

See also: [`renew`](@ref), [`FixImpurity`](@ref).
"""
function refresh(ai::Array{Impurity,1})
    # Renew internal dicts which encapsulate configuration parameters
    renew()

    # Renew quantum impurity problems
    FixImpurity(ai)
end

#=
*Remarks* :

In order to terminate the `Zen` code, the following two conditions
should be fulfilled at the same time.

* The argument `force_exit` is true.

* The `case.stop` file exists (from `query_stop()` function).

We usually use this functon to stop the whole DFT + DMFT iterations.
=#

"""
    monitor(force_exit::Bool = false)

Determine whether we need to terminate the `Zen` code.

See also: [`query_stop`](@ref).
"""
function monitor(force_exit::Bool = false)
    if force_exit && query_stop()
        exit(-1)
    end
end

"""
    suspend(second::I64 = 1)

Suspend the current process to wait the DFT engine. This function is
useful for charge fully self-consistent DFT + DMFT calculations.

Now this function only supports the `vasp` code. We have to improve it
to support more DFT engines.

See also: [`dft_core`](@ref).
"""
function suspend(second::I64 = 1)
    # Check second
    if second ≤ 0
        second = 5
    end

    # Sleep some seconds at first
    sleep(second)

    # Enter an infinite loop until some conditions are fulfilled.
    engine = get_d("engine")
    print("Waiting for $engine. Ticking")
    while true
        # Sleep
        sleep(second)
        #
        # Print some hints
        print(".")
        #
        # Check the stop condifion.
        # Here, we check the vasp.lock file. If it is absent, then we
        # break this loop
        !vaspq_lock() && break
    end
    println("\n")
end

"""
    suicide(it::IterInfo)

Kill the DFT engine abnormally. Now it supports the `vasp` code only.

See also: [`dft_core`](@ref).
"""
function suicide(it::IterInfo)
    # Stop it! Only for self-consistent DFT + DMFT iterations.
    if it.sc == 2
        engine = get_d("engine")
        #
        println("Try to kill the $engine app. Please waiting...")
        dft_stop(_engine_)
    end

    # Print the footer
    if ( it.sc == 1 && it.I₁ < it.M₃ ) || ( it.sc == 2 && it.I₃ < it.M₃ )
        println("Good news. The self-consistent iteration is converged.\n")
    else
        println("Sorry, maximum number of iterations is reached.\n")
    end
end

#=
### *Service Functions* : *Layer 2*
=#

"""
    dft_core(it::IterInfo, lr::Logger, sc::Bool = false)

Simple driver for DFT engine. It performs three tasks: (1) Examine
the runtime environment for the DFT engine. (2) Launch the DFT engine.
(3) Backup the output files by DFT engine for next iterations.

If `sc = true`, this function will read in the correction for density
matrix, and then feed it back to the DFT engine to continue the DFT + DMFT
calculations.

Now only the `vasp` and `quantum espresso` codes are supported. If you
want to support more DFT engines, this function must be adapted.

See also: [`adaptor_core`](@ref), [`dmft_core`](@ref), [`solver_core`](@ref).
"""
function dft_core(it::IterInfo, lr::Logger, sc::Bool = false)
    # Determine the chosen engine
    engine = get_d("engine")

    # Print the log
    prompt("DFT", cntr_it(it))
    prompt(lr.log, engine)

    # Enter dft directory
    cd("dft")

    # Activate the chosen DFT engine
    if !sc
        dft_call(_engine_, it)
    else
        # Save DFT energy for the current DFT + DMFT iteration
        it.et.dft = dft_resume(_engine_)
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    dmft_core(it::IterInfo, lr::Logger, task::I64)

Simple driver for DMFT engine. It performs three tasks: (1) Examine
the runtime environment for the DMFT engine. (2) Launch the DMFT engine.
(3) Backup the output files by DMFT engine for next iterations.

The argument `task` is used to specify running mode of the DMFT code.
Its value can be 1 or 2.

See also: [`adaptor_core`](@ref), [`dft_core`](@ref), [`solver_core`](@ref).
"""
function dmft_core(it::IterInfo, lr::Logger, task::I64)
    # Examine the argument `task`
    @assert task == 1 || task == 2

    # Print the log
    prompt("DMFT", cntr_it(it))
    prompt(lr.log, "dmft$task")

    # Enter dmft1 or dmft2 directory
    cd("dmft$task")

    # Activate the chosen DMFT engine
    # Solve the DMFT self-consistent equation
    # Generate DMFT correction for DFT charge density
    begin
        dmft_init(it, task)
        dmft_exec(it, task)
        dmft_save(it, task)
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    solver_core(it::IterInfo,
                lr::Logger,
                ai::Array{Impurity,1},
                force::Bool = false)

Simple driver for quantum impurity solvers. It performs three tasks: (1)
Examine the runtime environment for quantum impurity solver. (2) Launch
the quantum impurity solver. (3) Backup output files by quantum impurity
solver for next iterations.

If `force = true`, then we will try to solve all of the quantum impurity
problems explicitly, irrespective of their equivalences.

Now only the `ct_hyb1`, `ct_hyb2`, `hub1`, and `norg` quantum impurity
solvers are supported. If you want to support the other quantum impurity
solvers, this function must be adapted.

See also: [`adaptor_core`](@ref), [`dft_core`](@ref), [`dmft_core`](@ref).
"""
function solver_core(it::IterInfo,
                     lr::Logger,
                     ai::Array{Impurity,1},
                     force::Bool = false)
    # Sanity check
    @assert length(ai) == get_i("nsite")

    # Analyze the symmetry of quantum impurity problems
    println(blue("Analyze the quantum impurity problems..."))
    #
    # Print number of impurities
    println("  > Number of quantum impurity problems: ", length(ai))
    #
    # Determine the equivalence of quantum impurity problems
    equiv = abs.(get_i("equiv"))
    println("  > Equivalence of quantum impurity problems (abs): ", equiv)
    unique!(equiv)
    println("  > Equivalence of quantum impurity problems (uniq): ", equiv)
    #
    # Figure out which quantum impurity problem should be solved
    to_be_solved = fill(false, length(ai))
    for i in eachindex(equiv)
        ind = findfirst(x -> x.equiv == equiv[i], ai)
        isa(ind, Nothing) && continue
        to_be_solved[ind] = true
    end
    println("  > Quantum impurity problems (keep): ", findall(to_be_solved))
    println("  > Quantum impurity problems (skip): ", findall(.!to_be_solved))
    #
    # Reset DMFT energy
    it.et.dmft = 0.0
    #
    println(blue("Now we are ready to solve them..."))

    # Loop over each impurity site
    for i = 1:get_i("nsite")
        # Extract and show the Impurity strcut
        imp = ai[i]
        CatImpurity(imp)

        # The present quantum impurity problem need to be solved
        if to_be_solved[i] || force

            # Print the header
            println(green("It is interesting. Let us play with it.\n"))

            # Determine the chosen solver
            engine = get_s("engine")

            # Print the log
            prompt("Solvers", cntr_it(it))
            prompt(lr.log, engine)

            # Enter impurity.i directory
            cd("impurity.$i")

            # Activate the chosen quantum impurity solver
            solver_call(_solver_, it, imp)

            # Enter the parent directory
            cd("..")

        else

            # Print the header
            println(red("Mmm, it is not my job..."))

            # Well, the current quantum impurity problem is not solved.
            # We have to find out its sister which has been solved before.
            found = -1
            for j = 1:i
                # Impurity 𝑖 and Impurity 𝑗 are related by some kinds of
                # symmetry. Impurity 𝑗 has been solved before. So we can
                # copy its solution to Impurity 𝑖.
                if abs(imp.equiv) == abs(ai[j].equiv) && to_be_solved[j]
                    found = j
                end
            end
            # Sanity check
            @assert found > 0
            println(green("Maybe we can learn sth. from Impurity [$found]\n"))

            # Enter impurity.i directory
            cd("impurity.$i")

            # Next, we would like to copy solution from Impurity `found`
            # to the current Impurity 𝑖.
            solver_copy(_solver_, it, ai[found], imp)

            # Enter the parent directory
            cd("..")

        end

        # Well, now we would like to extract the DMFT energy.
        edmft = GetEdmft(imp)
        it.et.dmft = it.et.dmft + edmft
        println("  > DMFT interaction energy: [$i] -> $edmft eV")
    end # END OF I LOOP

    # Monitor the status
    monitor(true)
end

"""
    adaptor_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1})

Simple driver for the adaptor. It performs three tasks: (1) Initialize
the adaptor, to check whether the essential files exist. (2) Parse the
Kohn-Sham data output by the DFT engine, try to preprocess them, and
then transform them into IR format. (3) Backup the files by adaptor.

For the first task, both the `vasp` and `qe` adaptors are supported.
While for the second task, both the `plo` and `wannier` adaptors are
supported. If you want to support more adaptors, please adapt this
function by yourself.

See also: [`dft_core`](@ref), [`dmft_core`](@ref), [`solver_core`](@ref).
"""
function adaptor_core(it::IterInfo, lr::Logger, ai::Array{Impurity,1})
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
    engine = nameof(_engine_)
    prompt("Adaptor", cntr_it(it))
    prompt(lr.log, "adaptor::$engine")
    @time_call adaptor_call(_engine_, DFTData)

    #
    # A2: Process the original Kohn-Sham data
    #
    # Well, now we have the Kohn-Sham data. But they can not be used
    # directly. We have to check and process them carefully. Please
    # pay attention to that the DFTData dict will be modified in
    # the `plo_adaptor() / wannier_adaptor()` function.
    #
    # The `plo_adaptor() / wannier_adaptor()` function also have the
    # ability to calculate some selected physical quantities (such as
    # overlap matrix and density of states) to check the correctness
    # of the Kohn-Sham data. These features should be automatically
    # activated if you are in the REPL mode and there is a `case.test`
    # file in the present directory (i.e, the `dft` folder).
    #
    projtype = nameof(_adaptor_)
    prompt("Adaptor", cntr_it(it))
    prompt(lr.log, "adaptor::$projtype")
    @time_call adaptor_call(_adaptor_, DFTData, ai)

    #
    # A3: Output the processed Kohn-Sham data
    #
    # Ok, now the Kohn-Sham data are ready. We would like to write them
    # to some specified files with the IR format. Then these files will
    # be saved immediately.
    #
    prompt("Adaptor", cntr_it(it))
    prompt(lr.log, "adaptor::ir")
    @time_call ir_adaptor(DFTData)
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
    sigma_core(it::IterInfo,
               lr::Logger,
               ai::Array{Impurity,1},
               task::String = "reset")

Simple driver for functions for processing the self-energy functions
Σ and hybridization functions Δ (and local impurity levels ϵ).

Now it supports four tasks: `reset`, `dcount`, `split`, `gather`. It
won't change the current directory.

See also: [`mixer_core`](@ref).
"""
function sigma_core(it::IterInfo,
                    lr::Logger,
                    ai::Array{Impurity,1},
                    task::String = "reset")
    # Check the given task
    @assert task in ("reset", "dcount", "split", "gather")

    # Print the log
    prompt("Sigma", cntr_it(it))
    prompt(lr.log, "sigma::$task")

    # Launch suitable subroutine
    sigma_call(str_to_struct(task, "Mode"), it, ai)

    # Monitor the status
    monitor(true)
end

"""
    mixer_core(it::IterInfo,
               lr::Logger,
               ai::Array{Impurity,1},
               task::String = "Σ")

Simple driver for the mixer. It will try to mix the self-energy functions
Σ (or hybridization functions Δ, local impurity levels ϵ, density matrix
Γ) and generate a new one.

Now it supports four tasks: `Σ`, `Δ`, `E`, `Γ`. It won't change the
current directory.

See also: [`sigma_core`](@ref).
"""
function mixer_core(it::IterInfo,
                    lr::Logger,
                    ai::Array{Impurity,1},
                    task::String = "Σ")
    # Check the given task
    @assert task in ("Σ", "Δ", "E", "Γ")

    # Print the log
    prompt("Mixer", cntr_it(it))
    prompt(lr.log, "mixer::$task")

    # Check iteration number to see whether we have enough data to be mixed
    if it.sc == 1
        if it.I₁ ≤ 1
            return
        end
    else
        if task in ("Σ", "Δ", "E")
            if it.I₃ == 1 && it.I₁ == 1
                return
            end
        else
            @assert task == "Γ"
            if it.I₃ == 1 && it.I₂ == 1
                return
            end
        end
    end

    # Launch suitable subroutine
    mixer_call(str_to_struct(task, "Mixer"), it, ai)

    # Monitor the status
    monitor(true)
end

"""
    energy_core(it::IterInfo)

Simple driver for treating the total DFT + DMFT energy. It will print
the decomposition of total energy, and try to calculate the energy
difference between two successive DFT + DMFT iterations.

See also: [`Energy`](@ref), [`IterInfo`](@ref).
"""
function energy_core(it::IterInfo)
    # Try to print the header
    if it.sc == 1
        println("The DFT + DMFT Energy At Cycle [$(it.I₁) / $(it.M₃)]")
    else
        println("The DFT + DMFT Energy At Cycle [$(it.I₃) / $(it.M₃)]")
    end
    #
    println(repeat("==", 36))

    # For the first iteration
    if it.I₃ == 1
        println("  > E[DFT]   : $(it.et.dft) eV")
        println("  > E[DMFT]  : $(it.et.dmft) eV")
        println("  > E[CORR]  : $(it.et.corr) eV")
        println("  > E[DC]    : $(it.et.dc) eV")
        println("  > E[TOTAL] : $(it.et.total) eV")
    #
    # For the second and later iterations
    else
        # Calculate error bar
        err_dft = abs((it.et.dft - it.ep.dft) / it.et.dft) * 100
        err_dmft = abs((it.et.dmft - it.ep.dmft) / it.et.dmft) * 100
        err_corr = abs((it.et.corr - it.ep.corr) / it.et.corr) * 100
        err_dc = abs((it.et.dc - it.ep.dc) / it.et.dc) * 100
        err_total = abs((it.et.total - it.ep.total) / it.et.total) * 100
        #
        # Print energy and error bar
        println("  > E[DFT]   : $(it.et.dft) eV (err: $err_dft %)")
        println("  > E[DMFT]  : $(it.et.dmft) eV (err: $err_dmft %)")
        println("  > E[CORR]  : $(it.et.corr) eV (err: $err_corr %)")
        println("  > E[DC]    : $(it.et.dc) eV (err: $err_dc %)")
        println("  > E[TOTAL] : $(it.et.total) eV (err: $err_total %)")
        #
        # Calculate and show the difference
        dist = abs(it.et.total - it.ep.total)
        it.ce = ( dist < get_m("ec") )
        println("  > Calculated ΔE(TOTAL) = $dist ( convergence is $(it.ce) )")
    end

    # Try to print the footer
    println(repeat("==", 36), "\n")

    # Update it.ep with it.et
    it.ep = deepcopy(it.et)
end

#=
### *Service Functions* : *Layer 3*
=#

#=
*Remarks* :

The working directories include `dft`, `dmft1`, `dmft2`, and `impurity.i`.
If they exist already, it would be better to remove them at first.
=#

"""
    build_trees()

Prepare the working directories at advance.

See also: [`clear_trees`](@ref).
"""
function build_trees()
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
    clear_trees()

Remove the working directories finally.

See also: [`build_trees`](@ref).
"""
function clear_trees()
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

#=
*Service Functions*: *Layer 4* (*For IterInfo Struct*)
=#

"""
    incr_it(it::IterInfo)

Modify the internal counters in `IterInfo` struct. This function is
used in the one-shot DFT + DMFT calculations only.

See also: [`IterInfo`](@ref), [`zero_it`](@ref).
"""
function incr_it(it::IterInfo)
    @assert it.sc == 1
    #
    it.I₃ = 1
    it.I₁ = it.I₁ + 1
    it.I₄ = it.I₄ + 1
    #
    @assert it.I₁ ≤ it.M₃
end

"""
    incr_it(it::IterInfo, c::I64, v::I64)

Modify the internal counters in `IterInfo` struct. This function is
used in the fully charge self-consistent DFT + DMFT calculations only.
The argument `c` denotes the counter `it.I`, and `v` is the value.

See also: [`IterInfo`](@ref), [`zero_it`](@ref).
"""
function incr_it(it::IterInfo, c::I64, v::I64)
    @assert it.sc == 2
    @assert c in (1, 2, 3)
    @assert v ≥ 1

    if c == 1
        @assert v ≤ it.M₁
        it.I₁ = v
        it.I₄ = it.I₄ + 1
    elseif c == 2
        @assert v ≤ it.M₂
        it.I₂ = v
        it.I₄ = it.I₄ + 1
    elseif c == 3
        @assert v ≤ it.M₃
        it.I₃ = v
    end
end

"""
    zero_it(it::IterInfo)

Reset some (not all) counters in the `IterInfo` struct.

See also: [`IterInfo`](@ref), [`incr_it`](@ref).
"""
function zero_it(it::IterInfo)
    it.I₁ = 0
    it.I₂ = 0
end

"""
    prev_it(it::IterInfo)

Return the iteration information for previous DFT + DMFT step. This
function is suitable for one-shot calculation mode.

See also: [`mixer_core`](@ref), [`incr_it`](@ref).
"""
function prev_it(it::IterInfo)
    @assert it.sc == 1
    @assert it.I₁ ≥ 2
    return it.I₃, it.I₁ - 1
end

"""
    prev_it(it::IterInfo, c::I64)

Return the iteration information for previous DFT + DMFT step. This
function is suitable for fully self-consistent calculation mode. Here
the argument `c` denotes the counter `it.I`

See also: [`mixer_core`](@ref), [`incr_it`](@ref).
"""
function prev_it(it::IterInfo, c::I64)
    @assert it.sc == 2
    @assert c in (1, 2)

    if c == 1
        list = [(i3,i1) for i1 = 1:it.M₁, i3 = 1:it.M₃]
        newlist = reshape(list, it.M₁ * it.M₃)
        ind = findfirst(x -> x == (it.I₃, it.I₁), newlist)
        @assert ind ≥ 2
        return newlist[ind - 1]
    else
        # Special treatment for vasp code
        trueM₂ = ( is_vasp() ? 1 : it.M₂ )
        list = [(i3,i2) for i2 = 1:trueM₂, i3 = 1:it.M₃]
        newlist = reshape(list, trueM₂ * it.M₃)
        ind = findfirst(x -> x == (it.I₃, it.I₂), newlist)
        @assert ind ≥ 2
        return newlist[ind - 1]
    end
end

"""
    cntr_it(it::IterInfo)

Return the counters in the IterInfo struct as a format string.

See also: [`IterInfo`](@ref)
"""
function cntr_it(it::IterInfo)
    " [Cycle -> $(it.sc) : $(it.I₄) : $(it.I₃) : $(it.I₁) : $(it.I₂)] "
end

"""
    show_it(it::IterInfo, lr::Logger)

Try to record the iteration information in the `case.cycle` file.

See also: [`IterInfo`](@ref), [`Logger`](@ref).
"""
function show_it(it::IterInfo, lr::Logger)
    # Extract parameter `nsite`
    nsite = get_i("nsite")
    @assert nsite == length(it.nf)

    # Write the header
    if it.I₄ == 0
        # Write labels
        print(lr.cycle, "#    #    #    #    μ0          μ1          μ2          ")
        for t = 1:nsite
            print(lr.cycle, "Vdc$t        ")
        end
        print(lr.cycle, "Nd1         Nd2         ")
        for t = 1:nsite
            print(lr.cycle, "Ni[$t]       ")
        end
        print(lr.cycle, "Etot        ")
        println(lr.cycle, "C(C)    C(E)    C(S)")
        # Write separator
        println(lr.cycle, repeat('-', 3*8 + 4*5 + 6*12 + 24*nsite))
    # Write iteration information
    else
        @printf(lr.cycle, "%-5i", it.I₄)
        @printf(lr.cycle, "%-5i", it.I₃)
        @printf(lr.cycle, "%-5i", it.I₁)
        @printf(lr.cycle, "%-5i", it.I₂)
        #
        if it.μ₀ < 0.0
            @printf(lr.cycle, "%-12.7f", it.μ₀)
        else
            @printf(lr.cycle, "+%-11.7f", it.μ₀)
        end
        if it.μ₁ < 0.0
            @printf(lr.cycle, "%-12.7f", it.μ₁)
        else
            @printf(lr.cycle, "+%-11.7f", it.μ₁)
        end
        if it.μ₂ < 0.0
            @printf(lr.cycle, "%-12.7f", it.μ₂)
        else
            @printf(lr.cycle, "+%-11.7f", it.μ₂)
        end
        #
        for t = 1:nsite
            @printf(lr.cycle, "%-12.7f", it.dc[t])
        end
        #
        @printf(lr.cycle, "%-12.7f", it.n₁)
        @printf(lr.cycle, "%-12.7f", it.n₂)
        for t = 1:nsite
            @printf(lr.cycle, "%-12.7f", it.nf[t])
        end
        #
        @printf(lr.cycle, "%-12.7f", it.et.total)
        #
        @printf(lr.cycle, "%-8s", it.cc ? "true" : "false")
        @printf(lr.cycle, "%-8s", it.ce ? "true" : "false")
        @printf(lr.cycle, "%-8s", it.cs ? "true" : "false")
        #
        println(lr.cycle)
    end

    # Flush the IOStream
    flush(lr.cycle)
end

"""
    show_it(mode::String, iter::I64, max_iter::I64)

Try to record the iteration information in the terminal.
"""
function show_it(mode::String, iter::I64, max_iter::I64)
    @assert mode in ("dmft1", "dmft2", "dft")
    @assert iter ≥ 1
    @assert max_iter ≥ iter
    #
    print("Mode : [ $mode ], ")
    print("Requested Iteration : [ $max_iter ], ")
    println("Finished Iteration: [ $iter ]. \n")
end

"""
    conv_it(it::IterInfo)

Check whether the convergence flags are achieved.

See also: [`IterInfo`](@ref).
"""
function conv_it(it::IterInfo)
    if it.sc == 1
        conv = it.cs
    else
        conv = it.cc && it.ce && it.cs
    end
    conv
end
