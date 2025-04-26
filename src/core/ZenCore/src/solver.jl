#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/26
#

#=
### *Multiple Dispatchers*
=#

"""
    solver_call(::NULLSolver, it::IterInfo, imp::Impurity)
    solver_call(::CTHYB₁Solver, it::IterInfo, imp::Impurity)
    solver_call(::CTHYB₂Solver, it::IterInfo, imp::Impurity)
    solver_call(::HIASolver, it::IterInfo, imp::Impurity)
    solver_call(::NORGSolver, it::IterInfo, imp::Impurity)
    solver_call(::ATOMSolver, it::IterInfo, imp::Impurity)

Try to solve the quantum impurity problems by using various quantum
impurity solvers. It acts as a dispatcher. Now it supports `CTHYB₁`
(`ct_hyb1`), `CTHYB₂` (`ct_hyb2`), `HIA` (`hia`), `NORG` (`norg`),
and `ATOM` (`atomic`) solvers.

See also: [`_solver_`](@ref).
"""
function solver_call(::NULLSolver, it::IterInfo, imp::Impurity)
    sorry()
end
#
function solver_call(::CTHYB₁Solver, it::IterInfo, imp::Impurity)
    # For CTHYB₁ quantum impurity solver
    s_qmc1_init(it, imp)
    s_qmc1_exec(it)
    s_qmc1_save(it, imp)
end
#
function solver_call(::CTHYB₂Solver, it::IterInfo, imp::Impurity)
    # For CTHYB₂ quantum impurity solver
    s_qmc2_init(it)
    s_qmc2_exec(it)
    s_qmc2_save(it)
end
#
function solver_call(::HIASolver, it::IterInfo, imp::Impurity)
    # For HIA quantum impurity solver
    s_hub1_init(it)
    s_hub1_exec(it)
    s_hub1_save(it)
end
#
function solver_call(::NORGSolver, it::IterInfo, imp::Impurity)
    # For NORG quantum impurity solver
    s_norg_init(it, imp)
    s_norg_exec(it)
    s_norg_save(it, imp)
end
#
function solver_call(::ATOMSolver, it::IterInfo, imp::Impurity)
    # For atomic eigenvalues solver
    sorry()
end

"""
    solver_copy(::NULLSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    solver_copy(::CTHYB₁Solver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    solver_copy(::CTHYB₂Solver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    solver_copy(::HIASolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    solver_copy(::NORGSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    solver_copy(::ATOMSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Try to solve a quantum impurity problem by copying solution from another
equivalent quantum impurity problem. It acts as a dispatcher.

See also: [`_solver_`](@ref).
"""
solver_copy(::NULLSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    sorry()
#
solver_copy(::CTHYB₁Solver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    # For CTHYB₁ quantum impurity solver
    s_qmc1_copy(it, imp₁, imp₂)
#
solver_copy(::CTHYB₂Solver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    # For CTHYB₂ quantum impurity solver
    s_qmc2_copy(it, imp₁, imp₂)
#
solver_copy(::HIASolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    # For HIA quantum impurity solver
    s_hub1_copy(it, imp₁, imp₂)
#
solver_copy(::NORGSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    # For NORG quantum impurity solver
    s_norg_copy(it, imp₁, imp₂)
#
solver_copy(::ATOMSolver, it::IterInfo, imp₁::Impurity, imp₂::Impurity) =
    # For atomic eigenvalues solver
    sorry()

"""
    solver_sigma(::NULLSolver, imp::Impurity)
    solver_sigma(::CTHYB₁Solver, imp::Impurity)
    solver_sigma(::CTHYB₂Solver, imp::Impurity)
    solver_sigma(::HIASolver, imp::Impurity)
    solver_sigma(::NORGSolver, imp::Impurity)
    solver_sigma(::ATOMSolver, imp::Impurity)

Try to extract self-energy function from the output data of quantum
impurity solver. It acts as a dispatcher.

See also: [`_solver_`](@ref).
"""
solver_sigma(::NULLSolver, imp::Impurity) = sorry()
solver_sigma(::CTHYB₁Solver, imp::Impurity) = ctqmc_sigma(imp)
solver_sigma(::CTHYB₂Solver, imp::Impurity) = ctqmc_sigma(imp)
solver_sigma(::HIASolver, imp::Impurity) = sorry()
solver_sigma(::NORGSolver, imp::Impurity) = norg_sigma(imp)
solver_sigma(::ATOMSolver, imp::Impurity) = sorry()

"""
    solver_nimpx(::NULLSolver, imp::Impurity)
    solver_nimpx(::CTHYB₁Solver, imp::Impurity)
    solver_nimpx(::CTHYB₂Solver, imp::Impurity)
    solver_nimpx(::HIASolver, imp::Impurity)
    solver_nimpx(::NORGSolver, imp::Impurity)
    solver_nimpx(::ATOMSolver, imp::Impurity)

Try to extract impurity occupancy from the output data of quantum
impurity solver. It acts as a dispatcher.

See also: [`_solver_`](@ref).
"""
solver_nimpx(::NULLSolver, imp::Impurity) = sorry()
solver_nimpx(::CTHYB₁Solver, imp::Impurity) = ctqmc_nimpx(imp)
solver_nimpx(::CTHYB₂Solver, imp::Impurity) = ctqmc_nimpx(imp)
solver_nimpx(::HIASolver, imp::Impurity) = sorry()
solver_nimpx(::NORGSolver, imp::Impurity) = norg_nimpx(imp)
solver_nimpx(::ATOMSolver, imp::Impurity) = sorry()

"""
    solver_edmft(::NULLSolver)
    solver_edmft(::CTHYB₁Solver)
    solver_edmft(::CTHYB₂Solver)
    solver_edmft(::HIASolver)
    solver_edmft(::NORGSolver)
    solver_edmft(::ATOMSolver)

Try to extract interaction energy from the output data of quantum
impurity solver. It acts as a dispatcher.

See also: [`_solver_`](@ref).
"""
solver_edmft(::NULLSolver) = sorry()
solver_edmft(::CTHYB₁Solver) = ctqmc_edmft()
solver_edmft(::CTHYB₂Solver) = ctqmc_edmft()
solver_edmft(::HIASolver) = sorry()
solver_edmft(::NORGSolver) = norg_edmft()
solver_edmft(::ATOMSolver) = sorry()

#=
### *CTHYB₁ Quantum Impurity Solver*
=#

"""
    s_qmc1_init(it::IterInfo, imp::Impurity)

Check runtime environment of the CTHYB₁ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_exec`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_init(it::IterInfo, imp::Impurity)
    # Print the header
    println("Engine : CTHYB₁")
    println("Try to solve the quantum impurity problem: [$(imp.index)]")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for solver")

    # Generate configuration file for quantum impurity solver
    ctqmc_setup(imp)
    println("  > File solver.ctqmc.in is ready")

    # Prepare hybridization functions
    #
    # Extract frequency mesh and hybridization function from `dmft.delta`
    fmesh, Delta = read_delta(imp)
    #
    # Write frequency mesh and hybridization function to `solver.hyb.in`
    ctqmc_delta(fmesh, Delta)
    println("  > File solver.hyb.in is ready")

    # Prepare local impurity levels
    #
    # Extract local impurity levels from `dmft.eimpx`
    Eimpx = read_eimpx(imp)
    #
    # Write local impurity levels to `solver.eimp.in`
    ctqmc_eimpx(Eimpx)
    println("  > File solver.eimp.in is ready")
end

"""
    s_qmc1_exec(it::IterInfo)

Launch the CTHYB₁ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_exec(it::IterInfo)
    # Print the header
    println("Detect the runtime environment for solver")

    # Determine mpi prefix (whether the solver is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "solver", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of quantum impurity solver
    solver_home = query_solver(_solver_)
    println("  > Home directory for solver: ", solver_home)

    # Select suitable solver program
    solver_exe = "$solver_home/ctqmc"
    @assert isfile(solver_exe)
    println("  > Executable program is available: ", basename(solver_exe))

    # Assemble command
    if isnothing(mpi_prefix)
        solver_cmd = solver_exe
    else
        solver_cmd = split("$mpi_prefix $solver_exe", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', solver_cmd))")

    # Print the header
    println("Launch the computational engine (quantum impurity solver)")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$solver_cmd`, stdout = "solver.out"))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to solver.out.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is being executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Analyze the solver.out file during the calculation
    #
    # `c` is a time counter
    c = 0
    #
    # Enter infinite loop
    while true
        # Sleep five seconds
        sleep(5)

        # Increase the counter
        c = c + 1

        # Parse solver.out file
        lines = readlines("solver.out")
        filter!(x -> contains(x, "iter:"), lines)

        # Figure out the number of monte carlo sweeps
        if length(lines) > 0
            arr = line_to_array(lines[end])
            c_sweep = parse(I64, arr[5])
            t_sweep = parse(I64, arr[7])
            R = c_sweep / t_sweep
        else # Nothing
            c_sweep = 0
            R = 0.0
        end

        # Print the log to screen
        @printf("  > Elapsed %4i seconds, current sweeps: %i (%4.2f)\r", 5*c, c_sweep, R)

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the solver task to finish
    wait(t)

    # Extract how many monte carlo sampling blocks are executed
    lines = readlines("solver.out")
    filter!(x -> contains(x, "iter:"), lines)
    println("  > Finished after $(length(lines)) Monte Carlo sampling blocks")

    # Extract perturbation expansion order information
    println("Report From CTHYB₁ Quantum Impurity Solver")
    lines = readlines("solver.out")
    start = findlast(x -> contains(x, ">>> iter:"), lines) + 1
    finish = start + 20
    println("  [")
    foreach(x -> println(x), lines[start:finish])
    println("  ]")
end

"""
    s_qmc1_save(it::IterInfo, imp::Impurity)

Backup output files of the CTHYB₁ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_exec`](@ref).
"""
function s_qmc1_save(it::IterInfo, imp::Impurity)
    # Print the header
    println("Finalize the computational task")

    # Determine which files are important
    #
    # Major output
    fout = ["solver.out"]
    #
    # Green's functions
    fgrn = ["solver.grn.dat", "solver.green.dat"]
    #
    # Hybridization functions
    fhyb = ["solver.hyb.dat", "solver.hybri.dat"]
    #
    # Self-energy functions
    fsgm = ["solver.sgm.dat"]
    #
    # Auxiliary output files
    faux = ["solver.nmat.dat", "solver.paux.dat", "solver.prob.dat", "solver.hist.dat"]

    # Next, we have to backup the above files.
    foreach( x ->
        begin
            file_src = x
            file_dst = "$x.$(it.I₃).$(it.I₁)"
            cp(file_src, file_dst, force = true)
        end,
        union(fout, fgrn, fhyb, fsgm, faux)
    )
    println("  > Save the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    ctqmc_nimpx(imp)
    println("  > Extract the impurity occupancy from solver.nmat.dat: $(imp.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp.index] = imp.occup
end

"""
    s_qmc1_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Duplicate output files of the CTHYB₁ quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
are already in directory `impurity.2`.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_exec`](@ref).
"""
function s_qmc1_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    # Print the header
    println("Transfer results from impurity $(imp₁.index) to $(imp₂.index)")

    # Determine which files are important
    #
    # Major output
    fout = ["solver.out"]
    #
    # Green's functions
    fgrn = ["solver.grn.dat", "solver.green.dat"]
    #
    # Hybridization functions
    fhyb = ["solver.hyb.dat", "solver.hybri.dat"]
    #
    # Self-energy functions
    fsgm = ["solver.sgm.dat"]
    #
    # Auxiliary output files
    faux = ["solver.nmat.dat", "solver.paux.dat", "solver.prob.dat", "solver.hist.dat"]

    # Determine the index for imp₁
    index = imp₁.index

    # Next, we have to backup the above files.
    foreach( x ->
        begin
            file_src = "../impurity.$index/$x"
            file_dst = "$x"
            cp(file_src, file_dst, force = true)
        end,
        union(fout, fgrn, fhyb, fsgm, faux)
    )
    println("  > Copy the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    ctqmc_nimpx(imp₂)
    println("  > Extract the impurity occupancy from solver.nmat.dat: $(imp₂.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp₂.index] = imp₂.occup
end

#=
### *CTHYB₂ Quantum Impurity Solver*
=#

"""
    s_qmc2_init(it::IterInfo)

Check runtime environment of the CTHYB₂ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_exec`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_init(it::IterInfo)
    sorry()
end

"""
    s_qmc2_exec(it::IterInfo)

Launch the CTHYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_exec(it::IterInfo)
    sorry()
end

"""
    s_qmc2_save(it::IterInfo)

Backup output files of the CTHYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_save(it::IterInfo)
    sorry()
end

"""
    s_qmc2_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Duplicate output files of the CTHYB₂ quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
are already in directory `impurity.2`.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    sorry()
end

#=
### *HIA Quantum Impurity Solver*
=#

"""
    s_hub1_init(it::IterInfo)

Check runtime environment of the HIA quantum impurity solver. Prepare
the necessary input files.

See also: [`s_hub1_exec`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_init(it::IterInfo)
    sorry()
end

"""
    s_hub1_exec(it::IterInfo)

Launch the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_exec(it::IterInfo)
    sorry()
end

"""
    s_hub1_save(it::IterInfo)

Backup output files of the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_exec`](@ref).
"""
function s_hub1_save(it::IterInfo)
    sorry()
end

"""
    s_hub1_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Duplicate output files of the HIA quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
are already in directory `impurity.2`.

See also: [`s_hub1_init`](@ref), [`s_hub1_exec`](@ref).
"""
function s_hub1_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    sorry()
end

#=
### *NORG Quantum Impurity Solver*
=#

"""
    s_norg_init(it::IterInfo, imp::Impurity)

Check runtime environment of the NORG quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the RUC Team.

See also: [`s_norg_exec`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_init(it::IterInfo, imp::Impurity)
    # Print the header
    println("Engine : NORG")
    println("Try to solve the quantum impurity problem: [$(imp.index)]")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for solver")

    # Generate configuration file for quantum impurity solver
    norg_setup(imp)
    println("  > File solver.norg.in is ready")

    # Prepare hybridization functions
    #
    # Extract frequency mesh and hybridization function from `dmft.delta`
    fmesh, Delta = read_delta(imp)
    #
    # Write frequency mesh and hybridization function to `solver.hyb.in`
    norg_delta(fmesh, Delta)
    println("  > File solver.hyb.in is ready")

    # Prepare local impurity levels
    #
    # Extract local impurity levels from `dmft.eimpx`
    Eimpx = read_eimpx(imp)
    #
    # Write local impurity levels to `solver.eimp.in`
    norg_eimpx(Eimpx)
    println("  > File solver.eimp.in is ready")
end

"""
    s_norg_exec(it::IterInfo)

Launch the NORG quantum impurity solver.

This quantum impurity solver is from the RUC Team.

See also: [`s_norg_init`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_exec(it::IterInfo)
    # Print the header
    println("Detect the runtime environment for solver")

    # Determine mpi prefix (whether the solver is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "solver", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of quantum impurity solver
    solver_home = query_solver(_solver_)
    println("  > Home directory for solver: ", solver_home)

    # Select suitable solver program
    solver_exe = "$solver_home/norg"
    @assert isfile(solver_exe)
    println("  > Executable program is available: ", basename(solver_exe))

    # Assemble command
    if isnothing(mpi_prefix)
        solver_cmd = solver_exe
    else
        solver_cmd = split("$mpi_prefix $solver_exe", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', solver_cmd))")

    # Print the header
    println("Launch the computational engine (quantum impurity solver)")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$solver_cmd`, stdout = "solver.out"))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to solver.out.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is being executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Analyze the solver.out file during the calculation
    #
    # `c` is a time counter
    c = 0
    #
    # Enter infinite loop
    while true
        # Sleep five seconds
        sleep(5)

        # Increase the counter
        c = c + 1

        # Parse solver.out file
        lines = readlines("solver.out")
        filter!(x -> contains(x, "NORG begin"), lines)

        # Figure out the number of norg runs
        if length(lines) > 0
            nrun = length(lines)
        else # Nothing
            nrun = 0
        end

        # Print the log to screen
        @printf("  > Elapsed %4i seconds, current norg runs: %4i \r", 5*c, nrun)

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the solver task to finish
    wait(t)

    # Extract how many norg blocks are executed
    lines = readlines("solver.out")
    filter!(x -> contains(x, "NORG begin"), lines)
    nrun = length(lines)
    #
    lines = readlines("solver.out")
    filter!(x -> contains(x, "iter_norg_cnt"), lines)
    niter_norg = length(lines)
    println("  > Finished after $nrun norg runs ($niter_norg iterations)")

    # Extract impurity occupation information
    println("Report From NORG Quantum Impurity Solver")
    lines = readlines("solver.out")
    start = findlast(x -> contains(x, "iter_norg_cnt"), lines) + 1
    finish = findlast(x -> contains(x, "groundE_pre"), lines) - 1
    println("  [")
    foreach(x -> println(x), lines[start:finish])
    println("  ]")
end

"""
    s_norg_save(it::IterInfo, imp::Impurity)

Backup output files of the NORG quantum impurity solver.

This quantum impurity solver is from the RUC Team.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_save(it::IterInfo, imp::Impurity)
    # Print the header
    println("Finalize the computational task")

    # Determine which files are important
    #
    # Major output
    fout = ["solver.out"]
    #
    # Green's functions
    fgrn = ["gfimp.txt"]
    #
    # Bath Green's functions
    fhyb = ["g0imp.txt"]
    #
    # Self-energy functions
    fsgm = ["seimp.txt"]
    #
    # Auxiliary output files
    faux = ["hop.txt", "ose.txt", "nmat.txt"]

    # Next, we have to backup the above files.
    foreach( x ->
        begin
            file_src = x
            file_dst = "$x.$(it.I₃).$(it.I₁)"
            cp(file_src, file_dst, force = true)
        end,
        union(fout, fgrn, fhyb, fsgm, faux)
    )
    println("  > Save the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    norg_nimpx(imp)
    println("  > Extract the impurity occupancy from nmat.txt: $(imp.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp.index] = imp.occup
end

"""
    s_norg_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Duplicate output files of the NORG quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
are already in directory `impurity.2`.

This quantum impurity solver is from the RUC Team.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_copy(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    # Print the header
    println("Transfer results from impurity $(imp₁.index) to $(imp₂.index)")

    # Determine which files are important
    #
    # Major output
    fout = ["solver.out"]
    #
    # Green's functions
    fgrn = ["gfimp.txt"]
    #
    # Bath Green's functions
    fhyb = ["g0imp.txt"]
    #
    # Self-energy functions
    fsgm = ["seimp.txt"]
    #
    # Auxiliary output files
    faux = ["hop.txt", "ose.txt", "nmat.txt"]

    # Determine the index for imp₁
    index = imp₁.index

    # Next, we have to backup the above files.
    foreach( x ->
        begin
            file_src = "../impurity.$index/$x"
            file_dst = "$x"
            cp(file_src, file_dst, force = true)
        end,
        union(fout, fgrn, fhyb, fsgm, faux)
    )
    println("  > Copy the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    norg_nimpx(imp₂)
    println("  > Extract the impurity occupancy from nmat.txt: $(imp₂.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp₂.index] = imp₂.occup
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_setup(imp::Impurity)

Generate default configuration file (`solver.ctqmc.in`) for the CTHYB
quantum impurity solvers automatically (according to the information
encoded in the `Impurity` struct).

See also: [`Impurity`](@ref), [`ctqmc_atomx`](@ref).
"""
function ctqmc_setup(imp::Impurity)
    # Filename of configuration file
    fctqmc = "solver.ctqmc.in"

    # Extract parameters from `Impurity` struct
    #
    # mune is fixed to zero, because the chemical potential is adsorbed
    # in the local impurity level.
    nband = imp.nband
    norbs = 2 * nband
    ncfgs = 2^norbs
    Uc    = imp.upara
    Jz    = imp.jpara
    mune  = 0.0
    beta  = imp.beta

    # Get number of Matsubara frequency points
    nmesh = get_m("nmesh")    

    # Create the configuration file
    open(fctqmc, "w") do fout
        # Print some warning message
        println(fout, "# Generated by $__LIBNAME__ at ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
        println(fout, "# PLEASE DO NOT MODIFY IT MANUALLY!")
        println(fout)

        # Print the standard parameters
        println(fout, "# Standard parameters: Running mode")
        println(fout, "isscf = 1")
        println(fout)
        #
        println(fout, "# Standard parameters: System's size")
        println(fout, "nband = $nband")
        println(fout, "norbs = $norbs")
        println(fout, "ncfgs = $ncfgs")
        println(fout)
        #
        println(fout, "# Standard parameters: Interaction")
        println(fout, "Uc    = $Uc")
        println(fout, "Jz    = $Jz")
        println(fout)
        #
        println(fout, "# Standard parameters: Others")
        println(fout, "mune  = $mune")
        println(fout, "beta  = $beta")
        println(fout, "mfreq = $nmesh")
        println(fout)

        # Print the user-supplied parameters
        #
        # The correctness of these auxiliary parameters should be checked
        # in `config.jl/chk_dict()`.
        println(fout, "# Auxiliary parameters: By users")
        foreach(x -> println(fout, x), get_s("params"))
    end # END OF IOSTREAM
end

"""
    ctqmc_atomx(imp::Impurity)

Generate configuration file for the atomic problem solver.

See also: [`Impurity`](@ref), [`ctqmc_setup`](@ref).
"""
function ctqmc_atomx(imp::Impurity)
    sorry()
end

"""
    norg_setup(imp::Impurity)

Generate default configuration file (`solver.norg.in`) for the NORG
quantum impurity solvers automatically (according to the information
encoded in the `Impurity` struct).

See also: [`Impurity`](@ref), [`ctqmc_atomx`](@ref).
"""
function norg_setup(imp::Impurity)
    # Filename of configuration file
    fnorg = "solver.norg.in"

    # Extract parameters from `Impurity` struct
    #
    # mune is fixed to zero, because the chemical potential is adsorbed
    # in the local impurity level.
    nband = imp.nband
    norbs = 2 * nband
    Uc    = imp.upara
    Jz    = imp.jpara
    mune  = 0.0
    beta  = imp.beta

    # Get number of Matsubara frequency points
    nmesh = get_m("nmesh")

    # Create the configuration file
    open(fnorg, "w") do fout
        # Print some warning message
        println(fout, "# Generated by $__LIBNAME__ at ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
        println(fout, "# PLEASE DO NOT MODIFY IT MANUALLY!")
        println(fout)

        # Print the standard parameters
        println(fout, "# Standard parameters: Running mode")
        println(fout, "mode = norm")
        println(fout)
        #
        println(fout, "# Standard parameters: System's size")
        println(fout, "nband = $nband")
        println(fout, "norbs = $norbs")
        println(fout)
        #
        println(fout, "# Standard parameters: Interaction")
        println(fout, "Uc    = $Uc")
        println(fout, "Jz    = $Jz")
        println(fout)
        #
        println(fout, "# Standard parameters: Others")
        println(fout, "mune  = $mune")
        println(fout, "beta  = $beta")
        println(fout, "nmesh = $nmesh")
        println(fout)

        # Print the user-supplied parameters
        #
        # The correctness of these auxiliary parameters should be checked
        # in `config.jl/chk_dict()`.
        println(fout, "# Auxiliary parameters: By users")
        foreach(x -> println(fout, x), get_s("params"))
    end # END OF IOSTREAM
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_delta(fmesh::Array{F64,1}, Delta::Array{C64,4})

Write the hybridization functions to the `solver.hyb.in` file, which is
suitable for the CTHYB quantum impurity solver.

See also: [`ctqmc_eimpx`](@ref).
"""
function ctqmc_delta(fmesh::Array{F64,1}, Delta::Array{C64,4})
    # Extract key parameters from `Delta`
    _, nband, nmesh, nspin = size(Delta)

    # Write the hybridization functions to `solver.hyb.in`
    open("solver.hyb.in", "w") do fout
        # If nspin is 2, everything is OK. However, if nspin is 1, the
        # spin-down part is missed.
        for s = 1:nspin
            for p = 1:nband
                orb = (s - 1) * nband + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n",
                        orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end # END OF P LOOP
        end # END OF S LOOP

        # Supplement the spin-down part if it is missed. Here, we just
        # assume that the spin-up and spin-down parts are degenerated.
        if nspin == 1
            for s = 1:nspin
                for p = 1:nband
                    orb = (2 - 1) * nband + p
                    for m = 1:nmesh
                        z = Delta[p,p,m,s]
                        @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n",
                            orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                    end
                    println(fout)
                    println(fout)
                end # END OF P LOOP
            end # END OF S LOOP
        end
    end # END OF IOSTREAM
end

"""
    norg_delta(fmesh::Array{F64,1}, Delta::Array{C64,4})

Write the hybridization functions to the `solver.hyb.in` file, which is
suitable for the NORG quantum impurity solver.

See also: [`norg_eimpx`](@ref).
"""
norg_delta(fmesh::Array{F64,1}, Delta::Array{C64,4}) = ctqmc_delta(fmesh, Delta)

"""
    ctqmc_eimpx(Eimpx::Array{C64,3})

Write the local impurity levels to the `solver.eimp.in` file, which is
suitable for the CTHYB quantum impurity solver.

See also: [`ctqmc_delta`](@ref).
"""
function ctqmc_eimpx(Eimpx::Array{C64,3})
    # Extract key parameters from `Eimpx`
    _, nband, nspin = size(Eimpx)

    # Analyze the orbital degeneracy
    symm = GetSymmetry(Eimpx)

    # Write the local impurity levels to `solver.eimp.in`
    open("solver.eimp.in", "w") do fout
        # If nspin is 2, everything is OK. However, if nspin is 1, the
        # spin-down part is missed.
        for s = 1:nspin
            for p = 1:nband
                orb = (s - 1) * nband + p
                z = Eimpx[p, p, s]
                @printf(fout, "%4i%16.8f%4i\n", orb, real(z), symm[p,s])
            end # END OF P LOOP
        end # END OF S LOOP

        # Supplement the spin-down part if it is missed. Here, we just
        # assume that the spin-up and spin-down parts are degenerated.
        if nspin == 1
            for s = 1:nspin
                for p = 1:nband
                    orb = (2 - 1) * nband + p
                    z = Eimpx[p, p, s]
                    @printf(fout, "%4i%16.8f%4i\n", orb, real(z), symm[p,s])
                end # END OF P LOOP
            end # END OF S LOOP
        end
    end # END OF IOSTREAM
end

"""
    norg_eimpx(Eimpx::Array{C64,3})

Write the local impurity levels to the `solver.eimp.in` file, which is
suitable for the NORG quantum impurity solver.

See also: [`norg_delta`](@ref).
"""
norg_eimpx(Eimpx::Array{C64,3}) = ctqmc_eimpx(Eimpx)

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_sigma(imp::Impurity)

Parse the `solver.sgm.dat` file, which is generated by the CTHYB quantum
impurity solver, to extract the bare self-energy functions.

In the `sigma_gather()` function, these data will be combined to generate
the `sigma.bare` file, which is essential for the DMFT engine.

See also: [`GetSigma`](@ref), [`solver_sigma`](@ref).
"""
function ctqmc_sigma(imp::Impurity)
    # File name for self-energy functions
    fsgm = "solver.sgm.dat"

    # To make sure the data file is present
    @assert isfile(fsgm)

    # Extract parameters from Impurity struct
    nband = imp.nband
    norbs = nband * 2
    nspin = 2 # In the CTHYB impurity solver, nspin is fixed to 2.

    # We don't know how many frequency points are used a priori.
    # Here, we used a trick to determine `nmesh`.
    nline = countlines(fsgm)
    nmesh = convert(I64, nline / norbs - 2)

    # Create an array for frequency mesh
    fmesh = zeros(F64, nmesh)

    # Create an array for self-energy functions
    Σ = zeros(C64, nband, nband, nmesh, nspin)

    # Parse the data file for the self-energy functions
    open(fsgm, "r") do fin

        # Go through each spin orientation
        for s = 1:nspin
            # Go through each band
            for b = 1:nband
                # Go through each frequency point
                for m = 1:nmesh
                    # Split the line
                    arr = line_to_array(fin)

                    # Extract frequency
                    fmesh[m] = parse(F64, arr[2])

                    # Extract self-energy functions
                    _re = parse(F64, arr[3])
                    _im = parse(F64, arr[4])
                    Σ[b,b,m,s] = _re + _im * im
                end # END OF M LOOP

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF B LOOP
        end # END OF S LOOP

    end # END OF IOSTREAM

    # Return the desired arrays
    return fmesh, Σ
end

"""
    norg_sigma(imp::Impurity)

Parse the `seimp.txt` file, which is generated by the NORG quantum
impurity solver, to extract the bare self-energy functions.

In the `sigma_gather()` function, these data will be combined to generate
the `sigma.bare` file, which is essential for the DMFT engine.

See also: [`GetSigma`](@ref), [`solver_sigma`](@ref).
"""
function norg_sigma(imp::Impurity)
    # File name for self-energy functions
    fsgm = "seimp.txt"

    # To make sure the data file is present
    @assert isfile(fsgm)

    # Extract parameters from Impurity struct
    nband = imp.nband

    # Determine number of spin orientations
    lines = readlines(fsgm)
    filter!(x -> contains(x, "Re11"), lines)
    nspin = length(lines)

    # We don't know how many frequency points are used a priori.
    # Here, we used a trick to determine `nmesh`.
    nline = countlines(fsgm)
    if nspin == 1
        nmesh = nline - 3
    else
        nmesh = convert(I64, nline / 2 - 3)
    end

    # Create an array for frequency mesh
    fmesh = zeros(F64, nmesh)

    # Create an array for self-energy functions
    Σ = zeros(C64, nband, nband, nmesh, nspin)

    # Parse the data file for the self-energy functions
    open(fsgm, "r") do fin

        # Go through each spin orientation
        for s = 1:nspin
            # Skip one line
            readline(fin)

            # Go through each frequency point
            for m = 1:nmesh
                # Split the line
                arr = line_to_array(fin)
                val = parse.(F64, arr)

                # Extract frequency
                fmesh[m] = val[1]

                # Extract self-energy functions
                # Starting positions for real and imaginary parts
                rshift = 1
                ishift = convert(I64, nband * nband + 2)
                #
                # Go through each band
                for p = 1:nband
                    # Go through each band
                    for q = 1:nband
                        rshift = rshift + 1
                        ishift = ishift + 1
                        _re = val[rshift]
                        _im = val[ishift]
                        Σ[p,q,m,s] = _re + _im * im
                    end # END OF Q LOOP
                end # END OF P LOOP
            end # END OF M LOOP

            # Skip two lines
            readline(fin)
            readline(fin)
        end # END OF S LOOP

    end # END OF IOSTREAM

    # Return the desired arrays
    return fmesh, Σ
end

"""
    ctqmc_nimpx(imp::Impurity)

Parse the `solver.nmat.dat` file, which is generated by the CTHYB
quantum impurity solver, to extract the impurity occupancy. Then the
fields `nup`, `ndown`, and `occup` in Impurity struct will be updated.

See also: [`GetNimpx`](@ref), [`solver_nimpx`](@ref).
"""
function ctqmc_nimpx(imp::Impurity)
    # File name for impurity occupancy
    fnmat = "solver.nmat.dat"

    # To make sure the data file is present
    if !isfile(fnmat)
        return
    end

    # Open the solver.nmat.dat file for reading
    fin = open(fnmat, "r")

    # Parse the data file to extract total impurity occupancy
    readuntil(fin, "sup")
    nup   = parse(F64, line_to_array(fin)[1])
    #
    readuntil(fin, "sdn")
    ndown = parse(F64, line_to_array(fin)[1])
    #
    readuntil(fin, "sum")
    occup = parse(F64, line_to_array(fin)[1])

    # Close the solver.nmat.dat file
    close(fin)

    # Update the Impurity struct
    imp.nup = nup
    imp.ndown = ndown
    imp.occup = occup
end

"""
    norg_nimpx(imp::Impurity)

Parse the `nmat.txt` file, which is generated by the NORG
quantum impurity solver, to extract the impurity occupancy. Then the
fields `nup`, `ndown`, and `occup` in Impurity struct will be updated.

See also: [`GetNimpx`](@ref), [`solver_nimpx`](@ref).
"""
function norg_nimpx(imp::Impurity)
    # File name for impurity occupancy
    fnmat = "nmat.txt"

    # To make sure the data file is present
    if !isfile(fnmat)
        return
    end

    # Open the nmat.txt file for reading
    fin = open(fnmat, "r")

    # Parse the data file to extract total impurity occupancy
    readuntil(fin, "sup")
    nup   = parse(F64, line_to_array(fin)[1])
    #
    readuntil(fin, "sdn")
    ndown = parse(F64, line_to_array(fin)[1])
    #
    readuntil(fin, "sum")
    occup = parse(F64, line_to_array(fin)[1])

    # Close the nmat.txt file
    close(fin)

    # Update the Impurity struct
    imp.nup = nup
    imp.ndown = ndown
    imp.occup = occup
end

"""
    ctqmc_edmft()

Parse the `solver.paux.dat` file, which is generated by the CTHYB quantum
impurity solver, to extract the interaction energy.

See also: [`GetEdmft`](@ref), [`solver_edmft`](@ref).
"""
function ctqmc_edmft()
    # File name for DMFT energy
    fene = "solver.paux.dat"

    # To make sure the data file is present
    if !isfile(fene)
        return 0.0
    end

    # Parse the data file to extract potential energy
    lines = readlines(fene)
    filter!(x -> contains(x, "epot:"), lines)
    @assert length(lines) == 1
    epot = parse(F64, line_to_array(lines[1])[2])

    # Return the desired value
    return epot
end

"""
    norg_edmft()

Parse the `solver.out` file, which is generated by the NORG quantum
impurity solver, to extract the interaction energy.

See also: [`GetEdmft`](@ref), [`solver_edmft`](@ref).
"""
function norg_edmft()
    # File name for DMFT energy
    fene = "solver.out"

    # To make sure the data file is present
    if !isfile(fene)
        return 0.0
    end

    # Parse the data file to extract potential energy
    lines = readlines(fene)
    filter!(x -> contains(x, "norg ground state"), lines)
    @assert length(lines) == 1
    epot = parse(F64, line_to_array(lines[1])[5])

    # Return the desired value
    return epot
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    GetSigma(imp::Impurity)

Extract self-energy functions from the output files of various quantum
impurity solvers. The data will be combined in the `sigma_gather()`
function. Then they will be fed back to the DMFT engine. The working
directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`solver_sigma`](@ref).
"""
function GetSigma(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    #
    # Since this function is called by sigma_gather(). we have to change
    # the current directory.
    cd("impurity.$index")

    # Activate the corresponding `solver_sigma()` functions for various
    # quantum impurity solvers
    fmesh, sigma = solver_sigma(_solver_, imp)

    # Enter the parent directory
    cd("..")

    # Return the desired array
    return (fmesh, sigma)
end

"""
    GetNimpx(imp::Impurity)

Extract impurity occupancy from the output files of various quantum
impurity solvers. Then the field `occup` in the Impurity struct will
be updated, which will then be used to evaluate the double counting
term for self-energy functions. The working directory of this function
must be the root folder.

The argument `imp` may be modified in this function.

See also: [`Impurity`](@ref), [`solver_nimpx`](@ref).
"""
function GetNimpx(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    #
    # Since this function is called by sigma_dcount(). we have to change
    # the current directory.
    cd("impurity.$index")

    # Activate the corresponding solver_nimpx() functions for various
    # quantum impurity solvers
    solver_nimpx(_solver_, imp)

    # Enter the parent directory
    cd("..")
end

"""
    GetEdmft(imp::Impurity)

Extract interaction energy (i.e potential energy) from the output files
of various quantum impurity solvers. The input Impurity struct won't be
modified. The working directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`solver_edmft`](@ref).
"""
function GetEdmft(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    cd("impurity.$index")

    # Activate the corresponding solver_edmft() functions for various
    # quantum impurity solvers
    epot = solver_edmft(_solver_)

    # Enter the parent directory
    cd("..")

    # Return the desired value
    return epot
end

#=
### *Service Functions* : *Quantum Impurity Problems*
=#

"""
    GetSymmetry(Eimpx::Array{C64,3})

Analyze the symmetry according to the diagonal elements of the matrix of
the local impurity levels.

See also: [`ctqmc_eimpx`](@ref).
"""
function GetSymmetry(Eimpx::Array{C64,3})
    # Extract some key parameters from Eimpx
    _, nband, nspin = size(Eimpx)

    # Extract the diagonal elements of Eimpx
    # We only need the real parts
    eimp = zeros(F64, nband, nspin)
    for s = 1:nspin
        for b = 1:nband
            eimp[b,s] = round(real(Eimpx[b,b,s]), digits = 2)
        end
    end # END OF S LOOP

    # Create an array for symmetry
    symm = zeros(I64, nband, nspin)

    # Go through each spin orientation
    for s = 1:nspin
        # Get the unique energy levels
        E = unique(eimp[:,s])

        # Scan all the impurity levels. Try to find out its corresponding
        # value in `E`.
        for b = 1:nband
            symm[b,s] = findfirst(x -> x == eimp[b,s], E)
        end
    end # END OF S LOOP

    # Adjust the symmetry vectors for spin down
    #
    # Sometimes the symmetry vectors for spin up and spin down are the
    # same, we have to avoid this case.
    if nspin == 2
        max_ind = maximum(symm[:,1])
        @. symm[:,2] = symm[:,2] + max_ind
    end

    # Return the desired array
    return symm
end

"""
    GetImpurity()

Return an array of Impurity struct, which encapsulates useful information
about the quantum impurity problems.

See also: [`Impurity`](@ref).
"""
function GetImpurity()
    # Initialize an array of Impurity struct
    AI = Impurity[]

    # Go through each quantum impurity problem
    for i = 1:get_i("nsite")

        # Extract essential parameters
        str = get_i("atoms")[i]
        index = i
        atoms = String(line_to_array(str)[1])
        sites = parse(I64, line_to_array(str)[3])
        equiv = get_i("equiv")[i]
        shell = get_i("shell")[i]
        ising = get_i("ising")[i]
        occup = get_i("occup")[i]
        upara = get_i("upara")[i]
        jpara = get_i("jpara")[i]
        lpara = get_i("lpara")[i]
        beta  = get_m("beta")

        # Call the constructor
        Im = Impurity(index, atoms, sites,
                      equiv, shell, ising,
                      occup, upara, jpara, lpara, beta)

        # Save Im in AI
        push!(AI, Im)

    end # END OF I LOOP

    # Return the desired array
    return AI
end

"""
    FixImpurity(ai::Array{Impurity,1})

Update the quantum impurity problems encapsulated in `ai` according
to the configuration parameters.

See also: [`Impurity`](@ref).
"""
function FixImpurity(ai::Array{Impurity,1})
    # Sanity check
    @assert length(ai) == get_i("nsite")

    # Go through each quantum impurity problem
    for i = 1:get_i("nsite")
        ai[i].ising = get_i("ising")[i]
        ai[i].occup = get_i("occup")[i]
        ai[i].upara = get_i("upara")[i]
        ai[i].jpara = get_i("jpara")[i]
        ai[i].lpara = get_i("lpara")[i]
    end # END OF I LOOP
end

"""
    CatImpurity(imp::Impurity)

Display the Impurity struct that need to be solved.

See also: [`Impurity`](@ref).
"""
function CatImpurity(imp::Impurity)
    println()
    println(blue("[ Impurity $(imp.index) ]"))
    println(repeat("=", 20))
    println("  atoms : ", imp.atoms)
    println("  sites : ", imp.sites)
    println("  equiv : ", imp.equiv)
    println("  shell : ", imp.shell)
    println("  ising : ", imp.ising)
    println("  occup : ", imp.occup)
    println("  nup   : ", imp.nup)
    println("  ndown : ", imp.ndown)
    println("  upara : ", imp.upara)
    println("  jpara : ", imp.jpara)
    println("  lpara : ", imp.lpara)
    println("  beta  : ", imp.beta)
    println("  nband : ", imp.nband)
    println(repeat("=", 20))
end
