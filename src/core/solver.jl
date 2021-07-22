#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/07/21
#

#=
### *CT-HYB₁ Quantum Impurity Solver*
=#

"""
    s_qmc1_init(it::IterInfo, imp::Impurity)

Check runtime environment of the CT-HYB₁ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_exec`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_init(it::IterInfo, imp::Impurity)
    # Print the header
    println("Engine : CT-HYB₁")
    println("Try to solve the quantum impurity problem: ", imp.index)
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

Launch the CT-HYB₁ quantum impurity solver.

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
    solver_home = query_solver("ct_hyb1")
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

        # Figure out the task that is doing
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

    # Wait for the dmft task to finish
    wait(t)

    # Extract how many monte carlo sampling blocks are executed
    lines = readlines("solver.out")
    filter!(x -> contains(x, "iter:"), lines)
    println("  > Finished after $(length(lines)) Monte Carlo sampling blocks")

    # Extract perturbation expansion order information
    println("Statistics about diagrammatic quantum Monte Carlo algorithm")
    println("  > Order / Count / Percent / Error bar")
    lines = readlines("solver.hist.dat")
    filter!(!endswith("0.000000"), lines)
    filter!(!startswith("#"), lines)
    foreach(x -> println(x), lines)
end

"""
    s_qmc1_save(it::IterInfo, imp::Impurity)

Backup output files of the CT-HYB₁ quantum impurity solver.

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
    union(fout, fgrn, fhyb, fsgm, faux) )
    println("  > Save the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    ctqmc_nimpx(imp)
    println("  > Extract the impurity occupancy from solver.nmat.dat: $(imp.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp.index] = imp.occup
end

"""
    s_qmc1_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Backup output files of the CT-HYB₁ quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
already in directory `impurity.2`.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_exec`](@ref).
"""
function s_qmc1_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
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
    union(fout, fgrn, fhyb, fsgm, faux) )
    println("  > Copy the key output files")

    # Update the `occup` field in `imp` (Impurity struct)
    ctqmc_nimpx(imp₂)
    println("  > Extract the impurity occupancy from solver.nmat.dat: $(imp₂.occup)")

    # Update the `it` (IterInfo) struct
    it.nf[imp₂.index] = imp₂.occup
end

#=
### *CT-HYB₂ Quantum Impurity Solver*
=#

"""
    s_qmc2_init(it::IterInfo)

Check runtime environment of the CT-HYB₂ quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_exec`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_init(it::IterInfo)
    sorry()
end

"""
    s_qmc2_exec(it::IterInfo)

Launch the CT-HYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB₂")

    # Print the footer for a better visualization
    println()
end

"""
    s_qmc2_save(it::IterInfo)

Backup output files of the CT-HYB₂ quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_save(it::IterInfo)
    sorry()
end

"""
    s_qmc2_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Backup output files of the CT-HYB₂ quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
already in directory `impurity.2`.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    sorry()
end

#=
### *HUB-I Quantum Impurity Solver*
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
    # Print the header
    println("Engine : HIA")

    # Print the footer for a better visualization
    println()
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
    s_hub1_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Backup output files of the HIA quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
already in directory `impurity.2`.

See also: [`s_hub1_init`](@ref), [`s_hub1_exec`](@ref).
"""
function s_hub1_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    sorry()
end

#=
### *NORG Quantum Impurity Solver*
=#

"""
    s_norg_init(it::IterInfo)

Check runtime environment of the NORG quantum impurity solver. Prepare
the necessary input files.

See also: [`s_norg_exec`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_init(it::IterInfo)
    sorry()
end

"""
    s_norg_exec(it::IterInfo)

Launch the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_exec(it::IterInfo)
    # Print the header
    println("Engine : NORG")

    # Print the footer for a better visualization
    println()
end

"""
    s_norg_save(it::IterInfo)

Backup output files of the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_save(it::IterInfo)
    sorry()
end

"""
    s_norg_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)

Backup output files of the NORG quantum impurity solver. We just copy
selected output files from impurity.1 to impurity.2. Be careful, now we
already in directory `impurity.2`.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_save(it::IterInfo, imp₁::Impurity, imp₂::Impurity)
    sorry()
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_setup(imp::Impurity)

Generate default configuration file (`solver.ctqmc.in`) for the CT-QMC
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

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_delta(fmesh::Array{F64,1}, Delta::Array{C64,4})

Write the hybridization functions to the `solver.hyb.in` file, which is
suitable for the CT-QMC quantum impurity solver.

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
    ctqmc_eimpx(Eimpx::Array{C64,3})

Write the local impurity levels to the `solver.eimp.in` file, which is
suitable for the CT-QMC quantum impurity solver.

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

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    ctqmc_sigma(imp::Impurity)

Parse the `solver.sgm.dat` file to extract the bare self-energy functions.

In the `sigma_gather()` function, these data will be combined to generate
the `sigma.bare` file, which is essential for the DMFT engine.

See also: [`Impurity`](@ref), [`GetSigma`](@ref).
"""
function ctqmc_sigma(imp::Impurity)
    # File name for self-energy functions
    fsgm = "solver.sgm.dat"

    # To make sure the data file is present
    @assert isfile(fsgm)

    # Extract parameters from Impurity struct
    nband = imp.nband
    norbs = nband * 2
    nspin = 2 # In the CT-QMC impurity solver, nspin is fixed to 2.

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
    ctqmc_nimpx(imp::Impurity)

Parse the `solver.nmat.dat` file to extract the impurity occupancy. Then
the field `occup` in Impurity struct will be updated.

In this function, only the total impurity occupancy of the current site
is return. However, sometimes we need to known the spin-up and spin-down
components. Later, we will expand the Impurity struct and this function
to fulfill this requirement.

See also: [`Impurity`](@ref), [`GetNimpx`](@ref).
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
    ctqmc_energy()

Parse the `solver.paux.dat` file to extract the interaction energy.

See also: [`GetEnergy`](@ref).
"""
function ctqmc_energy()
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

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    GetSigma(imp::Impurity)

Extract self-energy functions from the output files of various quantum
impurity solvers. The data will be combined in the `sigma_gather()`
function. Then they will be fed back to the DMFT engine. The working
directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`ctqmc_sigma`](@ref).
"""
function GetSigma(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    #
    # Since this function is called by sigma_gather(). we have to change
    # the current directory.
    cd("impurity.$index")

    # Determine the chosen solver
    engine = get_s("engine")

    # Activate the corresponding `solver_sigma()` functions for various
    # quantum impurity solvers
    fmesh = nothing
    sigma = nothing
    @cswitch engine begin
        @case "ct_hyb1"
            fmesh, sigma = ctqmc_sigma(imp)
            break

        @case "ct_hyb2"
            fmesh, sigma = ctqmc_sigma(imp)
            break

        @case "hub1"
            sorry()
            break

        @case "norg"
            sorry()
            break

        @default
            sorry()
            break
    end

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

See also: [`Impurity`](@ref), [`ctqmc_nimpx`](@ref).
"""
function GetNimpx(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    #
    # Since this function is called by sigma_dcount(). we have to change
    # the current directory.
    cd("impurity.$index")

    # Determine the chosen solver
    engine = get_s("engine")

    # Activate the corresponding solver_nimpx() functions for various
    # quantum impurity solvers
    @cswitch engine begin
        @case "ct_hyb1"
            ctqmc_nimpx(imp)
            break

        @case "ct_hyb2"
            ctqmc_nimpx(imp)
            break

        @case "hub1"
            sorry()
            break

        @case "norg"
            sorry()
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    GetEnergy(imp::Impurity)

Extract interaction energy (i.e potential energy) from the output files
of various quantum impurity solvers. The input Impurity struct won't be
modified. The working directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`ctqmc_energy`](@ref).
"""
function GetEnergy(imp::Impurity)
    # Get the index for current quantum impurity problem
    index = imp.index

    # Change the directory
    cd("impurity.$index")

    # Determine the chosen solver
    engine = get_s("engine")

    # Activate the corresponding solver_energy() functions for various
    # quantum impurity solvers
    @cswitch engine begin
        @case "ct_hyb1"
            epot = ctqmc_energy()
            break

        @case "ct_hyb2"
            epot = ctqmc_energy()
            break

        @case "hub1"
            sorry()
            break

        @case "norg"
            sorry()
            break

        @default
            sorry()
            break
    end

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
            symm[b,s] = findfirst(isequal(eimp[b]), E)
        end
    end # END OF S LOOP

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
    CatImpurity(imp::Impurity)

Display the Impurity struct that need to be solved.
"""
function CatImpurity(imp::Impurity)
    println()
    println(blue("[ Impurity $(imp.index) ]"))
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
end
