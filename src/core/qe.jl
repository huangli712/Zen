#
# Project : Pansy
# Source  : qe.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/10/17
#

#=
### *Multiple Dispatchers*
=#

"""
    dft_call(::QEEngine, it::IterInfo)

Try to carry out full DFT calculation with the `qe` code. It is only a
dispatcher. Similar function is defined in `vasp.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_call(::QEEngine, it::IterInfo)
    qe_init(it)
    qe_exec(it, true)  # Self-consistent calculation
    qe_exec(it, false) # Non-self-consistent calculation
    qe_save(it)
end

"""
    dft_stop(::QEEngine)

Try to terminate DFT calculation and kill running process of the DFT
backend. It supports the `qe` code. It is only a dispatcher. Similar
function is defined in `vasp.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_stop(::QEEngine)
    sorry()
end

"""
    dft_resume(::QEEngine)

Try to wake up the DFT backend and resume the DFT calculation. It only
supports the `qe` code. It is only a dispatcher. Similar function is
defined in `vasp.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_resume(::QEEngine)
    sorry()
end

"""
    adaptor_call(::QEEngine, D::Dict{Symbol,Any})

It is a dispatcher for the DFT-DMFT adaptor. It calls `qe_adaptor()`
function to deal with the outputs of the DFT backend (such as qe) and
generate key dataset for the next level adaptor (`WANNIERAdaptor`). Note
that similar function is also defined in `vasp.jl`.

See also: [`qe_adaptor`](@ref).
"""
function adaptor_call(::QEEngine, D::Dict{Symbol,Any})
    qeq_files()
    qe_adaptor(D)
end

#=
### *Driver Functions*
=#

"""
    qe_adaptor(D::Dict{Symbol,Any})

Adaptor support for the `quantum espresso` (`pwscf` code). It will parse
the output files of the `quantum espresso` code, extract the Kohn-Sham
dataset, and then fulfill the `DFTData` dict (i.e `D`).

The following output files of the `quantum espresso` code are needed:

* `scf.out`
* `nscf.out`

Note in the input file of the `quantum espresso` code, the `verbosity`
parameter must be set to 'high'.

See also: [`wannier_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function qe_adaptor(D::Dict{Symbol,Any})
    # Q01: Print the header
    println("Adaptor : QUANTUM ESPRESSO")
    println("Try to extract the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # Q02: Read in lattice structure
    D[:latt] = qeio_lattice(pwd(), false)

    # Q03: Read in kmesh and the corresponding weights
    D[:kmesh], D[:weight] = qeio_kmesh(pwd())

    # Q04: Read in band structure and the corresponding occupancies
    D[:enk], D[:occupy] = qeio_eigen(pwd())

    # Q05: Read in fermi level
    D[:fermi] = qeio_fermi(pwd(), false)

    # Q06: Generate MLWFs for the QE + WANNIER mode
    is_wannier() && qe_to_wan(D)

    # Q06: Generate projected local orbitals for the QE + PLO mode
    is_plo() && qe_to_plo(D)
end

"""
    qe_to_wan(D::Dict{Symbol,Any})

Try to call the `wannier90` and `pw2wannier90` codes to generate the
maximally-localized wannier functions. the `DFTData` dict (i.e `D`)
will not be modified.

See also: [`wannier_adaptor`](@ref).
"""
function qe_to_wan(D::Dict{Symbol,Any})
    # Check the validity of the original dict
    key_list = [:latt, :kmesh, :enk, :fermi]
    for k in key_list
        @assert haskey(D, k)
    end

    # Extract key parameters
    case = get_c("case") # Prefix for quantum espresso
    sp = get_d("lspins") # Is it a spin-polarized system?

    # Now this feature require quantum espresso as a DFT engine
    @assert is_qe() && is_wannier()

    # W01: Execute the wannier90 code to generate w90.nnkp
    if sp # For spin-polarized system
        # Spin up
        wannier_init(D, "up")
        wannier_exec("up", op = "-pp")
        wannier_save("up", op = "-pp")
        #
        # Spin down
        wannier_init(D, "dn")
        wannier_exec("dn", op = "-pp")
        wannier_save("dn", op = "-pp")
    else # For spin-unpolarized system
        wannier_init(D)
        wannier_exec(op = "-pp")
        wannier_save(op = "-pp")
    end

    # W02: Execute the pw2wannier90 code to generate necessary files
    # for the wannier90 code
    if sp # For spin-polarized system
        # Spin up
        pw2wan_init(case, "up")
        pw2wan_exec(case, "up")
        pw2wan_save("up")
        #
        # Spin down
        pw2wan_init(case, "dn")
        pw2wan_exec(case, "dn")
        pw2wan_save("dn")
    else # For spin-unpolarized system
        pw2wan_init(case)
        pw2wan_exec(case)
        pw2wan_save()
    end

    # W03: Execute the wannier90 code again to generate wannier functions
    if sp # For spin-polarized system
        # Spin up
        wannier_init(D, "up")
        wannier_exec("up")
        wannier_save("up")
        #
        # Spin down
        wannier_init(D, "dn")
        wannier_exec("dn")
        wannier_save("dn")
    else # For spin-unpolarized system
        wannier_init(D)
        wannier_exec()
        wannier_save()
    end
end

"""
    qe_to_plo(D::Dict{Symbol,Any})

Postprocess outputs of the `quantum espresso` (`pwscf` code), call the
`wannier90` and `pw2wannier90` codes to generate the projected local
orbitals (which are not maximally-localized wannier functions). The key
data are fed into the `DFTData` dict (i.e `D`).

Most of the functions used in the `qe_to_plo()` function are implemented
in another file (`wannier.jl`).

See also: [`plo_adaptor`](@ref), [`qe_adaptor`](@ref).
"""
function qe_to_plo(D::Dict{Symbol,Any})
    # Check the validity of the original dict
    key_list = [:latt, :kmesh, :enk, :fermi]
    for k in key_list
        @assert haskey(D, k)
    end

    # Extract key parameters
    case = get_c("case") # Prefix for quantum espresso
    sp = get_d("lspins") # Is it a spin-polarized system?

    # Now this feature require quantum espresso as a DFT engine
    @assert is_qe() && is_plo()

    # P01: Execute the wannier90 code to generate w90.nnkp
    if sp # For spin-polarized system
        # Spin up
        wannier_init(D, "up")
        wannier_exec("up", op = "-pp")
        wannier_save("up", op = "-pp")
        #
        # Spin down
        wannier_init(D, "dn")
        wannier_exec("dn", op = "-pp")
        wannier_save("dn", op = "-pp")
    else # For spin-unpolarized system
        wannier_init(D)
        wannier_exec(op = "-pp")
        wannier_save(op = "-pp")
    end

    # P02: Execute the pw2wannier90 code to generate necessary files for
    # the wannier90 code
    if sp # For spin-polarized system
        # Spin up
        pw2wan_init(case, "up")
        pw2wan_exec(case, "up")
        pw2wan_save("up")
        #
        # Spin down
        pw2wan_init(case, "dn")
        pw2wan_exec(case, "dn")
        pw2wan_save("dn")
    else # For spin-unpolarized system
        pw2wan_init(case)
        pw2wan_exec(case)
        pw2wan_save()
    end

    # This is enough. We do not really need the wannier functions. The
    # case.amn file is what we need.

    # P03: Read accurate band eigenvalues from w90.eig
    #
    # D[:enk] will be updated
    # Be careful, the eigenvalues will be calibrated in plo_fermi().
    # This is different from what we see in wannier_adaptor().
    if sp # For spin-polarized system
        # Spin up
        eigs_up = w90_read_eigs("up")
        nband, nkpt = size(eigs_up)
        eigs_up = reshape(eigs_up, (nband, nkpt, 1))
        #
        # Spin down
        eigs_dn = w90_read_eigs("dn")
        nband, nkpt = size(eigs_dn)
        eigs_dn = reshape(eigs_dn, (nband, nkpt, 1))
        #
        # Sanity check
        @assert size(eigs_up) == size(eigs_dn)
        #
        # Concatenate eigs_up and eigs_dn
        D[:enk] = cat(eigs_up, eigs_dn, dims = 3)
        #
        # Sanity check
        @assert size(D[:enk]) == (nband, nkpt, 2)
    else # For spin-unpolarized system
        eigs = w90_read_eigs()
        nband, nkpt = size(eigs)
        eigs = reshape(eigs, (nband, nkpt, 1))
        D[:enk] = deepcopy(eigs)
        @assert size(D[:enk]) == (nband, nkpt, 1)
    end

    # P04: Read projected local orbitals from w90.amn
    #
    # D[:chipsi] will be created
    if sp # For spin-polarized system
        # Spin up
        Aup = w90_read_amat("up")
        nproj, nband, nkpt = size(Aup)
        Aup = reshape(Aup, (nproj, nband, nkpt, 1))
        #
        # Spin down
        Adn = w90_read_amat("dn")
        nproj, nband, nkpt = size(Adn)
        Adn = reshape(Adn, (nproj, nband, nkpt, 1))
        #
        # Sanity check
        @assert size(Aup) == size(Adn)
        #
        # Concatenate Aup and Adn
        D[:chipsi] = cat(Aup, Adn, dims = 4)
        #
        # Sanity check
        @assert size(D[:chipsi]) == (nproj, nband, nkpt, nspin)
    else # For spin-unpolarized system
        Amn = w90_read_amat()
        nproj, nband, nkpt = size(Amn)
        Amn = reshape(Amn, (nproj, nband, nkpt, 1))
        D[:chipsi] = deepcopy(Amn)
        @assert size(D[:chipsi]) == (nproj, nband, nkpt, 1)
    end

    # P05: Setup the PrTrait and PrGroup structs
    #
    # D[:PT] and D[:PG] will be created
    latt =D[:latt]
    if sp # For spin-polarized system
        # Spin up
        PT_up, PG_up = w90_make_group(latt, "up")
        #
        # Spin down
        PT_dn, PG_dn = w90_make_group(latt, "dn")
        #
        # Merge PT_up and PT_dn
        @assert PT_up == PT_dn
        D[:PT] = deepcopy(PT_up)
        #
        # Merge PG_up and PG_dn
        @assert PG_up == PG_dn
        D[:PG] = deepcopy(PG_up)
    else # For spin-unpolarized system
        PT, PG = w90_make_group(latt)
        D[:PT] = deepcopy(PT)
        D[:PG] = deepcopy(PG)
    end
end

"""
    qe_init(it::IterInfo)

Check the runtime environment of `quantum espresso` (`pwscf`), prepare
necessary input files.

See also: [`qe_exec`](@ref), [`qe_save`](@ref).
"""
function qe_init(it::IterInfo)
    # Print the header
    println("Engine : QUANTUM ESPRESSO")
    println("Try to perform ab initio electronic structure calculation")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for quantum espresso")

    # Prepare essential input files
    # Copy QE.INP (It is used as a template.)
    cp("../QE.INP", joinpath(pwd(), "QE.INP"), force = true)
    println("  > File QE.INP is ready")
    #
    # Create the real input file, case.scf and case.nscf.
    ControlNL, AtomicSpeciesBlock = qec_input(it)
    case = get_c("case")
    println("  > File $case.scf is ready")
    println("  > File $case.nscf is ready")
    #
    # Check the pseudopotentials
    pseudo_dir = strip(ControlNL["pseudo_dir"],''')
    upf = map(x -> joinpath(pseudo_dir, x.upf), AtomicSpeciesBlock.data)
    for f in upf
        @assert isfile(f)
        println("  > File $(basename(f)) is ready")
    end
end

"""
    qe_exec(it::IterInfo, scf::Bool = true)

Execute the `quantum espresso` (`pwscf`) program, monitor the convergence
progress, and output the relevant information. The argument `scf` controls
which input file should be used. If `scf == true`, then the input file is
`case.scf`, or else it is `case.nscf`.

In order to execute this function correctly, you have to setup the
following environment variables:

* QE_HOME

and make sure the file `MPI.toml` is available.

See also: [`qe_init`](@ref), [`qe_save`](@ref).
"""
function qe_exec(it::IterInfo, scf::Bool = true)
    # Print the header
    println("Detect the runtime environment for quantum espresso")

    # Determine mpi prefix (whether the code is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of quantum espresso
    qe_home = query_dft(_engine_)
    println("  > Home directory for quantum espresso: ", qe_home)

    # Select suitable quantum espresso program
    # We use the same code pw.x for with or without spin-orbit coupling
    if get_d("lspinorb")
        qe_exe = "$qe_home/pw.x"
    else
        qe_exe = "$qe_home/pw.x"
    end
    @assert isfile(qe_exe)
    println("  > Executable program is available: ", basename(qe_exe))

    # Assemble command
    if isnothing(mpi_prefix)
        qe_cmd = qe_exe
    else
        qe_cmd = split("$mpi_prefix $qe_exe", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', qe_cmd))")

    # Determine suitable input and output files
    case = get_c("case")
    if scf
        finp = "$case.scf"
        fout = "scf.out"
    else
        finp = "$case.nscf"
        fout = "nscf.out"
    end
    println("  > Self-consistent DFT calculation: $scf")
    println("  > Applying input file: $finp")
    println("  > Applying output file: $fout")

    # Print the header
    println("Launch the computational engine quantum espresso")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$qe_cmd`, stdin = finp, stdout = fout))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to `fout`.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # To ensure that the task is executed
    while true
        sleep(2)
        istaskstarted(t) && break
    end

    # Analyze the `fout` file during the calculation
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

        # For self-consistent DFT calculation mode
        if scf

            # Parse the `fout` file
            iters = readlines(fout)
            filter!(x -> contains(x, "iteration #"), iters)
            ethrs = readlines(fout)
            filter!(x -> contains(x, "ethr ="), ethrs)

            # Figure out the number of iterations (`ni`) and deltaE (`dE`)
            if length(ethrs) > 0
                arr = line_to_array(iters[end])
                ni = parse(I64, arr[3])
                arr = line_to_array(ethrs[end])
                dE = strip(arr[3],',')
            else # The first iteration has not been finished
                ni = 0
                dE = "unknown"
            end

            # Print the log to screen
            @printf("  > Elapsed %4i seconds, %3i iterations (dE = %12s)\r", 5*c, ni, dE)

        # For non-self-consistent DFT calculation mode
        else

            # Parse the `fout` file
            lines = readlines(fout)
            filter!(x -> contains(x, "Computing kpt #"), lines)

            # Figure out how many k-points are finished
            if length(lines) > 0
                arr = line_to_array(lines[end])
                ikpt = parse(I64, arr[4])
                nkpt = parse(I64, arr[6])
            else # The first k-point has not been finished
                ikpt = 0
                nkpt = 0
            end

            # Print the log to screen
            @printf("  > Elapsed %4i seconds, %4i of %4i k-points\r", 5*c, ikpt, nkpt)

        end

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the quantum espresso task to finish
    wait(t)

    # Extract how many iterations are executed
    if scf
        iters = readlines(fout)
        filter!(x -> contains(x, "iteration #"), iters)
        println("  > Converged after $(length(iters)) iterations")
    # Extract how many k-points are finished
    else
        lines = readlines(fout)
        filter!(x -> contains(x, "Computing kpt #"), lines)
        println("  > Calculated eigenvalues for $(length(lines)) k-points")
    end
end

"""
    qe_save(it::IterInfo)

Backup the output files of `quantum espresso` (`pwscf`) if necessary.
Furthermore, the DFT fermi level in `IterInfo` struct is also updated
(i.e. `IterInfo.μ₀`).

See also: [`qe_init`](@ref), [`qe_exec`](@ref).
"""
function qe_save(it::IterInfo)
    # Print the header
    println("Finalize the computational task")

    # Store the data files
    #
    # Create list of files
    case = get_c("case")
    fl = ["$case.scf", "$case.nscf", "scf.out", "nscf.out"]
    #
    # Go through the file list, backup the files one by one.
    for i in eachindex(fl)
        f = fl[i]
        cp(f, "$f.$(it.I₃)", force = true)
    end
    println("  > Save the key output files")

    # Anyway, the DFT fermi level is extracted from scf.out, and its
    # value will be saved at IterInfo.μ₀.
    it.μ₀ = qeio_fermi(pwd())
    println("  > Extract the fermi level from scf.out: $(it.μ₀) eV")

    # We also try to read the DFT band energy from scf.out, and its
    # value will be saved at IterInfo.et.
    it.et.dft = qeio_energy(pwd())
    println("  > Extract the DFT band energy from scf.out: $(it.et.dft) eV")
end

#=
### *Service Functions* : *Group A*
=#

"""
    qec_input(it::IterInfo)

It will parse the `QE.INP` file at first.  The `QE.INP` is a standard,
but mini input file for `quantum espresso` (`pwscf`). It only includes
three namelists (namely `control`, `system`, and `electrons`) and three
cards (namely `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, and `K_POINTS`).
If you want to support more input entries, please make your own
modifications here.

Then this function will try to customize these namelists and cards
according to the setup in `case.toml`.

At last, it will try to generate the input files for `quantum espresso`
(`pwscf`). They are `case.scf` and `case.nscf`. As shown by their names,
one file is for the self-consistent calculation, while another one is
for the non-self-consistent calculation.

The return values of this function are namelist (`control`) and card
(`ATOMIC_SPECIES`), which will be used to check the pseudopotentials
within the `qe_init()` function.

See also: [`QENamelist`](@ref), [`QECard`](@ref).
"""
function qec_input(it::IterInfo)
    # Check the file status
    finput = "QE.INP"
    @assert isfile(finput)

    # Parse three namelists, control, system, and electrons.
    lines = readlines(finput)
    ControlNL = parse(QENamelist, lines, "control")
    SystemNL = parse(QENamelist, lines, "system")
    ElectronsNL = parse(QENamelist, lines, "electrons")

    # Parse three cards, ATOMIC_SPECIES, ATOMIC_POSITIONS, and K_POINTS.
    line = read(finput, String)
    AtomicSpeciesBlock = parse(AtomicSpeciesCard, line)
    AtomicPositionsBlock = parse(AtomicPositionsCard, line)
    KPointsBlock = parse(KPointsCard, line)

    # Customize the namelists and cards according to case.toml
    #
    # For smearing
    smear = get_d("smear")
    smear == "mp2"   && begin
        SystemNL["occupations"] = "'smearing'"
        SystemNL["smearing"] = "'m-p'"
    end
    #
    smear == "mp1"   && begin
        SystemNL["occupations"] = "'smearing'"
        SystemNL["smearing"] = "'m-p'"
    end
    #
    smear == "gauss" && begin
        SystemNL["occupations"] = "'smearing'"
        SystemNL["smearing"] = "'gauss'"
    end
    #
    smear == "tetra" && begin
        SystemNL["occupations"] = "'tetrahedra'"
        delete!(SystemNL, "smearing")
    end

    # For kmesh density
    #
    # Note that if kmesh == "file", the original setup in QE.INP is
    # kept. In other words, KPointsBlock will not be changed.
    kmesh = get_d("kmesh")
    if isa(KPointsBlock,AutoKmeshCard)
        shift = copy(KPointsBlock.data.shift)
        kmesh == "accurate" && begin
            KPointsBlock = AutoKmeshCard([10, 10, 10], shift)
        end
        kmesh == "medium"   && begin
            KPointsBlock = AutoKmeshCard([08, 08, 08], shift)
        end
        kmesh == "coarse"   && begin
            KPointsBlock = AutoKmeshCard([06, 06, 06], shift)
        end
    end

    # For magnetic moment
    magmom = get_d("magmom")
    if !isa(magmom, Missing)
        moment = parse.(F64, line_to_array(magmom))
        @assert length(moment) == parse(I64, SystemNL["ntyp"])
        for i in eachindex(moment)
            SystemNL["starting_magnetization($i)"] = moment[i]
        end
    end

    # For symmetry
    lsymm = get_d("lsymm")
    if lsymm
        SystemNL["nosym"] = ".false."
    else # Ignore the symmetry completely
        SystemNL["nosym"] = ".true."
    end

    # For spin polarizations
    lspins = get_d("lspins")
    if lspins
        SystemNL["nspin"] = 2
        @assert haskey(SystemNL, "starting_magnetization(1)")
    else
        SystemNL["nspin"] = 1
    end

    # For spin-orbit coupling
    lspinorb = get_d("lspinorb")
    if lspinorb
        delete!(SystemNL, "nspin")
        SystemNL["noncolin"] = ".true."
        SystemNL["lspinorb"] = ".true."
    else
        SystemNL["noncolin"] = ".false."
        SystemNL["lspinorb"] = ".false."
    end

    # Special treatment for verbosity
    ControlNL["verbosity"] = "'high'"

    # Special treatment for pseudo_dir
    pseudo_dir = strip(ControlNL["pseudo_dir"],''') # Get rid of `
    pseudo_dir = joinpath("..", pseudo_dir)
    ControlNL["pseudo_dir"] = "'$pseudo_dir'" # Add ' back

    # Build input files for quantum espresso
    #
    # Get case's name
    case = get_c("case")
    #
    # Setup filenames
    fscf = "$case.scf"
    fnscf = "$case.nscf"
    #
    # For case.scf
    ControlNL["calculation"] = "'scf'"
    open(fscf, "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
        write(fout, AtomicPositionsBlock)
        write(fout, KPointsBlock)
    end
    #
    # For case.nscf
    ControlNL["calculation"] = "'nscf'"
    delete!(ControlNL, "restart_mode")
    #
    # We have to specify the k-points explicitly during the
    # non-self-consistent calculations.
    begin
        # The tetrahedron algorithm requires an uniform, automatically
        # generted k-mesh. This can not be fulfilled in the calculations.
        # So we have to use the smearing algorithm to calculate the
        # occupations.
        SystemNL["occupations"] = "'smearing'"
        kmesh == "accurate" && begin
            KPointsBlock = SpecialPointsCard(10)
        end
        #
        kmesh == "medium"   && begin
            KPointsBlock = SpecialPointsCard(08)
        end
        #
        kmesh == "coarse"   && begin
            KPointsBlock = SpecialPointsCard(06)
        end
        #
        kmesh == "file"     && begin
            @assert isa(KPointsBlock, SpecialPointsCard)
        end
    end
    #
    open(fnscf, "w") do fout
        write(fout, ControlNL)
        write(fout, SystemNL)
        write(fout, ElectronsNL)
        write(fout, AtomicSpeciesBlock)
        write(fout, AtomicPositionsBlock)
        write(fout, KPointsBlock)
    end

    # Return the namelist (control) and the card (ATOMIC_SPECIES), which
    # will be used to check whether the pseudopotential files are ready.
    return ControlNL, AtomicSpeciesBlock
end

#=
### *Service Functions* : *Group B*
=#

"""
    qeq_files(f::String)

Check the essential output files by `quantum espresso` (`pwscf`). Here
`f` means only the directory that contains the desired files.

See also: [`adaptor_core`](@ref).
"""
function qeq_files(f::String)
    fl = ["scf.out", "nscf.out"]
    for i in eachindex(fl)
        @assert isfile( joinpath(f, fl[i]) )
    end
end

"""
    qeq_files()

Check the essential output files by `quantum espresso` (`pwscf`) in the
current directory.

See also: [`adaptor_core`](@ref).
"""
qeq_files() = qeq_files(pwd())

#=
### *Service Functions* : *Group C*
=#

"""
    qeio_energy(f::String)

Reading quantum espresso's `scf.out` file, return DFT total energy, which
will be used to determine the DFT + DMFT energy. Here `f` means only the
directory that contains `scf.out`.
"""
function qeio_energy(f::String)
    # Try to figure out whether the scf.out file is valid
    lines = readlines(joinpath(f, "scf.out"))
    filter!(x -> contains(x, "!    total energy"), lines)
    @assert length(lines) == 1

    # Extract the total energy
    etot = parse(F64, line_to_array(lines[end])[5])

    # Return the desired value
    return etot
end

"""
    qeio_energy()

Reading quantum espresso's `scf.out` file, return DFT total energy, which
will be used to determine the DFT + DMFT energy.
"""
qeio_energy() = qeio_energy(pwd())

"""
    qeio_lattice(f::String, silent::Bool = true)

Reading quantum espresso's `scf.out` file, return crystallography
information. Here `f` means only the directory that contains `scf.out`.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
function qeio_lattice(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse lattice")
    !silent && println("  > Open and read scf.out")

    # Read in all lines in `scf.out`
    lines = readlines(joinpath(f, "scf.out"))

    # Sorry, `scf.out` does not contain any information about case.
    _case = get_c("case")

    # Get the scaling factor
    # Shall we convert it into SI?
    ind = findfirst(x -> contains(x, "lattice parameter"), lines)
    @assert ind > 0
    scale = parse(F64, line_to_array(lines[ind])[5])

    # Get the total number of atoms
    ind = findfirst(x -> contains(x, "number of atoms/cell"), lines)
    @assert ind > 0
    natom = parse(I64, line_to_array(lines[ind])[5])

    # Get the number of sorts of atoms
    ind = findfirst(x -> contains(x, "number of atomic types"), lines)
    @assert ind > 0
    nsort = parse(I64, line_to_array(lines[ind])[6])

    # Now all the parameters are ready, we would like to create
    # `Lattice` struct here.
    latt = Lattice(_case, scale, nsort, natom)

    # Get the lattice vectors
    ind = findfirst(x -> contains(x, "crystal axes:"), lines)
    @assert ind > 0
    latt.lvect[1, :] = parse.(F64, line_to_array(lines[ind+1])[4:6])
    latt.lvect[2, :] = parse.(F64, line_to_array(lines[ind+2])[4:6])
    latt.lvect[3, :] = parse.(F64, line_to_array(lines[ind+3])[4:6])

    # Get the symbol list
    ind = findfirst(x -> contains(x, "atomic species   valence"), lines)
    @assert ind > 0
    for i = 1:nsort
        latt.sorts[i, 1] = string( line_to_array(lines[ind+i])[1] )
        latt.sorts[i, 2] = 0 # Init it with 0
    end

    # Get the atom list and the coordinates of atoms
    ind = findfirst(x -> contains(x, "Cartesian axes"), lines)
    @assert ind > 0
    for i = 1:natom
        latt.atoms[i] = line_to_array(lines[ind+2+i])[2]
        latt.coord[i, :] = parse.(F64, line_to_array(lines[ind+2+i])[7:9])
    end

    # Well, now we try to update latt.sorts further
    for i = 1:natom
        for j = 1:nsort
            if latt.atoms[i] == latt.sorts[j, 1]
                latt.sorts[j, 2] = latt.sorts[j, 2] + 1
            end
        end
    end

    # Print some useful information to check
    !silent && println("  > System: ", latt._case)
    !silent && println("  > Atoms: ", latt.atoms)

    # Return the desired struct
    return latt
end

"""
    qeio_lattice()

Reading quantum espresso's `scf.out` file, return crystallography
information.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
qeio_lattice() = qeio_lattice(pwd())

"""
    qeio_kmesh(f::String)

Reading quantum espresso's `nscf.out` file, return `kmesh` and `weight`.
Here `f` means only the directory that contains `nscf.out`.

Note in `scf.out`, the 𝑘-mesh is not uniform. So we have to read k-mesh
from the `nscf.out`. In addition, the verbosity parameter must be set to
'high' in the input file.

See also: [`qeio_tetra`](@ref), [`irio_kmesh`](@ref).
"""
function qeio_kmesh(f::String)
    # Print the header
    println("Parse kmesh and weight")
    println("  > Open and read nscf.out")

    # Read in all lines in `nscf.out`
    lines = readlines(joinpath(f, "nscf.out"))

    # Extract number of 𝑘-points
    ind = findfirst(x -> contains(x, "number of k points="), lines)
    @assert ind > 0
    nkpt = parse(I64, line_to_array(lines[ind])[5])

    # Create arrays
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt)

    # Read in the 𝑘-points and their weights
    for i = 1:nkpt
        k1k2 = line_to_array(lines[ind+1+i])[5:6]
        kmesh[i, 1:2] = parse.(F64, k1k2)
        #
        k3 = line_to_array(lines[ind+1+i])[7]
        k3 = strip(k3, [',', ')']) # Get rid of some chars
        kmesh[i, 3] = parse(F64, k3)
        #
        w  = line_to_array(lines[ind+1+i])[10]
        weight[i] = parse(F64, w)
    end

    # Normalize the weight
    # We have to make sure the sum of weights is equal to nkpt.
    wsum = sum(weight)
    @. weight = weight / wsum * nkpt

    # Print some useful information to check
    println("  > Number of k-points: ", nkpt)
    println("  > Total sum of weights: ", sum(weight))
    println("  > Shape of Array kmesh: ", size(kmesh))
    println("  > Shape of Array weight: ", size(weight))

    # Return the desired arrays
    return kmesh, weight
end

"""
    qeio_kmesh()

Reading quantum espresso's `nscf.out` file, return `kmesh` and `weight`.

See also: [`qeio_tetra`](@ref), [`irio_kmesh`](@ref).
"""
qeio_kmesh() = qeio_kmesh(pwd())

"""
    qeio_eigen(f::String)

Reading quantum espresso's `nscf.out` file, return energy band structure
information. Here `f` means only the directory that contains `nscf.out`.

Note that in `scf.out`, the eigenvalues may be not defined on the uniform
𝑘-mesh. So we have to read eigenvalues from the `nscf.out` file.

Note that the eigenvalues read from `nscf.out` is somewhat coarse. They
should be updated by the values read from `case.eig`.

See also: [`irio_eigen`](@ref).
"""
function qeio_eigen(f::String)
    # Print the header
    println("Parse enk and occupy")
    println("  > Open and read nscf.out")

    # Read in all lines in `nscf.out`
    lines = readlines(joinpath(f, "nscf.out"))

    # Extract number of 𝑘-points
    ind = findfirst(x -> contains(x, "number of k points="), lines)
    @assert ind > 0
    nkpt = parse(I64, line_to_array(lines[ind])[5])

    # Extract number of bands
    ind = findfirst(x -> contains(x, "number of Kohn-Sham states"), lines)
    @assert ind > 0
    nband = parse(I64, line_to_array(lines[ind])[5])

    # Determine number of spins
    #
    # The `nscf.out` file contains `SPIN UP` and `SPIN DOWN` blocks if
    # the system is spin-polarized.
    ind = findall(x -> contains(x, " SPIN "), lines)
    @assert length(ind) == 2 || length(ind) == 0
    if length(ind) == 2
        nspin = 2
    else
        nspin = 1
    end

    # Create arrays
    enk = zeros(F64, nband, nkpt, nspin)
    occupy = zeros(F64, nband, nkpt, nspin)

    # Read in the energy bands and the corresponding occupations
    #
    # Determine the start of data block
    start = findfirst(x -> contains(x, "End of band structure calculation"), lines)
    @assert start > 0
    #
    # Determine how many lines are there for each k-point
    elements_per_line = 8
    nrow = div(nband, elements_per_line)
    nrem = rem(nband, elements_per_line)
    #
    # Go through each spin
    for s = 1:nspin
        # Skip `SPIN UP` and `SPIN DOWN` lines.
        if nspin == 2
            start = start + 3
        end
        #
        # Go through each k-point
        for k = 1:nkpt
            # Read eigenvalues
            start = start + 3
            if nrow > 1 # nband > elements_per_line
                for r = 1:nrow
                    start = start + 1
                    bs = (r - 1) * elements_per_line + 1
                    be = (r - 1) * elements_per_line + elements_per_line
                    enk[bs:be,k,s] = parse.(F64, line_to_array(lines[start]))
                end # END OF R LOOP
                if nrem > 0
                    start = start + 1
                    bs = nrow * elements_per_line + 1
                    be = nband
                    @assert nrem == be - bs + 1
                    enk[bs:be,k,s] = parse.(F64, line_to_array(lines[start]))
                end
            else
                @assert nrow == 1
                start = start + 1
                enk[:,k,s] = parse.(F64, line_to_array(lines[start]))
            end
            #
            # Read occupations
            start = start + 2
            if nrow > 1 # nband > elements_per_line
                for r = 1:nrow
                    start = start + 1
                    bs = (r - 1) * elements_per_line + 1
                    be = (r - 1) * elements_per_line + elements_per_line
                    occupy[bs:be,k,s] = parse.(F64, line_to_array(lines[start]))
                end # END OF R LOOP
                if nrem > 0
                    start = start + 1
                    bs = nrow * elements_per_line + 1
                    be = nband
                    @assert nrem == be - bs + 1
                    occupy[bs:be,k,s] = parse.(F64, line_to_array(lines[start]))
                end
            else
                @assert nrow == 1
                start = start + 1
                occupy[:,k,s] = parse.(F64, line_to_array(lines[start]))
            end
        end # END OF K LOOP
    end # END OF S LOOP

    # Print some useful information to check
    println("  > Number of DFT bands: ", nband)
    println("  > Number of k-points: ", nkpt)
    println("  > Number of spins: ", nspin)
    println("  > Shape of Array enk: ", size(enk))
    println("  > Shape of Array occupy: ", size(occupy))

    # return the desired arrays
    return enk, occupy
end

"""
    qeio_eigen()

Reading quantum espresso's `nscf.out` file, return energy band structure
information.

See also: [`irio_eigen`](@ref).
"""
qeio_eigen() = qeio_eigen(pwd())

"""
    qeio_fermi(f::String, silent::Bool = true)

Reading quantum espresso's `nscf.out` file, return the fermi level. Here
`f` means only the directory that contains `scf.out`.

See also: [`irio_fermi`](@ref).
"""
function qeio_fermi(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse fermi level")
    !silent && println("  > Open and read nscf.out")

    # Try to figure out whether the scf.out file is valid
    lines = readlines(joinpath(f, "nscf.out"))
    filter!(x -> contains(x, "Fermi energy"), lines)
    @assert length(lines) == 1

    # Extract the fermi level
    fermi = parse(F64, line_to_array(lines[end])[5])

    # Print some useful information to check
    !silent && println("  > Fermi level: $fermi eV")

    # Return the desired data
    return fermi
end

"""
    qeio_fermi()

Reading quantum espresso's `scf.out` file, return the fermi level.

See also: [`irio_fermi`](@ref).
"""
qeio_fermi() = qeio_fermi(pwd())

#=
*Remarks* :

The following codes are internally used to parse the input files of
`quantum espresso` (`pwscf`). Please do not change them or export them.
=#

#=
### *Customized Structs : Abstract Types*
=#

"""
    QEInputEntry

An abstract type representing an input component of `quantum espresso`
(`pwscf`). Note that all other input types (such as `QENamelist` and
`QECard`) should subtype `QEInputEntry`.  It is used to build a internal
type system.

See also: [`QECard`](@ref), [`QENamelist`](@ref).
"""
abstract type QEInputEntry end

"""
    QECard

It represents abstract cards in the input file of `quantum espresso`
(`pwscf`).  It is used to build the internal type system. The input
file of `quantum espresso` (`pwscf`) consists of various cards and
namelists, represented by `QECard` and `QENamelist`, respectively.

See also: [`QENamelist`](@ref).
"""
abstract type QECard <: QEInputEntry end

"""
    KPointsCard

Represent an abstract card (`K-POINTS`) in the input file of
`quantum espresso` (`pwscf`).
"""
abstract type KPointsCard <: QECard end

#=
### *Customized Structs : K-Grids*
=#

"""
    ReciprocalPoint

Represent a special point of the 3D Brillouin zone. Each of them has
a weight `w`.

### Members

* coord  -> Coordinates, i.e., ``k_x``, ``k_y``, and ``k_z``.
* weight -> Weight for the ``k``-point.

See also: [`MonkhorstPackGrid`](@ref).
"""
struct ReciprocalPoint
    coord  :: Vector{F64}
    weight :: F64

    # Inner constructor
    function ReciprocalPoint(coord, weight)
        @assert length(coord) == 3
        @assert weight > 0.0
        return new(coord, weight)
    end
end

"""
    ReciprocalPoint(x::F64, y::F64, z::F64, w::F64)

Constructor for `ReciprocalPoint`.

See also: [`MonkhorstPackGrid`](@ref).
"""
function ReciprocalPoint(x::F64, y::F64, z::F64, w::F64)
    return ReciprocalPoint([x, y, z], w)
end

"""
    MonkhorstPackGrid

Represent the Monkhorst-Pack grid.

### Members

* mesh  -> A length-three vector specifying the ``k``-point grid
           (``nk_1 × nk_2 × nk_3``) as in Monkhorst-Pack grids.
* shift -> A length-three vector specifying whether the grid is displaced
           by half a grid step in the corresponding directions.

See also: [`ReciprocalPoint`](@ref).
"""
struct MonkhorstPackGrid
    mesh  :: Vector{I64}
    shift :: Vector{Bool}

    # Inner constructor
    function MonkhorstPackGrid(mesh, shift)
        @assert length(mesh) == 3
        @assert length(shift) == 3
        @assert all(mesh .≥ 1)
        #
        if eltype(shift) != Bool
            shift = Bool.(shift)
        end
        #
        return new(mesh, shift)
    end
end

"""
    MonkhorstPackGrid(k₁::I64, k₂::I64, k₃::I64, s₁::I64, s₂::I64, s₃::I64)

Constructor for `MonkhorstPackGrid`.

See also: [`ReciprocalPoint`](@ref).
"""
function MonkhorstPackGrid(k₁::I64, k₂::I64, k₃::I64, s₁::I64, s₂::I64, s₃::I64)
    k = [k₁, k₂, k₃]
    s = [s₁, s₂, s₃]
    return MonkhorstPackGrid(k, s)
end

#=
### *Customized Structs : Basic Entries*
=#

"""
    AtomicSpecies

Represent each line of the `ATOMIC_SPECIES` card in the input file of
`quantum espresso` (`pwscf`).

### Members

* atom -> Label of the atom. Maximum total length cannot exceed
          3 characters.
* mass -> Mass of the atomic species in atomic unit. Used only when
          performing molecular dynamics (MD) run or structural
          optimization runs using damped MD.
* upf  -> File containing pseudopotential for this species.

### Examples

```julia-repl
julia> AtomicSpecies("C1", 12, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
AtomicSpecies("C1", 12.0, "C.pbe-n-kjpaw_psl.1.0.0.UPF")
```

See also: [`AtomicSpeciesCard`](@ref).
"""
struct AtomicSpecies
    atom :: String
    mass :: F64
    upf  :: String

    # Inner constructor
    function AtomicSpecies(atom::Union{AbstractChar,AbstractString}, mass, upf)
        @assert length(atom) ≤ 3
        return new(string(atom), mass, upf)
    end
end

"""
    AtomicPosition

Represent each line of the `ATOMIC_POSITIONS` card in the input file of
`quantum espresso` (`pwscf`).

### Members

* atom   -> Label of the atom as specified in `AtomicSpecies`. It
            accepts at most 3 characters.
* pos    -> Atomic positions. A three-element vector of floats.
* if_pos -> Component `i` of the force for this atom is multiplied
            by `if_pos(i)`, which must be either `0` or `1`.  Used
            to keep selected atoms and/or selected components fixed
            in MD dynamics or structural optimization run.

### Examples

```julia-repl
julia> AtomicPosition('O', [0, 0, 0])
AtomicPosition("O", [0.0, 0.0, 0.0], Bool[1, 1, 1])
```

See also: [`AtomicPositionsCard`](@ref).
"""
struct AtomicPosition
    atom   :: String
    pos    :: Vector{F64}
    if_pos :: Vector{Bool}

    # Inner constructor
    function AtomicPosition(atom::Union{AbstractChar,AbstractString}, pos, if_pos)
        @assert length(atom) ≤ 3
        @assert length(pos) == 3
        @assert length(if_pos) == 3
        return new(string(atom), pos, if_pos)
    end
end

"""
    AtomicSpecies(x::AtomicPosition, mass, upf)

Constructor for `AtomicSpecies`.

### Examples

```julia-repl
julia> AtomicSpecies(
           AtomicPosition('S', [0.500000000, 0.288675130, 1.974192764]),
           32.066,
           "S.pz-n-rrkjus_psl.0.1.UPF",
       )
AtomicSpecies("S", 32.066, "S.pz-n-rrkjus_psl.0.1.UPF")
```

See also: [`AtomicSpeciesCard`](@ref).
"""
AtomicSpecies(x::AtomicPosition, mass, upf) = AtomicSpecies(x.atom, mass, upf)

"""
    AtomicPosition(atom, pos)

Constructors for `AtomicPosition`.

See also: [`AtomicPositionsCard`](@ref).
"""
AtomicPosition(atom, pos) = AtomicPosition(atom, pos, trues(3))

"""
    AtomicPosition(x::AtomicSpecies, pos, if_pos)

Constructors for `AtomicPosition`.

### Examples

```julia-repl
julia> AtomicPosition(
           AtomicSpecies('S', 32.066, "S.pz-n-rrkjus_psl.0.1.UPF"),
           [0.500000000, 0.288675130, 1.974192764],
       )
AtomicPosition("S", [0.5, 0.28867513, 1.974192764], Bool[1, 1, 1])
```

See also: [`AtomicPositionsCard`](@ref).
"""
AtomicPosition(x::AtomicSpecies, pos, if_pos) = AtomicPosition(x.atom, pos, if_pos)

#=
### *Customized Structs : Input Blocks*
=#

"""
    QENamelist

Represent a namelist in the input file of `quantum espresso` (`pwscf`),
a basic Fortran data structure.

### Members

* name -> Name of the namelist. It should be `control`, `system`, or
          `electrons`. If you want to support more namelists, please
          make your own modifications.
* data -> A dict containing pairs of key and value.

See also: [`QECard`](@ref).
"""
mutable struct QENamelist <: QEInputEntry
    name :: AbstractString
    data :: Dict{AbstractString,Any}
end

"""
    Base.haskey(qnl::QENamelist, key::AbstractString)

Examine the existence of an entry (specified by `key`) in the namelist
object (`qnl`).

See also: [`QENamelist`](@ref).
"""
Base.haskey(qnl::QENamelist, key::AbstractString) = haskey(qnl.data, key)

"""
    Base.getindex(qnl::QENamelist, key::AbstractString)

Return an entry (specified by `key`) in the namelist object (`qnl`).

See also: [`QENamelist`](@ref).
"""
Base.getindex(qnl::QENamelist, key::AbstractString) = qnl.data[key]

"""
    Base.setindex!(qnl::QENamelist, value, key::AbstractString)

Modify an entry (specified by `key`) in the namelist object (`qnl`).

See also: [`QENamelist`](@ref).
"""
function Base.setindex!(qnl::QENamelist, value, key::AbstractString)
    qnl.data[key] = value
end

"""
    Base.delete!(qnl::QENamelist, key::AbstractString)

Remove an entry (specified by `key`) in the namelist object (`qnl`).

See also: [`QENamelist`](@ref).
"""
function Base.delete!(qnl::QENamelist, key::AbstractString)
    delete!(qnl.data, key)
end

"""
    AtomicSpeciesCard

Represent the `ATOMIC_SPECIES` card in the input file of
`quantum espresso` (`pwscf`).

### Members

* data -> A vector containing `AtomicSpecies`.

See also: [`AtomicSpecies`](@ref).
"""
struct AtomicSpeciesCard <: QECard
    data :: Vector{AtomicSpecies}
end

"""
    AtomicPositionsCard

Represent the `ATOMIC_POSITIONS` card in the input file of
`quantum espresso` (`pwscf`).

### Members

* data   -> A vector containing `AtomicPosition`s.
* option -> The scheme about how to define atomic positions.

See also: [`AtomicPosition`](@ref).
"""
struct AtomicPositionsCard <: QECard
    data   :: Vector{AtomicPosition}
    option :: String

    # Inner constructor
    function AtomicPositionsCard(data, option = "alat")
        @assert option in ("alat",
                           "bohr",
                           "angstrom",
                           "crystal",
                           "crystal_sg")
        return new(data, option)
    end
end

"""
    AutoKmeshCard

Represent the `K_POINTS` card in the input file of `quantum espresso`
(`pwscf`) (be compatible with the `automatic` mode only).

See also: [`KPointsCard`](@ref).
"""
struct AutoKmeshCard <: KPointsCard
    data :: MonkhorstPackGrid
end

"""
    AutoKmeshCard(mesh::Vector{I64}, shift::Vector{I64})

Constructor for `AutoKmeshCard`.

See also: [`KPointsCard`](@ref).
"""
function AutoKmeshCard(mesh::Vector{I64}, shift::Vector{Bool})
    return AutoKmeshCard(MonkhorstPackGrid(mesh, shift))
end

"""
    GammaPointCard

Represent the `K_POINTS` card in the input file of `quantum espresso`
(`pwscf`) (be compatible with the `gamma` mode).

See also: [`KPointsCard`](@ref).
"""
struct GammaPointCard <: KPointsCard end

"""
    SpecialKPointsCard

Represent the `K_POINTS` card in the input file of `quantum espresso`
(`pwscf`) (be compatible with the `tpiba` or `crystal` mode).

### Members

* data   -> A vector containing `ReciprocalPoint`s.
* option -> The way about how to define ``k``-mesh.

See also: [`KPointsCard`](@ref).
"""
struct SpecialPointsCard <: KPointsCard
    data   :: Vector{ReciprocalPoint}
    option :: String

    # Inner constructor
    function SpecialPointsCard(data, option = "tpiba")
        @assert option in ("tpiba",
                           "crystal",
                           "tpiba_b",
                           "crystal_b",
                           "tpiba_c",
                           "crystal_c")
        return new(data, option)
    end
end

"""
    SpecialPointsCard(nkx::I64, nky::I64, nkz::I64, option::String)

Constructor for `SpecialPointsCard`. This function is insprired by the
`kmesh.pl` tool as included in the `wannier90` package.

See also: [`KPointsCard`](@ref).
"""
function SpecialPointsCard(nkx::I64, nky::I64, nkz::I64, option::String)
    # Sanity check
    @assert nkx ≥ 1
    @assert nky ≥ 1
    @assert nkz ≥ 1

    # Calculate total number of k-points
    nkpt = nkx * nky * nkz

    # Calculate weight per k-point
    w = 1.0 / nkpt

    # Generate uniform k-mesh. Note that the k-point is represented
    # by ReciprocalPoint struct.
    data = Vector{ReciprocalPoint}(undef,nkpt)
    c = 0 # Counter
    #
    for x = 0:nkx - 1
        for y = 0:nky - 1
            for z = 0:nkz - 1
                c = c + 1
                kx = float(x) / nkx
                ky = float(y) / nky
                kz = float(z) / nkz
                data[c] = ReciprocalPoint(kx, ky, kz, w)
            end
        end
    end
    #
    @assert c == nkpt

    # Call the default constructor
    return SpecialPointsCard(data, option)
end

"""
    SpecialPointsCard(nkx::I64, option::String = "crystal")

Constructor for `SpecialPointsCard`.

See also: [`KPointsCard`](@ref).
"""
function SpecialPointsCard(nkx::I64, option::String = "crystal")
    return SpecialPointsCard(nkx, nkx, nkx, option)
end

#=
### *Predefined Regex*

*Remarks* :

Note that these regular expressions are followed by various combinations
of the `i`, `m`, and `x` flags. These flags have the following meanings:

* `i` : Do case-insensitive pattern matching.
* `m` : Treat string as multiple lines.
* `s` : Treat string as single line.
* `x` : Tells the regular expression parser to ignore most whitespace
        that is neither backslashed nor within a character class.

See: https://docs.julialang.org/en/v1/manual/strings/#Regular-Expressions
=#

#=
### *Remarks* : For ATOMIC_SPECIES Card
=#

const ATOMIC_SPECIES_BLOCK = r"""
^ [ \t]* ATOMIC_SPECIES [ \t]* \R+
(?P<block>
    (?:
        ^ [ \t]*
            \S+ [ \t]+
            (?:[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+
            \S+ [ \t]*
        \R?
    )+
)
"""imx

const ATOMIC_SPECIES_ITEM = r"""
^ [ \t]*
    (?P<name>\S+) [ \t]+
    (?P<mass>[-+]?[0-9]*\.[0-9]+|[0-9]+\.?[0-9]*) [ \t]+
    (?P<pseudo>\S+) [ \t]*
\R?
"""mx

#=
*Remarks* : For ATOMIC_POSITIONS Card
=#

const ATOMIC_POSITIONS_BLOCK = r"""
^ \s* ATOMIC_POSITIONS \s*                   # Atomic positions start with that string
[{(]? \s* (?P<units>\S+?)? \s* [)}]? \s* $\R # The units are after the string in optional brackets
(?P<block>                                   # This is the block of positions
    (
        (
            \s*                              # White space in front of the element spec is ok
            (
                [A-Za-z]+[A-Za-z0-9]{0,2}    # Element spec
                (
                    \s+                      # White space in front of the number
                    [-|+]?                   # Plus or minus in front of the number (optional)
                    (
                        (
                            \d*              # optional decimal in the beginning .0001 is ok, for example
                            [\.]             # There has to be a dot followed by
                            \d+              # at least one decimal
                        )
                        |                    # OR
                        (
                            \d+              # at least one decimal, followed by
                            [\.]?            # an optional dot ( both 1 and 1. are fine)
                            \d*              # And optional number of decimals (1.00001)
                        )                    # followed by optional decimals
                    )
                    ([E|e|d|D][+|-]?\d+)?    # optional exponents E+03, e-05
                ){3}                         # I expect three float values
                ((\s+[0-1]){3}\s*)?          # Followed by optional ifpos
                \s*                          # Followed by optional white space
                |
                \#.*                         # If a line is commented out, that is also ok
                |
                \!.*                         # Comments also with excl. mark in fortran
            )
            |                                # OR
            \s*                              # A line only containing white space
         )
        \R                                   # line break at the end
    )+                                       # A positions block should be one or more lines
)
"""imx

const ATOMIC_POSITIONS_ITEM = r"""
^                                       # Linestart
[ \t]*                                  # Optional white space
(?P<name>[A-Za-z]+[A-Za-z0-9]{0,2})\s+  # get the symbol, max 3 chars, starting with a char
(?P<x>                                  # Get x
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<y>                                  # Get y
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]+
(?P<z>                                  # Get z
    [\-|\+]?(\d*[\.]\d+ | \d+[\.]?\d*)
    ([E|e|d|D][+|-]?\d+)?
)
[ \t]*
(?P<fx>[01]?)                           # Get fx
[ \t]*
(?P<fy>[01]?)                           # Get fx
[ \t]*
(?P<fz>[01]?)                           # Get fx
"""mx

#=
*Remarks* : For K_POINTS Card
=#

const K_POINTS_AUTOMATIC_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* automatic [ \t]* [)}]? [ \t]* \R
^ [ \t]*
    (\d+) [ \t]+
    (\d+) [ \t]+
    (\d+) [ \t]+
    (\d+) [ \t]+
    (\d+) [ \t]+
    (\d+) [ \t]*
\R?
"""imx

const K_POINTS_GAMMA_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* gamma [ \t]* [)}]? [ \t]* \R*
"""imx

const K_POINTS_SPECIAL_BLOCK = r"""
^ [ \t]* K_POINTS [ \t]* [{(]? [ \t]* (?P<type>\S+?)? [ \t]* [)}]? [ \t]* \R+
    ^ [ \t]* \S+ [ \t]* \R+  # nks
        (?P<block>
            (?:
                ^ [ \t]*
                    \S+ [ \t]+
                    \S+ [ \t]+
                    \S+ [ \t]+
                    \S+ [ \t]*
                \R*
            )+
        )
"""imx

const K_POINTS_SPECIAL_ITEM = r"""
^ [ \t]*
    (\S+) [ \t]+
    (\S+) [ \t]+
    (\S+) [ \t]+
    (\S+) [ \t]*
\R?
"""mx

#=
### *Input File Parsers*
=#

"""
    Base.parse(::Type{QENamelist}, strs::Vector{String}, name::String)

Parse the `QENamelist` object. The name of the namelist is specified
by argument `name`.

See also: [`QENamelist`](@ref).
"""
function Base.parse(::Type{QENamelist}, strs::Vector{String}, name::String)
    # Try to parse `strs` to extract the data
    #
    # Prepare necessary data structures
    group_data = []
    group_meet = false
    group_name = "unknown"
    #
    # Go through each line
    for l in eachindex(strs)
        # Get rid of the blanks and `,`
        strip_line = strip(strip(strs[l]), ',')
        # Meet a namelist
        if startswith(strip_line, "&")
            group_name = lowercase(split(strip_line, "&")[2])
            # It is the namelist what we try to locate
            if group_name == name
                group_meet = true
            end
        end

        # End of namelist
        if startswith(strip_line, "/")
            group_meet = false
        end

        # Store the data in group_data
        if group_meet
            push!(group_data, strip_line)
        end
    end
    #
    # The first element is not useful. We have to remove it.
    if isempty(group_data)
        return nothing # Fail to figure out a namelist
    else
        popfirst!(group_data)
    end

    # Try to build a QENamelist object
    #
    # Prepare a dictionary, it will contain the namelist's data.
    NLData = Dict{AbstractString,Any}()
    #
    # Go throught each element in group_data
    for i in eachindex(group_data)
        # There are multiple entries in the element (line)
        if count(",", group_data[i]) > 0
            pairs = split(group_data[i], ",")
            for j in eachindex(pairs)
                key, value = map(x -> strip(x), split(pairs[j], "="))
                NLData[key] = value
            end
        # There is only one entry in the element (line)
        else
            key, value = map(x -> strip(x), split(group_data[i], "="))
            NLData[key] = value
        end
    end

    # Return the desired object
    return QENamelist(name, NLData)
end

"""
    Base.parse(::Type{T}, str::AbstractString)

Parse the `QECard` object. Now we support the following cards:

* `ATOMIC_SPECIES` (`AtomicSpeciesCard`)
* `ATOMIC_POSITIONS` (`AtomicPositionsCard`)
* `K_POINTS` (`AutoKmeshCard`, `GammaPointCard`, `SpecialKPointsCard`)

See also: [`QECard`](@ref).
"""
function Base.parse(::Type{T}, str::AbstractString) where {T<:QECard}
    x = tryparse(T, str)
    if x === nothing
        error("cannot find card `$T`!")
    else
        return x
    end
end

"""
    Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)

Try to parse the `AtomicSpeciesCard` object.

See also: [`AtomicSpeciesCard`](@ref).
"""
function Base.tryparse(::Type{AtomicSpeciesCard}, str::AbstractString)
    m = match(ATOMIC_SPECIES_BLOCK, str)

    # Function `match` only searches for the first match of the regular
    # expression, so it could be a `nothing`.
    if m !== nothing
        content = only(m.captures)
        return AtomicSpeciesCard(
            map(eachmatch(ATOMIC_SPECIES_ITEM, content)) do matched
                captured = matched.captures
                atom, mass, upf = captured[1], parse(F64, captured[2]), captured[3]
                AtomicSpecies(atom, mass, upf)
            end,
        )
    end
end

"""
    Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)

Try to parse the `AtomicPositionsCard` object.

See also: [`AtomicPositionsCard`](@ref).
"""
function Base.tryparse(::Type{AtomicPositionsCard}, str::AbstractString)
    m = match(ATOMIC_POSITIONS_BLOCK, str)

    # Function `match` only searches for the first match of the regular
    # expression, so it could be a `nothing`.
    if m !== nothing
        if string(m.captures[1]) === nothing
            @info "No option is specified, 'alat' is assumed."
            option = "alat"
        else
            option = string(m.captures[1])
        end
        content = m.captures[2]
        return AtomicPositionsCard(
            map(eachmatch(ATOMIC_POSITIONS_ITEM, content)) do matched
                # The `matched` cannot be a `nothing` since we have
                # tested by the block regular expression.
                captured = matched.captures
                #
                # The `if_pos` field is optionally given by users. If
                # they do not give, we provide the default values `1`.
                if_pos = map(x -> isempty(x) ? 1 : parse(I64, x), captured[11:13])
                #
                # The `atom` and `pos` fields are mandatory. So we do
                # not need special treatment.
                atom, pos = captured[1],
                    map(x -> parse(F64, x), [captured[2], captured[5], captured[8]])
                AtomicPosition(atom, pos, if_pos)
            end,
            option,
        )
    end
end

"""
    Base.tryparse(::Type{KPointsCard}, str::AbstractString)

Try to parse the `KPointsCard` object.

See also: [`KPointsCard`](@ref).
"""
function Base.tryparse(::Type{KPointsCard}, str::AbstractString)
    for T in (AutoKmeshCard, GammaPointCard, SpecialPointsCard)
        x = tryparse(T, str)
        if x !== nothing
            return x
        end
    end
end

"""
    Base.tryparse(::Type{AutoKmeshCard}, str::AbstractString)

Try to parse the `AutoKmeshCard` object.

See also: [`AutoKmeshCard`](@ref).
"""
function Base.tryparse(::Type{AutoKmeshCard}, str::AbstractString)
    m = match(K_POINTS_AUTOMATIC_BLOCK, str)

    if m !== nothing
        data = map(x -> parse(I64, x), m.captures)
        return AutoKmeshCard(MonkhorstPackGrid(data[1:3], data[4:6]))
    end
end

"""
    Base.tryparse(::Type{GammaPointCard}, str::AbstractString)

Try to parse the `GammaPointCard` object.

See also: [`GammaPointCard`](@ref).
"""
function Base.tryparse(::Type{GammaPointCard}, str::AbstractString)
    m = match(K_POINTS_GAMMA_BLOCK, str)
    return m === nothing ? nothing : GammaPointCard()
end

"""
    Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)

Try to parse the `SpecialPointsCard` object.

See also: [`SpecialPointsCard`](@ref).
"""
function Base.tryparse(::Type{SpecialPointsCard}, str::AbstractString)
    m = match(K_POINTS_SPECIAL_BLOCK, str)

    if m !== nothing
        option = m.captures[1] === nothing ? "tpiba" : m.captures[1]
        return SpecialPointsCard(
            map(eachmatch(K_POINTS_SPECIAL_ITEM, m.captures[2])) do matched
                ReciprocalPoint(map(x -> parse(F64, x), matched.captures)...)
            end,
            option,
        )
    end
end

#=
### *Input File Writers*
=#

"""
    Base.write(io::IO, x::QENamelist)

Write the `QENamelist` object to `IOStream`.

See also: [`QENamelist`](@ref).
"""
function Base.write(io::IO, x::QENamelist)
    println(io, " &$(x.name)")
    for key in keys(x.data)
        println(io, "    $key = ", x.data[key], ",")
    end
    println(io, " /")
end

"""
    Base.write(io::IO, x::AtomicSpeciesCard)

Write the `AtomicSpeciesCard` object to `IOStream`.

See also: [`AtomicSpeciesCard`](@ref).
"""
function Base.write(io::IO, x::AtomicSpeciesCard)
    println(io, "ATOMIC_SPECIES")
    for i = 1:length(x.data)
        AS = x.data[i]
        println(io, " $(AS.atom)  $(AS.mass)  $(AS.upf)")
    end
end

"""
    Base.write(io::IO, x::AtomicPositionsCard)

Write the `AtomicPositionsCard` object to `IOStream`.

See also: [`AtomicPositionsCard`](@ref).
"""
function Base.write(io::IO, x::AtomicPositionsCard)
    println(io, "ATOMIC_POSITIONS {$(x.option)}")
    for i =  1:length(x.data)
        AP = x.data[i]
        print(io, " $(AP.atom) ")
        @printf(io, "%6.3f %6.3f %6.3f\n", AP.pos...)
    end
end

"""
    Base.write(io::IO, x::AutoKmeshCard)

Write the `AutoKmeshCard` object to `IOStream`.

See also: [`AutoKmeshCard`](@ref).
"""
function Base.write(io::IO, x::AutoKmeshCard)
    println(io, "K_POINTS {automatic}")
    MPG = x.data
    @printf(io, "%3i%3i%3i%2i%2i%2i\n", MPG.mesh..., MPG.shift...)
end

"""
    Base.write(io::IO, x::GammaPointCard)

Write the `GammaPointCard` object to `IOStream`.

See also: [`GammaPointCard`](@ref).
"""
function Base.write(io::IO, x::GammaPointCard)
    println(io, "K_POINTS {gamma}")
end

"""
    Base.write(io::IO, x::SpecialPointsCard)

Write the `SpecialPointsCard` object to `IOStream`.

See also: [`SpecialPointsCard`](@ref).
"""
function Base.write(io::IO, x::SpecialPointsCard)
    println(io, "K_POINTS {$(x.option)}")
    nks = length(x.data)
    println(io, "  $nks")
    for i = 1:nks
        RP = x.data[i]
        @printf(io, " %11.7f%11.7f%11.7f%16.12f\n", RP.coord..., RP.weight)
    end
end
