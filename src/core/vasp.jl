#
# Project : Pansy
# Source  : vasp.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/10/17
#

#=
### *Multiple Dispatchers*
=#

"""
    dft_call(::VASPEngine, it::IterInfo)

Try to carry out full DFT calculation with the `vasp` code. It is only a
dispatcher. Similar function is defined in `qe.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_call(::VASPEngine, it::IterInfo)
    vasp_init(it)
    vasp_exec(it)
    vasp_save(it)
end

"""
    dft_stop(::VASPEngine)

Try to terminate DFT calculation and kill running process of the DFT
backend. It supports the `vasp` code. It is only a dispatcher. Similar
function is defined in `qe.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_stop(::VASPEngine)
    vasp_stop()
end

"""
    dft_resume(::VASPEngine)

Try to wake up the DFT backend and resume the DFT calculation. It only
supports the `vasp` code. It is only a dispatcher. Similar function is
defined in `qe.jl` as well.

See also: [`_engine_`](@ref).
"""
function dft_resume(::VASPEngine)
    # Reactivate the DFT engine
    @time_call vasp_back()

    # Wait the DFT engine to finish its job and sleep
    suspend(2)

    # Get the DFT energy
    edft = vaspio_energy()

    return edft
end

"""
    adaptor_call(::VASPEngine, D::Dict{Symbol,Any})

It is a dispatcher for the DFT-DMFT adaptor. It calls `vasp_adaptor()`
function to deal with the outputs of the DFT backend (such as vasp) and
generate key dataset for the next level adaptor (`PLOAdaptor`). Note
that similar function is also defined in `qe.jl`.

See also: [`vasp_adaptor`](@ref).
"""
function adaptor_call(::VASPEngine, D::Dict{Symbol,Any})
    vaspq_files()
    vasp_adaptor(D)
end

#=
### *Driver Functions*
=#

"""
    vasp_adaptor(D::Dict{Symbol,Any})

Adaptor support for vasp code. It will parse the output files of vasp
code, extract the Kohn-Sham dataset,  and then fulfill the `DFTData`
dict (i.e `D`).

The following vasp's output files must be presented:

* `POSCAR`
* `IBZKPT`
* `EIGENVAL`
* `LOCPROJ`
* `DOSCAR`
* `OSZICAR`

See also: [`plo_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function vasp_adaptor(D::Dict{Symbol,Any})
    # V01: Print the header
    println("Adaptor : VASP")
    println("Try to extract the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # V02: Read in lattice structure
    D[:latt] = vaspio_lattice(pwd(), false)

    # V03: Read in kmesh and the corresponding weights
    D[:kmesh], D[:weight] = vaspio_kmesh(pwd())

    # V04: Read in band structure and the corresponding occupancies
    D[:enk], D[:occupy] = vaspio_eigen(pwd())

    # V05: Read in raw projectors, traits, and groups
    D[:PT], D[:PG], D[:chipsi] = vaspio_projs(pwd())

    # V06: Read in fermi level
    D[:fermi] = vaspio_fermi(pwd(), false)

    # V07: Read in tetrahedron data if they are available
    if get_d("smear") == "tetra"
        D[:volt], D[:itet] = vaspio_tetra(pwd())
    end
end

"""
    vasp_init(it::IterInfo)

Check the runtime environment of vasp, prepare necessary input files.

See also: [`vasp_exec`](@ref), [`vasp_save`](@ref).
"""
function vasp_init(it::IterInfo)
    # Print the header
    println("Engine : VASP")
    println("Try to perform ab initio electronic structure calculation")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for vasp")

    # Prepare essential input files
    #
    # Copy POTCAR and POSCAR
    cp("../POTCAR", joinpath(pwd(), "POTCAR"), force = true)
    cp("../POSCAR", joinpath(pwd(), "POSCAR"), force = true)
    println("  > File POTCAR is ready")
    println("  > File POSCAR is ready")
    #
    # How about INCAR
    if it.I₃ == 0
        # Generate INCAR automatically
        vaspc_incar(it.μ₀, it.sc)
    else
        # Maybe we need to update INCAR file here
        vaspc_incar(it.μ₁, it.sc)
    end
    println("  > File INCAR is ready")
    #
    # Well, perhaps we need to generate the KPOINTS file by ourselves.
    get_d("kmesh") == "file" && vaspc_kpoints()
    println("  > File KPOINTS is ready")

    # Check essential input files
    flist = ["INCAR", "POSCAR", "POTCAR"]
    get_d("kmesh") == "file" && push!(flist, "KPOINTS")
    for i in eachindex(flist)
        filename = flist[i]
        if !isfile(filename)
            error("Please make sure the file $filename is available")
        end
    end
end

"""
    vasp_exec(it::IterInfo)

Execute the vasp program, monitor the convergence progress, and output
the relevant information. Especially, if it.sc == 2 (self-consistent
mode), this function will launch the vasp code, then return immediately.

In order to execute this function correctly, you have to setup the
following environment variables:

* VASP_HOME

and make sure the file `MPI.toml` is available.

See also: [`vasp_init`](@ref), [`vasp_save`](@ref).
"""
function vasp_exec(it::IterInfo)
    # Print the header
    println("Detect the runtime environment for vasp")

    # Determine mpi prefix (whether the vasp is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of vasp
    vasp_home = query_dft(_engine_)
    println("  > Home directory for vasp: ", vasp_home)

    # Select suitable vasp program
    if get_d("lspinorb")
        vasp_exe = "$vasp_home/vasp_ncl"
    else
        vasp_exe = "$vasp_home/vasp_std"
    end
    @assert isfile(vasp_exe)
    println("  > Executable program is available: ", basename(vasp_exe))

    # Assemble command
    if isnothing(mpi_prefix)
        vasp_cmd = vasp_exe
    else
        if it.sc == 2
            vasp_cmd = split("mpiexec -n 1 $vasp_exe", " ")
        else
            vasp_cmd = split("$mpi_prefix $vasp_exe", " ")
        end
    end
    println("  > Assemble command: $(prod(x -> x * ' ', vasp_cmd))")

    # Print the header
    println("Launch the computational engine vasp")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to vasp.out.
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

    # Special treatment for self-consistent mode
    if it.sc == 2
        println("Let's Rock & Roll\n")
        return
    end

    # Analyze the vasp.out file during the calculation
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

        # Parse vasp.out file
        iters = readlines("vasp.out")
        filter!(x -> contains(x, "DAV:"), iters)

        # Figure out the number of iterations (`ni`) and deltaE (`dE`)
        if length(iters) > 0
            arr = line_to_array(iters[end])
            ni = parse(I64, arr[2])
            dE = arr[4]
        else # The first iteration has not been finished
            ni = 0
            dE = "unknown"
        end

        # Print the log to screen
        @printf("  > Elapsed %4i seconds, %3i iterations (dE = %12s)\r", 5*c, ni, dE)

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the vasp task to finish
    wait(t)

    # Extract how many iterations are executed
    iters = readlines("vasp.out")
    filter!(x -> contains(x, "DAV:"), iters)
    println("  > Converged after $(length(iters)) iterations")
end

"""
    vasp_save(it::IterInfo)

Backup the output files of vasp if necessary. Furthermore, the DFT fermi
level and the DFT band energy in `IterInfo` struct will also be updated
(i.e `IterInfo.μ₀` and `IterInfo.et.dft`).

See also: [`vasp_init`](@ref), [`vasp_exec`](@ref).
"""
function vasp_save(it::IterInfo)
    # Special treatment for self-consistent mode
    it.sc == 2 && return

    # Print the header
    println("Finalize the computational task")

    # Store the data files
    #
    # Create list of files
    fl = ["INCAR", "vasp.out"]
    #
    # Go through the file list, backup the files one by one.
    for i in eachindex(fl)
        f = fl[i]
        cp(f, "$f.$(it.I₃)", force = true)
    end
    println("  > Save the key output files")

    # Anyway, the DFT fermi level is extracted from DOSCAR, and its
    # value will be saved at IterInfo.μ₀.
    it.μ₀ = vaspio_fermi(pwd())
    println("  > Extract the fermi level from DOSCAR: $(it.μ₀) eV")

    # We also try to read the DFT band energy from OSZICAR, and its
    # value will be saved at IterInfo.et.
    it.et.dft = vaspio_energy(pwd())
    println("  > Extract the DFT band energy from OSZICAR: $(it.et.dft) eV")
end

"""
    vasp_back()

Reactivate the vasp engine to continue the charge self-consistent
DFT + DMFT calculation. It will prepare the file `GAMMA`, and try
to create a lock file (`vasp.lock`). Then the vasp engine will wake
up and continue to work automatically.

See also: [`vasp_stop`](@ref).
"""
function vasp_back()
    # Read in the correction for density matrix which is produced by
    # the dmft1 code (Dyson).
    println("Read correction for density matrix")
    _, kwin, gcorr = read_gcorr("../dmft2/dmft.gcorr")

    # Write the GAMMA file for vasp
    println("Write correction for density matrix")
    vaspc_gcorr(kwin, gcorr)

    # Create a vasp.lock file to wake up the vasp
    println("Reactivate the vasp engine (vasp.lock)")
    vaspc_lock()
end

"""
    vasp_stop()

Stop the vasp engine by creating a STOPCAR file in the working folder.
It will not return until the vasp engine is completely termined.

See also: [`vasp_back`](@ref).
"""
function vasp_stop()
    # Create the STOPCAR file in the dft folder.
    vaspc_stopcar()

    # Maybe vasp.lock is necessary.
    vaspc_lock("dft")

    # Sleep until the STOPCAR is deleted automatically, which means
    # that the vasp process is termined completely.
    while true
        # Sleep two seconds
        sleep(2)

        # Check the status of STOPCAR
        !vaspq_stopcar() && break
    end
end

#=
### *Service Functions* : *Group A*
=#

#=
*Remarks* :

The parameter `NBANDS` is quite important. Sometimes its default value
is too small. The adaptor will fail to generate reasonable projectors.
At this case, you will see an error thrown by the `try_diag()` function.
The solution is quite simple, i.e., increasing `NBANDS` a bit.

The current algorithm is suitable for paramagnetic systems. But it has
not been tested for `magnetically ordered materials`.
=#

"""
    vaspc_incar(fermi::F64, sc_mode::I64)

Generate an `INCAR` file. It will be used only when the DFT engine
is vasp.

See also: [`vaspc_kpoints`](@ref).
"""
function vaspc_incar(fermi::F64, sc_mode::I64)
    # Open the iostream
    ios = open("INCAR", "w")

    # Standard part
    case = get_c("case")
    write(ios, "System   = $case \n")
    write(ios, "PREF     = Accurate \n")
    write(ios, "EDIFF    = 1E-8 \n")
    write(ios, "ALGO     = Normal \n")
    write(ios, "LASPH    = .TRUE. \n")
    write(ios, "LMAXMIX  = 6 \n")
    write(ios, "NCORE    = 4 \n")
    write(ios, "LWAVE    = .TRUE. \n")

    # Customize your INCAR according to the case.toml
    #
    # For smearing
    smear = get_d("smear")
    smear == "mp2"   && write(ios, "ISMEAR   = 2 \n")
    smear == "mp1"   && write(ios, "ISMEAR   = 1 \n")
    smear == "gauss" && write(ios, "ISMEAR   = 0 \n")
    smear == "tetra" && write(ios, "ISMEAR   =-5 \n")

    # For kmesh density
    #
    # If kmesh == "file", then vaspc_kpoints() will be used to generate
    # the KPOINTS file.
    kmesh = get_d("kmesh")
    kmesh == "accurate" && write(ios, "KSPACING = 0.1 \n")
    kmesh == "medium"   && write(ios, "KSPACING = 0.2 \n")
    kmesh == "coarse"   && write(ios, "KSPACING = 0.4 \n")

    # For magnetic moment
    magmom = get_d("magmom")
    if !isa(magmom, Missing)
        write(ios, "MAGMOM   = $magmom \n")
    end

    # For symmetry
    lsymm = get_d("lsymm")
    if lsymm
        write(ios, "ISYM     = 2 \n")
    else # Ignore the symmetry completely
        write(ios, "ISYM     =-1 \n")
    end

    # For spin polarizations
    lspins = get_d("lspins")
    if lspins
        write(ios, "ISPIN    = 2 \n")
    else
        write(ios, "ISPIN    = 1 \n")
    end

    # For spin-orbit coupling
    lspinorb = get_d("lspinorb")
    if lspinorb
        write(ios, "LSORBIT  = .TRUE. \n")
    else
        write(ios, "LSORBIT  = .FALSE. \n")
    end

    # For optimized projectors
    ewidth = 4.0 # A magic number
    write(ios, "LORBIT   = 14 \n")
    #
    emin = fermi - ewidth
    write(ios, "EMIN     = $emin \n")
    #
    emax = fermi + ewidth
    write(ios, "EMAX     = $emax \n")

    # For local orbitals and projectors
    lproj = get_d("lproj")
    sproj = get_d("sproj")
    if !isa(lproj, Missing) && !isa(sproj, Missing)
        if lproj
            for p in eachindex(sproj)
                str = sproj[p]
                write(ios, "LOCPROJ  = $str \n")
            end
        end
    end

    # For number of bands
    nbands = vaspio_nband(pwd())
    #
    # Special treatment for Ce. For Ce, the default nbands is 10, which
    # is too small. Another solution is to modify vaspio_nband().
    if nbands ≤ 10
        nbands = nbands + 8
    end
    write(ios, "NBANDS   = $nbands \n")

    # Special treatment for sc_mode == 2
    #
    # We would like to activate ICHARG = 5 to perform charge fully
    # self-consistent DFT + DMFT calculations.
    if sc_mode == 2
        write(ios, "\n")
        write(ios, "ICHARG   = 5 \n")
        write(ios, "NELM     = 1000 \n")
        write(ios, "NELMIN   = 1000 \n")
        write(ios, "NELMDL   = -8 \n")
    end

    # Close the iostream
    close(ios)
end

"""
    vaspc_kpoints(mp_scheme::Bool = true, n::I64 = 9)

Generate a valid `KPOINTS` file for vasp.

See also: [`vaspc_incar`](@ref).
"""
function vaspc_kpoints(mp_scheme::Bool = true, n::I64 = 9)
    # If the `KPOINTS` file is available, we do nothing or else we will
    # try to create a new one.
    if !isfile("KPOINTS")
        # Open the iostream
        ios = open("KPOINTS", "w")

        # Write the body
        write(ios, "Automatic K-mesh Generation \n")
        write(ios, "0 \n")
        if mp_scheme
            write(ios, "Monkhorst-Pack \n")
        else
            write(ios, "Gamma \n")
        end
        write(ios, " $n  $n  $n \n")
        write(ios, " 0  0  0 \n")

        # Close the iostream
        close(ios)
    end
end

"""
    vaspc_gcorr(kwin::Array{I64,3}, gcorr::Array{C64,4})

Generate a valid `GAMMA` file for vasp. The vasp will need this file
when `ICHARG = 5`.

See also: [`write_gcorr`](@ref), [`read_gcorr`](@ref).
"""
function vaspc_gcorr(kwin::Array{I64,3}, gcorr::Array{C64,4})
    # Extract the dimensional parameters
    _, xbnd, nkpt, nspin = size(gcorr)
    @assert nspin in (1,2) # Current limitation

    # Determine filename for correction for density matrix
    fgcorr = "GAMMA"

    # Write the data
    open(fgcorr, "w") do fout
        # Print the header. It seems it is in old format.
        # Please check the fileio.F/READGAMMA_HEAD() subroutine in
        # vasp's source code for more details.
        @printf(fout, " %i  -1  ! Number of k-points, default number of bands \n", nkpt)

        # Go through each spin
        for s = 1:nspin
            # Go through each 𝑘-point
            for k = 1:nkpt
                # Determine the band window
                bs = kwin[k,s,1] # Start of band window
                be = kwin[k,s,2] # End of band window
                cbnd = be - bs + 1

                # Sanity check
                @assert cbnd ≤ xbnd

                # Write the band window
                @printf(fout, " %i  %i  %i\n", k, bs, be)

                # Go through each band
                for p = 1:cbnd
                    for q = 1:cbnd
                        z = gcorr[p,q,k,1]
                        @printf(fout, " %.14f  %.14f", real(z), imag(z))
                    end
                    println(fout) # Create a new line
                end
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  > Write gcorr matrix into: dft/$fgcorr")
end

"""
    vaspc_stopcar()

Create the `STOPCAR` file in the dft directory to stop the vasp engine.
Vasp will stop at the next electronic step, i.e. `WAVECAR` and `CHGCAR`
might contain non-converged results.
"""
function vaspc_stopcar()
    # Create STOPCAR
    fstop = "dft/STOPCAR"
    #
    open(fstop, "w") do fout
        println(fout, "LABORT = .TRUE.")
    end
    #
    println("  > Create STOPCAR for vasp: $fstop")
end

"""
    vaspc_lock(dir::String = ".")

Create the `vasp.lock` file. This file is relevant for `ICHARG = 5`.
The vasp program runs only when the `vasp.lock` file is present in the
current directory. Its default working directory is just `dft`, but we
can specify it via argument `dir`.

See also: [`vaspq_lock`](@ref).
"""
function vaspc_lock(dir::String = ".")
    touch(joinpath(dir, "vasp.lock"))
end

#=
### *Service Functions* : *Group B*
=#

"""
    vaspq_stopcar()

Return whether the `STOPCAR` file is available. Its working directory
might be `root` or `dft`.

See also: [`vaspc_stopcar`](@ref).
"""
function vaspq_stopcar()
    return isfile("dft/STOPCAR") || isfile("STOPCAR")
end

"""
    vaspq_lock()

Return whether the `vasp.lock` file is available. Its working directory
might be `root` or `dft`.

See also: [`vaspc_lock`](@ref).
"""
function vaspq_lock()
    return isfile("dft/vasp.lock") || isfile("vasp.lock")
end

"""
    vaspq_files(f::String)

Check the essential output files by vasp. Here `f` means only the
directory that contains the desired files.

See also: [`adaptor_core`](@ref).
"""
function vaspq_files(f::String)
    # Define file list
    fl = ["POSCAR",
          "IBZKPT",
          "EIGENVAL",
          "LOCPROJ",
          "DOSCAR",
          "CHGCAR",
          "OSZICAR"]

    # Check them one by one
    for i in eachindex(fl)
        @assert isfile( joinpath(f, fl[i]) )
    end
end

"""
    vaspq_files()

Check the essential output files by vasp in the current directory.

See also: [`adaptor_core`](@ref).
"""
vaspq_files() = vaspq_files(pwd())

#=
### *Service Functions* : *Group C*
=#

#=
*Remarks* :

In vasp, the `NBANDS` parameter is determined automatically. According
to the wiki of vasp, it is equal to `nelect / 2 + latt.natom / 2` for
paramagnetic case by default. However, it may be insufficient for
generating projectors. For example, for SrVO``_3``, the default `NBANDS`
parameter is only 20. It is too small to determine the five V``_{3d}``
projectors. It would be better to increase it to 30.

Here, we increase `NBANDS` to `1.6 * NBANDS`. The `1.6` is a magic
number, you can adjust it by yourself. For magnetic ordered systems,
perhaps we nedd a larger factor.
=#

"""
    vaspio_nband(f::String)

Reading vasp's `POSCAR` and `POTCAR` files, evaluating number of bands. It
will be used to create the `INCAR` file. Here `f` means only the directory
that contains `POSCAR` and `POTCAR`.

See also: [`vaspc_incar`](@ref), [`vaspio_valence`](@ref), [`vaspio_lattice`](@ref).
"""
function vaspio_nband(f::String)
    # Extract crystallography information from `POSCAR`
    latt = vaspio_lattice(f)

    # Extract number of valence electrons from `POTCAR`
    zval = vaspio_valence(f)

    # Sanity check
    nsort, _ = size(latt.sorts)
    @assert nsort == length(zval)

    # Evaluate number of valence electrons in total
    nelect = sum(@. latt.sorts[:, 2] * zval)

    # Evaluate number of bands
    if get_d("lspins")
        nband = floor(I64, nelect * 1.2)
    else
        nband = floor(I64, (nelect / 2 + latt.natom / 2) * 1.6)
    end

    # Return the desired nband
    return nband
end

"""
    vaspio_nband()

Reading vasp's `POSCAR` and `POTCAR` files, evaluating number of bands. It
will be used to create the `INCAR` file.

See also: [`vaspc_incar`](@ref), [`vaspio_valence`](@ref), [`vaspio_lattice`](@ref).
"""
vaspio_nband() = vaspio_nband(pwd())

"""
    vaspio_valence(f::String)

Reading vasp's `POTCAR` file, return `ZVAL`. Here `f` means only the
directory that contains `POTCAR`. The information about `ZVAL` will
be used to determine `NBANDS` in the `INCAR` file.

See also: [`vaspc_incar`](@ref), [`vaspio_nband`](@ref).
"""
function vaspio_valence(f::String)
    # Open the iostream
    fin = open(joinpath(f, "POTCAR"), "r")

    # Read the POTCAR
    strs = readlines(fin)

    # Extract the lines that contains the `ZVAL`
    filter!(x -> contains(x, "ZVAL"), strs)

    # Extract ZVAL, convert it into float, than save it.
    zval = map(x -> parse(F64, line_to_array(x)[6]), strs)

    # Close the iostream
    close(fin)

    # Return the desired arrays
    return zval
end

"""
    vaspio_valence()

Reading vasp's `POTCAR` file, return `ZVAL`. The information about
`ZVAL` will be used to determine `NBANDS` in the `INCAR` file.

See also: [`vaspc_incar`](@ref), [`vaspio_nband`](@ref).
"""
vaspio_valence() = vaspio_valence(pwd())

"""
    vaspio_energy(f::String)

Reading vasp's `OSZICAR` file, return DFT total energy, which will be
used to determine the total DFT + DMFT energy. Here `f` means only the
directory that contains `OSZICAR`.
"""
function vaspio_energy(f::String)
    # Open the iostream
    fin = open(joinpath(f, "OSZICAR"), "r")

    # Read the OSZICAR
    strs = readlines(fin)

    # Extract ZVAL, convert it into float, than save it.
    etot = parse(F64, line_to_array(strs[end])[3])

    # Close the iostream
    close(fin)

    # Return the desired value
    return etot
end

"""
    vaspio_energy()

Reading vasp's `OSZICAR` file, return DFT total energy, which will be
used to determine the total DFT + DMFT energy.
"""
vaspio_energy() = vaspio_energy(pwd())

"""
    vaspio_procar(f::String)

Reading vasp's `PROCAR` file, extract orbital weight information. Here `f`
means only the directory that contains `PROCAR`.

This function is not invoked directly during the DFT + DMFT iteration. It
is designed for users merely. They can use it to judge which orbitals are
the most relevant, and then apply the obtained information to customize
their case.toml configuration file (specifically, the `window` parameter
in the `dft` block).

!!! note

    This function is only called by tools/analyze.jl or used in REPL.
"""
function vaspio_procar(f::String)
    # Open the iostream
    fin = open(joinpath(f, "PROCAR"), "r")

    # We have to determine the dimensional parameters at first
    #
    # (1) Determine key parameters: nkpt, nband, and natom.
    readline(fin)
    arr = line_to_array(fin)
    nkpt = parse(I64, arr[4])
    nband = parse(I64, arr[8])
    natom = parse(I64, arr[12])
    seekstart(fin) # Rewind the stream
    #
    # (2) Determine key parameters: norbs.
    readuntil(fin, "ion ")
    arr = line_to_array(fin)
    norbs = length(arr) - 1
    @assert norbs == 9 || norbs == 16
    seekstart(fin) # Rewind the stream
    #
    # (3) Determine key parameters: nspin.
    nspin = 0
    for line in eachline(fin)
        if contains(line, "k-points")
            nspin = nspin + 1
        end
    end
    seekstart(fin) # Rewind the stream
    #
    # (4) Determine whether spin-orbit coupling is activated.
    readuntil(fin, "ion ")
    soc = false
    if natom == 1
        lc = 0
        for line in eachline(fin)
            lc = lc + 1
            if contains(line, "ion ")
                break
            end
        end

        if lc == 3
            soc = false
        elseif ls == 6
            soc = true
        else
            error("Something wrong in PROCAR")
        end
    else # natom > 1
        lc = 0
        for line in eachline(fin)
            if contains(line, "tot ")
                lc = lc + 1
            end
            if contains(line, "ion ")
                break
            end
        end

        if lc == 1
            soc = false
        elseif lc == 4
            soc = true
        else
            error("Something wrong in PROCAR")
        end
    end
    seekstart(fin) # Rewind the stream
    #
    # (5) Additional check. `soc` must be compatible with `nspin`.
    if soc
        @assert nspin == 1
    end
    #
    # (6) Debug
    print("Parameters: ")
    @show norbs, natom, nband, nkpt, nspin, soc

    # Build `str`, which is used to tell the parser how to skip the
    # blank lines for various cases.
    #
    # (1) Init the string
    fstr = ""
    #
    # (2) Single atom vs. multiple atoms
    if natom == 1
        fstr = fstr * "1"
    else
        fstr = fstr * "2"
    end
    #
    # (3) d system vs. f system
    if norbs == 9
        fstr = fstr * "d"
    else
        fstr = fstr * "f"
    end
    #
    # (4) SOC or not
    if soc
        fstr = fstr * "s"
    else
        fstr = fstr * "n"
    end
    #
    # (5) Eight valid cases
    @assert fstr in ["1ds", "1dn", "1fs", "1fn", "2ds", "2dn", "2fs", "2fn"]
    println("PROCAR format: $fstr")

    # Create arrays
    # The `worb` is used to save the raw data, while `oab` is used to
    # save the processed data. The `enk` denotes the eigenvalues and
    # `occ` denotes the occupations.
    worb = zeros(F64, norbs, natom, nband, nkpt, nspin)
    oab = zeros(F64, norbs, natom, nband, nspin)
    enk = zeros(F64, nband, nkpt, nspin)
    occ = zeros(F64, nband, nkpt, nspin)

    # Start to parse the `PROCAR` file
    println("Parsing PROCAR...")

    # Go through each spin orientation
    for s = 1:nspin
        # Blank line
        readline(fin)

        # Go through each k-point
        for k = 1:nkpt
            # Blank lines
            readline(fin)
            readline(fin)

            # Check k-index
            arr = line_to_array(fin)
            @assert k == parse(I64, arr[2])
            println("Finishing spin $s k-point $k")

            # Go through each band for the given k-point and spin
            for b = 1:nband
                # Blank line
                readline(fin)

                # Check band index
                arr = line_to_array(fin)
                @assert b == parse(I64, arr[2])

                # Parse eigenvalues and occupations
                enk[b, k, s] = parse(F64, arr[5])
                occ[b, k, s] = parse(F64, arr[8])

                # Blank lines
                readline(fin)
                readline(fin)

                # Parse weights
                # Go through each atom and orbital
                for a = 1:natom
                    arr = line_to_array(fin)
                    @assert parse(I64, arr[1]) == a
                    @assert norbs == length(arr) - 2
                    worb[:, a, b, k, s] = parse.(F64, arr[2:end-1])
                end

                # Skip useless lines
                @cswitch fstr begin
                    @case "1ds"
                        foreach(x -> readline(fin), 1:1:3)
                        break

                    @case "1dn"
                        break

                    @case "1fs"
                        foreach(x -> readline(fin), 1:1:3)
                        break

                    @case "1fn"
                        break

                    @case "2ds"
                        foreach(x -> readline(fin), 1:1:(3*(1+natom)+1))
                        break

                    @case "2dn"
                        readline(fin)
                        break

                    @case "2fs"
                        foreach(x -> readline(fin), 1:1:(3*(2+natom)+2))
                        break

                    @case "2fn"
                        readline(fin)
                        readline(fin)
                        break

                    @default
                        sorry()
                        break
                end # @cswitch

                # Skip blank line which starts with "ion ..."
                readline(fin)

                # Parse phase factors
                # Go through each atom and orbital
                for a = 1:natom
                    readline(fin)
                end

                # Skip useless lines
                @cswitch fstr begin
                    @case "1ds"
                        break

                    @case "1dn"
                        break

                    @case "1fs"
                        break

                    @case "1fn"
                        break

                    @case "2ds"
                        readline(fin)
                        break

                    @case "2dn"
                        readline(fin)
                        break

                    @case "2fs"
                        readline(fin)
                        readline(fin)
                        break

                    @case "2fn"
                        readline(fin)
                        readline(fin)
                        break

                    @default
                        sorry()
                        break
                end # @cswitch
            end # END OF B LOOP
        end # END OF K LOOP
    end # END OF S LOOP

    # Close the iostream
    close(fin)

    # Try to build `oab` by k-summation
    for s = 1:nspin
        for b = 1:nband
            for a = 1:natom
                for o = 1:norbs
                    oab[o, a, b, s] = sum(worb[o, a, b, :, s]) / float(nkpt)
                end
            end
        end # END OF B LOOP
    end # END OF S LOOP

    # Return the desired arrays
    return oab, enk, occ
end

"""
    vaspio_procar()

Reading vasp's `PROCAR` file, extract orbital weight information.

This function is not invoked directly during the DFT + DMFT iteration. It
is designed for users merely. They can use it to judge which orbitals are
the most relevant, and then apply the obtained information to customize
their case.toml configuration file (specifically, the `window` parameter
in the `dft` block).

!!! note

    This function is only called by tools/analyze.jl or used in REPL.
"""
vaspio_procar() = vaspio_procar(pwd())

"""
    vaspio_lattice(f::String, silent::Bool = true)

Reading vasp's `POSCAR` file, return crystallography information. Here `f`
means only the directory that contains `POSCAR`.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
function vaspio_lattice(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse lattice")
    !silent && println("  > Open and read POSCAR")

    # Open the iostream
    fin = open(joinpath(f, "POSCAR"), "r")

    # Get the case
    _case = string(strip(readline(fin)))

    # Get the scaling factor
    scale = parse(F64, readline(fin))

    # Get the lattice vectors
    lvect = zeros(F64, 3, 3)
    lvect[1, :] = parse.(F64, line_to_array(fin))
    lvect[2, :] = parse.(F64, line_to_array(fin))
    lvect[3, :] = parse.(F64, line_to_array(fin))

    # Get the symbol list
    symbols = line_to_array(fin)

    # Get the number of sorts of atoms
    nsort = length(symbols)

    # Get the number list
    numbers = parse.(I64, line_to_array(fin))

    # Get the total number of atoms
    natom = sum(numbers)

    # Now all the parameters are ready, we would like to create
    # a `Lattice` struct here.
    latt = Lattice(_case, scale, nsort, natom)

    # Update latt using the available data
    latt.lvect = lvect
    for i = 1:nsort
        latt.sorts[i, 1] = string(symbols[i])
        latt.sorts[i, 2] = numbers[i]
    end

    # Get the atom list
    k = 0
    for i = 1:nsort
        for j = 1:numbers[i]
            k = k + 1
            latt.atoms[k] = symbols[i]
        end
    end
    # Sanity check
    @assert k == natom

    # Get the coordinates of atoms
    readline(fin)
    for i = 1:natom
        latt.coord[i, :] = parse.(F64, line_to_array(fin)[1:3])
    end

    # Close the iostream
    close(fin)

    # Print some useful information to check
    !silent && println("  > System: ", latt._case)
    !silent && println("  > Atoms: ", latt.atoms)

    # Return the desired struct
    return latt
end

"""
    vaspio_lattice()

Reading vasp's `POSCAR` file, return crystallography information.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
vaspio_lattice() = vaspio_lattice(pwd())

"""
    vaspio_kmesh(f::String)

Reading vasp's `IBZKPT` file, return `kmesh` and `weight`. Here `f` means
only the directory that contains `IBZKPT`.

See also: [`vaspio_tetra`](@ref), [`irio_kmesh`](@ref).
"""
function vaspio_kmesh(f::String)
    # Print the header
    println("Parse kmesh and weight")
    println("  > Open and read IBZKPT")

    # Open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # Extract number of 𝑘-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # Create arrays
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt)

    # Read in the 𝑘-points and their weights
    for i = 1:nkpt
        arr = parse.(F64, line_to_array(fin))
        kmesh[i, 1:3] = arr[1:3]
        weight[i] = arr[4]
    end

    # Close the iostream
    close(fin)

    # Print some useful information to check
    println("  > Number of k-points: ", nkpt)
    println("  > Total sum of weights: ", sum(weight))
    println("  > Shape of Array kmesh: ", size(kmesh))
    println("  > Shape of Array weight: ", size(weight))

    # Return the desired arrays
    return kmesh, weight
end

"""
    vaspio_kmesh()

Reading vasp's `IBZKPT` file, return `kmesh` and `weight`.

See also: [`vaspio_tetra`](@ref), [`irio_kmesh`](@ref).
"""
vaspio_kmesh() = vaspio_kmesh(pwd())

"""
    vaspio_tetra(f::String)

Reading vasp's `IBZKPT` file, return tetrahedra information. Here `f`
means only the directory that contains `IBZKPT`.

See also: [`vaspio_kmesh`](@ref), [`irio_tetra`](@ref).
"""
function vaspio_tetra(f::String)
    # Print the header
    println("Parse tetrahedron")
    println("  > Open and read IBZKPT")

    # Open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # Extract number of 𝑘-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # Read in the 𝑘-points and their weights
    #
    # Skip nkpt lines
    for i = 1:nkpt
        readline(fin)
    end

    # Read in the tetrahedra information
    #
    # Skip one empty line
    readline(fin)

    # Extract total number of tetrahedra and volume of a tetrahedron
    arr = line_to_array(fin)
    ntet = parse(I64, arr[1])
    volt = parse(F64, arr[2])

    # Create arrays
    itet = zeros(I64, ntet, 5)

    # Parse the input tetrahedra information
    for t = 1:ntet
        itet[t, :] = parse.(I64, line_to_array(fin))
    end

    # Close the iostream
    close(fin)

    # Print some useful information to check
    println("  > Number of tetrahedra: ", ntet)
    println("  > Volume of one tetrahedron: ", volt)
    println("  > Shape of Array itet: ", size(itet))

    # Return the desired arrays
    return volt, itet
end

"""
    vaspio_tetra()

Reading vasp's `IBZKPT` file, return tetrahedra information.

See also: [`vaspio_kmesh`](@ref), [`irio_tetra`](@ref).
"""
vaspio_tetra() = vaspio_tetra(pwd())

"""
    vaspio_eigen(f::String)

Reading vasp's `EIGENVAL` file, return energy band information. Here `f`
means only the directory that contains `EIGENVAL`.

Sometimes the `EIGENVAL` file does not contain any useful data, then we
turn to the `LOCPROJ` file to obtain the energy band information.

See also: [`irio_eigen`](@ref).
"""
function vaspio_eigen(f::String)
    # Print the header
    println("Parse enk and occupy")

    # Check whether the `EIGENVAL` file contains valid data
    lines = readlines(joinpath(f, "EIGENVAL"))

    # Read EIGENVAL
    if length(lines) ≥ 10
        println("  > Open and read EIGENVAL")

        # Open the iostream
        fin = open(joinpath(f, "EIGENVAL"), "r")

        # Determine number of spins
        nspin = parse(I64, line_to_array(fin)[end])
        @assert nspin == 1 || nspin == 2

        # Skip for lines
        for i = 1:4
            readline(fin)
        end

        # Read in some key parameters: nelect, nkpt, nbands
        _, nkpt, nband = parse.(I64, line_to_array(fin))

        # Create arrays
        enk = zeros(F64, nband, nkpt, nspin)
        occupy = zeros(F64, nband, nkpt, nspin)

        # Read in the energy bands and the corresponding occupations
        for i = 1:nkpt
            readline(fin)
            readline(fin)
            for j = 1:nband
                arr = line_to_array(fin)
                for s = 1:nspin
                    enk[j, i, s] = parse(F64, arr[s+1])
                    occupy[j, i, s] = parse(F64, arr[s+1+nspin])
                end # END OF S LOOP
            end # END OF J LOOP
        end # END OF I LOOP

        # close the iostream
        close(fin)

        # Print some useful information to check
        println("  > Number of DFT bands: ", nband)
        println("  > Number of k-points: ", nkpt)
        println("  > Number of spins: ", nspin)
        println("  > Shape of Array enk: ", size(enk))
        println("  > Shape of Array occupy: ", size(occupy))

        # return the desired arrays
        return enk, occupy

    # Read LOCPROJ
    else
        println("  > Open and read LOCPROJ")

        # Open the iostream
        fin = open(joinpath(f, "LOCPROJ"), "r")

        # Extract number of spins (nspin), number of k-points (nkpt),
        # number of bands (nband), and number of projectors (nproj).
        nspin, nkpt, nband, nproj = parse.(I64, line_to_array(fin)[1:4])
        @assert nspin == 1 || nspin == 2

        #@show nspin, nkpt, nband, nproj
        for i = 1:nproj
            readline(fin)
        end

        # Create arrays
        enk = zeros(F64, nband, nkpt, nspin)
        occupy = zeros(F64, nband, nkpt, nspin)

        # Read in the energy bands and the corresponding occupations
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    readline(fin)
                    arr = line_to_array(fin)
                    _s, _k, _b = parse.(I64, arr[2:4])
                    @assert _s == s
                    @assert _k == k
                    @assert _b == b

                    enk[b,k,s] = parse(F64, arr[5])
                    occupy[b,k,s] = parse(F64,arr[6])
                    for p = 1:nproj
                        readline(fin)
                    end
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP

        # Close the iostream
        close(fin)

        # Print some useful information to check
        println("  > Number of DFT bands: ", nband)
        println("  > Number of k-points: ", nkpt)
        println("  > Number of spins: ", nspin)
        println("  > Shape of Array enk: ", size(enk))
        println("  > Shape of Array occupy: ", size(occupy))

        # return the desired arrays
        return enk, occupy
    end
end

"""
    vaspio_eigen()

Reading vasp's `EIGENVAL` file, return energy band information.

See also: [`irio_eigen`](@ref).
"""
vaspio_eigen() = vaspio_eigen(pwd())

"""
    vaspio_projs(f::String)

Reading vasp's `LOCPROJ` file, return raw projector matrix. Here `f` means
only the directory that contains `LOCPROJ`.

See also: [`irio_projs`](@ref).
"""
function vaspio_projs(f::String)
    # Print the header
    println("Parse projector")
    println("  > Open and read LOCPROJ")

    # Open the iostream
    fin = open(joinpath(f, "LOCPROJ"), "r")

    # Extract number of spins (nspin), number of k-points (nkpt),
    # number of bands (nband), and number of projectors (nproj).
    nspin, nkpt, nband, nproj = parse.(I64, line_to_array(fin)[1:4])
    @assert nspin == 1 || nspin == 2

    # Extract raw information about projectors
    sites = zeros(I64, nproj)
    descs = fill("", nproj)
    for i = 1:nproj
        arr = line_to_array(fin)
        sites[i] = parse(I64, arr[2])
        # Get rid of the ":" char
        descs[i] = replace(arr[end], ":" => "")
    end

    # Try to build PrTrait struct. The raw information about projectors
    # should be encapsulated in it.
    PT = PrTrait[]
    for i = 1:nproj
        push!(PT, PrTrait(sites[i], descs[i]))
    end

    # Try to split these projectors into groups.
    #
    # At first, we collect the tuple (site, l) for all projectors.
    site_l = Tuple[]
    for i in eachindex(PT)
        push!(site_l, (PT[i].site, PT[i].l))
    end
    #
    # Second, we figure out the unique (site, l) pairs. Here, we use the
    # union! function, which is much more efficient than the unique! one.
    union!(site_l)
    #
    # Third, we create a array of PrGroup struct (except for site and
    # l, most of its member variables need to be corrected). Note
    # that for a given PrGroup, the projectors indexed by PrGroup.Pr
    # should share the same site and l.
    PG = PrGroup[]
    for i in eachindex(site_l)
        push!(PG, PrGroup(site_l[i]...))
    end
    #
    # Fourth, for each PrGroup, we scan all of the projectors to find
    # out those with correct site and l; record their indices; and
    # save them at PrGroup.Pr array.
    for i in eachindex(PG)
        site, l = PG[i].site, PG[i].l
        PG[i].Pr = findall(x -> (x.site, x.l) == (site, l), PT)
    end
    #
    # Finally, check correctness
    @assert nproj == sum(x -> length(x.Pr), PG)

    # Create arrays
    chipsi = zeros(C64, nproj, nband, nkpt, nspin)

    # Read in raw projector data
    for s = 1:nspin
        for k = 1:nkpt
            for b = 1:nband
                # Skip two empty lines
                readline(fin)
                readline(fin)

                # Parse the data line by line
                for p = 1:nproj
                    chipsi[p, b, k, s] = line_to_cmplx(fin)
                end
            end # END OF B LOOP
        end # END OF K LOOP
    end # END OF S LOOP

    # Close the iostream
    close(fin)

    # Print some useful information to check
    println("  > Number of DFT bands: ", nband)
    println("  > Number of k-points: ", nkpt)
    println("  > Number of spins: ", nspin)
    println("  > Number of projectors: ", nproj)
    for i in eachindex(PT)
        if nproj ≥ 10
            @printf("    [ Trait %2i ]", i)
        else
            @printf("    [ Trait %1i ]", i)
        end
        print("  site -> ", PT[i].site)
        print("  l -> ", PT[i].l)
        print("  m -> ", PT[i].m)
        println("  desc -> ", PT[i].desc)
    end
    #
    println("  > Number of groups: ", length(PG))
    for i in eachindex(PG)
        print("    [ Group $i ]")
        print("  site -> ", PG[i].site)
        print("  l -> ", PG[i].l)
        print("  corr -> ", PG[i].corr)
        print("  shell -> ", PG[i].shell)
        println("  Pr -> ", PG[i].Pr)
    end
    println("  > Shape of Array chipsi: ", size(chipsi))

    # Return the desired arrays
    # Note: PG should be further setup at plo_group() function.
    return PT, PG, chipsi
end

"""
    vaspio_projs()

Reading vasp's `LOCPROJ` file, return raw projector matrix.

See also: [`irio_projs`](@ref).
"""
vaspio_projs() = vaspio_projs(pwd())

"""
    vaspio_fermi(f::String, silent::Bool = true)

Reading vasp's `DOSCAR` file, return the fermi level. Here `f` means
only the directory that contains `DOSCAR`.

Sometimes the `DOSCAR` file does not contain the necessary data. Thus
we have to turn to the `LOCPROJ` file.

See also: [`irio_fermi`](@ref).
"""
function vaspio_fermi(f::String, silent::Bool = true)
    # Print the header
    !silent && println("Parse fermi level")

    # Try to figure out whether the DOSCAR file is valid
    lines = readlines(joinpath(f, "DOSCAR"))

    # Read DOSCAR
    if length(lines) ≥ 6
        # Print the header
        !silent && println("  > Open and read DOSCAR")

        # Open the iostream
        fin = open(joinpath(f, "DOSCAR"), "r")

        # Skip five empty lines
        for i = 1:5
            readline(fin)
        end

        # Extract the fermi level
        fermi = parse(F64, line_to_array(fin)[4])

        # Close the iostream
        close(fin)

        # Print some useful information to check
        !silent && println("  > Fermi level: $fermi eV")

        # Return the desired data
        return fermi

    # Read LOCPROJ
    else
        # Print the header
        !silent && println("  > Open and read LOCPROJ")

        # Open the iostream
        fin = open(joinpath(f, "LOCPROJ"), "r")

        # Extract the fermi level
        fermi = parse(F64, line_to_array(fin)[5])

        # Close the iostream
        close(fin)

        # Print some useful information to check
        !silent && println("  > Fermi level: $fermi eV")

        # Return the desired data
        return fermi
    end
end

"""
    vaspio_fermi()

Reading vasp's `DOSCAR` file, return the fermi level.

See also: [`irio_fermi`](@ref).
"""
vaspio_fermi() = vaspio_fermi(pwd())

"""
    vaspio_charge(f::String)

Reading vasp's `CHGCAR` file, return the charge density. Here `f` means
only the directory that contains `CHGCAR`.

See also: [`irio_charge`](@ref).
"""
function vaspio_charge(f::String) end

"""
    vaspio_charge()

Reading vasp's `CHGCAR` file, return the charge density.

See also: [`irio_charge`](@ref).
"""
vaspio_charge() = vaspio_charge(pwd())
