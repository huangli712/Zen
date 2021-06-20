#
# Project : Pansy
# Source  : dmft.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/21
#

#
# Driver Functions
#

"""
    dmft_init(it::IterInfo, task::I64)

Initialize the dynamical mean-field theory engine. Prepare the necessary
files, and generate the configuration file.

See also: [`dmft_exec`](@ref), [`dmft_save`](@ref).
"""
function dmft_init(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Print the header
    println("Engine : DMFT$(subscript(task))")
    println("Try to solve the dynamical mean-field self-consistent equation")
    println("Current directory: ", pwd())
    println("Prepare necessary input files for dmft")

    # Well, determine which files are necessary. They are defined in
    # `fsig`, `fir1`, `fir2`, and `fdmft`.
    #
    # Self-energy functions and their double counting terms
    fsig = ["sigma.bare", "sigma.dc"]
    #
    # Parameter sets within the IR format
    fir1 = ["params.ir", "maps.ir", "groups.ir", "windows.ir"]
    #
    # Kohn-Sham data within the IR format
    fir2 = ["lattice.ir", "kmesh.ir", "eigen.ir", "projs.ir"]
    #
    # Tetrahedron data are available
    if get_d("smear") === "tetra"
        push!(fir2, "tetra.ir")
    end
    #
    # Main control file for the DMFT engine
    fdmft = ("dmft.in")

    # Next, we have to copy Kohn-Sham data from `dft` to `dmft1`.
    foreach( x ->
        begin
            file_src = joinpath("../dft", x)
            file_dst = x
            cp(file_src, file_dst, force = true)
        end,
    union(fir1, fir2) )

    # Extract key parameters, which should be written into the `dmft.in`.
    axis = get_m("axis")
    beta = get_m("beta")
    mc = ( get_m("mc") isa Missing ? 0.001 : get_m("mc") )
    lfermi = ( get_m("lfermi") isa Missing ? true : get_m("lfermi") )
    ltetra = ( get_d("smear") === "tetra" )

    # Generate essential input files, such as dmft.in, dynamically.
    # If the `dmft.in` file exists already, it will be overwritten.
    open("dmft.in", "w") do fout
        println(fout, "task = $task")
        println(fout, "axis = $axis")
        println(fout, "beta = $beta")
        println(fout, "mc   = $mc")
        println(fout, "lfermi = $lfermi")
        println(fout, "ltetra = $ltetra")
    end

    # Check essential input files
    flist = (fdmft, fsig..., fir1..., fir2...)
    for i in eachindex(flist)
        filename = flist[i]
        if !isfile(filename)
            error("Please make sure the file $filename is available")
        end
        println("  > $filename is ready")
    end
end

"""
    dmft_exec(it::IterInfo, task::I64)

Execute the dynamical mean-field theory engine.

See also: [`dmft_init`](@ref), [`dmft_save`](@ref).
"""
function dmft_exec(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Print the header
    println("Detect the runtime environment for dmft")

    # Determine mpi prefix (whether the dmft is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dmft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  > Using $numproc processors (MPI)")

    # Get the home directory of DMFT engine
    dmft_home = query_dmft()
    println("  > Home directory for dmft: ", dmft_home)

    # Select suitable dmft program
    dmft_exe = "$dmft_home/dmft"
    @assert isfile(dmft_exe)

    # Assemble command
    if isnothing(mpi_prefix)
        dmft_cmd = dmft_exe
    else
        dmft_cmd = split("$mpi_prefix $dmft_exe", " ")
    end
    println("  > Assemble command: $(prod(x -> x * ' ', dmft_cmd))")

    # Print the header
    println("Launch the computational engine (dmft)")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$dmft_cmd`, stdout = "dmft.out"))
    end
    println("  > Create a task")

    # Launch it, the terminal output is redirected to dmft.out.
    # Note that the task runs asynchronously. It will not block
    # the execution.
    schedule(t)
    println("  > Add the task to the scheduler's queue")
    println("  > Waiting ...")

    # Analyze the dmft.out file during the calculation
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

        # Parse dmft.out file
        lines = readlines("dmft.out")
        filter!(x -> contains(x, "Task"), lines)

        # Figure out the task that is doing
        if length(lines) > 0
            arr = line_to_array(lines[end])
            job = arr[5]
        else # Nothing
            job = "unknown"
        end

        # Print the log to screen
        @printf("  > Elapsed %4i seconds, current task: %s\r", 5*c, job)

        # Break the loop
        istaskdone(t) && break
    end
    #
    # Keep the last output
    println()

    # Wait for the dmft task to finish
    wait(t)

    # Extract how many iterations are executed
    lines = readlines("dmft.out")
    filter!(x -> contains(x, "Task"), lines)
    println("  > Finished after $(length(lines)) tasks")
end

"""
    dmft_save(it::IterInfo, task::I64)

Backup the files outputed by the dynamical mean-field theory engine.

See also: [`dmft_init`](@ref), [`dmft_exec`](@ref).
"""
function dmft_save(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Print the header
    println("Finalize the computational task")

    # Create a list of files that need to be backup
    fdmf1 = ["dmft.out"]
    fdmf2 = ["dmft.fermi"]
    fdmf3 = ["dmft.eimps", "dmft.eimpx"]
    fdmf4 = ["dmft.green", "dmft.delta", "dmft.weiss"]
    fdmf5 = ["dmft.gamma"]

    # Be careful, the final file list depends on the task
    if task == 1
        file_list = union(fdmf1, fdmf2, fdmf3, fdmf4)
    else
        file_list = union(fdmf1, fdmf2, fdmf5)
    end

    # Store the data files
    for i in eachindex(file_list)
        f = file_list[i]
        if task == 1
            cp(f, "$f.$(it.I₃).$(it.I₁)", force = true)
        else
            cp(f, "$f.$(it.I₃).$(it.I₂)", force = true)
        end
    end
    println("  > Save the key output files")

    # Extract the fermi level, and use it to update the IterInfo struct.
    fermi = read_fermi()
    task == 1 ? it.μ₁ = fermi : it.μ₂ = fermi
    println("  > Extract the fermi level from dmft.fermi: $fermi eV")

    # Print the footer for a better visualization
    println()
end

#
# Service Functions: For I/O Operations
#

"""
    read_fermi()

Parse the dmft1/dmft.fermi file to extract the chemical potential.

See also: [`dmft_save`](@ref).
"""
function read_fermi()
    fname = "dmft.fermi"
    fermi = 0.0

    if isfile(fname)
        str = readline("dmft.fermi")
        fermi = parse(F64, line_to_array(str)[3])
    end

    return fermi
end

#
# Service Functions: For I/O Operations
#

"""
    read_delta(imp::Impurity)

Extract hybridization functions from `impurity.i/dmft.delta` file, which
is generated by `sigma_split()`. These data are essential for quantum
impurity solvers. The working directory of this function must be equal
to `impurity.i`.

The frequency mesh is also extracted in this function.

See also: [`Impurity`](@ref), [`read_eimpx`](@ref).
"""
function read_delta(imp::Impurity)
    # Get index of the quantum impurity problem
    index = imp.index

    # Get number of orbitals for the quantum impurity problem
    nband = imp.nband

    # Examine the current directory
    dirname = basename(pwd())
    dirvect = split(dirname, ".")
    @assert dirvect[1] == "impurity"
    @assert parse(I64, dirvect[2]) == index

    # Declare empty arrays for frequency mesh and hybridization functions
    fmesh = nothing
    Delta = nothing

    # Parse the `impurity.i/dmft.delta` file
    open("dmft.delta", "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite ≥ 1
        @assert qdim ≥ nband

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, nband, nband, nmesh, nspin)

        # Go through each spin orientation
        for s = 1:nspin
            # Analyze the important parameters
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            _d = parse(I64, line_to_array(strs)[7])
            @assert _t == index
            @assert _s == s
            @assert _d == nband

            # Go through each frequency point
            for m = 1:nmesh
                # Extract frequency
                fmesh[m] = parse(F64, line_to_array(fin)[3])

                # Parse the hybridization functions
                for q = 1:nband
                    for p = 1:nband
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Delta[p,q,m,s] = _re + _im * im
                    end
                end
            end # END OF M LOOP

            # Skip two lines
            readline(fin)
            readline(fin)
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Return the desired arrays
    return fmesh, Delta
end

"""
    read_delta(ai::Array{Impurity,1}, fhyb::String = "dmft1/dmft.delta")

Read the `dmft1/dmft.delta` file, extract the hybridization functions Δ
and the corresponding frequency mesh ω. The working directory of this
function must be the root folder.

See also: [`sigma_split`](@ref), [`read_eimpx`](@ref).
"""
function read_delta(ai::Array{Impurity,1}, fhyb::String = "dmft1/dmft.delta")
    # Declare the frequency mesh and hybridization function
    fmesh = []
    Delta = []

    # Make sure the existence of hybridization functions
    @assert isfile(fhyb)

    # Parse `fhyb`, extract the hybridization functions
    open(fhyb, "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite == length(ai)

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, qdim, qdim, nmesh, nspin, nsite)

        # Read the data
        # Go through each quantum impurity site and spin
        for t = 1:nsite
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                _d = parse(I64, line_to_array(strs)[7])
                @assert _t == t
                @assert _s == s
                @assert _d == ai[t].nband

                # Go through each frequency point
                for m = 1:nmesh
                    # Parse frequency mesh
                    fmesh[m] = parse(F64, line_to_array(fin)[3])
                    # Parse hybridization functions
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            _re, _im = parse.(F64, line_to_array(fin)[3:4])
                            Delta[p,q,m,s,t] = _re + _im * im
                        end
                    end
                end # END OF M LOOP

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF S LOOP
        end # END OF T LOOP
    end # END OF IOSTREAM
    println("  > Read hybridization functions from: $fhyb")
    println("  > Shape of Array fmesh: ", size(fmesh))
    println("  > Shape of Array Delta: ", size(Delta))

    # Return the desired arrays
    return fmesh, Delta
end

#
# Service Functions: For I/O Operations
#

"""
    read_eimpx(imp::Impurity)

Extract local impurity levels from `impurity.i/dmft.eimpx` file, which
is generated by `sigma_split()`. These data are essential for quantum
impurity solvers. The working directory for this function must be equal
to `impurity.i`.

See also: [`Impurity`](@ref), [`read_delta`](@ref).
"""
function read_eimpx(imp::Impurity)
    # Get index of the quantum impurity problem
    index = imp.index

    # Get number of orbitals for the quantum impurity problem
    nband = imp.nband

    # Examine the current directory
    dirname = basename(pwd())
    dirvect = split(dirname, ".")
    @assert dirvect[1] == "impurity"
    @assert parse(I64, dirvect[2]) == index

    # Declare an empty array for local impurity levels
    Eimpx = nothing

    # Parse the `dmft.eimpx` file
    open("dmft.eimpx", "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite ≥ 1
        @assert qdim ≥ nband

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, nband, nband, nspin)

        # Go through each spin orientation
        for s = 1:nspin
            # Analyze the important parameters
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            _d = parse(I64, line_to_array(strs)[7])
            @assert _t == index
            @assert _s == s
            @assert _d == nband

            # Parse local impurity levels
            for q = 1:nband
                for p = 1:nband
                    _re, _im = parse.(F64, line_to_array(fin)[3:4])
                    Eimpx[p,q,s] = _re + _im * im
                end
            end

            # Skip two lines
            readline(fin)
            readline(fin)
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Return the desired array
    return Eimpx
end

"""
    read_eimpx(ai::Array{Impurity,1}, flev::String = "dmft1/dmft.eimpx")

Read the `dmft1/dmft.eimpx` file, extract the local impurity levels εᵢ. The
working directory of this function must be the root folder.

See also: [`sigma_split`](@ref), [`read_delta`](@ref).
"""
function read_eimpx(ai::Array{Impurity,1}, flev::String = "dmft1/dmft.eimpx")
    # Declare the local impurity levels
    Eimpx = []

    # Make sure the existence of local impurity levels
    @assert isfile(flev)

    # Parse `flev`, extract the local impurity levels
    open(flev, "r") do fin
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim  = parse(I64, line_to_array(fin)[4])
        @assert nsite == length(ai)

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, qdim, qdim, nspin, nsite)

        # Read the data
        # Go through each quantum impurity site and spin
        for t = 1:nsite
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                _d = parse(I64, line_to_array(strs)[7])
                @assert _t == t
                @assert _s == s
                @assert _d == ai[t].nband

                # Parse local impurity levels
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Eimpx[p,q,s,t] = _re + _im * im
                    end
                end

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF S LOOP
        end # END OF T LOOP
    end # END OF IOSTREAM

    println("  > Read local impurity levels from: $flev")
    println("  > Shape of Array Eimpx: ", size(Eimpx))

    # Return the desired array
    return Eimpx
end

#
# Service Functions: For I/O Operations
#

"""
    read_gamma()

Read the `dmft2/dmft.gamma` file. It contains the correction for density
matrix which is from electronic correlation.

This function also return the 𝑘-mesh, which is useful for mixing the Γ
matrix with the `Kerker` algorithm.

See also: [`write_gamma`](@ref).
"""
function read_gamma()
    # Declare the arrays for 𝑘-mesh and correction for density matrix
    kmesh = nothing
    kwin = nothing
    gamma = nothing

    # Make sure the data file is available
    fgamma = "dmft2/dmft.gamma"
    @assert isfile(fgamma)

    # Parse `fgamma`, extract 𝑘-mesh and correction for density matrix
    open(fgamma, "r") do fin
        # Get the dimensional parameters
        nkpt  = parse(I64, line_to_array(fin)[4])
        nspin = parse(I64, line_to_array(fin)[3])
        qbnd  = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create arrays
        kmesh = zeros(F64, nkpt, 3)
        kwin = zeros(I64, nkpt, nspin, 2)
        gamma = zeros(C64, qbnd, qbnd, nkpt, nspin)

        # Read the data
        # Go through each spin and 𝑘-point
        for s = 1:nspin
            for k = 1:nkpt
                # Parse indices and dimensional parameter
                #
                # For spin
                strs = readline(fin)
                _s = parse(I64, line_to_array(strs)[3])
                @assert _s == s
                #
                # For 𝑘-point
                strs = readline(fin)
                _k = parse(I64, line_to_array(strs)[3])
                @assert _k == k
                kmesh[k,1:3] = parse.(F64, line_to_array(strs)[4:6])
                #
                # For band window
                strs = readline(fin)
                cbnd = parse(I64, line_to_array(strs)[3])
                @assert cbnd ≤ qbnd
                bs = parse(I64, line_to_array(strs)[5])
                be = parse(I64, line_to_array(strs)[7])
                @assert bs ≤ be
                @assert be - bs + 1 == cbnd
                kwin[k,s,1] = bs
                kwin[k,s,2] = be

                # Parse Γ matrix
                for q = 1:cbnd
                    for p = 1:cbnd
                        strs = readline(fin)
                        _re = parse(F64, line_to_array(strs)[3])
                        _im = parse(F64, line_to_array(strs)[4])
                        gamma[p,q,k,s] = _re + _im*im
                    end
                end

                # Skip two lines
                readline(fin)
                readline(fin)
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM
    println("  Read gamma matrix from: $fgamma")

    # Return the desired arrays
    return kmesh, kwin, gamma
end

#
# Service Functions: For I/O Operations
#

"""
    write_delta(fmesh::Array{F64,1}, Delta::Array{C64,5}, ai::Array{Impurity,1})

Split hybridization functions and the corresponding frequency mesh into
the `impurity.i/dmft.delta` file, which is essential for the quantum
impurity solver. The working directory of this function must be the
root folder.

See also: [`Impurity`](@ref), [`read_delta`](@ref), [`write_eimpx`](@ref).
"""
function write_delta(fmesh::Array{F64,1}, Delta::Array{C64,5}, ai::Array{Impurity,1})
   # Extract the dimensional parameters
    _, qdim, nmesh, nspin, nsite = size(Delta)

    # Go through each quantum impurity problem
    for t = 1:nsite
        # Determine filename for hybridization functions
        fhyb = "impurity.$t/dmft.delta"

        # Write the data
        open(fhyb, "w") do fout
            # Write dimensional parameters
            @printf(fout, "# nsite: %4i\n", nsite)
            @printf(fout, "# nspin: %4i\n", nspin)
            @printf(fout, "# nmesh: %4i\n", nmesh)
            @printf(fout, "# qdim : %4i\n", qdim)

            # Write separators
            println(fout)
            println(fout)

            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through each frequency point
                for m = 1:nmesh
                    @printf(fout, "w:%6i%16.8f\n", m, fmesh[m])
                    # Go through the orbital space
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            z = Delta[p,q,m,s,t]
                            @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                        end
                    end
                end # END OF M LOOP

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end # END OF IOSTREAM

        # Print message to the screen
        println("  > Split hybridization functions for site $t into: $fhyb")
        println("  > Shape of Array fmesh: ", size(fmesh))
        println("  > Shape of Array Delta: ", size(Delta[:,:,:,:,t]))
    end # END OF T LOOP
end

"""
    write_delta(fmesh::Array{F64,1}, Delta::Array{C64,5}, ai::Array{Impurity,1}, fhyb::String)

Write hybridization functions into the `fhyb` file. This function is usually
called by `mixer_delta()` to update the `dmft1/dmft.delta` file. The working
directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`read_delta`](@ref), [`write_eimpx`](@ref).
"""
function write_delta(fmesh::Array{F64,1}, Delta::Array{C64,5}, ai::Array{Impurity,1}, fhyb::String)
    # Extract the dimensional parameters
    _, qdim, nmesh, nspin, nsite = size(Delta)

    # Determine filename for hybridization functions
    @assert fhyb == "dmft1/dmft.delta"

    # Write the data
    open(fhyb, "w") do fout
        # Write dimensional parameters
        @printf(fout, "# nsite: %4i\n", nsite)
        @printf(fout, "# nspin: %4i\n", nspin)
        @printf(fout, "# nmesh: %4i\n", nmesh)
        @printf(fout, "# qdim : %4i\n", qdim)

        # Write separators
        println(fout)
        println(fout)

        # Go through each quantum impurity problem
        for t = 1:nsite
            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through each frequency point
                for m = 1:nmesh
                    @printf(fout, "w:%6i%16.8f\n", m, fmesh[m])
                    # Go through the orbital space
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            z = Delta[p,q,m,s,t]
                            @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                        end
                    end
                end # END OF M LOOP

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end # END OF T LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  > Write hybridization functions into: $fhyb")
    println("  > Shape of Array fmesh: ", size(fmesh))
    println("  > Shape of Array Delta: ", size(Delta))
end

#
# Service Functions: For I/O Operations
#

"""
    write_eimpx(Eimpx::Array{C64,4}, ai::Array{Impurity,1})

Split local impurity levels into the `impurity.i/dmft.eimpx` file, which
is essential for the quantum impurity solver. The working directory of
this function must be the root folder.

See also: [`Impurity`](@ref), [`read_eimpx`](@ref), [`write_delta`](@ref).
"""
function write_eimpx(Eimpx::Array{C64,4}, ai::Array{Impurity,1})
    # Extract the dimensional parameters
    _, qdim, nspin, nsite = size(Eimpx)

    # Go through each quantum impurity problems
    for t = 1:nsite
        # Determine filename for local impurity levels
        flev = "impurity.$t/dmft.eimpx"

        # Write the data
        open(flev, "w") do fout
            # Write dimensional parameters
            @printf(fout, "# nsite: %4i\n", nsite)
            @printf(fout, "# nspin: %4i\n", nspin)
            @printf(fout, "# qdim : %4i\n", qdim)

            # Write separators
            println(fout)
            println(fout)

            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through the orbital space
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        z = Eimpx[p,q,s,t]
                        @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                    end
                end

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end # END OF IOSTREAM

        # Print message to the screen
        println("  > Split local impurity levels for site $t into: $flev")
        println("  > Shape of Array Eimpx: ", size(Eimpx[:,:,:,t]))
    end # END OF T LOOP
end

"""
    write_eimpx(Eimpx::Array{C64,4}, ai::Array{Impurity,1}, flev::String)

Write local impurity levels into the `flev` file. This function is usually
called by `mixer_eimpx()` to update the `dmft1/dmft.eimpx` file. The
working directory of this function must be the root folder.

See also: [`Impurity`](@ref), [`read_eimpx`](@ref), [`write_delta`](@ref).
"""
function write_eimpx(Eimpx::Array{C64,4}, ai::Array{Impurity,1}, flev::String)
    # Extract the dimensional parameters
    _, qdim, nspin, nsite = size(Eimpx)

    # Determine filename for local impurity levels
    @assert flev == "dmft1/dmft.eimpx"

    # Write the data
    open(flev, "w") do fout
        # Write dimensional parameters
        @printf(fout, "# nsite: %4i\n", nsite)
        @printf(fout, "# nspin: %4i\n", nspin)
        @printf(fout, "# qdim : %4i\n", qdim)

        # Write separators
        println(fout)
        println(fout)

        # Go through each quantum impurity problem
        for t = 1:nsite
            # Go through each spin
            for s = 1:nspin
                # Write key parameters
                @printf(fout, "# site:%4i  spin:%4i  dims:%4i\n", t, s, ai[t].nband)

                # Go through the orbital space
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        z = Eimpx[p,q,s,t]
                        @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                    end
                end

                # Write separators
                println(fout)
                println(fout)
            end # END OF S LOOP
        end # END OF T LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  > Write local impurity levels into: $flev")
    println("  > Shape of Array Eimpx: ", size(Eimpx))
end

#
# Service Functions: For I/O Operations
#

"""
    write_gamma(kmesh::Array{F64,2}, kwin::Array{I64,3}, gamma::Array{C64,4})

Write correction for density matrix to `dmft2/dmft.gamma` file. This
function is usually called by `mixer_gamma()` function. The working
directory of this function must be the root folder.
"""
function write_gamma(kmesh::Array{F64,2}, kwin::Array{I64,3}, gamma::Array{C64,4})
    # Extract the dimensional parameters
    _, qbnd, nkpt, nspin = size(gamma)

    # Determine filename for correction for density matrix
    fgamma = "dmft2/dmft.gamma"

    # Write the data
    open(fgamma, "w") do fout
        # Write dimensional parameters
        @printf(fout, "# nkpt : %4i\n", nkpt)
        @printf(fout, "# nspin: %4i\n", nspin)
        @printf(fout, "# qbnd : %4i\n", qbnd)

        # Write separators
        println(fout)
        println(fout)

        # Go through each spin
        for s = 1:nspin
            # Go through 𝑘-point
            for k = 1:nkpt                
                # Write key parameters
                #
                # For spin
                @printf(fout, "# spin:%4i\n", s)
                #
                # For 𝑘-point
                @printf(fout, "# kpt:%4i  %16.12f%16.12f%16.12f\n", k, kmesh[k,1:3]...)
                #
                # For band window
                bs = kwin[k,s,1]
                be = kwin[k,s,2]
                cbnd = be - bs + 1
                @printf(fout, "# cbnd:%4i  bs:%4i  be:%4i\n", cbnd, bs, be)

                # Go through the orbital space
                for q = 1:cbnd
                    for p = 1:cbnd
                        z = gamma[p,q,k,s]
                        @printf(fout, "%4i%4i%16.8f%16.8f\n", p, q, real(z), imag(z))
                    end
                end

                # Write separators
                println(fout)
                println(fout)
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  Write gamma matrix into: $fgamma")
end
