#
# Project : Pansy
# Source  : dmft.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/25
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

    # Extract key parameters
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
    println("Engine : DMFT$(subscript(task))")

    # Get the home directory of DMFT engine
    dmft_home = query_dmft()

    # Determine mpi prefix (whether the dmft is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dmft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  Para : Using $numproc processors")

    # Select suitable dmft program
    dmft_exe = "$dmft_home/dmft"
    @assert isfile(dmft_exe)
    println("  Exec : $dmft_exe")

    # Assemble command
    if isnothing(mpi_prefix)
        dmft_cmd = dmft_exe
    else
        dmft_cmd = split("$mpi_prefix $dmft_exe", " ")
    end

    # Launch it, the terminal output is redirected to dmft.out
    run(pipeline(`$dmft_cmd`, stdout = "dmft.out"))

    # Print the footer for a better visualization
    println()
end

"""
    dmft_save(it::IterInfo, task::I64)

Backup the files outputed by the dynamical mean-field theory engine.

See also: [`dmft_init`](@ref), [`dmft_exec`](@ref).
"""
function dmft_save(it::IterInfo, task::I64)
    # Check the task
    @assert task in (1, 2)

    # Create a list of files that need to be backup
    fdmf1 = ["dmft.out"]
    fdmf2 = ["dmft.eimps", "dmft.eimpx", "dmft.fermi"]
    fdmf3 = ["dmft.grn_l", "dmft.hyb_l", "dmft.wss_l"]

    # Be careful, the final file list depends on the task
    if task == 1
        file_list = union(fdmf1, fdmf2, fdmf3)
    else
        file_list = fdmf1
    end

    # Store the data files
    for i in eachindex(file_list)
        f = file_list[i]
        cp(f, "$f.$(it.dmft_cycle).$(it.dmft1_iter)", force = true)
    end
end
