#
# Project : Pansy
# Source  : wannier.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/10/09
#

#=
### *Driver Functions*
=#

"""
    wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})

Adaptor support. It will read and parse the outputs of the `wannier90`
code, convert the data into IR format. The data contained in `D` dict
will be updated.

Be careful, now this adaptor only supports `quantum espresso` (`pwscf`).

See also: [`qe_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function wannier_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : WANNIER")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # Check the validity of the original dict
    key_list = [:latt, :kmesh, :weight, :enk, :occupy, :fermi]
    for k in key_list
        @assert haskey(D, k)
    end

    # Extract key parameters
    sp = get_d("lspins") # Is it a spin-polarized system

    # W01: Read energy window (outer window) from w90.wout
    if sp # For spin-polarized system
        # Spin up
        ewin_up = w90_read_wout("up")
        #
        # Spin down
        ewin_dn = w90_read_wout("dn")
        #
        # Calibrate the energy window by the fermi level
        ewin_up = ewin_up .- D[:fermi]
        ewin_dn = ewin_dn .- D[:fermi]
    else # For spin-unpolarized system
        ewin = w90_read_wout()
        ewin = ewin .- D[:fermi]
    end

    # W02: Read accurate band eigenvalues from w90.eig
    #
    # D[:enk] will be updated
    if sp # For spin-polarized system
        # Spin up
        eigs_up = w90_read_eigs("up") .- D[:fermi]
        nband, nkpt = size(eigs_up)
        eigs_up = reshape(eigs_up, (nband, nkpt, 1))
        #
        # Spin down
        eigs_dn = w90_read_eigs("dn") .- D[:fermi]
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
        eigs = w90_read_eigs() .- D[:fermi]
        nband, nkpt = size(eigs)
        eigs = reshape(eigs, (nband, nkpt, 1))
        D[:enk] = deepcopy(eigs)
        @assert size(D[:enk]) == (nband, nkpt, 1)
    end
    #
    # Calibrate the eigenvalues to force the fermi level to be zero
    # Be careful, the original eigenvalues from qeio_eigen() have
    # not been calibrated.
    #
    # The calibration has been done above.
    #
    # @. D[:enk] = D[:enk] - D[:fermi]

    # W03: Deduce band window from energy window
    if sp # For spin-polarized system
        # Spin up
        bwin_up = w90_find_bwin(ewin_up, eigs_up)
        #
        # Spin down
        bwin_dn = w90_find_bwin(ewin_dn, eigs_dn)
    else # For spin-unpolarized system
        bwin = w90_find_bwin(ewin, eigs)
    end

    # W04: Read transform matrix from w90_u.mat
    if sp # For spin-polarized system
        # Spin up
        umat_up = w90_read_umat("up")
        #
        # Spin down
        umat_dn = w90_read_umat("dn")
    else # For spin-unpolarized system
        umat = w90_read_umat()
    end

    # W05: Read disentanglement matrix from w90_u_dis.mat
    if sp # For spin-polarized system
        # Spin up
        udis_up = w90_read_udis(bwin_up, "up")
        #
        # Spin down
        udis_dn = w90_read_udis(bwin_dn, "dn")
    else # For spin-unpolarized system
        udis = w90_read_udis(bwin)
    end

    # W06: Build projection matrix
    #
    # D[:chipsi] will be created
    if sp # For spin-polarized system
        # Spin up
        proj_up = w90_make_chipsi(umat_up, udis_up)
        nproj, nband, nkpt = size(proj_up)
        proj_up = reshape(proj_up, (nproj, nband, nkpt, 1))
        #
        # Spin down
        proj_dn = w90_make_chipsi(umat_dn, udis_dn)
        nproj, nband, nkpt = size(proj_dn)
        proj_dn = reshape(proj_dn, (nproj, nband, nkpt, 1))
        #
        # Sanity check
        @assert size(proj_up) == size(proj_dn)
        #
        # Concatenate proj_up and proj_dn
        D[:chipsi] = cat(proj_up, proj_dn, dims = 4)
        #
        # Sanity check
        @assert size(D[:chipsi]) == (nproj, nband, nkpt, 2)
    else # For spin-unpolarized system
        proj = w90_make_chipsi(umat, udis)
        nproj, nband, nkpt = size(proj)
        D[:chipsi] = reshape(proj, (nproj, nband, nkpt, 1))
    end

    # W07: Setup the PrTrait and PrGroup structs
    #
    # D[:PT] and D[:PG] will be created
    latt = D[:latt]
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

    # W08: Setup the band window for projections
    #
    # If you do not want to filter the projections, please use another
    # version of w90_make_window(), i.e, w90_make_window(PG, eigs). It
    # is quite clear that the current version is much more efficient.
    #
    # D[:PW] will be created
    if sp # For spin-polarized system
        # Spin up
        PW_up = w90_make_window(PG_up, ewin_up, bwin_up)
        #
        # Spin down
        PW_dn = w90_make_window(PG_dn, ewin_dn, bwin_dn)
        #
        # Merge PW_up and PW_dn
        D[:PW] = w90_make_window(PW_up, PW_dn)
    else # For spin-unpolarized system
        PW = w90_make_window(PG, ewin, bwin)
        D[:PW] = deepcopy(PW)
    end

    # W09: Create connections/mappings between projectors (or band
    # windows) and quantum impurity problems
    #
    # D[:MAP] will be created
    D[:MAP] = w90_make_map(D[:PG], ai)

    # W10: Setup the PrGroup strcut further
    #
    # D[:PG] will be updated
    w90_make_group(D[:MAP], D[:PG])

    # W11: Transform the projectors
    #
    # D[:Rchipsi] will be created
    D[:Rchipsi] = w90_make_chipsi(D[:PG], D[:chipsi])

    # W12: Filter the projectors
    #
    # D[:Fchipsi] will be created
    D[:Fchipsi] = w90_make_chipsi(D[:PW], D[:Rchipsi])

    # Are the projectors correct?
    #
    # We will try to calculate some physical quantitites, which
    # will be written to external files or terminal for reference.
    isinteractive() &&
    isfile(query_case()*".test") &&
    wannier_monitor(D)
end

#=
### *Service Functions* : *Group A*
=#

"""
    wannier_init(D::Dict{Symbol,Any}, sp::String = "")

Try to generate the `w90.win` file, which is the essential input for
the `wannier90` code. Here, we always use `w90` as the seedname. If
the system is spin polarized (the argument `sp` is `up` or `dn`), then
the seednames will be `w90up` and `w90dn`, respectively.

See also: [`wannier_exec`](@ref), [`wannier_save`](@ref).
"""
function wannier_init(D::Dict{Symbol,Any}, sp::String = "")
    # Print the header
    println("Generate input files for wannier90")

    # Extract necessary data from D
    # These data are read in qe_adaptor()
    latt  = D[:latt]
    kmesh = D[:kmesh]
    enk   = D[:enk]
    fermi = D[:fermi]

    # Extract the nband parameter
    nband, _, _ = size(enk)

    # Try to prepare control parameters
    w90c = w90_make_ctrl(latt, nband, fermi)

    # Try to define projections
    proj = w90_make_proj()

    # Try to write w90.win
    #
    # Setup filename correctly. It depends on the spin.
    fwin = "w90" * sp * ".win"
    #
    # Write seedname.win
    open(fwin, "w") do fout
        w90_write_win(fout, w90c)
        w90_write_win(fout, proj)
        w90_write_win(fout, latt)
        w90_write_win(fout, kmesh)
    end
    #
    println("  > File $fwin is created")
end

"""
    wannier_exec(sp::String = ""; op::String = "")

Execute the wannier90 program, monitor the convergence progress, and
output the relevant information. The argument `sp` denotes the spin
component, while `op` specifies the running mode for wannier90. If
`op == -pp`, the wannier90 code will try to generate the `w90.nnkp`.

See also: [`wannier_init`](@ref), [`wannier_save`](@ref).
"""
function wannier_exec(sp::String = ""; op::String = "")
    # Print the header
    println("Detect the runtime environment for wannier90")

    # Get the home directory of wannier90
    #
    # We can not guarantee that the wannier90 code is always installed
    # within the directory of quantum espresso.
    wannier90_home = query_dft("wannier90")
    println("  > Home directory for wannier90: ", wannier90_home)

    # Select suitable wannier90 program
    wannier90_exe = "$wannier90_home/wannier90.x"
    @assert isfile(wannier90_exe)
    println("  > Executable program is available: ", basename(wannier90_exe))

    # Assemble command
    seedname = "w90" * sp
    #
    if op == "-pp" # As a preprocessor to generate w90.nnkp
        wannier90_cmd = split("$wannier90_exe $op $seedname", " ")
    else # Standard run to generate wannier function
        wannier90_cmd = split("$wannier90_exe $seedname", " ")
    end
    #
    println("  > Assemble command: $(prod(x -> x * ' ', wannier90_cmd))")

    # Determine suitable output file
    finp = "w90" * sp * ".win"
    fout = "w90" * sp * ".wout"
    println("  > Applying input file: $finp")
    println("  > Applying output file: $fout")

    # Print the header
    println("Launch the computational engine wannier90")

    # Create a task, but do not run it immediately
    t = @task begin
        run(pipeline(`$wannier90_cmd`, stdout = fout))
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

    # Wait for the wannier90 task to finish
    wait(t)

    # Extract some relevant information
    if op != "-pp"
        # For disentanglement
        iters = readlines(fout)
        filter!(x -> contains(x, " <-- DIS"), iters)
        arr = line_to_array(iters[end])
        niter = parse(I64, arr[1])
        delta = parse(F64, arr[4])
        println("  > Disentangled after $niter iterations (Δ = $delta)")
        #
        # For wannierization
        iters = readlines(fout)
        filter!(x -> contains(x, " <-- CONV"), iters)
        arr = line_to_array(iters[end])
        niter = parse(I64, arr[1])
        delta = parse(F64, arr[2])
        println("  > Wannierized after $niter iterations (Δ = $delta)")
    end
end

"""
    wannier_save(sp::String = ""; op::String = "")

Backup and check the output files of wannier90 if necessary. Actually,
there are no files that need to be stored. We just check whether they
have been created correctly.

See also: [`wannier_init`](@ref), [`wannier_exec`](@ref).
"""
function wannier_save(sp::String = ""; op::String = "")
    # Print the header
    println("Finalize the computational task")

    # Check the output files of wannier90
    if op == "-pp"
        fnnkp = "w90" * sp * ".nnkp"
        if isfile(fnnkp)
            println("  > File $fnnkp is created")
        else
            error("File $fnnkp is not created")
        end
    else
        fhr = "w90" * sp * "_hr.dat"
        fu = "w90" * sp * "_u.mat"
        fdis = "w90" * sp * "_u_dis.mat"
        #
        flist = (fhr, fu, fdis)
        for i in eachindex(flist)
            filename = flist[i]
            if isfile(filename)
                println("  > File $filename is created")
            else
                error("File $filename is not created")
            end
        end
    end
end

"""
    wannier_monitor(D::Dict{Symbol,Any})

Try to check and examine whether the obtained wannier functions are
correct and reasonable.

See also: [`wannier_adaptor`](@ref), [`plo_monitor`](@ref).
"""
function wannier_monitor(D::Dict{Symbol,Any})
    # Calculate and output overlap matrix
    ovlp = calc_ovlp(D[:PW], D[:Fchipsi], D[:weight])
    view_ovlp(D[:PG], ovlp)

    # Calculate and output density matrix
    dm = calc_dm(D[:PW], D[:Fchipsi], D[:weight], D[:occupy])
    view_dm(D[:PG], dm)

    # Calculate and output local hamiltonian
    hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
    view_hamk(D[:PG], hamk)
end

#=
### *Service Functions* : *Group B*
=#

"""
    w90_make_ctrl(latt:Lattice, nband::I64, fermi::F64)

Try to make the control parameters for the `w90.win` file. The `latt`
object represent the crystallography information, and `nband` is the
number of Kohn-Sham states outputed by the DFT code, `fermi` is the
fermi level. This function is called by `wannier_init()`.

See also: [`w90_make_proj`](@ref).
"""
function w90_make_ctrl(latt::Lattice, nband::I64, fermi::F64)
    # Create a dict to store the configurations, which will be returned.
    w90c = Dict{String,Any}()

    # Generate a dict, which defines a mapping from orbital definition
    # to number of orbitals.
    #
    # Perhaps you have to extend it to support your cases.
    orb_dict = Dict{String,I64}(
                 "s"   => 1,
                 "l=0" => 1,
                 "p"   => 3,
                 "l=1" => 3,
                 "d"   => 5,
                 "l=2" => 5,
                 "f"   => 7,
                 "l=3" => 7,
             )

    # Get number of wannier functions, `num_wann`.
    #
    # Step 1, get the string for projection.
    sproj = get_d("sproj")
    @assert length(sproj) ≥ 2
    # The first element of sproj should specify the type of wannier
    # function. Now only the MLWF and SAWF modes are supported.
    @assert sproj[1] in ("mlwf", "sawf")
    #
    # Step 2, calculate num_wann.
    num_wann = 0
    # Go throught each element of sproj
    for i = 2:length(sproj)
        # Extract atomic symbol and orbital specification
        str_atm, str_orb = strip.(split(sproj[i], ":"))
        # We have to remove all white spaces in str_orb. It is necessary.
        str_orb = replace(str_orb, " " => "")
        # Extract the corresponding number of orbitals
        num_orb = orb_dict[str_orb]
        # Find out how many atoms are there in the lattice.
        num_atm = 0
        for j = 1:latt.natom
            if latt.atoms[j] == str_atm
                num_atm = num_atm + 1
            end
        end
        @assert num_atm ≥ 1
        # Update num_wann
        num_wann = num_wann + num_orb * num_atm
    end
    #
    # Step 3, store num_wann in the dict.
    w90c["num_wann"] = num_wann

    # Get number of bands, `num_bands`.
    @assert nband ≥ num_wann
    w90c["num_bands"] = nband

    # Deal with the disentanglement setup
    #
    # Step 1, get the setup for disentanglement.
    window = get_d("window")
    @assert window isa Vector{F64}
    @assert length(window) ≥ 2
    #
    # Step 2, calibrate the energy window by the fermi level
    window = window .+ fermi
    #
    # Step 3, determine the disentanglement parameters and store them.
    if length(window) == 2
        w90c["dis_win_min"]  = window[1]
        w90c["dis_win_max"]  = window[2]
    elseif length(window) == 4
        w90c["dis_win_min"]  = window[1]
        w90c["dis_win_max"]  = window[2]
        w90c["dis_froz_min"] = window[3]
        w90c["dis_froz_max"] = window[4]
    else
        error("Wrong window's definition")
    end

    # Some additional but necessary parameters for wannier90
    #
    # Acticate the convergence window
    w90c["conv_window"] = 3
    #
    # We should write the hamiltonian
    w90c["write_hr"] = ".true."
    #
    # We should write the transform matrix
    w90c["write_u_matrices"] = ".true."
    #
    # Support symmetry-adapted wannier functions
    if sproj[1] == "sawf"
        w90c["site_symmetry"] = ".true."
    end
    #
    # Special treatment for spinors
    if get_d("lspinorb")
        w90c["num_wann"] = num_wann * 2
        w90c["spinors"] = ".true."
    end

    # Return the required object.
    return w90c
end

"""
    w90_make_proj()

Try to make the projection block for the `w90.win` file. We will not
check the validness of these projections here. This function is called
by the `wannier_init()` function.

See also: [`w90_make_ctrl`](@ref).
"""
function w90_make_proj()
    # Create a string array to store the definitions for projections
    proj = String[]

    # Get the string for projection
    sproj = get_d("sproj")

    # Check the length of sproj
    @assert length(sproj) ≥ 2

    # The first element of sproj should specify the type of wannier
    # function. Now only the MLWF and SAWF modes are supported.
    @assert sproj[1] in ("mlwf", "sawf")

    # Transfer sproj to proj
    for i = 2:length(sproj)
        push!(proj, sproj[i])
    end

    # Return the desired array
    return proj
end

#=
### *Service Functions* : *Group C*
=#

"""
    w90_make_map(PG::Array{PrGroup,1}, ai::Array{Impurity,1})

Create connections / mappings between projectors (or band windows) and
quantum impurity problems. Return a `Mapping` struct. This function is
just a copy of the `plo_map()` function.

See also: [`PrGroup`](@ref), [`PrWindow`](@ref), [`Mapping`](@ref).
"""
function w90_make_map(PG::Array{PrGroup,1}, ai::Array{Impurity,1})
    # Print the header
    println("Establish mapping")

    # Extract key parameters
    #
    # Here, `nsite` is the number of quantum impurity problems, `ngrp` is
    # the number of groups for projectors, and `nwnd` is the number of
    # band windows.
    #
    # In this code, we have to ensure `nwnd` is always equal to `ngrp`.
    nsite = get_i("nsite")
    ngrp = length(PG)
    nwnd = ngrp

    # Additional check for nsite
    @assert nsite == length(ai)

    # The lshell creates a mapping from shell (string) to l (integer).
    # It is used to parse Impurity.shell to extract the `l` parameter.
    lshell = Dict{String,I64}(
                 "s"     => 0,
                 "p"     => 1,
                 "d"     => 2,
                 "f"     => 3,
                 "d_t2g" => 2, # Only a subset of d orbitals
                 "d_eg"  => 2, # Only a subset of d orbitals
             )

    # Loop over each site (the quantum impurity problem) to gather some
    # relevant information, such as `site` and `l`. We use an array of
    # Tuple (site_l) to record them.
    site_l = Tuple[]
    for i = 1:nsite
        # Determine site
        sites = ai[i].sites

        # Determine l and its specification
        shell = ai[i].shell
        l = get(lshell, shell, nothing)

        # Push the data into site_l
        push!(site_l, (sites, l, shell))
    end
    println("  > Figure out the traits of quantum impurity problems")

    # Create the Mapping struct
    Map = Mapping(nsite, ngrp, nwnd)

    # Determine Map.i_grp (imp -> grp) and Map.g_imp (grp -> imp)
    for i = 1:nsite
        SL = site_l[i]
        for g = 1:ngrp
            if (PG[g].site, PG[g].l) == (SL[1], SL[2])
                Map.i_grp[i] = g
                Map.g_imp[g] = i
            end
        end
    end

    # Examine Map.i_grp
    #
    # For a given quantum impurity problem, we can always find out the
    # corresponding group of projectors.
    @assert all(x -> (0 < x ≤ ngrp), Map.i_grp)
    #
    println("  > Map quantum impurity problems to groups  (i_grp)")

    # Examine Map.g_imp
    #
    # For a given group of projectors, if we fail to find out the
    # corresponding quantum impurity problem, it must be non-correlated.
    @assert all(x -> (0 ≤ x ≤ nsite), Map.g_imp)
    #
    println("  > Map groups to quantum impurity problems  (g_imp)")

    # Setup Map.i_wnd and Map.w_imp
    #
    # They are actually copies of i_grp and g_imp
    Map.i_wnd[:] = Map.i_grp[:]
    #
    println("  > Map quantum impurity problems to windows (i_wnd)")
    #
    Map.w_imp[:] = Map.g_imp[:]
    #
    println("  > Map windows to quantum impurity problems (w_imp)")

    # Return the desired struct
    return Map
end

"""
    w90_make_group(latt::Lattice, sp::String = "")

Try to read the `w90.nnkp` file, parse the `projections` block. Finally,
it will return arrays of PrTrait and PrGroup objects, which contain the
definitions of projectors. The argument `latt` is essential. It includes
the atomic coordinates for all lattice sites, which are quite useful to
distinguish these projectors. And the argument `sp` is optional. It is
used when the system is spin-polarized.

See also: [`PrTrait`](@ref), [`PrGroup`](@ref).
"""
function w90_make_group(latt::Lattice, sp::String = "")
    # Print the header
    println("Build traits and groups")

    # Read and parse the `w90.nnkp` file
    #
    # Build the filename
    fnnkp = "w90" * sp * ".nnkp"
    println("  > Open and read $fnnkp")
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    #
    # Read it and figure out the projections block
    lines = readlines(fnnkp)
    ind = findfirst(x -> contains(x, "begin projections"), lines)
    @assert ind > 0
    start = ind + 1
    #
    # Create some arrays to store the projectors
    nproj = parse(I64, lines[start]) # Number of projectors
    coord = zeros(F64, nproj, 3)     # Atomic coordinates
    l_vec = zeros(I64, nproj) # Quantum number 𝑙
    m_vec = zeros(I64, nproj) # Quantum number 𝑚
    #
    # Parse the projectors
    for i = 1:nproj
        start = start + 1
        arr = line_to_array(lines[start])
        coord[i,:] = parse.(F64, arr[1:3])
        l_vec[i] = parse(I64, arr[4])
        m_vec[i] = parse(I64, arr[5])
        start = start + 1
    end

    # Try to build PrTrait struct. The raw information about projectors
    # should be encapsulated in it.
    #
    # Define all possible specifications for projectors. Be careful, the
    # following orbital sequence is different from that defined in the
    # constructor of PrTrait struct.
    spec = ("s",
            "pz", "px", "py",
            "dz2", "dxz", "dyz", "dx2-y2", "dxy",
            "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)"
    )
    #
    # Generate PrTrait struct one by one
    PT = PrTrait[]
    for i = 1:nproj
        # Try to figure out which site it is. We just compare the atomic
        # coordinates in coord with those saved in latt.coord.
        site = -1
        for j = 1:latt.natom
            if latt.coord[j,1:3] == coord[i,1:3]
                site = j # That is it
                break
            end
        end
        #
        # Sanity check
        @assert site > 0
        #
        # Get 𝑙, 𝑚, and desc.
        l = l_vec[i]
        m = m_vec[i]
        desc = spec[m + l*l]
        #
        # Save the PrTrait struct
        # Here we call the default constructor of PrTrait
        push!(PT, PrTrait(site, l, m, desc))
    end

    # The following codes are borrowed from vaspio_projs().

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

    # Print some useful information to check
    println("  > Number of projectors: ", nproj)
    for i in eachindex(PT)
        println("    [ PrTrait $i ]")
        println("      site -> ", PT[i].site)
        println("      l -> ", PT[i].l)
        println("      m -> ", PT[i].m)
        println("      desc -> ", PT[i].desc)
    end
    #
    println("  > Number of groups: ", length(PG))
    for i in eachindex(PG)
        println("    [ PrGroup $i ]")
        println("      site -> ", PG[i].site)
        println("      l -> ", PG[i].l)
        println("      corr -> ", PG[i].corr)
        println("      shell -> ", PG[i].shell)
        println("      Pr -> ", PG[i].Pr)
    end

    # Return the desired arrays
    # Note: PG will be further fixed at another w90_make_group() function.
    return PT, PG
end

"""
    w90_make_group(MAP::Mapping, PG::Array{PrGroup,1})

Use the information contained in the `Mapping` struct to further setup
the `PrGroup` struct. This function is just a copy of the `plo_group()`
function.

See also: [`PIMP`](@ref), [`Mapping`](@ref), [`PrGroup`](@ref).
"""
function w90_make_group(MAP::Mapping, PG::Array{PrGroup,1})
    # Print the header
    println("Complete groups")

    # Scan the groups of projectors, setup them one by one.
    for g in eachindex(PG)
        # Examine PrGroup, check number of projectors
        @assert 2 * PG[g].l + 1 == length(PG[g].Pr)

        # Extract the index of quantum impurity problem for the current
        # group of projectors
        s = MAP.g_imp[g]

        # Yes, this group of projectors has a corresponding quantum
        # impurity problem. We have to further modify it. If `s` is
        # 0, it means that the group of projectors is non-correlated.
        if s != 0
            # Setup corr property
            PG[g].corr = true
            println("  > Turn group $g (site: $(PG[g].site)) into correlated")

            # Setup shell property
            # Later it will be used to generate `Tr`
            PG[g].shell = get_i("shell")[s]
            println("  > Treat group $g (site: $(PG[g].site)) as $(PG[g].shell) orbitals")
        end

        # Setup Tr array further
        @cswitch PG[g].shell begin
            @case "s"
                PG[g].Tr = Matrix{ComplexF64}(I, 1, 1)
                break

            @case "p"
                PG[g].Tr = Matrix{ComplexF64}(I, 3, 3)
                break

            @case "d"
                PG[g].Tr = Matrix{ComplexF64}(I, 5, 5)
                break

            @case "f"
                PG[g].Tr = Matrix{ComplexF64}(I, 7, 7)
                break

            @case "d_t2g"
                PG[g].Tr = zeros(C64, 3, 5)
                # For vasp
                is_vasp() && begin
                    PG[g].Tr[1, 1] = 1.0 + 0.0im
                    PG[g].Tr[2, 2] = 1.0 + 0.0im
                    PG[g].Tr[3, 4] = 1.0 + 0.0im
                end
                # For quantum espresso + wannier90
                is_qe() && begin
                    PG[g].Tr[1, 1] = 1.0 + 0.0im # WRONG
                    PG[g].Tr[2, 2] = 1.0 + 0.0im # WRONG
                    PG[g].Tr[3, 4] = 1.0 + 0.0im # WRONG
                end
                break

            @case "d_eg" # TO_BE_CHECK
                PG[g].Tr = zeros(C64, 2, 5)
                # For vasp
                is_vasp() && begin
                    PG[g].Tr[1, 3] = 1.0 + 0.0im
                    PG[g].Tr[2, 5] = 1.0 + 0.0im
                end
                # For quantum espresso + wannier90
                is_qe() && begin
                    PG[g].Tr[1, 3] = 1.0 + 0.0im # WRONG
                    PG[g].Tr[2, 5] = 1.0 + 0.0im # WRONG
                end
                break

            @default
                sorry()
                break
        end
        println("  > Build transformation matrix for group $g (site: $(PG[g].site))")
    end # END OF G LOOP
end

"""
    w90_make_window(PG::Array{PrGroup,1}, enk::Array{F64,3})

Make band window to filter the projections. Actually, all of the Kohn-Sham
eigenvalues are retained, so the band window is always `[1, nband]`. This
function will return an array of `PrWindow` struct.

See also: [`PrWindow`](@ref).
"""
function w90_make_window(PG::Array{PrGroup,1}, enk::Array{F64,3})
    # Print the header
    println("Generate windows")

    # Extract the key parameters
    nband, nkpt, _ = size(enk)

    # Initialize an array of PrWindow struct
    PW = PrWindow[]

    # Scan the groups of projectors, setup PrWindow for them.
    for p in eachindex(PG)
        # Setup global band window
        bwin = (1, nband)
        #
        # Setup momentum-dependent band window
        # Be careful, here we assume nspin = 1.
        kwin = zeros(I64, nkpt, 1, 2)
        fill!(view(kwin, :, :, 1), 1)
        fill!(view(kwin, :, :, 2), nband)
        #
        # Create the `PrWindow` struct, and push it into the PW array.
        push!(PW, PrWindow(kwin, bwin))
        #
        # Print some useful information
        println("  > Create window $p: $bwin <--> ($(PW[p].bmin), $(PW[p].bmax))")
    end # END OF P LOOP

    # Return the desired array
    return PW
end

"""
    w90_make_window(PG::Array{PrGroup,1},
                    ewin::Tuple{F64,F64},
                    bwin::Array{I64,2})

Make band window to filter the projections. Actually, only those relevant
bands (which are restricted by the energy window `ewin` or the band window
`bwin`) are retained. This function will return an array of `PrWindow`
struct.

See also: [`PrWindow`](@ref).
"""
function w90_make_window(PG::Array{PrGroup,1},
                         ewin::Tuple{F64,F64},
                         bwin::Array{I64,2})
    # Print the header
    println("Generate windows")

    # Extract the key parameters
    nkpt, _ = size(bwin)

    # Initialize an array of PrWindow struct
    PW = PrWindow[]

    # Scan the groups of projectors, setup PrWindow for them.
    for p in eachindex(PG)
        # Setup momentum-dependent band window
        # Be careful, here we assume nspin = 1.
        kwin = zeros(I64, nkpt, 1, 2)
        #
        # We copy bwin to kwin directly
        @. kwin[:,1,:] = bwin
        #
        # Create the `PrWindow` struct, and push it into the PW array.
        push!(PW, PrWindow(kwin, ewin))
        #
        # Print some useful information
        println("  > Create window $p: $ewin <--> ($(PW[p].bmin), $(PW[p].bmax))")
    end # END OF P LOOP

    # Return the desired array
    return PW
end

"""
    w90_make_window(PWup::Array{PrWindow,1}, PWdn::Array{PrWindow,1})

Try to merge two arrays of `PrWindow` struct and generate a new one.
Actually, the new array is similar to the olds. We only modify one of
its members, `kwin`.

See also: [`PrWindow`](@ref).
"""
function w90_make_window(PWup::Array{PrWindow,1}, PWdn::Array{PrWindow,1})
    # Sanity check
    #
    # The two arrays must be the same. The comparison of two PrWindow
    # structs are defined in types.jl.
    @assert PWup == PWdn

    # Create an empty array for PrWindow struct
    PW = PrWindow[]

    # Go through old array of PrWindow struct
    for p in eachindex(PWup)
        # Extract `kwin`
        kwin1 = PWup[p].kwin
        kwin2 = PWdn[p].kwin
        #
        # Check `kwin`, its dimension for spin orientation must be 1
        @assert size(kwin1, 2) == 1
        @assert size(kwin2, 2) == 1
        #
        # Merge two `kwin`
        kwin = cat(kwin1, kwin2, dims = 2)
        #
        # Extract `bwin`
        bwin  = PWup[p].bwin
        #
        # Create new PrWindow and push it into PW
        push!(PW, PrWindow(kwin, bwin))
    end

    # Return the desired array
    return PW
end

"""
    w90_make_chipsi(umat::Array{C64,3}, udis::Array{C64,3})

Try to merge the transform matrix `umat` with the disentanglement matrix
`udis` to construct the final projection matrix `chipsi`, which actually
is the overlap matrix between the wannier functions and the Kohn-Sham
wave functions.

See also: [`w90_read_umat`](@ref), [`w90_read_udis`](@ref).
"""
function w90_make_chipsi(umat::Array{C64,3}, udis::Array{C64,3})
    # Print the header
    println("Generate projection matrix")

    # Extract key parameters
    nproj, _, nkpt = size(umat) # (nproj, nproj, nkpt)
    nband, _nproj, _nkpt = size(udis)
    #
    # Sanity check
    @assert nproj == _nproj
    @assert nkpt == _nkpt

    # Create arrays
    utmp = zeros(C64, nband, nproj)
    chipsi = zeros(C64, nproj, nband, nkpt)

    # Go through each k-point
    for k = 1:nkpt
        # Einstein summation
        for j = 1:nproj
            for i = 1:nband
                utmp[i,j] = zero(C64)
                for m = 1:nproj
                    utmp[i,j] = utmp[i,j] + udis[i,m,k] * umat[m,j,k]
                end
            end # END OF I LOOP
        end # END OF J LOOP
        # Calculate conjugate transpose
        chipsi[:,:,k] = utmp'
    end # END OF K LOOP

    # Print some useful information
    println("  > Number of k-points: ", nkpt)
    println("  > Number of DFT bands: ", nband)
    println("  > Number of wannier functions: ", nproj)
    println("  > Shape of Array chipsi: ", size(chipsi))

    # Return the desired array
    return chipsi
end

"""
    w90_make_chipsi(PG::Array{PrGroup,1}, chipsi::Array{C64,4})

Perform global rotations or transformations for the projectors. In
this function, the projectors will be classified into different
groups, and then they will be rotated group by group. This function
is just a copy of the `plo_rotate()` function.

See also: [`PrGroup`](@ref), [`plo_rotate`](@ref).
"""
function w90_make_chipsi(PG::Array{PrGroup,1}, chipsi::Array{C64,4})
    # Print the header
    println("Rotate projectors")

    # Extract some key parameters from raw projector matrix
    nproj, nband, nkpt, nspin = size(chipsi)
    @assert nproj ≥ 1

    # Initialize new array. It stores the rotated projectors.
    # Now it is empty, but we will allocate memory for it later.
    Rchipsi = Array{C64,4}[]

    # Go through each PrGroup and perform the rotation
    for i in eachindex(PG)
        # Determine the range of original projectors
        p1 = PG[i].Pr[1]
        p2 = PG[i].Pr[end]

        # Determine the number of projectors after rotation
        ndim = size(PG[i].Tr)[1]
        @assert size(PG[i].Tr)[2] == (p2 - p1 + 1)

        # Create a temporary array R
        R = zeros(C64, ndim, nband, nkpt, nspin)
        @assert nband ≥ ndim

        # Rotate chipsi by Tr, the results are stored at R.
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    R[:, b, k, s] = PG[i].Tr * chipsi[p1:p2, b, k, s]
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP

        # Push R into Rchipsi to save it
        push!(Rchipsi, R)

        # Print some useful information
        println("  > Rotate group $i (site: $(PG[i].site)): number of local orbitals -> $ndim")
    end # END OF I LOOP

    # Return the desired array
    return Rchipsi
end

"""
    w90_make_chipsi(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}}

Filter the projector matrix according to band window. This function is
just a copy of the `plo_filter()` function.

See also: [`PrWindow`](@ref), [`plo_filter`](@ref).
"""
function w90_make_chipsi(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    # Print the header
    println("Filter projectors")

    # Initialize new array. It stores the filtered projectors.
    # Now it is empty, but we will allocate memory for it later.
    Fchipsi = Array{C64,4}[]

    # Go through each PrWindow
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nband, nkpt, nspin = size(chipsi[p])

        # Create a temporary array F
        F = zeros(C64, ndim, PW[p].nbnd, nkpt, nspin)
        @assert nband ≥ PW[p].nbnd ≥ ndim

        # Go through each spin and k-point
        for s = 1:nspin
            for k = 1:nkpt
                # Select projectors which live in the given band window
                # `ib1` and `ib2` are the boundaries.
                ib1 = PW[p].kwin[k, s, 1]
                ib2 = PW[p].kwin[k, s, 2]
                @assert ib1 ≤ ib2

                # `ib3` are total number of bands for given `s` and `k`
                ib3 = ib2 - ib1 + 1

                # Sanity check
                @assert ib3 ≤ PW[p].nbnd

                # We just copy data from chipsi[p] to F
                F[:, 1:ib3, k, s] = chipsi[p][:, ib1:ib2, k, s]
            end # END OF K LOOP
        end # END OF S LOOP

        # Push F into Fchipsi to save it
        push!(Fchipsi, F)

        # Print some useful information
        println("  > Apply window $p: maximum number of bands -> $(PW[p].nbnd)")
    end # END OF P LOOP

    # Return the desired array
    return Fchipsi
end

"""
    w90_make_hamr()
"""
function w90_make_hamr()
end

"""
    w90_make_hamk()
"""
function w90_make_hamk()
end

"""
    w90_find_bwin(ewin::Tuple{F64,F64}, enk::Array{F64,3})

During the disentanglement procedure, we can define an outer energy
window to restrict the Kohn-Sham eigenvalues. This function will return
the corresponding band window, which will be used to displace the
disentanglement matrix. Here, `ewin` is the outer energy window which
is extracted from `w90.wout`, and `enk` is the Kohn-Sham eigenvalues.

This function works for spin-unpolarized case only.

See also: [`w90_read_udis`](@ref), [`w90_read_wout`](@ref).
"""
function w90_find_bwin(ewin::Tuple{F64,F64}, enk::Array{F64,3})
    # Print the header
    println("Extract band window for disentanglement")

    # Extract key parameters
    nband, nkpt, nspin = size(enk)
    @assert nspin == 1

    # Setup energy window
    emin, emax = ewin

    # Create an array for momentum-dependent band window
    bwin = zeros(I64, nkpt, 2)

    # Go through each k-point to figure out the band window
    for k = 1:nkpt
        bmin = findfirst(x -> x > emin, enk[:,k,1])
        bmax = findfirst(x -> x > emax, enk[:,k,1]) - 1
        bwin[k,1] = bmin
        bwin[k,2] = bmax
        @assert nband ≥ bmax > bmin ≥ 1
    end

    # Print some useful information
    println("  > Number of k-points: ", nkpt)
    println("  > Minimum band index: ", minimum(bwin[:,1]))
    println("  > Maixmum band index: ", maximum(bwin[:,2]))
    println("  > Shape of Array bwin: ", size(bwin))

    # Return the desired array
    return bwin
end

#=
### *Service Functions* : *Group D*
=#

"""
    w90_read_amat(sp::String = "")

Try to read and parse the `w90.amn` file to get the ``A_{mn}`` matrix,
which represents the projection of the Bloch states onto trail localized
orbitals. The argument `sp` denotes the spin component.

Note that this function has not been used so far.
"""
function w90_read_amat(sp::String = "")
    # Print the header
    println("Parse trail projection")

    # Build the filename
    famn = "w90" * sp * ".amn"
    println("  > Open and read $famn")

    # Read the `w90.amn` file
    lines = readlines(famn)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nband, nkpt, nproj = parse.(I64, line_to_array(lines[2]))

    # Create an array for ``A_{mn}``
    Amn = zeros(C64, nproj, nband, nkpt)

    # We use `start` to record the line index
    start = 2

    # Parse the data and fill in the Amn array
    for k = 1:nkpt
        for p = 1:nproj
            for b = 1:nband
                # Increase the counter
                start = start + 1
                # Convert string to array
                arr = line_to_array(lines[start])
                # Determine indices for bands, projectors, and 𝑘-points
                _b, _p, _k = parse.(I64, arr[1:3])
                # Sanity check
                @assert _b == b
                @assert _p == p
                @assert _k == k
                # Fill in the array
                _re, _im = parse.(F64, arr[4:5])
                Amn[p,b,k] = _re + im*_im
            end # END OF B LOOP
        end # END OF P LOOP
    end  # END OF K LOOP
    #
    @assert start - 2 == nband * nproj * nkpt

    # Print some useful information to check
    println("  > Number of wannier functions: ", nproj)
    println("  > Number of DFT bands: ", nband)
    println("  > Number of k-points: ", nkpt)
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    println("  > Shape of Array Amn: ", size(Amn))

    # Return the desired array
    return Amn
end

"""
    w90_read_eigs(sp::String = "")

Try to read and parse the `w90.eig` file to get the band eigenvalues.
Note that the eigenvalues from `nscf.out` are not accurate enough. We
will use the data extracted from `w90.eig` to update them. The argument
`sp` denotes the spin component.

See also: [`qeio_eigen`](@ref).
"""
function w90_read_eigs(sp::String = "")
    # Print the header
    println("Refine band eigenvalues")

    # Build the filename
    feig = "w90" * sp * ".eig"
    println("  > Open and read $feig")

    # Read the `w90.eig` file
    lines = readlines(feig)
    @assert length(lines) ≥ 1

    # Determine the key parameters from the last line!
    nband, nkpt = parse.(I64, line_to_array(lines[end])[1:2])

    # Create an array for ``E_{nk}``
    eigs = zeros(F64, nband, nkpt)

    # We use `start` to record the line index
    start = 0

    # Parse the data and fill in the eigs array
    for k = 1:nkpt
        for b = 1:nband
            # Increase the counter
            start = start + 1
            # Convert string to array
            arr = line_to_array(lines[start])
            # Determine indices for bands and 𝑘-points
            _b, _k = parse.(I64, arr[1:2])
            # Sanity check
            @assert _b == b
            @assert _k == k
            # Fill in the array
            _v = parse(F64, arr[3])
            eigs[b,k] = _v
        end # END OF B LOOP
    end # END OF K LOOP
    #
    @assert start == nkpt * nband

    # Print some useful information to check
    println("  > Number of DFT bands: ", nband)
    println("  > Number of k-points: ", nkpt)
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    println("  > Shape of Array enk: ", size(eigs))

    # Return the desired array
    return eigs
end

"""
    w90_read_hmat(sp::String = "")

Try to read and parse the `w90_hr.dat` file, return the hamiltonian
matrix in WF representation, the Wigner-Seitz grid points, and their
weights (degeneracies). The argument `sp` denotes the spin component.
The data return by this function can be used to validate the projection
matrix further.

See also: [`wannier_monitor`](@ref).
"""
function w90_read_hmat(sp::String = "")
    # Print the header
    println("Parse hamiltonian in WF basis")

    # Build the filename
    fhr = "w90" * sp * "_hr.dat"
    println("  > Open and read $fhr")

    # Read the `w90_hr.dat` file
    lines = readlines(fhr)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nproj = parse(I64, lines[2]) # Number of wannier functions
    nrpt = parse(I64, lines[3])  # Number of Wigner-Seitz grid points

    # Create arrays
    hamr = zeros(C64, nproj, nproj, nrpt) # Hamiltonian matrix
    rvec = zeros(I64, nrpt, 3) # Wigner-Seitz grids
    rdeg = zeros(I64, nrpt)    # Degeneracies or weights

    # We use `start` to record the line index
    start = 3

    # Read the degeneracies
    #
    # Determine how many lines are there for this block
    nrow = div(nrpt, 15)
    nrem = rem(nrpt, 15)
    @assert nrow ≥ 1
    #
    # Read the first nrow lines
    # There are 15 elements per line.
    for r = 1:nrow
        start = start + 1
        rs = (r - 1) * 15 + 1
        re = (r - 1) * 15 + 15
        rdeg[rs:re] = parse.(I64, line_to_array(lines[start]))
    end
    #
    # Read the final line
    if nrem > 0
        start = start + 1
        rs = nrow * 15 + 1
        re = nrpt
        @assert nrem == re - rs + 1
        rdeg[rs:re] = parse.(I64, line_to_array(lines[start]))
    end

    # Try to read the hamiltonian
    for r = 1:nrpt
        for i = 1:nproj
            for j = 1:nproj
                # Increase the counter
                start = start + 1
                # Convert string to line
                arr = line_to_array(lines[start])
                # Get Wigner-Seitz grid
                rvec[r,:] = parse.(I64, arr[1:3])
                # Get indices of wannier functions
                _j, _i = parse.(I64, arr[4:5])
                # Sanity check
                @assert _j == j
                @assert _i == i
                # Fill in the hamiltonian
                _re, _im = parse.(F64, arr[6:7])
                hamr[j, i, r] = _re + im * _im
            end # END OF J LOOP
        end # END OF I LOOP
    end # END OF R LOOP

    # Print some useful information to check
    println("  > Number of wannier functions: ", nproj)
    println("  > Number of Wigner-Seitz points: ", nrpt)
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    println("  > Shape of Array rdeg: ", size(rdeg))
    println("  > Shape of Array rvec: ", size(rvec))
    println("  > Shape of Array hamr: ", size(hamr))

    # Return the desired arrays
    return rdeg, rvec, hamr
end

"""
    w90_read_umat(sp::String = "")

Try to read and parse the `w90_u.mat` file, return the u-matrix, which
gives the unitary rotations from the optimal subspace to the optimally
smooth states. The argument `sp` denotes the spin component.

See also: [`w90_read_udis`](@ref).
"""
function w90_read_umat(sp::String = "")
    # Print the header
    println("Parse umat")

    # Build the filename
    fu = "w90" * sp * "_u.mat"
    println("  > Open and read $fu")

    # Read the `w90_u.mat` file
    lines = readlines(fu)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nkpt, nproj = parse.(I64, line_to_array(lines[2])[1:2])

    # Create array
    umat = zeros(C64, nproj, nproj, nkpt)

    # We use `start` to record the line index
    start = 2

    # Try to read the u-matrix
    for k = 1:nkpt
        # Increase the counter
        start = start + 2 # Skip one empty line and one 𝑘-point line
        #
        # Parse the data
        for j = 1:nproj
            for i = 1:nproj
                # Increase the counter
                start = start + 1
                # Parse the line and fill in the array
                _re, _im = parse.(F64, line_to_array(lines[start]))
                umat[i, j, k] = _re + im * _im
            end # END OF I LOOP
        end # END OF J LOOP
    end # END OF K LOOP

    # Print some useful information to check
    println("  > Number of wannier functions: ", nproj)
    println("  > Number of k-points: ", nkpt)
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    println("  > Shape of Array umat: ", size(umat))

    # Return the desired array
    return umat
end

"""
    w90_read_udis(bwin::Array{I64,2}, sp::String = "")

Try to read and parse the `w90_u_dis.mat` file. Return the udis-matrix,
which gives the `nproj` dimension optimal subspace from the original
bloch states. Actually, it is the transform matrix for disentanglement.
The argument `sp` denotes the spin component, the band window `bwin` is
from `w90_make_window()` and `w90_read_wout()`.

See also: [`w90_read_umat`](@ref).
"""
function w90_read_udis(bwin::Array{I64,2}, sp::String = "")
    # Print the header
    println("Parse udis")

    # Build the filename
    fdis = "w90" * sp * "_u_dis.mat"
    println("  > Open and read $fdis")

    # Read the `w90_u_dis.mat` file
    lines = readlines(fdis)
    @assert length(lines) ≥ 3

    # Determine the key parameters
    nkpt, nproj, nband = parse.(I64, line_to_array(lines[2]))

    # Create array
    utmp = zeros(C64, nband, nproj)
    udis = zeros(C64, nband, nproj, nkpt)

    # We use `start` to record the line index
    start = 2

    # Try to read the udis-matrix
    for k = 1:nkpt
        # Increase the counter
        start = start + 2 # Skip one empty line and one 𝑘-point line
        #
        # Parse the data
        for j = 1:nproj
            for i = 1:nband
                # Increase the counter
                start = start + 1
                # Parse the line and fill in the array
                _re, _im = parse.(F64, line_to_array(lines[start]))
                utmp[i, j] = _re + im * _im
            end # END OF I LOOP
        end # END OF J LOOP
        #
        # Shift the raw disentanglement matrix by the band window
        bs = bwin[k,1]
        be = bwin[k,2]
        nb = be - bs + 1
        udis[bs:be, :, k] = utmp[1:nb, :]
    end # END OF K LOOP

    # Print some useful information to check
    println("  > Number of DFT bands: ", nband)
    println("  > Number of wannier functions: ", nproj)
    println("  > Number of k-points: ", nkpt)
    println("  > Spin orientation: ", sp == "" ? "none" : sp)
    println("  > Shape of Array udis: ", size(udis))

    # Return the desired array
    return udis
end

"""
    w90_read_wout(sp::String = "")

Try to read and parse the `w90.wout` file. Return the energy window for
disentanglement procedure. The argument `sp` denotes the spin component.

See also: [`w90_make_window`](@ref).
"""
function w90_read_wout(sp::String = "")
    # Print the header
    println("Extract energy window for disentanglement")

    # Build the filename
    fout = "w90" * sp * ".wout"
    println("  > Open and read $fout")
    println("  > Spin orientation: ", sp == "" ? "none" : sp)

    # Read the `w90.wout` file
    lines = readlines(fout)
    @assert length(lines) ≥ 3

    # Locate the line
    filter!(x -> contains(x, "Outer:"), lines)
    @assert length(lines) == 1 # It should contain only 1 line

    # Extract the energy window
    arr = line_to_array(lines[1])
    wmin = parse(F64, arr[3])
    wmax = parse(F64, arr[5])
    println("  > Outer energy window: ($wmin, $wmax) eV")

    # Return the desired energy window
    return (wmin, wmax)
end

#=
### *Service Functions* : *Group E*
=#

"""
    w90_write_win(io::IOStream, w90c::Dict{String,Any})

Write control parameters into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, w90c::Dict{String,Any})
    for key in keys(w90c)
        val = w90c[key]
        println(io, key, " = ", val)
    end
    # Add an empty line
    println(io)
end

"""
    w90_write_win(io::IOStream, proj::Array{String,1})

Write projection block into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, proj::Array{String,1})
    println(io, "begin projections")
    #
    for i = 1:length(proj)
        println(io, proj[i])
    end
    #
    println(io, "end projections\n")
end

"""
    w90_write_win(io::IOStream, latt::Lattice)

Write crystallography information into `w90.win`.

See also: [`Lattice`](@ref), [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, latt::Lattice)
    # Extract parameters
    natom = latt.natom

    # Convert atomiclength (bohr) to angstrom
    # 1 bohr = 0.529177249 angstrom
    lvect = latt.lvect * (latt.scale * 0.529177249)

    # Print atoms_frac block
    println(io, "begin atoms_frac")
    #
    for i = 1:natom
        @printf(io, "%4s%12.8f%12.8f%12.8f\n", latt.atoms[i], latt.coord[i,:]...)
    end
    #
    println(io, "end atoms_frac\n")

    # Print unit_cell_cart block
    println(io, "begin unit_cell_cart")
    #
    for i = 1:3
        @printf(io, "%12.8f%12.8f%12.8f\n", lvect[i,:]...)
    end
    #
    println(io, "end unit_cell_cart\n")
end

"""
    w90_write_win(io::IOStream, kmesh::Array{F64,2})

Write 𝑘-mesh block into `w90.win`.

See also: [`wannier_init`](@ref).
"""
function w90_write_win(io::IOStream, kmesh::Array{F64,2})
    # Extract parameters
    nkpt, ndir = size(kmesh)

    # Sanity check
    @assert ndir == 3

    # Write the block for k-points
    ndiv = ceil(I64, nkpt^(1/3))
    #
    @printf(io, "mp_grid :%3i%3i%3i\n", ndiv, ndiv, ndiv)
    println(io, "begin kpoints")
    #
    for k = 1:nkpt
        @printf(io, "%12.8f%12.8f%12.8f\n", kmesh[k,:]...)
    end
    #
    println(io, "end kpoints\n")
end

#=
### *Service Functions* : *Group X*
=#

"""
    pw2wan_init(case::String, sp::String = "")

Check the runtime environment of `pw2wannier90`, prepare necessary input
files (it is `case.pw2wan`). The argument `case` means the prefix for
`quantum espresso`, and `sp` determines the spin component which can be
empty string, `up`, or `dn`.

See also: [`pw2wan_exec`](@ref), [`pw2wan_save`](@ref).
"""
function pw2wan_init(case::String, sp::String = "")
    # Print the header
    println("Generate input files for pw2wannier90")

    # Try to create a QENamelist object.
    #
    # Setup name of namelist. It is always fixed.
    name = "inputpp"
    #
    # Create an empty dict
    NLData = Dict{AbstractString,Any}()
    #
    # Update the dict
    NLData["outdir"] = "'./'"
    NLData["prefix"] = "'$case'"
    if sp == ""
        NLData["seedname"] = "'w90'"
        NLData["spin_component"] = "'none'"
    elseif sp == "up"
        NLData["seedname"] = "'w90up'"
        NLData["spin_component"] = "'up'"
    elseif sp == "dn"
        NLData["seedname"] = "'w90dn'"
        NLData["spin_component"] = "'down'"
    end
    NLData["write_mmn"] = ".true."
    NLData["write_amn"] = ".true."
    NLData["write_dmn"] = ".true."
    NLData["write_unk"] = ".false."
    #
    # Create a QENamelist
    QNL = QENamelist(name, NLData)

    # Try to write case.pw2wan
    fwan = case * sp * ".pw2wan"
    open(fwan, "w") do fout
        write(fout, QNL) # This write function is defined in qe.jl
    end
    #
    println("  > File $fwan is created")
end

"""
    pw2wan_exec(case::String, sp::String = "")

Execute the pw2wannier90 program, monitor the convergence progress, and
output the relevant information. Here, `case` means the prefix for input
files, and `sp` is `up` or `dn`.

See also: [`pw2wan_init`](@ref), [`pw2wan_save`](@ref).
"""
function pw2wan_exec(case::String, sp::String = "")
    # Print the header
    println("Detect the runtime environment for pw2wannier90")

    # Get the home directory of pw2wannier90
    # It is actually the same with that of quantum espresso.
    qe_home = query_dft("qe")
    println("  > Home directory for pw2wannier90: ", qe_home)

    # Select suitable pw2wannier90 program
    qe_exe = "$qe_home/pw2wannier90.x"
    @assert isfile(qe_exe)
    println("  > Executable program is available: ", basename(qe_exe))

    # Assemble command
    qe_cmd = split("$qe_exe", " ")
    println("  > Assemble command: $(prod(x -> x * ' ', qe_cmd))")

    # Determine suitable input and output files
    finp = case * sp * ".pw2wan"
    fout = "pw2wan.out"
    println("  > Applying input file: $finp")
    println("  > Applying output file: $fout")

    # Print the header
    println("Launch the computational engine pw2wannier90")

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

    # Wait for the pw2wannier90 task to finish
    wait(t)
end

"""
    pw2wan_save(sp::String = "")

Backup and check the output files of the `pw2wannier90` code if necessary.
The argument `sp` specifies the spin component. It could be empty string,
`up`, or `dn`.

See also: [`pw2wan_init`](@ref), [`pw2wan_exec`](@ref).
"""
function pw2wan_save(sp::String = "")
    # Print the header
    println("Finalize the computational task")

    # Check the output files of pw2wannier90
    #
    # Determine filenames
    seedname = "w90" * sp
    famn = seedname * ".amn"
    fmmn = seedname * ".mmn"
    feig = seedname * ".eig"
    fdmn = seedname * ".dmn"
    fsym = seedname * ".sym"
    #
    # Determine list of files
    flist = [famn, fmmn, feig]
    #
    # If the symmetry-adapted wannier function mode is chosen, more
    # output files are created by pw2wannier90.
    sproj = get_d("sproj")
    if sproj[1] == "sawf"
        push!(flist, fdmn)
        push!(flist, fsym)
    end
    #
    # Check the availability of the output files
    for i in eachindex(flist)
        filename = flist[i]
        if isfile(filename)
            println("  > File $filename is created")
        else
            error("File $filename is not created")
        end
    end
end
