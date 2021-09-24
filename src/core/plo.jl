#
# Project : Pansy
# Source  : plo.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/09/24
#

#=
### *Driver Functions*
=#

"""
    plo_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})

Adaptor support. It will preprocess the raw projector matrix. The dict
`D` contains all of the necessary Kohn-Sham data, which will be modified
in this function.

If it is in the REPL mode and the file `case.test` is present in the
current directory, this function will try to calculate some physical
observables, such as the density matrix, overlap matrix, partial density
of states, and local hamiltonian, which will be written to external
files or terminal for further reference.

See also: [`vasp_adaptor`](@ref), [`ir_adaptor`](@ref).
"""
function plo_adaptor(D::Dict{Symbol,Any}, ai::Array{Impurity,1})
    # Print the header
    println("Adaptor : PLO")
    println("Try to process the Kohn-Sham dataset")
    println("Current directory: ", pwd())

    # Check the validity of the original dict
    key_list = [:enk, :fermi, :chipsi, :PG]
    for k in key_list
        @assert haskey(D, k)
    end

    # How about the original projectors? We will try to calculate some
    # physical observables to check the quality of the projectors.
    isinteractive() &&
    isfile(query_case()*".test") &&
    plo_monitor(D)

    # P01: Create connections/mappings between projectors (or band
    # windows) and quantum impurity problems
    #
    # D[:MAP] will be created
    D[:MAP] = plo_map(D[:PG], ai)

    # P02: Adjust the band structure
    #
    # D[:enk] will be updated
    plo_fermi(D[:enk], D[:fermi])

    # P03: Setup the PrGroup strcut further
    #
    # D[:PG] will be updated
    plo_group(D[:MAP], D[:PG])

    # P04: Setup the band / energy window for projectors
    #
    # D[:PW] will be created
    D[:PW] = plo_window(D[:PG], D[:enk])

    # P05: Transform the projectors
    #
    # D[:Rchipsi] will be created
    D[:Rchipsi] = plo_rotate(D[:PG], D[:chipsi])

    # P06: Filter the projectors
    #
    # D[:Fchipsi] will be created
    D[:Fchipsi] = plo_filter(D[:PW], D[:Rchipsi])

    # P07: Orthogonalize and normalize the projectors
    #
    # D[:Fchipsi] will be updated. It contains the final data
    # for projector matrix.
    plo_orthog(D[:PW], D[:Fchipsi])

    # Are the projectors correct?
    #
    # We will try to calculate some physical quantitites, which
    # will be written to external files or terminal for reference.
    #
    # These physical quantities include density matrix, overlap
    # matrix, local hamiltonian, full hamiltonian, and partial
    # density of states. Of course, it is time-comsuming to do
    # these things. So it is a good idea to turn off this feature
    # if everything is on the way.
    isinteractive() &&
    isfile(query_case()*".test") &&
    plo_monitor(D)
end

#=
### *Service Functions* : *Group A*
=#

"""
    plo_map(PG::Array{PrGroup,1}, ai::Array{Impurity,1})

Create connections / mappings between projectors (or band windows) and
quantum impurity problems. Return a `Mapping` struct.

See also: [`PrGroup`](@ref), [`PrWindow`](@ref), [`Mapping`](@ref).
"""
function plo_map(PG::Array{PrGroup,1}, ai::Array{Impurity,1})
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
    println("  > Create quantum impurity problems -> groups (i_grp)")

    # Examine Map.g_imp
    #
    # For a given group of projectors, if we fail to find out the
    # corresponding quantum impurity problem, it must be non-correlated.
    @assert all(x -> (0 ≤ x ≤ nsite), Map.g_imp)
    #
    println("  > Create groups -> quantum impurity problems (g_imp)")

    # Setup Map.i_wnd and Map.w_imp
    #
    # They are actually copies of i_grp and g_imp
    Map.i_wnd[:] = Map.i_grp[:]
    #
    println("  > Create quantum impurity problems -> windows (i_wnd)")
    #
    Map.w_imp[:] = Map.g_imp[:]
    #
    println("  > Create windows -> quantum impurity problems (w_imp)")

    # Return the desired struct
    return Map
end

"""
    plo_fermi(enk::Array{F64,3}, fermi::F64)

Calibrate the band structure to enforce the fermi level to be zero.

See also: [`vaspio_fermi`](@ref), [`irio_fermi`](@ref).
"""
function plo_fermi(enk::Array{F64,3}, fermi::F64)
    println("Calibrate eigenvalues")
    @. enk = enk - fermi
    println("  > Reset DFT fermi level to zero")
end

#=
*Remarks* :

Until now, the `PG` array was only created in `vasp.jl/vaspio_projs()`.

In this function, `corr`, `shell`,  and `Tr` which are members of
`PrGroup` struct will be modified according to users' configuration,
in other words, the `case.toml` file.
=#

"""
    plo_group(MAP::Mapping, PG::Array{PrGroup,1})

Use the information contained in the `Mapping` struct to further complete
the `PrGroup` struct.

See also: [`PIMP`](@ref), [`Mapping`](@ref), [`PrGroup`](@ref).
"""
function plo_group(MAP::Mapping, PG::Array{PrGroup,1})
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
                PG[g].Tr[1, 1] = 1.0 + 0.0im
                PG[g].Tr[2, 2] = 1.0 + 0.0im
                PG[g].Tr[3, 4] = 1.0 + 0.0im
                break

            @case "d_eg" # TO_BE_CHECK
                PG[g].Tr = zeros(C64, 2, 5)
                PG[g].Tr[1, 3] = 1.0 + 0.0im
                PG[g].Tr[2, 5] = 1.0 + 0.0im
                break

            @default
                sorry()
                break
        end
        println("  > Build transformation matrix for group $g (site: $(PG[g].site))")
    end # END OF G LOOP
end

#=
*Remarks* :

Here, `window` means energy window or band window. When nwin is 1, it
means that all `PrGroup` share the same window. When nwin is equal to
length(PG), it means that each `PrGroup` should have its own window.

If nwin is neither `1` nor `length(PG)`, there must be something wrong.
=#

"""
    plo_window(PG::Array{PrGroup,1}, enk::Array{F64,3})

Calibrate the band window to filter the Kohn-Sham eigenvalues.

See also: [`PrWindow`](@ref), [`get_win1`](@ref), [`get_win2`](@ref).
"""
function plo_window(PG::Array{PrGroup,1}, enk::Array{F64,3})
    # Print the header
    println("Generate windows")

    # Preprocess the input. Get how many windows there are.
    window = get_d("window")
    nwin = convert(I64, length(window) / 2)
    println("  > Number of recognized windows: $nwin")

    # Sanity check
    @assert nwin == 1 || nwin == length(PG)

    # Initialize an array of PrWindow struct
    PW = PrWindow[]

    # Initialize an array of Tuple to record the window for correlated
    # groups of projectors
    CW = Tuple[]

    # Scan the groups of projectors, setup PrWindow for them.
    for p in eachindex(PG)
        # Determine bwin. Don't forget it is a Tuple. bwin = (emin, emax).
        if nwin == 1
            # All `PrGroup` shares the same window
            bwin = (window[1], window[2])
        else
            # Each `PrGroup` has it own window
            bwin = (window[2*p-1], window[2*p])
        end

        # Record the current window if the corresponding group of
        # projectors is correlated. Later it will be used to analyze
        # the correctness of band window.
        if PG[p].corr
            push!(CW, bwin)
        end

        # Examine `bwin` further. Its elements should obey the order. This
        # window must be defined by band indices (they are integers) or
        # energies (two float numbers).
        @assert bwin[2] > bwin[1]
        @assert typeof(bwin[1]) == typeof(bwin[2])
        @assert bwin[1] isa Integer || bwin[1] isa AbstractFloat

        # The `bwin` is only the global window. But we actually need a
        # momentum-dependent and spin-dependent window. This is `kwin`.
        if bwin[1] isa Integer
            kwin = get_win1(enk, bwin)
        else
            kwin = get_win2(enk, bwin)
        end

        # Create the `PrWindow` struct, and push it into the PW array.
        push!(PW, PrWindow(kwin, bwin))

        # Print some useful information
        println("  > Create window $p: $bwin <--> ($(PW[p].bmin), $(PW[p].bmax))")
    end # END OF P LOOP

    # Well, now CW contains all the windows for correlated groups of
    # projectors. In Zen, we assume that all of the correlated groups of
    # projectors must share the same energy / band windows. In other
    # words, the real windows of them must be the same (the indices for
    # these windows can be different). For example, please consider the
    # following case. We have three PrGroups and (of course) PrWindows:
    #
    # PrGroup 1 <-> PrWindow 1
    # PrGroup 2 <-> PrWindow 2
    # PrGroup 3 <-> PrWindow 3
    #
    # where PrGroup 1 and PrGroup 2 are correlated. Then, PrWindow 1
    # and PrWindow 2 must be the same (though their indices are not the
    # same). In other words, their `kwin` must be the same.
    #
    # In order to achieve this goal, the unique elements in CW must be
    # 1. This is the underlying idea for the following codes.
    unique!(CW)
    @assert length(CW) == 1

    # Chech PrWindow for correlated groups again
    #
    # Get the first window for correlated group
    c1 = findfirst(x -> x.corr, PG)
    PW₁ = PW[c1]
    #
    # Then search next window for correlated group
    c2 = c1
    while true
        c2 = findnext(x -> x.corr, PG, c2 + 1)
        # Find nothing, break the cycle
        isa(c2, Nothing) && break
        # Find new window, then we have to compare it with PW₁ and
        # make sure they are the same window.
        PW₂ = PW[c2]
        @assert PW₁ == PW₂
    end
    println("  > Verify windows for correlated groups")

    # Return the desired array
    return PW
end

#=
*Remarks* :

`PG[i].Tr` must be a matrix. Its size must be `(ndim, p2 - p1 + 1)`.
=#

"""
    plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})

Perform global rotations or transformations for the projectors. In
this function, the projectors will be classified into different
groups, and then they will be rotated group by group.

See also: [`PrGroup`](@ref), [`plo_filter`](@ref), [`plo_orthog`](@ref).
"""
function plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})
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

#=
*Theory* :

The projector matrix (`chipsi`) is defined as follows:
```math
P^{\mathbf{R}}_{m\nu}(\mathbf{k})
    =
    \langle \chi^{\mathbf{R}}_{\mathbf{k}m} | \Psi_{\mathbf{k}\nu} \rangle.
```
Here, ``\mathbf{R}`` denotes the correlated atom within the primitive
unit cell, around which the local orbital ``|\chi^{\mathbf{R}}_{m}\rangle``
is centered, and ``m = 1, \cdots, M`` is an orbital index within the
correlated subset. ``|\Psi_{\mathcal{k}\nu}\rangle`` denotes the Bloch
wave functions, and ``\nu`` is an band index.
=#

"""
    plo_filter(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}}

Filter the projector matrix according to band window.

See also: [`PrWindow`](@ref), [`plo_rotate`](@ref), [`plo_orthog`](@ref).
"""
function plo_filter(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
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
    plo_orthog(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})

Orthogonalize and normalize the projectors.

See also: [`PrWindow`](@ref), [`plo_rotate`](@ref), [`plo_filter`](@ref).
"""
function plo_orthog(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    # Print the header
    println("Normalize projectors")

    # Preprocess the input. Get how many windows there are.
    window = get_d("window")
    nwin = convert(I64, length(window) / 2)

    # Sanity check
    @assert nwin == 1 || nwin == length(PW)

    # Choose suitable service functions
    if nwin == 1
        # All the PrGroups share the same energy / band window, we should
        # orthogonalize and normalize the projectors as a whole.
        #
        println("  > Try to orthogonalize the projectors as a whole")
        #
        try_blk1(PW, chipsi)
    else
        # Each PrGroup has its own energy / band window, we have to
        # orthogonalize and normalize the projectors group by group.
        #
        println("  > Try to orthogonalize the projectors group by group")
        #
        try_blk2(PW, chipsi)
    end

    # Print some useful information
    for p in eachindex(chipsi)
        println("  > Final shape of Array chipsi (group $p): $(size(chipsi[p]))")
    end
end

"""
    plo_monitor(D::Dict{Symbol,Any})

Generate some key physical quantities by using the projectors and the
Kohn-Sham band structures. It is used for debug only.

See also: [`plo_adaptor`](@ref).
"""
function plo_monitor(D::Dict{Symbol,Any})
    if haskey(D, :MAP)
        # Calculate and output overlap matrix
        ovlp = calc_ovlp(D[:PW], D[:Fchipsi], D[:weight])
        view_ovlp(D[:PG], ovlp)

        # Calculate and output density matrix
        dm = calc_dm(D[:PW], D[:Fchipsi], D[:weight], D[:occupy])
        view_dm(D[:PG], dm)

        # Calculate and output local hamiltonian
        hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
        view_hamk(D[:PG], hamk)

        # Calculate and output full hamiltonian
        hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])
        view_hamk(hamk)

        # Calculate and output density of states
        if get_d("smear") == "tetra"
            mesh, dos = calc_dos(D[:PW], D[:Fchipsi], D[:itet], D[:enk])
            view_dos(mesh, dos)
        end
    else
        # Calculate and output overlap matrix
        ovlp = calc_ovlp(D[:chipsi], D[:weight])
        view_ovlp(ovlp)

        # Calculate and output density matrix
        dm = calc_dm(D[:chipsi], D[:weight], D[:occupy])
        view_dm(dm)

        # Calculate and output local hamiltonian
        #hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:weight], D[:enk])
        #view_hamk(D[:PG], hamk)
    end
end

#=
### *Service Functions* : *Group B*
=#

"""
    get_win1(enk::Array{F64,3}, bwin::Tuple{I64,I64})

Return momentum- and spin-dependent band window (case 1). The users
provide only the band indices for the window.

See also: [`plo_window`](@ref).
"""
function get_win1(enk::Array{F64,3}, bwin::Tuple{I64,I64})
    # Here, `bwin` defines the global band window.
    bmin, bmax = bwin

    # Extract some key parameters
    nband, nkpt, nspin = size(enk)
    @assert nband ≥ bmin && nband ≥ bmax

    # Create array `kwin`, which is used to record the band window
    # for each k-point and each spin.
    kwin = zeros(I64, nkpt, nspin, 2)

    # Fill `kwin` with global band boundaries
    fill!(view(kwin, :, :, 1), bmin)
    fill!(view(kwin, :, :, 2), bmax)

    # Return the desired array
    return kwin
end

"""
    get_win2(enk::Array{F64,3}, bwin::Tuple{F64,F64})

Return momentum- and spin-dependent band window (case 2). The users
provide only the maximum and minimum energies for the window.

See also: [`plo_window`](@ref).
"""
function get_win2(enk::Array{F64,3}, bwin::Tuple{F64,F64})
    # Here, `bwin` defines the global energy window.
    emin, emax = bwin

    # Sanity check. We should make sure there is an overlap between
    # [emin, emax] and the band structure.
    if emax < minimum(enk) || emin > maximum(enk)
        error("Energy window does not overlap with the band structure")
    end

    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Create array `kwin`, which is used to record the band window
    # for each k-point and each spin.
    kwin = zeros(I64, nkpt, nspin, 2)

    # Scan the band structure to determine `kwin`
    for s = 1:nspin
        for k = 1:nkpt
            # For lower boundary
            ib1 = 1
            while enk[ib1, k, s] < emin
                ib1 = ib1 + 1
            end

            # For upper boundary
            ib2 = nband
            while enk[ib2, k, s] > emax
                ib2 = ib2 - 1
            end

            # Check the boundaries
            @assert ib1 < ib2

            # Save the boundaries. The ib1 and ib2 mean the lower and
            # upper boundaries, respectively.
            kwin[k, s, 1] = ib1
            kwin[k, s, 2] = ib2
        end # END OF K LOOP
    end # END OF S LOOP

    # Return the desired array
    return kwin
end

#=
*Remarks* :

We assume that the energy / band windows for all of the projectors are
the same. In other words, `PW` only has an unique PrWindow object.
=#

"""
    try_blk1(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})

Try to orthogonalize and normalize the projectors as a whole.

See also: [`PrWindow`](@ref), [`try_blk2`](@ref), [`plo_orthog`](@ref).
"""
function try_blk1(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    # Extract some key parameters
    nkpt = size(chipsi[1], 3)
    nspin = size(chipsi[1], 4)

    # Determine number of projectors contained in each group.
    # The `ndims` is a array.
    dims = map(x -> size(x, 1), chipsi)

    # The `block` is used to store the first index and the last
    # index for each group of projectors.
    block = Tuple[]
    start = 0
    for i in eachindex(dims)
        push!(block, (start + 1, start + dims[i]))
        start = start + dims[i]
    end

    # Create a temporary array
    max_proj = sum(dims)
    max_band = PW[1].nbnd
    M = zeros(C64, max_proj, max_band)

    # Loop over spins and k-points
    for s = 1:nspin
        for k = 1:nkpt
            # Determine band indices
            ib1 = PW[1].kwin[k, s, 1]
            ib2 = PW[1].kwin[k, s, 2]

            # Determine band window
            ib3 = ib2 - ib1 + 1

            # Sanity check
            @assert max_band ≥ ib3
            @assert ib3 ≥ max_proj

            # Try to combine all of the groups of projectors
            for p in eachindex(PW)
                M[block[p][1]:block[p][2], 1:ib3] = chipsi[p][:, 1:ib3, k, s]
            end

            # Orthogonalize and normalize it
            try_diag(view(M, :, 1:ib3))

            # Copy the results back to the sources
            for p in eachindex(PW)
                chipsi[p][:, 1:ib3, k, s] = M[block[p][1]:block[p][2], 1:ib3]
            end
        end # END OF K LOOP
    end # END OF S LOOP
end

"""
    try_blk2(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})

Try to orthogonalize the projectors group by group.

See also: [`PrWindow`](@ref), [`try_blk1`](@ref), [`plo_orthog`](@ref).
"""
function try_blk2(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    # Go through each PrWindow / PrGroup
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd == PW[p].nbnd

        # Loop over spins and k-points
        for s = 1:nspin
            for k = 1:nkpt
                # Determine band indices
                ib1 = PW[p].kwin[k, s, 1]
                ib2 = PW[p].kwin[k, s, 2]

                # Determine band window
                ib3 = ib2 - ib1 + 1

                # Sanity check
                @assert ib3 ≤ PW[p].nbnd
                @assert ib3 ≥ ndim

                # Make a view for the desired subarray
                M = view(chipsi[p], 1:ndim, 1:ib3, k, s)

                # Orthogonalize it (chipsi[p] is update at the same time)
                try_diag(M)
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF P LOOP
end

#=
*Theory* :

First, we try to calcuate the overlap matrix:
```math
\mathcal{O}_{mm'}(\mathbf{k})
    =
    \langle \chi_{\mathbf{k}m} | \chi_{\mathbf{k}m'} \rangle
    =
    \sum_{\nu \in \mathcal{W}} P_{m\nu}(\mathbf{k}) P^{*}_{\nu m'}(\mathbf{k}),
```
and its inverse square root:
```math
\mathcal{S}_{mm'}(\mathbf{k})
    =
    \left\{ \mathcal{O}(\mathbf{k})^{-\frac{1}{2}} \right\}_{mm'}.
```
Then the renormalized projector reads:
```math
\bar{P}_{m\nu}(\mathbf{k})
    =
    \sum_{m'} \mathcal{S}_{mm'} (\mathbf{k}) P_{m'\nu}(\mathbf{k}).
```
=#

"""
    try_diag(M::AbstractArray{C64,2})

Orthogonalize the given matrix.

See also: [`try_blk1`](@ref), [`try_blk2`](@ref).
"""
function try_diag(M::AbstractArray{C64,2})
    # Calculate overlap matrix, it must be a hermitian matrix.
    ovlp = M * M'
    @assert ishermitian(ovlp)

    # Diagonalize the overlap matrix, return eigenvalues and eigenvectors.
    vals, vecs = eigen(Hermitian(ovlp))
    @assert all(vals .> 0)

    # Calculate the renormalization factor
    sqrt_vals = map(x -> 1.0 / sqrt(x), vals)
    S = vecs * Diagonal(sqrt_vals) * vecs'

    # Renormalize the input matrix
    copy!(M, S * M)
end

#=
### *Service Functions* : *Group C*
=#

"""
    calc_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. For raw projectors only.

See also: [`view_ovlp`](@ref).
"""
function calc_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)
    @assert nband ≥ nproj

    # Create overlap array
    ovlp = zeros(F64, nproj, nproj, nspin)

    # Build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            A = view(chipsi, :, :, k, s)
            ovlp[:, :, s] = ovlp[:, :, s] + real(A * A') * wght
        end
    end # END OF S LOOP

    # Return the desired array
    return ovlp
end

"""
    calc_ovlp(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. For normalized projectors only.

See also: [`view_ovlp`](@ref), [`PrWindow`](@ref).
"""
function calc_ovlp(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1})
    # Create an empty array. Next we will fill it.
    ovlp = Array{F64,3}[]

    # Go through each PrWindow / PrGroup
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd == PW[p].nbnd

        # Create a temporary array
        V = zeros(F64, ndim, ndim, nspin)

        # Build overlap array
        for s = 1:nspin
            for k = 1:nkpt
                wght = weight[k] / nkpt
                A = view(chipsi[p], :, :, k, s)
                V[:, :, s] = V[:, :, s] + real(A * A') * wght
            end
        end # END OF S LOOP

        # Push V into ovlp to save it
        push!(ovlp, V)
    end # END OF P LOOP

    # Return the desired array
    return ovlp
end

"""
    calc_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. For raw projectors only.

See also: [`view_dm`](@ref).
"""
function calc_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)
    @assert nband ≥ nproj

    # Evaluate spin factor
    sf = (nspin == 1 ? 2 : 1)

    # Create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # Build density matrix array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[:, k, s]
            A = view(chipsi, :, :, k, s)
            dm[:, :, s] = dm[:, :, s] + real(A * Diagonal(occs) * A') * wght
        end
    end # END OF S LOOP

    # Return the desired array
    return dm
end

"""
    calc_dm(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. For normalized projectors only.

See also: [`view_dm`](@ref), [`PrWindow`](@ref).
"""
function calc_dm(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1}, occupy::Array{F64,3})
    # Create an empty array. Next we will fill it.
    dm = Array{F64,3}[]

    # Go through each PrWindow / PrGroup
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd == PW[p].nbnd

        # Evaluate spin factor
        sf = (nspin == 1 ? 2 : 1)

        # Create a temporary array
        M = zeros(F64, ndim, ndim, nspin)

        # Build density matrix array
        for s = 1:nspin
            for k = 1:nkpt
                wght = weight[k] / nkpt * sf
                occs = occupy[PW[p].bmin:PW[p].bmax, k, s]
                A = view(chipsi[p], :, :, k, s)
                M[:, :, s] = M[:, :, s] + real(A * Diagonal(occs) * A') * wght
            end
        end # END OF S LOOP

        # Push M into dm to save it
        push!(dm, M)
    end # END OF P LOOP

    # Return the desired array
    return dm
end

"""
    calc_hamk(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1}, enk::Array{F64,3})

Try to build the local hamiltonian. For normalized projectors only.

See also: [`view_hamk`](@ref), [`PrWindow`](@ref).
"""
function calc_hamk(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, weight::Array{F64,1}, enk::Array{F64,3})
    # Create an empty array. Next we will fill it.
    hamk = Array{C64,3}[]

    # Go through each PrWindow / PrGroup
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd == PW[p].nbnd

        # Create a temporary array
        H = zeros(C64, ndim, ndim, nspin)

        # Build hamiltonian array
        for s = 1:nspin
            for k = 1:nkpt
                wght = weight[k] / nkpt
                eigs = enk[PW[p].bmin:PW[p].bmax, k, s]
                A = view(chipsi[p], :, :, k, s)
                H[:, :, s] = H[:, :, s] + (A * Diagonal(eigs) * A') * wght
            end
        end # END OF S LOOP

        # Push H into hamk to save it
        push!(hamk, H)
    end # END OF P LOOP

    # Return the desired array
    return hamk
end

#=
*Remarks* :

We assume that the energy / band windows for all of the projectors are
the same. In other words, `PW` only has an unique PrWindow object.
=#

"""
    calc_hamk(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, enk::Array{F64,3})

Try to build the full hamiltonian. For normalized projectors only.

See also: [`view_hamk`](@ref), [`PrWindow`](@ref).
"""
function calc_hamk(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, enk::Array{F64,3})
    # Extract some key parameters
    nkpt = size(chipsi[1], 3)
    nspin = size(chipsi[1], 4)

    # Determine number of projectors contained in each group.
    # The `ndims` is a array.
    dims = map(x -> size(x, 1), chipsi)

    # The `block` is used to store the first index and the last
    # index for each group of projectors.
    block = Tuple[]
    start = 0
    for i in eachindex(dims)
        push!(block, (start + 1, start + dims[i]))
        start = start + dims[i]
    end

    # Create a temporary array
    max_proj = sum(dims)
    max_band = PW[1].nbnd
    M = zeros(C64, max_proj, max_band)

    # Create a array for the hamiltonian
    H = zeros(C64, max_proj, max_proj, nkpt, nspin)

    # Loop over spins and k-points
    for s = 1:nspin
        for k = 1:nkpt
            # Determine band indices
            ib1 = PW[1].kwin[k, s, 1]
            ib2 = PW[1].kwin[k, s, 2]

            # Determine band window
            ib3 = ib2 - ib1 + 1

            # Sanity check
            @assert max_band ≥ ib3

            # Try to combine all of the groups of projectors
            for p in eachindex(PW)
                M[block[p][1]:block[p][2], 1:ib3] = chipsi[p][:, 1:ib3, k, s]
            end

            # Build hamiltonian array
            eigs = enk[ib1:ib2, k, s]
            A = view(M, :, 1:ib3)
            H[:, :, k, s] = H[:, :, k, s] + (A * Diagonal(eigs) * A')
        end # END OF K LOOP
    end # END OF S LOOP

    # Return the desired array
    return H
end

"""
    calc_dos(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, itet::Array{I64,2}, enk::Array{F64,3})

Try to calculate the partial density of states using the analytical
tetrahedron method.

See also: [`view_dos`](@ref), [`PrWindow`](@ref).
"""
function calc_dos(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}, itet::Array{I64,2}, enk::Array{F64,3})
    # Create array of mesh
    MA = Array{F64,1}[]

    # Create array of density of states
    DA = Array{F64,3}[]

    # Go through each PrWindow / PrGroup
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd == PW[p].nbnd

        # Create the mesh. It depends on PrWindow.bwin.
        #
        # Extract the band window / energy window
        emin, emax = PW[p].bwin
        #
        # If it is the band window, then we will create a mesh from
        # minimum(enk[emin,:,:]) to maximum(enk[emax,:,:]).
        if emin isa Integer
            M = collect(minimum(enk[emin,:,:]):0.01:maximum(enk[emax,:,:]))
        #
        # If it is the energy window, then we will create a mesh from
        # emin to emax.
        else
            @assert emin isa AbstractFloat
            M = collect(emin:0.01:emax)
        end
        #
        # Evaluate number of mesh points
        npnts = length(M)

        # Create a temporary array for density of states
        D = zeros(F64, ndim, nspin, npnts)

        # Go through each mesh point
        for i in 1:npnts
            # Obtain the integration weights for density of states by
            # using the analytical tetrahedron method.
            W = bzint(M[i], itet, enk[PW[p].bmin:PW[p].bmax, :, :])

            # Perform brillouin zone summation
            for s = 1:nspin, k = 1:nkpt, b = 1:nbnd, q = 1:ndim
                D[q, s, i] = D[q, s, i] + W[b, k, s] * abs( chipsi[p][q, b, k, s] )^2
            end
        end

        # Push M into MA to save it
        push!(MA, M)

        # Push D into DA to save it
        push!(DA, D)
    end # END OF P LOOP

    # Return the desired array
    return MA, DA
end

#=
### *Service Functions* : *Group D*
=#

"""
    view_ovlp(ovlp::Array{F64,3})

Output the overlap matrix to screen. For raw projectors only.

See also: [`calc_ovlp`](@ref).
"""
function view_ovlp(ovlp::Array{F64,3})
    # Print the header
    println("<- Overlap Matrix ->")

    # Extract some key parameters
    _, nproj, nspin = size(ovlp)

    # Output the data
    for s = 1:nspin
        println("Spin: $s")
        for p = 1:nproj
            foreach(x -> @printf("%12.7f", x), ovlp[p, :, s])
            println()
        end
    end # END OF S LOOP
end

"""
    view_ovlp(PG::Array{PrGroup,1}, ovlp::Array{Array{F64,3},1})

Output the overlap matrix to screen. For normalized projectors only.

See also: [`calc_ovlp`](@ref), [`PrGroup`](@ref).
"""
function view_ovlp(PG::Array{PrGroup,1}, ovlp::Array{Array{F64,3},1})
    # Print the header
    println("<- Overlap Matrix ->")

    # Go through each PrGroup
    for p in eachindex(PG)
        println("Site -> $(PG[p].site) L -> $(PG[p].l) Shell -> $(PG[p].shell)")

        # Extract some key parameters
        _, ndim, nspin = size(ovlp[p])

        # Output the data
        for s = 1:nspin
            println("Spin: $s")
            for q = 1:ndim
                foreach(x -> @printf("%12.7f", x), ovlp[p][q, 1:ndim, s])
                println()
            end
        end # END OF S LOOP
    end # END OF P LOOP
end

"""
    view_dm(dm::Array{F64,3})

Output the density matrix to screen. For raw projectors only.

See also: [`calc_dm`](@ref).
"""
function view_dm(dm::Array{F64,3})
    # Print the header
    println("<- Density Matrix ->")

    # Extract some key parameters
    _, nproj, nspin = size(dm)

    # Output the data
    for s = 1:nspin
        println("Spin: $s")
        for p = 1:nproj
            foreach(x -> @printf("%12.7f", x), dm[p, :, s])
            println()
        end
    end # END OF S LOOP
end

"""
    view_dm(PG::Array{PrGroup,1}, dm::Array{Array{F64,3},1})

Output the density matrix to screen. For normalized projectors only.

See also: [`calc_dm`](@ref), [`PrGroup`](@ref).
"""
function view_dm(PG::Array{PrGroup,1}, dm::Array{Array{F64,3},1})
    # Print the header
    println("<- Density Matrix ->")

    # Go through each PrGroup
    for p in eachindex(PG)
        println("Site -> $(PG[p].site) L -> $(PG[p].l) Shell -> $(PG[p].shell)")

        # Extract some key parameters
        _, ndim, nspin = size(dm[p])

        # Output the data
        for s = 1:nspin
            println("Spin: $s")
            for q = 1:ndim
                foreach(x -> @printf("%12.7f", x), dm[p][q, 1:ndim, s])
                println()
            end
        end # END OF S LOOP
    end # END OF P LOOP
end

"""
    view_hamk(PG::Array{PrGroup,1}, hamk::Array{Array{C64,3},1})

Output the local hamiltonian to screen. For normalized projectors only.

See also: [`calc_hamk`](@ref), [`PrGroup`](@ref).
"""
function view_hamk(PG::Array{PrGroup,1}, hamk::Array{Array{C64,3},1})
    # Print the header
    println("<- Local Hamiltonian ->")

    # Go through each PrGroup
    for p in eachindex(PG)
        println("Site -> $(PG[p].site) L -> $(PG[p].l) Shell -> $(PG[p].shell)")

        # Extract some key parameters
        _, ndim, nspin = size(hamk[p])

        # Output the data
        for s = 1:nspin
            println("Spin: $s")

            # Real parts
            println("Re:")
            for q = 1:ndim
                foreach(x -> @printf("%12.7f", x), real(hamk[p][q, 1:ndim, s]))
                println()
            end

            # Imag parts
            println("Im:")
            for q = 1:ndim
                foreach(x -> @printf("%12.7f", x), imag(hamk[p][q, 1:ndim, s]))
                println()
            end
        end # END OF S LOOP
    end # END OF P LOOP
end

#=
*Remarks* :

The data file `hamk.chk` is used to debug. It should not be read by the
DMFT engine. That is the reason why we name this function as `view_hamk`
and put it in plo.jl.
=#

"""
    view_hamk(hamk::Array{C64,4})

Output the full hamiltonian to `hamk.chk`. For normalized projectors only.

See also: [`calc_hamk`](@ref).
"""
function view_hamk(hamk::Array{C64,4})
    # Extract some key parameters
    nproj, _, nkpt, nspin = size(hamk)

    # Output the data
    open("hamk.chk", "w") do fout
        # Write the header
        println(fout, "# File: hamk.chk")
        println(fout, "# Data: hamk[nproj,nproj,nkpt,nspin]")
        println(fout)
        println(fout, "nproj -> $nproj")
        println(fout, "nkpt  -> $nkpt")
        println(fout, "nspin -> $nspin")
        println(fout)

        # Write the body
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nproj
                    for p = 1:nproj
                        z = hamk[p, b, k, s]
                        @printf(fout, "%16.12f %16.12f\n", real(z), imag(z))
                    end
                end
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM
end

"""
    view_dos(mesh::Array{Array{F64,1},1}, dos::Array{Array{F64,3},1})

Output the density of states to `dos.chk`. For normalized projectors only.

See also: [`calc_dos`](@ref).
"""
function view_dos(mesh::Array{Array{F64,1},1}, dos::Array{Array{F64,3},1})
    # Go through each PrGroup
    for p in eachindex(dos)
        # Extract some key parameters
        ndim, nspin, npnts = size(dos[p])
        @assert npnts == length(mesh[p])

        # Output the data
        open("dos.chk.$p", "w") do fout
            # Write the header
            println(fout, "# File: dos.chk")
            println(fout, "# Data: mesh[npnts] and dos[ndim,nspin,npnts]")
            println(fout)
            println(fout, "npnts -> $npnts")
            println(fout, "ndim  -> $ndim")
            println(fout, "nspin -> $nspin")
            println(fout)

            # Write the body
            for m = 1:npnts
                @printf(fout, "%12.7f", mesh[p][m])
                for s = 1:nspin
                    foreach(x -> @printf(fout, "%12.7f", x), dos[p][:, s, m])
                end
                println(fout)
            end
        end # END OF IOSTREAM
    end # END OF P LOOP
end
