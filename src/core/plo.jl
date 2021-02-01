#
# project : pansy
# source  : plo.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/02/02
#

#
# Driver Functions
#

"""
    plo_adaptor(D::Dict{Symbol,Any}, debug::Bool = false)

Adaptor support. It will postprocess the raw projector matrix. The dict
D contains all of the necessary Kohn-Sham data, which will be modified
in this function.
"""
function plo_adaptor(D::Dict{Symbol,Any}, debug::Bool = false)
    # S01: Print the header
    println("  < PLO Adaptor >")

    # S02: Check the validity of the original dict
    key_list = [:PG, :enk, :fermi, :chipsi]
    for k in key_list
        @assert haskey(D, k)
    end

    # S03: Setup the PrGroup strcut further
    println("    Complete Groups")
    plo_group(D[:PG])

    # S04: Create the PrUnion struct
    println("    Generate Unions")
    D[:PU] = plo_union(D[:PG])

    # S05: Adjust the band structure
    println("    Calibrate Fermi Level")
    plo_fermi(D[:enk], D[:fermi])

    # S06: Setup the band / energy window for projectors
    println("    Calibrate Band Window")
    D[:PW] = plo_window(D[:PG], D[:enk])

    # S07: Transform the projector matrix
    println("    Rotate Projectors")
    D[:chipsi_r] = plo_rotate(D[:PG], D[:chipsi])

    # S08: Filter the projector matrix
    println("    Filter Projectors")
    D[:chipsi_f] = plo_filter(D[:PW], D[:chipsi_r])

    # S09: To make sure the projectors orthogonalize with each other
    println("    Normalize Projectors")
    plo_orthog(D[:PW], D[:chipsi_f])
    exit(-1)

    # S10: Write the density matrix and overlap matrix for checking
    if debug
        println("DEBUG!")
    end
end

#
# Service Functions (Group A)
#

"""
    plo_group(PG::Array{PrGroup,1})

Use the information contained in the PIMP dict to further complete
the PrGroup struct.
"""
function plo_group(PG::Array{PrGroup,1})

#
# Remarks:
#
# 1. Until now, the PG array was only created in vasp.jl/vaspio_projs().
#
# 2. In this function, `corr`, `shell`,  and `Tr` which are members of
#    PrGroup struct will be modified according to users' configuration
#    in other words, the case.toml file.
#

    # Additional check for the parameters contained in PIMP dict
    @assert get_i("nsite") === length(get_i("atoms"))
    @assert get_i("nsite") === length(get_i("shell"))

    # The lshell creates a mapping from shell (string) to l (integer).
    # It is used to parse get_i("shell") to extract the `l` parameter.
    lshell = Dict{String,I64}(
                 "s"     => 0,
                 "p"     => 1,
                 "d"     => 2,
                 "f"     => 3,
                 "d_t2g" => 2, # Only a subset of d orbitals
                 "d_eg"  => 2, # Only a subset of d orbitals
             )

    # Loop over each site (the quantum impurity problem) to gather some
    # relevant information, such as `site` and `l`. We use a Array of
    # Tuple (site_l) to record them.
    site_l = Tuple[]
    for i = 1:get_i("nsite")
        # Determine site
        str = get_i("atoms")[i]
        site = parse(I64, line_to_array(str)[3])

        # Determine l and its specification
        str = get_i("shell")[i]
        l = get(lshell, str, nothing)

        # Push the data into site_l
        push!(site_l, (site, l, str))
    end

    # Scan the groups of projectors, setup them one by one.
    for g in eachindex(PG)
        # Examine PrGroup, check number of projectors
        @assert 2 * PG[g].l + 1 === length(PG[g].Pr)

        # Loop over each site (quantum impurity problem)
        for i in eachindex(site_l)
            SL = site_l[i]
            # Well, find out the required PrGroup
            if (PG[g].site, PG[g].l) === (SL[1], SL[2])
                # Setup corr property
                PG[g].corr = true

                # Setup shell property. Later it will be used to generate `Tr`
                PG[g].shell = SL[3]
            end
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
    end
end

"""
    plo_union(PG::Array{PrGroup,1})

Generate the PrUnion struct, which is used to describe the projections.
"""
function plo_union(PG::Array{PrGroup,1})
    # Initialize an array of PrUnion struct
    PU = PrUnion[]

    # Go through each PrGroup
    for i in eachindex(PG)
        # Determine how many projectors should be included in this group
        # according to the rotation matrix (transformation matrix)
        ndim = size(PG[i].Tr)[1]

        # Create a PrUnion struct and push it into the PU array
        push!(PU, PrUnion(PG[i].site, PG[i].l, ndim, PG[i].corr, PG[i].shell))
    end

    # Setup Pr property for each PrUnion struct
    for i in eachindex(PU)
        for j = 1:PU[i].ndim
            PU[i].Pr[j] = j
        end
    end

    # Return the desired arrays
    return PU
end

"""
    plo_fermi(enk::Array{F64,3}, fermi::F64)

Calibrate the band structure to enforce the fermi level to be zero.
"""
function plo_fermi(enk::Array{F64,3}, fermi::F64)
    @. enk = enk - fermi
end

"""
    plo_window(PG::Array{PrGroup,1}, enk::Array{F64,3})

Calibrate the band window to filter the Kohn-Sham eigenvalues.
"""
function plo_window(PG::Array{PrGroup,1}, enk::Array{F64,3})

#
# Remarks:
#
# Here, `window` means energy window or band window. When nwin is 1, it
# means that all PrGroup share the same window. When nwin is equal to
# length(PG), it means that each PrGroup should have its own window.
#

    # Preprocess the input. Get how many windows there are.
    window = get_d("window")
    nwin = convert(I64, length(window) / 2)

    # Sanity check
    @assert nwin === 1 || nwin === length(PG)

    # Initialize an array of PrWindow struct
    PW = PrWindow[]

    # Scan the groups of projectors, setup PrWindow for them.
    for p in eachindex(PG)
        # Determine bwin. Don't forget it is a Tuple. bwin = (emin, emax).
        if nwin === 1
            # All PrGroup shares the same window
            bwin = (window[1], window[2])
        else
            # Each PrGroup has it own window
            bwin = (window[2*p-1], window[2*p])
        end

        # Examine bwin further. Its elements should obey the order. This
        # window must be defined by band indices (they are integers) or
        # energies (two float numbers).
        @assert bwin[2] > bwin[1]
        @assert typeof(bwin[1]) === typeof(bwin[2])
        @assert bwin[1] isa Integer || bwin[1] isa AbstractFloat

        # The `bwin` is only the global window. But we actually need a
        # momentum-dependent and spin-dependent window. This is `kwin`.
        if bwin[1] isa Integer
            kwin = get_win1(enk, bwin)
        else
            kwin = get_win2(enk, bwin)
        end

        # Create the PrWindow struct, and push it into the PW array.
        push!(PW, PrWindow(kwin, bwin))
    end

    # Return the desired arrays
    return PW
end

"""
    plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})

Perform global rotations or transformations for the projectors.
"""
function plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})
    # Extract some key parameters from raw projector matrix
    nproj, nband, nkpt, nspin = size(chipsi)

    # Initialize new array. It stores the rotated projectors.
    # Now it is empty, but we will allocate memory for it later.
    chipsi_r = Array{C64,4}[]

#
# Remarks:
#
# PG[i].Tr must be a matrix. Its size must be (ndim, p2 - p1 + 1).
#

    # Go through each PrGroup and perform the rotation 
    for i in eachindex(PG)
        # Determine the range of original projectors
        p1 = PG[i].Pr[1]
        p2 = PG[i].Pr[end]

        # Determine the number of projectors after rotation
        ndim = size(PG[i].Tr)[1]

        # Create a temporary array M
        M = zeros(C64, ndim, nband, nkpt, nspin)

        # Rotate chipsi by Tr, the results are stored at M.
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    M[:, b, k, s] = PG[i].Tr * chipsi[p1:p2, b, k, s]
                end
            end
        end

        # Push M into chipsi_r to save it
        push!(chipsi_r, M)
    end

    # Return the desired arrays
    return chipsi_r
end

"""
    plo_filter(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1}}

Filter the projector matrix according to band window.
"""
function plo_filter(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    # Initialize new array. It stores the filtered projectors.
    # Now it is empty, but we will allocate memory for it later.
    chipsi_f = Array{C64,4}[]

    # Go through each PrWindow / PrGroup / PrUnion
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nband, nkpt, nspin = size(chipsi[p])

        # Create a temporary array M
        M = zeros(C64, ndim, PW[p].nbnd, nkpt, nspin)

        # Go through each spin and k-point
        for s = 1:nspin
            for k = 1:nkpt
                # Select projectors which live in the given band window
                # `ib1` and `ib2` are the boundaries.
                ib1 = PW[p].kwin[k, s, 1]
                ib2 = PW[p].kwin[k, s, 2]

                # `ib3` are total number of bands for given `s` and `k`
                ib3 = ib2 - ib1 + 1

                # Sanity check
                @assert ib3 <= PW[p].nbnd

                # We just copy data from chipsi to M
                M[:, 1:ib3, k, s] = chipsi[p][:, ib1:ib2, k, s]
            end
        end

        # Push M into chipsi_f
        push!(chipsi_f, M) 
    end

    # Return the desired arrays
    return chipsi_f
end

"""
    plo_orthog(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})

Try to orthogonalize the projectors group by group.
"""
function plo_orthog(PW::Array{PrWindow,1}, chipsi::Array{Array{C64,4},1})
    for p in eachindex(PW)
        # Extract some key parameters
        ndim, nbnd, nkpt, nspin = size(chipsi[p])
        @assert nbnd === PW[p].nbnd

        # Loop over spins and kpoints
        for s = 1:nspin
            for k = 1:nkpt
                # Determine band index and band window
                ib1 = PW[p].kwin[k, s, 1]
                ib2 = PW[p].kwin[k, s, 2]
                ib3 = ib2 - ib1 + 1

                # Make a view for the desired subarray
                M = view(chipsi[p], 1:ndim, 1:ib3, k, s)

                # Orthogonalize it (chipsi is update at the same time)
                try_diag(M)
            end
        end
    end
end




function get_win1(enk::Array{F64,3}, bwin::Tuple{I64,I64})
    bmin, bmax = bwin

    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Create arrays
    # The kwin is used to record the band window for each kpt and each spin
    kwin = zeros(I64, nkpt, nspin, 2)

    fill!(kwin[:, :, 1], bmin)
    fill!(kwin[:, :, 2], bmax)
end

function get_win2(enk::Array{F64,3}, bwin::Tuple{F64,F64})
    emin, emax = bwin

    # Sanity check. Here we should make sure there is an overlap between
    # [emin, emax] and band structure.
    if emax < minimum(enk) || emin > maximum(enk)
        error("Energy window does not overlap with the band structure")
    end

    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Create arrays
    # The kwin is used to record the band window for each kpt and each spin
    kwin = zeros(I64, nkpt, nspin, 2)

    # Scan the band structure to determine kwin
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
            @assert ib1 <= ib2

            # Save the boundaries
            # The ib1 and ib2 mean the lower and upper boundaries, respectively.
            kwin[k, s, 1] = ib1
            kwin[k, s, 2] = ib2
        end
    end

    return kwin
end

"""
    plo_diag(M::AbstractArray{C64,2})

Orthogonalize the given matrix.
"""
function plo_diag(M::AbstractArray{C64,2})
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

#
# Service Functions (Group B)
#

"""
    plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. A general version.
"""
function plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Create overlap array
    ovlp = zeros(F64, nproj, nproj, nspin)

    # Build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            A = view(chipsi, :, :, k, s)
            ovlp[:, :, s] = ovlp[:, :, s] + real(A * A') * wght
        end
    end

    # Return the desired array
    return ovlp
end

"""
    plo_ovlp(PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. It should be block-diagonal.
"""
function plo_ovlp(PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Create overlap array
    ovlp = zeros(F64, nproj, nproj, nspin)

    # Build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            for p in eachindex(PU)
                q1 = PU[p].Pr[1]
                q2 = PU[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                ovlp[q1:q2, q1:q2, s] = ovlp[q1:q2, q1:q2, s] + real(A * A') * wght
            end
        end
    end

    # Return the desired array
    return ovlp
end

"""
    plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. A general version.
"""
function plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Evaluate spin factor
    sf = (nspin === 1 ? 2 : 1)

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
    end

    # Return the desired array
    return dm
end

"""
    plo_dm(bmin::I64, bmax::I64, PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. It should be block-diagonal.
"""
function plo_dm(bmin::I64, bmax::I64, PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Sanity check
    @assert nband === bmax - bmin + 1

    # Evaluate spin factor
    sf = (nspin === 1 ? 2 : 1)

    # Create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # Build density matrix array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[bmin:bmax, k, s]
            for p in eachindex(PU)
                q1 = PU[p].Pr[1]
                q2 = PU[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                dm[q1:q2, q1:q2, s] = dm[q1:q2, q1:q2, s] + real(A * Diagonal(occs) * A') * wght
            end
        end
    end

    # Return the desired array
    return dm
end

"""
    plo_hamk(chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})

Try to build the local hamiltonian. A general version.
"""
function plo_hamk(chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Create hamiltonian array
    hamk = zeros(C64, nproj, nproj, nspin)

    # Build hamiltonian array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            eigs = enk[:, k, s]
            A = view(chipsi, :, :, k, s)
            hamk[:, :, s] = hamk[:, :, s] + (A * Diagonal(eigs) * A') * wght
        end
    end

    # Return the desired array
    return hamk
end

"""
    plo_hamk(bmin::I64, bmax::I64, PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})

Try to build the local hamiltonian. It should be block-diagonal.
"""
function plo_hamk(bmin::I64, bmax::I64, PU::Array{PrUnion,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Sanity check
    @assert nband === bmax - bmin + 1

    # Create hamiltonian array
    hamk = zeros(C64, nproj, nproj, nspin)

    # Build hamiltonian array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            eigs = enk[bmin:bmax, k, s]
            for p in eachindex(PU)
                q1 = PU[p].Pr[1]
                q2 = PU[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                hamk[q1:q2, q1:q2, s] = hamk[q1:q2, q1:q2, s] + (A * Diagonal(eigs) * A') * wght
            end
        end
    end

    # Return the desired array
    return hamk
end

"""
    plo_dos(itet::Array{I64,2}, enk::Array{F64,3})

Try to calculate the density of states.
"""
function plo_dos(bmin::I64, bmax::I64, PU::Array{PrUnion,1}, chipsi::Array{C64,4}, itet::Array{I64,2}, enk::Array{F64,3})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Sanity check
    @assert nband === bmax - bmin + 1

    mesh = collect(-4.0:0.01:4.0)
    pdos = zeros(F64, length(mesh), nproj)

    for i in eachindex(mesh)
        wght = tetra_dos(mesh[i], itet, enk[bmin:bmax, :, :])
        for s = 1:nspin
            for k = 1:nkpt
                for p in eachindex(PU)
                    q1 = PU[p].Pr[1]
                    q2 = PU[p].Pr[end]
                    for q = q1:q2
                        for b = 1:nband
                            pdos[i, q] = pdos[i, q] + wght[b, k, s] * real( chipsi[q, b, k, s] * conj(chipsi[q, b, k, s]) )
                        end
                    end
                end
            end
        end
        print(mesh[i], " ")
        foreach(x -> @printf("%12.7f", x), pdos[i, :])
        println()
    end
end

#
# Service Functions (Group C)
#

"""
    view_ovlp(ovlp::Array{F64,3})

Output the overlap matrix. A general version.
"""
function view_ovlp(ovlp::Array{F64,3})
    # Extract some key parameters
    _, nproj, nspin = size(ovlp)

    # Output the data
    println("<- Overlap Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p = 1:nproj
            foreach(x -> @printf("%12.7f", x), ovlp[p, :, s])
            println()
        end
    end
end

"""
    view_ovlp(PU::Array{PrUnion,1}, ovlp::Array{F64,3})

Output the overlap matrix. It should be block-diagonal.
"""
function view_ovlp(PU::Array{PrUnion,1}, ovlp::Array{F64,3})
    # Extract some key parameters
    _, nproj, nspin = size(ovlp)

    # Output the data
    println("<- Overlap Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PU)
            println("site -> $(PU[p].site) l -> $(PU[p].l) shell -> $(PU[p].shell)")
            q1 = PU[p].Pr[1]
            q2 = PU[p].Pr[end]
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), ovlp[q, q1:q2, s])
                println()
            end
        end
    end
end

"""
    view_dm(dm::Array{F64,3})

Output the density matrix. A general version.
"""
function view_dm(dm::Array{F64,3})
    # Extract some key parameters
    _, nproj, nspin = size(dm)

    # Output the data
    println("<- Density Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p = 1:nproj
            foreach(x -> @printf("%12.7f", x), dm[p, :, s])
            println()
        end
    end
end

"""
    view_dm(PU::Array{PrUnion,1}, dm::Array{F64,3})

Output the density matrix. It should be block-diagonal.
"""
function view_dm(PU::Array{PrUnion,1}, dm::Array{F64,3})
    # Extract some key parameters
    _, nproj, nspin = size(dm)

    # Output the data
    println("<- Density Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PU)
            println("site -> $(PU[p].site) l -> $(PU[p].l) shell -> $(PU[p].shell)")
            q1 = PU[p].Pr[1]
            q2 = PU[p].Pr[end]
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), dm[q, q1:q2, s])
                println()
            end
        end
    end
end

"""
    view_hamk(hamk::Array{C64,3})

Output the local hamiltonian. A general version.
"""
function view_hamk(hamk::Array{C64,3})
    # Extract some key parameters
    _, nproj, nspin = size(hamk)

    # Output the data
    println("<- Local Hamiltonian ->")
    for s = 1:nspin
        println("Spin: $s")
        println("Re:")
        for p in 1:nproj
            foreach(x -> @printf("%12.7f", x), real(hamk[p, :, s]))
            println()
        end
        println("Im:")
        for q in 1:nproj
            foreach(x -> @printf("%12.7f", x), imag(hamk[q, :, s]))
            println()
        end
    end
end

"""
    view_hamk(PU::Array{PrUnion,1}, hamk::Array{C64,3})

Output the local hamiltonian. It should be block-diagonal.
"""
function view_hamk(PU::Array{PrUnion,1}, hamk::Array{C64,3})
    # Extract some key parameters
    _, nproj, nspin = size(hamk)

    # Output the data
    println("<- Local Hamiltonian ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PU)
            println("site -> $(PU[p].site) l -> $(PU[p].l) shell -> $(PU[p].shell)")
            println("Re:")
            q1 = PU[p].Pr[1]
            q2 = PU[p].Pr[end]
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), real(hamk[q, q1:q2, s]))
                println()
            end
            println("Im:")
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), imag(hamk[q, q1:q2, s]))
                println()
            end
        end
    end
end

"""
    view_dos()
"""
function view_dos() end
