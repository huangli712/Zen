#
# project : pansy
# source  : plo.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/25
#

#
# Driver Functions
#

"""
    plo_adaptor(debug::Bool = false)

Adaptor support. It will postprocess the raw projector matrix.
"""
function plo_adaptor(debug::Bool = false)
    # S01: Print the header
    println("  < PLO Adaptor >")

    # S02: Setup the PrGroup strcut
    println("    Grouping")
    if haskey(KohnShamData, "PG")
        plo_group(KohnShamData["PG"])
    else
        error("The KohnShamData dict does not contain the key: PG")
    end

    # S03: Transform the projector matrix
    println("    Rotating")
    if haskey(KohnShamData, "PG") && haskey(KohnShamData, "chipsi")
        KohnShamData["PGT"], KohnShamData["chipsi"] = plo_rotate(KohnShamData["PG"], KohnShamData["chipsi"])
    else
        error("The KohnShamData dict does not contain the keys: PG and chipsi")
    end

    # S04: Adjust the band structure
    println("    Leveling")
    if haskey(KohnShamData, "fermi") && haskey(KohnShamData, "enk")
    else
        error("The KohnShamData dict does not contain the keys: fermi and enk")
    end

    # S05:
    println("    Filtering")

    # S06:
    println("    Orthogonalizing")

    # S07:
    if debug
        # TODO
        println("DEBUG!")
    end
end

"""
    plo_group(PG::Array{PrGroup,1})

Use the information contained in the PIMP dict to further complete
the PrGroup struct.
"""
function plo_group(PG::Array{PrGroup,1})
    # additional check for the parameters contained in PIMP dict
    @assert get_i("nsite") === length(get_i("atoms"))
    @assert get_i("nsite") === length(get_i("shell"))

    # lshell defines a mapping from shell (string) to l (integer)
    lshell = Dict{String,I64}(
                 "s"     => 0,
                 "p"     => 1,
                 "d"     => 2,
                 "f"     => 3,
                 "d_t2g" => 2, # only a subset of d orbitals
                 "d_eg"  => 2, # only a subset of d orbitals
             )

    # loop over each site (quantum impurity problem)
    for i = 1:get_i("nsite")
        # determine site
        str = get_i("atoms")[i]
        site = parse(I64, line_to_array(str)[3])

        # determine l
        str = get_i("shell")[i]
        l = haskey(lshell, str) ? lshell[str] : nothing

        # scan the groups of projectors
        for g in eachindex(PG)
            # examine PrGroup, check number of projectors
            @assert 2 * PG[g].l + 1 === length(PG[g].Pr)

            # well, find out the required PrGroup
            if (PG[g].site, PG[g].l) === (site, l)
                # setup corr property
                PG[g].corr = true

                # setup shell property
                PG[g].shell = str
            end

            # setup Tr array further
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

                @case "d_eg"
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
end

"""
    plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})

Perform global rotations or transformations for the projectors
"""
function plo_rotate(PG::Array{PrGroup,1}, chipsi::Array{C64,4})
    # create a array of PrGroupT struct
    PGT = PrGroupT[]
    for i in eachindex(PG)
        # determine how many projectors should be included in this group
        # according to the rotation matrix (transformation matrix)
        ndim = size(PG[i].Tr)[1]
        # create PrGroupT and push it into PGT array
        push!(PGT, PrGroupT(PG[i].site, PG[i].l, ndim, PG[i].corr, PG[i].shell))
    end

    # setup the Pr property of PrGroupT struct
    c = 0
    for i in eachindex(PGT)
        for j = 1:PGT[i].ndim
            c = c + 1
            PGT[i].Pr[j] = c
        end
    end
    # until now PrGroupT struct is ready

    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # tell us how many projectors there are after the rotation
    nproj_ = sum(x -> x.ndim, PGT)
    @assert nproj_ === sum(x -> size(x.Tr)[1], PG)
    @assert nproj_ <= nproj

    # create new arrays
    chipsi_r = zeros(C64, nproj_, nband, nkpt, nspin)

    # perform rotation or transformation
    # note: PG[i].Tr must be a matrix. its size must be (q2 - q1 + 1, p2 - p1 + 1)
    for s = 1:nspin
        for k = 1:nkpt
            for b = 1:nband
                for i in eachindex(PG)
                    p1 = PG[i].Pr[1]
                    p2 = PG[i].Pr[end]
                    q1 = PGT[i].Pr[1]
                    q2 = PGT[i].Pr[end]
                    chipsi_r[q1:q2, b, k, s] = PG[i].Tr * chipsi[p1:p2, b, k, s]
                end
            end
        end
    end

    # return the desired arrays
    return PGT, chipsi_r
end

"""
    plo_window(enk::Array{F64,3}, emax::F64, emin::F64, chipsi::Array{C64,4})

Extract the projectors within a given energy window
"""
function plo_window(enk::Array{F64,3}, emax::F64, emin::F64, chipsi::Array{C64,4})
    # sanity check
    # make sure there is an overlap between [emin, emax] and band structure
    if emax < minimum(enk) || emin > maximum(enk)
        error("Energy window does not overlap with the band structure")
    end

    # extract some key parameters
    nband, nkpt, nspin = size(enk)

    # create arrays
    # ib_window is used to record the band window for each kpt and each spin
    ib_window = zeros(I64, nkpt, nspin, 2)

    # scan the band structure to determine ib_window
    for s = 1:nspin
        for k = 1:nkpt
            # for lower boundary
            ib1 = 1
            while enk[ib1, k, s] < emin
                ib1 = ib1 + 1
            end

            # for upper boundary
            ib2 = nband
            while enk[ib2, k, s] > emax
                ib2 = ib2 - 1
            end

            # check the boundaries
            @assert ib1 <= ib2

            # save the boundaries
            # ib1 and ib2 mean the lower and upper boundaries, respectively
            ib_window[k, s, 1] = ib1
            ib_window[k, s, 2] = ib2
        end
    end

    # try to find out the global minimum and maximum band indices
    ib_min = minimum(ib_window[:, :, 1])
    ib_max = maximum(ib_window[:, :, 2])

    # try to find out the maximum number of selected bands in the window
    nbmax = ib_max - ib_min + 1

    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # create arrays
    # chipsi_w is used to store the required projectors
    chipsi_w = zeros(C64, nproj, nbmax, nkpt, nspin)

    # select projectors which live in the given band window
    # copy data from chipsi to chipsi_w
    for s = 1:nspin
        for k = 1:nkpt
            ib1 = ib_window[k, s, 1]
            ib2 = ib_window[k, s, 2]
            ib3 = ib2 - ib1 + 1
            @assert ib3 <= nbmax
            chipsi_w[:, 1:ib3, k, s] = chipsi[:, ib1:ib2, k, s]
        end
    end

    # return the desired arrays
    return ib_min, ib_max, ib_window, chipsi_w
end

"""
    plo_orthog(window::Array{I64,3}, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4})

Try to orthogonalize the projectors group by group (site_l by site_l)
"""
function plo_orthog(window::Array{I64,3}, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # loop over spins and kpoints
    for s = 1:nspin
        for k = 1:nkpt
            # determine band index and band window
            b1 = window[k, s, 1]
            b2 = window[k, s, 2]
            nb = b2 - b1 + 1
            for p in eachindex(PGT)
                # determine projector index
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]

                # make a view for the desired subarray
                M = view(chipsi, q1:q2, 1:nb, k, s)

                # orthogonalize it (chipsi is update at the same time)
                plo_diag(M)
            end
        end
    end
end

"""
    plo_diag(M::AbstractArray{C64,2})

Orthogonalize the given matrix
"""
function plo_diag(M::AbstractArray{C64,2})
    # calculate overlap matrix, it must be a hermitian matrix
    ovlp = M * M'
    @assert ishermitian(ovlp)

    # diagonalize the overlap matrix, return eigenvalues and eigenvectors
    vals, vecs = eigen(Hermitian(ovlp))
    @assert all(vals .> 0)

    # calculate the renormalization factor
    sqrt_vals = map(x -> 1.0 / sqrt(x), vals)
    S = vecs * Diagonal(sqrt_vals) * vecs'

    # renormalize the input matrix
    copy!(M, S * M)
end

"""
    plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. General version
"""
function plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # create overlap array
    ovlp = zeros(F64, nproj, nproj, nspin)

    # build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            A = view(chipsi, :, :, k, s)
            ovlp[:, :, s] = ovlp[:, :, s] + real(A * A') * wght
        end
    end

    # return the desired array
    return ovlp
end

"""
    plo_ovlp(PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap matrix out of projectors. It should be block-diagonal
"""
function plo_ovlp(PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # create overlap array
    ovlp = zeros(F64, nproj, nproj, nspin)

    # build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            for p in eachindex(PGT)
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                ovlp[q1:q2, q1:q2, s] = ovlp[q1:q2, q1:q2, s] + real(A * A') * wght
            end
        end
    end

    # return the desired array
    return ovlp
end

"""
    plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. General version
"""
function plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # evaluate spin factor
    sf = (nspin === 1 ? 2 : 1)

    # create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # build density matrix array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[:, k, s]
            A = view(chipsi, :, :, k, s)
            dm[:, :, s] = dm[:, :, s] + real(A * Diagonal(occs) * A') * wght
        end
    end

    # return the desired array
    return dm
end

"""
    plo_dm(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors. It should be block-diagonal
"""
function plo_dm(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # sanity check
    @assert nband === bmax - bmin + 1

    # evaluate spin factor
    sf = (nspin === 1 ? 2 : 1)

    # create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # build density matrix array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[bmin:bmax, k, s]
            for p in eachindex(PGT)
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                dm[q1:q2, q1:q2, s] = dm[q1:q2, q1:q2, s] + real(A * Diagonal(occs) * A') * wght
            end
        end
    end

    # return the desired array
    return dm
end

"""
    plo_hamk(chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})

Try to build the local hamiltonian. General version
"""
function plo_hamk(chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # create hamiltonian array
    hamk = zeros(C64, nproj, nproj, nspin)

    # build hamiltonian array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            eigs = enk[:, k, s]
            A = view(chipsi, :, :, k, s)
            hamk[:, :, s] = hamk[:, :, s] + (A * Diagonal(eigs) * A') * wght
        end
    end

    # return the desired array
    return hamk
end

"""
    plo_hamk(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})

Try to build the local hamiltonian. It should be block-diagonal
"""
function plo_hamk(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, enk::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # sanity check
    @assert nband === bmax - bmin + 1

    # create hamiltonian array
    hamk = zeros(C64, nproj, nproj, nspin)

    # build hamiltonian array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            eigs = enk[bmin:bmax, k, s]
            for p in eachindex(PGT)
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]
                A = view(chipsi, q1:q2, :, k, s)
                hamk[q1:q2, q1:q2, s] = hamk[q1:q2, q1:q2, s] + (A * Diagonal(eigs) * A') * wght
            end
        end
    end

    # return the desired array
    return hamk
end

"""
    plo_dos(itet::Array{I64,2}, enk::Array{F64,3})

Try to calculate the density of states
"""
function plo_dos(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, itet::Array{I64,2}, enk::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # sanity check
    @assert nband === bmax - bmin + 1

    mesh = collect(-4.0:0.01:4.0)
    pdos = zeros(F64, length(mesh), nproj)

    for i in eachindex(mesh)
        wght = tetra_dos(mesh[i], itet, enk[bmin:bmax, :, :])
        for s = 1:nspin
            for k = 1:nkpt
                for p in eachindex(PGT)
                    q1 = PGT[p].Pr[1]
                    q2 = PGT[p].Pr[end]
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

"""
    view_ovlp(ovlp::Array{F64,3})

Output the overlap matrix. General version
"""
function view_ovlp(ovlp::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(ovlp)

    # output the data
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
    view_ovlp(PGT::Array{PrGroupT,1}, ovlp::Array{F64,3})

Output the overlap matrix. It should be block-diagonal
"""
function view_ovlp(PGT::Array{PrGroupT,1}, ovlp::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(ovlp)

    # output the data
    println("<- Overlap Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PGT)
            println("site -> $(PGT[p].site) l -> $(PGT[p].l) shell -> $(PGT[p].shell)")
            q1 = PGT[p].Pr[1]
            q2 = PGT[p].Pr[end]
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), ovlp[q, q1:q2, s])
                println()
            end
        end
    end
end

"""
    view_dm(dm::Array{F64,3})

Output the density matrix. General version
"""
function view_dm(dm::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(dm)

    # output the data
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
    view_dm(PGT::Array{PrGroupT,1}, dm::Array{F64,3})

Output the density matrix. It should be block-diagonal
"""
function view_dm(PGT::Array{PrGroupT,1}, dm::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(dm)

    # output the data
    println("<- Density Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PGT)
            println("site -> $(PGT[p].site) l -> $(PGT[p].l) shell -> $(PGT[p].shell)")
            q1 = PGT[p].Pr[1]
            q2 = PGT[p].Pr[end]
            for q = q1:q2
                foreach(x -> @printf("%12.7f", x), dm[q, q1:q2, s])
                println()
            end
        end
    end
end

"""
    view_hamk(hamk::Array{C64,3})

Output the local hamiltonian. General version
"""
function view_hamk(hamk::Array{C64,3})
    # extract some key parameters
    _, nproj, nspin = size(hamk)

    # output the data
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
    view_hamk(PGT::Array{PrGroupT,1}, hamk::Array{C64,3})

Output the local hamiltonian. It should be block-diagonal
"""
function view_hamk(PGT::Array{PrGroupT,1}, hamk::Array{C64,3})
    # extract some key parameters
    _, nproj, nspin = size(hamk)

    # output the data
    println("<- Local Hamiltonian ->")
    for s = 1:nspin
        println("Spin: $s")
        for p in eachindex(PGT)
            println("site -> $(PGT[p].site) l -> $(PGT[p].l) shell -> $(PGT[p].shell)")
            println("Re:")
            q1 = PGT[p].Pr[1]
            q2 = PGT[p].Pr[end]
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
