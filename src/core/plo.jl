#
# project : pansy
# source  : plo.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/30
#

"""
    plo_group(PG::Array{PrGroup,1})

Use the information contained in the PIMP dict to further complete
the PrGroup struct
"""
function plo_group(PG::Array{PrGroup,1})
    # check the parameters contained in PIMP dict
    @assert _i("nsite") === length(_i("atoms"))
    @assert _i("nsite") === length(_i("shell"))

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
    for i = 1:_i("nsite")
        # determine site
        str = _i("atoms")[i]
        site = parse(I64, line_to_array(str)[3])

        # determine l
        str = _i("shell")[i]
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
    for spin = 1:nspin
        for kpt = 1:nkpt
            # for lower boundary
            ib1 = 1
            while enk[ib1, kpt, spin] < emin
                ib1 = ib1 + 1
            end

            # for upper boundary
            ib2 = nband
            while enk[ib2, kpt, spin] > emax
                ib2 = ib2 - 1
            end

            # check the boundaries
            @assert ib1 <= ib2

            # save the boundaries
            ib_window[kpt, spin, 1] = ib1
            ib_window[kpt, spin, 2] = ib2
        end
    end

    # try to find out the minimum and maximum band indices
    ib_min = minimum(ib_window[:, :, 1])
    ib_max = maximum(ib_window[:, :, 2])

    # try to find out the maximum number of selected bands
    nbmax = ib_max - ib_min + 1

    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # create arrays
    # chipsi_ is used to store the required projectors
    chipsi_ = zeros(C64, nproj, nbmax, nkpt, nspin)

    # select projectors which lie in the given window
    # copy data from chipsi to chipsi_
    for spin = 1:nspin
        for kpt = 1:nkpt
            ib1 = ib_window[kpt, spin, 1]
            ib2 = ib_window[kpt, spin, 2]
            ib3 = ib2 - ib1 + 1
            @assert ib3 <= nbmax
            chipsi_[:, 1:ib3, kpt, spin] = chipsi[:, ib1:ib2, kpt, spin]
        end
    end

    # return the desired arrays
    return ib_min, ib_max, ib_window, chipsi_
end

"""
    plo_orthog()

Try to orthogonalize the projectors group by group (site by site)
"""
function plo_orthog(window::Array{I64,3}, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)
    #@show nproj, nband, nkpt, nspin

    # create arrays
    TmpMat = zeros(C64, nproj, nband)

    for spin = 1:nspin
        for kpt = 1:nkpt
            fill!(TmpMat, 0.0 + 0.0im)
            b1 = window[kpt, spin, 1]
            b2 = window[kpt, spin, 2]
            nb = b2 - b1 + 1
            for p in eachindex(PGT)
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]
                TmpMat[q1:q2, 1:nb] = chipsi[q1:q2, 1:nb, kpt, spin]
                #println("here")
                #@show TmpMat[q1:q2, 1:nb]
                STmpMat = plo_diag(TmpMat[q1:q2, 1:nb])
                #@show STmpMat
                chipsi[q1:q2, 1:nb, kpt, spin] = STmpMat
            end
        end
    end
end

"""
    plo_diag(M::Array{C64,2})
"""
function plo_diag(M::Array{C64,2})
    ovlp = M * M'
    @assert ishermitian(ovlp)

    vals, vecs = eigen(Hermitian(ovlp))
    @assert all(vals .> 0)

    sqrt_vals = map(x -> 1.0 / sqrt(x), vals)
    S = vecs * Diagonal(sqrt_vals) * vecs'

    return S * M
end

"""
    plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap out of projectors
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
            A = chipsi[:, :, k, s]
            ovlp[:, :, s] = ovlp[:, :, s] + real(A * A') * wght
        end
    end

    # return the desired array
    return ovlp
end

"""
    plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})

Calculate the overlap out of projectors
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
                A = chipsi[q1:q2, :, k, s]
                ovlp[q1:q2, q1:q2, s] = ovlp[q1:q2, :, s] + real(A * A') * wght
            end
        end
    end
    @show ovlp[:, :, 1]

    # return the desired array
    return ovlp
end

"""
    plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors
"""
function plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # evaluate spin factor
    sf = (nspin == 1 ? 2 : 1)

    # create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[:, k, s]
            A = chipsi[:, :, k, s]
            dm[:, :, s] = dm[:, :, s] + real(A * Diagonal(occs) * A') * wght
        end
    end

    # return the desired array
    return dm
end

"""
    plo_dm(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})

Calculate the density matrix out of projectors
"""
function plo_dm(bmin::I64, bmax::I64, PGT::Array{PrGroupT,1}, chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)
    @assert nband === bmax - bmin + 1

    # evaluate spin factor
    sf = (nspin == 1 ? 2 : 1)

    # create density matrix array
    dm = zeros(F64, nproj, nproj, nspin)

    # build overlap array
    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt * sf
            occs = occupy[bmin:bmax, k, s]
            for p in eachindex(PGT)
                q1 = PGT[p].Pr[1]
                q2 = PGT[p].Pr[end]
                A = chipsi[q1:q2, :, k, s]
                dm[q1:q2, q1:q2, s] = dm[q1:q2, q1:q2, s] + real(A * Diagonal(occs) * A') * wght
            end
        end
    end
    @show dm[:, :, 1]

    # return the desired array
    return dm
end

"""
    view_ovlp(ovlp::Array{F64,3})

Output the overlap matrix, only for debug
"""
function view_ovlp(ovlp::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(ovlp)

    # output the data
    println("<- Overlap Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p1 = 1:nproj
            map(x -> @printf("%12.7f", x), ovlp[p1, :, s])
            println()
        end
    end
end

"""
    view_ovlp()
"""
function view_ovlp()
end

"""
    view_dm(dm::Array{F64,3})

Output the density matrix, only for debug
"""
function view_dm(dm::Array{F64,3})
    # extract some key parameters
    _, nproj, nspin = size(dm)

    # output the data
    println("<- Density Matrix ->")
    for s = 1:nspin
        println("Spin: $s")
        for p1 = 1:nproj
            map(x -> @printf("%12.7f", x), dm[p1, :, s])
            println()
        end
    end
end

"""
    view_dm()
"""
function view_dm()
end
