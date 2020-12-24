#
# project : pansy
# source  : plo.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/24
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

                # setup Tr array further
                @cswitch str begin
                    @case "s"
                        PG[g].Tr = Diagonal(fill(1.0, 1))
                        break

                    @case "p"
                        PG[g].Tr = Diagonal(fill(1.0, 3))
                        break

                    @case "d"
                        PG[g].Tr = Diagonal(fill(1.0, 5))
                        break

                    @case "f"
                        PG[g].Tr = Diagonal(fill(1.0, 7))
                        break

                    @case "d_t2g"
                        PG[g].Tr = zeros(F64, 3, 5)
                        PG[g].Tr[1, 1] = 1.0 
                        PG[g].Tr[2, 2] = 1.0 
                        PG[g].Tr[3, 4] = 1.0 
                        break

                    @case "d_eg"
                        PG[g].Tr = zeros(F64, 2, 5)
                        PG[g].Tr[1, 3] = 1.0
                        PG[g].Tr[2, 1] = 1.0
                        break

                    @default
                        sorry()
                        break
                end 
            end
        end
    end

    for i in eachindex(PG)
        @show PG[i]
    end
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
