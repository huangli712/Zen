#
# project : pansy
# source  : tetra.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/06
#

"""
    TetraWeight

Integration weights for analytical tetrahedron algorithm

.cw -> blochl corrections for dweight
.dw -> density of states weights at the four corners of a given tetrahedron
.tw -> integration weights at the four corners of a given tetrahedron
"""
struct TetraWeight
    cw :: F64
    dw :: Array{F64,1}
    tw :: Array{F64,1}
end

"""
    tetra_p_ek1()

Blochl algorithm, case 1, for fully unoccupied tetrahedron
"""
function tetra_p_ek1()
    # integration weights, apply equation (B1)
    tw = zeros(F64, 4)

    # density of states weights
    dw = zeros(F64, 4)

    # corrections for dweight
    cw = 0.0

    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek12(z::F64, e::Array{F64,1})

Blochl algorithm, case 2, for partially occupied tetrahedron
"""
function tetra_p_ek12(z::F64, e::Array{F64,1})
    # sainty check
    @assert length(e) === 4

    # setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze1 = z - e[1]
    e21 = e[2] - e[1]
    e31 = e[3] - e[1]
    e41 = e[4] - e[1]

    # intermediate variable, apply equation (B6)
    c = ze1 * ze1 * ze1 / ( 4.0 * e21 * e31 * e41 )
    dc = 3.0 * ze1 * ze1 / ( 4.0 * e21 * e31 * e41 )

    # integration weights
    tw = zeros(F64, 4)
    #
    # apply equation (B2)
    tw[1] = c * ( 4.0 - ze1 * ( 1.0 / e21 + 1.0 / e31 + 1.0 / e41 ) )
    #
    # apply equation (B3)
    tw[2] = c * ze1 / e21
    #
    # apply equation (B4)
    tw[3] = c * ze1 / e31
    #
    # apply equation (B5)
    tw[4] = c * ze1 / e41

    # density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = 4.0 * dc - ( dc * ze1 + c ) * ( 1.0 / e21 + 1.0 / e31 + 1.0 / e41 )
    dw[2] = dc * ze1 / e21 + c / e21
    dw[3] = dc * ze1 / e31 + c / e31
    dw[4] = dc * ze1 / e41 + c / e41

    # corrections for dweight
    cw = 6.0 * ze1 / ( e21 * e31 * e41 )

    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek23(z::F64, e::Array{F64,1})

Blochl algorithm, case 3, for partially occupied tetrahedron
"""
function tetra_p_ek23(z::F64, e::Array{F64,1})
    # sainty check
    @assert length(e) === 4

    # setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze1 = z - e[1]
    ze2 = z - e[2]
    ze3 = z - e[3]
    ze4 = z - e[4]
    e21 = e[2] - e[1]
    e31 = e[3] - e[1]
    e32 = e[3] - e[2]
    e41 = e[4] - e[1]
    e42 = e[4] - e[2]

    # intermediate variables
    # apply equation (B11)
    c1  = ze1 * ze1 / ( 4.0 * e41 * e31 )
    dc1 = ze1 / ( 2.0 * e41 * e31 )

    # apply equation (B12)
    c2  = - ze1 * ze2 * ze3 / ( 4.0 * e41 * e32 * e31 )
    dc2 = ( - ze2 * ze3 - ze1 * ze3 - ze1 * ze2 ) / ( 4.0 * e41 * e32 * e31 )

    # apply equation (B13)
    c3  = - ze2 * ze2 * ze4 / ( 4.0 * e42 * e32 * e41 )
    dc3 = ( - 2.0 * ze2 * ze4 - ze2 * ze2 ) / ( 4.0 * e42 * e32 * e41 )

    # integration weights
    tw = zeros(F64, 4)
    #
    # apply equation (B7)
    tw[1] = c1 - ( c1 + c2 ) * ze3 / e31 - ( c1 + c2 + c3 ) * ze4 / e41
    #
    # apply equation (B8)
    tw[2] = c1 + c2 + c3 - ( c2 + c3 ) * ze3 / e32 - c3 * ze4 / e42
    #
    # apply equation (B9)
    tw[3] = ( c1 + c2 ) * ze1 / e31 + ( c2 + c3 ) * ze2 / e32
    #
    # apply equation (B10)
    tw[4] = ( c1 + c2 + c3 ) * ze1 / e41 + c3 * ze2 / e42

    # density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = dc1 - ( ( dc1 + dc2 ) * ze3 + c1 + c2 ) / e31 - ( ( dc1 + dc2 + dc3 ) * ze4 + c1 + c2 + c3 ) / e41
    dw[2] = dc1 + dc2 + dc3 - ( ( dc2 + dc3 ) * ze3 + c2 + c3 ) / e32 - ( dc3 * ze4 + c3 ) / e42
    dw[3] = ( ( dc1 + dc2 ) * ze1 + c1 + c2 ) / e31 + ( ( dc2 + dc3 ) * ze2 + c2 + c3 ) / e32
    dw[4] = ( ( dc1 + dc2 + dc3 ) * ze1 + c1 + c2 + c3 ) / e41 + ( dc3 * ze2 + c3 ) / e42

    # corrections for dweight
    cw = 6.0 * ( 1.0 - ( e31 + e42 ) * ze2 / ( e32 * e42 ) ) / ( e31 * e41 )

    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek34(z::F64, e::Array{F64,1})

Blochl algorithm, case 4, for partially occupied tetrahedron
"""
function tetra_p_ek34(z::F64, e::Array{F64,1})
    # sainty check
    @assert length(e) === 4

    # setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze4 = z - e[4]
    e41 = e[4] - e[1]
    e42 = e[4] - e[2]
    e43 = e[4] - e[3]

    # intermediate variables, apply equation (B18)
    c = - ze4 * ze4 * ze4 / ( 4.0 * e41 * e42 * e43 )
    dc = - 3.0 * ze4 * ze4 / ( 4.0 * e41 * e42 * e43 )

    # integration weights
    tw = zeros(F64, 4)
    #
    # apply equation (B14)
    tw[1] = 0.25 + c * ze4 / e41
    #
    # apply equation (B15)
    tw[2] = 0.25 + c * ze4 / e42
    #
    # apply equation (B16)
    tw[3] = 0.25 + c * ze4 / e43
    #
    # apply equation (B17)
    tw[4] = 0.25 - c * ( 4.0 + ( 1.0 / e41 + 1.0 / e42 + 1.0 / e43 ) * ze4 )

    # density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = ( dc * ze4 + c ) / e41
    dw[2] = ( dc * ze4 + c ) / e42
    dw[3] = ( dc * ze4 + c ) / e43
    dw[4] = - 4.0 * dc - ( 1.0 / e41 + 1.0 / e42 + 1.0 / e43) * ( dc * ze4 + c )

    # corrections for dweight
    cw = 6.0 * ze4 / ( e41 * e42 * e43 )

    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek4()

Blochl algorithm, case 5, for fully occupied tetrahedron
"""
function tetra_p_ek4()
    # integration weights, apply equation (B19)
    tw = fill(0.25, 4)

    # density of states weights
    dw = zeros(F64, 4)

    # corrections for dweight
    cw = zero

    TetraWeight(cw, dw, tw)
end

function tetra_weight1()
end

"""
    tetra_weight2(z::F64, e::Array{F64,1})

Peter E. Blochl algorithm for (integrated) density of states and relevant
integration weights. Blochl corrections are taken into considersions as
well. See Phys. Rev. B, 49, 16223 (1994) for more details
"""
function tetra_weight2(z::F64, e::Array{F64,1})
    # sort the corner energy according to increasing values
    sort!(e)

    # remove possible degenerancies in e
    for i = 1:3
        if abs( e[i] - e[i+1] ) < eps(F64)
            e[i] = e[i] + eps(F64) / float(i)
        end
    end

    for i = 1:4
        if abs( e[i] - z ) < eps(F64) / 10.0
            e[i] = e[i] + eps(F64) / 10.0 / float(i)
        end
    end

    # find the case, to calculate TetraWeight (dw, tw, and cw)
    # case 1, fully unoccupied tetrahedron
    if z < e[1]
        TW = tetra_p_ek1()

    # case 2, partially occupied tetrahedron
    elseif z < e[2] && z > e[1]
        TW = tetra_p_ek12(z, e)

    # case 3, partially occupied tetrahedron
    elseif z < e[3] && z > e[2]
        TW = tetra_p_ek23(z, e)

    # case 4, partially occupied tetrahedron
    elseif z < e[4] && z > e[3]
        TW = tetra_p_ek34(z, e)

    # case 5, fully occupied tetrahedron
    elseif z > e[4]
        TW = tetra_p_ek4()

    end

    # add up Blochl corrections for density of states weights
    # apply equation (22)
    for i = 1:4
        do j = 1:4
            TW.dw[i] = TW.dw[i] + ( e[j] - e[i] ) * TW.cw * 0.025
        end
    end

    # add up Blochl corrections for integration weights
    # apply equation (22)
    for i = 1:4
        for j = 1:4
            TW.tw[i] = TW.tw[i] + ( e[j] - e[i] ) * sum(TW.dw) * 0.025
        end
    end

    return TW
end
