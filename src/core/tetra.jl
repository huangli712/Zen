#
# project : pansy
# source  : tetra.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/04
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
    # apply equation (B1)
    cw = 0.0
    dw = zeros(F64, 4)
    tw = zeros(F64, 4)

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
    # ze1: ze_{i} = e - e_{i}
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
    dw[1] = 4.0 * dc - ( dc * ze1 + c ) * ( 1.0 / e21 + 1.0 / e31 + 1.0 / e41 )
    dw[2] = dc * ze1 / e21 + c / e21
    dw[3] = dc * ze1 / e31 + c / e31
    dw[4] = dc * ze1 / e41 + c / e41

    # corrections for dweight
    cw = 6.0 * ze1 / ( e21 * e31 * e41 )

    TetraWeight(cw, dw, tw)
end

function tetra_p_ek23
end

function tetra_p_ek34
end

function tetra_p_ek4
end

function tetra_weight1()
end

function tetra_weight2()
end
