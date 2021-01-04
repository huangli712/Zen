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
    cw = 0.0
    dw = zeros(F64, 4)
    tw = zeros(F64, 4)

    TetraWeight(cw, dw, tw)
end

function tetra_p_ek12
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
