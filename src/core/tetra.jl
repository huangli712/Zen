#
# Project : Pansy
# Source  : tetra.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/03/02
#

#
# Remarks:
#
# This file provides some functions to do the brillouin zone integration.
# They are very useful for calculating the density of states. So far the
# following algorithms are supported:
#     (1) Gaussian broadening method.
#     (2) Fermi-Dirac broadening method.
#     (3) Analytical tetrahedron algorithm with Blochl corrections.
# Note that you have to modify the line 87-89 to choose suitable driver.
# Perhaps you also need to modify the `gamm` parameter to obtain more
# reasonable results. Now the default algorithm is (3).
#

#
# Customized Structs
#

"""
    TetraWeight

Struct. Integration weights for analytical tetrahedron algorithm.

## Members

.cw -> Blochl corrections for `dw`.\n
.dw -> Density of states weights at the four corners of a given tetrahedron.\n
.tw -> Integration weights at the four corners of a given tetrahedron.\n
"""
struct TetraWeight
    cw :: F64
    dw :: Array{F64,1}
    tw :: Array{F64,1}
end

#
# Driver Functions
#

"""
    bzint(z::F64, itet::Array{I64,2}, enk::Array{F64,3})

Compute tetrahedron integrated weights for Brillouin zone integration.
It is used to calculate (partial) density of states.

See also: [`gauss_weight`](@ref), [`fermi_weight`](@ref), [`tetra_weight`](@ref).
"""
function bzint(z::F64, itet::Array{I64,2}, enk::Array{F64,3})
    # Extract some key parameters
    ntet, ndim = size(itet)
    @assert ndim === 5

    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Energy at the corner of tetrahedron
    zc = zeros(F64, 4)

    # Mass of tetrahedron
    mtet = sum( itet[:, 1] )

    # Initialize weights. It shares the same size with `enk`.
    W = zeros(F64, nband, nkpt, nspin)

    # Go through each tetrahedron, spin, and band.
    for t = 1:ntet
        for s = 1:nspin
            for b = 1:nband

                # Store the four corner energies of one tetrahedron in `zc`
                for c = 1:4
                    k = itet[t, c + 1]
                    zc[c] = enk[b, k, s]
                end

                # Actually calculates weights for 4 corners of the tetrahedron
                #TW = gauss_weight(z, zc)
                #TW = fermi_weight(z, zc)
                TW = tetra_weight(z, zc)

                # Stores integration weights for irreducible k-points
                for c = 1:4
                    k = itet[t, c + 1]
                    W[b, k, s] = W[b, k, s] + TW.dw[c] * float(itet[t, 1])
                end
            end
        end
    end

    # Normalize the weights properly
    @. W = W / float(mtet)

    # Return the desired arrays
    return W
end

#
# Service Functions (Level 1)
#

"""
    gauss_weight(z::F64, e::Array{F64,1})

Gaussian broadening algorithm for (integrated) density of states and relevant
integration weights.

See also: [`TetraWeight`](@ref).
"""
function gauss_weight(z::F64, e::Array{F64,1})
    # Integration weights, apply equation (B1).
    tw = zeros(F64, 4)

    # Density of states weights
    dw = zeros(F64, 4)

    # Corrections for dweight
    cw = 0.0

    # Create TetraWeight struct
    TW = TetraWeight(cw, dw, tw)

    # Parameter for gaussian broadening
    gamm = 0.25

    # Further setup the integration weights
    for i = 1:4
        dummy = ( z - e[i] ) / gamm
        TW.tw[i] = 0.125 * ( 1.0 - erf(-dummy) )
        TW.dw[i] = 0.25 * exp(-dummy^2.0) / ( sqrt(pi) * gamm )
    end

    # Return the TetraWeight struct
    return TW
end

"""
    fermi_weight(z::F64, e::Array{F64,1})

Fermi-Dirac broadening algorithm for (integrated) density of states and
relevant integration weights.

See also: [`TetraWeight`](@ref).
"""
function fermi_weight(z::F64, e::Array{F64,1})
    # Integration weights, apply equation (B1)
    tw = zeros(F64, 4)

    # Density of states weights
    dw = zeros(F64, 4)

    # Corrections for dweight
    cw = 0.0

    # Create TetraWeight struct
    TW = TetraWeight(cw, dw, tw)

    # Parameter for gaussian broadening
    gamm = 0.25

    # Further setup the integration weights
    for i = 1:4
        dummy = ( z - e[i] ) / gamm
        TW.tw[i] = 1.0 / (1.0 +  exp(-dummy) )
        TW.dw[i] = 0.5 / (1.0 + cosh( dummy) )
    end

    # Return the TetraWeight struct
    return TW
end

"""
    tetra_weight(z::F64, e::Array{F64,1})

Peter E. Blochl algorithm for (integrated) density of states and relevant
integration weights. Blochl corrections are taken into considersions as
well. See Phys. Rev. B, 49, 16223 (1994) for more details.

See also: [`TetraWeight`](@ref).
"""
function tetra_weight(z::F64, e::Array{F64,1})
    # Sort the corner energy according to increasing values
    sort!(e)

    # Remove possible degenerancies in e
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

    # Find the case to calculate TetraWeight (including `dw`, `tw`, and `cw`)
    #
    # Case 1, fully unoccupied tetrahedron.
    if z < e[1]
        TW = tetra_p_ek1()
    #
    # Case 2, partially occupied tetrahedron.
    elseif z < e[2] && z > e[1]
        TW = tetra_p_ek12(z, e)

    # Case 3, partially occupied tetrahedron.
    elseif z < e[3] && z > e[2]
        TW = tetra_p_ek23(z, e)
    #
    # Case 4, partially occupied tetrahedron.
    elseif z < e[4] && z > e[3]
        TW = tetra_p_ek34(z, e)
    #
    # Case 5, fully occupied tetrahedron.
    elseif z > e[4]
        TW = tetra_p_ek4()
    #
    end

    # Add up Blochl corrections for density of states weights
    # Apply equation (22)
    for i = 1:4
        for j = 1:4
            TW.dw[i] = TW.dw[i] + ( e[j] - e[i] ) * TW.cw * 0.025
        end
    end

    # Add up Blochl corrections for integration weights
    # Apply equation (22)
    for i = 1:4
        for j = 1:4
            TW.tw[i] = TW.tw[i] + ( e[j] - e[i] ) * sum(TW.dw) * 0.025
        end
    end

    # Return the TetraWeight struct
    return TW
end

#
# Service Functions (Level 2)
#

"""
    tetra_p_ek1()

Blochl algorithm, case 1, for fully unoccupied tetrahedron.

See also: [`tetra_weight`](@ref).
"""
function tetra_p_ek1()
    # Integration weights, apply equation (B1)
    tw = zeros(F64, 4)

    # Density of states weights
    dw = zeros(F64, 4)

    # Corrections for dweight
    cw = 0.0

    # Return TetraWeight struct
    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek12(z::F64, e::Array{F64,1})

Blochl algorithm, case 2, for partially occupied tetrahedron.

See also: [`tetra_weight`](@ref).
"""
function tetra_p_ek12(z::F64, e::Array{F64,1})
    # Sainty check
    @assert length(e) === 4

    # Setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze1 = z - e[1]
    #
    e21 = e[2] - e[1]
    e31 = e[3] - e[1]
    e41 = e[4] - e[1]

    # Intermediate variable, apply equation (B6)
    c = ze1 * ze1 * ze1 / ( 4.0 * e21 * e31 * e41 )
    dc = 3.0 * ze1 * ze1 / ( 4.0 * e21 * e31 * e41 )

    # Integration weights
    tw = zeros(F64, 4)
    #
    # Apply equation (B2)
    tw[1] = c * ( 4.0 - ze1 * ( 1.0 / e21 + 1.0 / e31 + 1.0 / e41 ) )
    #
    # Apply equation (B3)
    tw[2] = c * ze1 / e21
    #
    # Apply equation (B4)
    tw[3] = c * ze1 / e31
    #
    # Apply equation (B5)
    tw[4] = c * ze1 / e41

    # Density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = 4.0 * dc - ( dc * ze1 + c ) * ( 1.0 / e21 + 1.0 / e31 + 1.0 / e41 )
    dw[2] = dc * ze1 / e21 + c / e21
    dw[3] = dc * ze1 / e31 + c / e31
    dw[4] = dc * ze1 / e41 + c / e41

    # Corrections for dweight
    cw = 6.0 * ze1 / ( e21 * e31 * e41 )

    # Return TetraWeight struct
    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek23(z::F64, e::Array{F64,1})

Blochl algorithm, case 3, for partially occupied tetrahedron.

See also: [`tetra_weight`](@ref).
"""
function tetra_p_ek23(z::F64, e::Array{F64,1})
    # Sainty check
    @assert length(e) === 4

    # Setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze1 = z - e[1]
    ze2 = z - e[2]
    ze3 = z - e[3]
    ze4 = z - e[4]
    #
    e21 = e[2] - e[1]
    e31 = e[3] - e[1]
    e32 = e[3] - e[2]
    e41 = e[4] - e[1]
    e42 = e[4] - e[2]

    # Intermediate variables
    # Apply equation (B11)
    c1  = ze1 * ze1 / ( 4.0 * e41 * e31 )
    dc1 = ze1 / ( 2.0 * e41 * e31 )
    #
    # Apply equation (B12)
    c2  = - ze1 * ze2 * ze3 / ( 4.0 * e41 * e32 * e31 )
    dc2 = ( - ze2 * ze3 - ze1 * ze3 - ze1 * ze2 ) / ( 4.0 * e41 * e32 * e31 )
    #
    # Apply equation (B13)
    c3  = - ze2 * ze2 * ze4 / ( 4.0 * e42 * e32 * e41 )
    dc3 = ( - 2.0 * ze2 * ze4 - ze2 * ze2 ) / ( 4.0 * e42 * e32 * e41 )

    # Integration weights
    tw = zeros(F64, 4)
    #
    # Apply equation (B7)
    tw[1] = c1 - ( c1 + c2 ) * ze3 / e31 - ( c1 + c2 + c3 ) * ze4 / e41
    #
    # Apply equation (B8)
    tw[2] = c1 + c2 + c3 - ( c2 + c3 ) * ze3 / e32 - c3 * ze4 / e42
    #
    # Apply equation (B9)
    tw[3] = ( c1 + c2 ) * ze1 / e31 + ( c2 + c3 ) * ze2 / e32
    #
    # Apply equation (B10)
    tw[4] = ( c1 + c2 + c3 ) * ze1 / e41 + c3 * ze2 / e42

    # Density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = dc1 - ( ( dc1 + dc2 ) * ze3 + c1 + c2 ) / e31 - ( ( dc1 + dc2 + dc3 ) * ze4 + c1 + c2 + c3 ) / e41
    dw[2] = dc1 + dc2 + dc3 - ( ( dc2 + dc3 ) * ze3 + c2 + c3 ) / e32 - ( dc3 * ze4 + c3 ) / e42
    dw[3] = ( ( dc1 + dc2 ) * ze1 + c1 + c2 ) / e31 + ( ( dc2 + dc3 ) * ze2 + c2 + c3 ) / e32
    dw[4] = ( ( dc1 + dc2 + dc3 ) * ze1 + c1 + c2 + c3 ) / e41 + ( dc3 * ze2 + c3 ) / e42

    # Corrections for dweight
    cw = 6.0 * ( 1.0 - ( e31 + e42 ) * ze2 / ( e32 * e42 ) ) / ( e31 * e41 )

    # Return TetraWeight struct
    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek34(z::F64, e::Array{F64,1})

Blochl algorithm, case 4, for partially occupied tetrahedron.

See also: [`tetra_weight`](@ref).
"""
function tetra_p_ek34(z::F64, e::Array{F64,1})
    # Sainty check
    @assert length(e) === 4

    # Setup common variables
    # zei: ze_{i} = e - e_{i}
    # eij: e_{ij} = e_{i} - e_{j}
    ze4 = z - e[4]
    #
    e41 = e[4] - e[1]
    e42 = e[4] - e[2]
    e43 = e[4] - e[3]

    # Intermediate variables, apply equation (B18)
    c = - ze4 * ze4 * ze4 / ( 4.0 * e41 * e42 * e43 )
    dc = - 3.0 * ze4 * ze4 / ( 4.0 * e41 * e42 * e43 )

    # Integration weights
    tw = zeros(F64, 4)
    #
    # Apply equation (B14)
    tw[1] = 0.25 + c * ze4 / e41
    #
    # Apply equation (B15)
    tw[2] = 0.25 + c * ze4 / e42
    #
    # Apply equation (B16)
    tw[3] = 0.25 + c * ze4 / e43
    #
    # Apply equation (B17)
    tw[4] = 0.25 - c * ( 4.0 + ( 1.0 / e41 + 1.0 / e42 + 1.0 / e43 ) * ze4 )

    # Density of states weights
    dw = zeros(F64, 4)
    #
    dw[1] = ( dc * ze4 + c ) / e41
    dw[2] = ( dc * ze4 + c ) / e42
    dw[3] = ( dc * ze4 + c ) / e43
    dw[4] = - 4.0 * dc - ( 1.0 / e41 + 1.0 / e42 + 1.0 / e43) * ( dc * ze4 + c )

    # Corrections for dweight
    cw = 6.0 * ze4 / ( e41 * e42 * e43 )

    # Return TetraWeight struct
    TetraWeight(cw, dw, tw)
end

"""
    tetra_p_ek4()

Blochl algorithm, case 5, for fully occupied tetrahedron.

See also: [`tetra_weight`](@ref).
"""
function tetra_p_ek4()
    # Integration weights, apply equation (B19)
    tw = fill(0.25, 4)

    # Density of states weights
    dw = zeros(F64, 4)

    # Corrections for dweight
    cw = 0.0

    # Return TetraWeight struct
    TetraWeight(cw, dw, tw)
end
