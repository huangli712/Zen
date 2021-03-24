#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/24
#

"""
    sigma_reset(lr::Logger)

Create initial self-energy functions and write them to `sigma.bare`. 

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset(lr::Logger)
    prompt(lr.log, "sigma::reset")
    prompt("Sigma : Reset")

    # Create frequency mesh
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    fmesh = zeros(F64, nmesh)
    if axis === 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / beta
        end
    else # Real axis
        sorry()
    end

    # Create self-energy functions
    # sigma = zeros(C64, nmesh, nsite)
end

"""
    sigma_dcount(lr::Logger)

Calculate double counting terms for self-energy functions and write
them to `sigma.dc`.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(lr::Logger)
end

function sigma_split(lr::Logger)
end

function sigma_gather(lr::Logger)
end
