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
"""
function sigma_reset(lr::Logger)
    prompt(lr.log, "sigma::reset")
    prompt("Sigma : Reset")

    # Create frequency mesh
    fmesh = zeros(F64, get_m("nmesh"))
    if get_m("axis") === 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / get_m("beta")
        end
    else # Real axis
        sorry()
    end
    @show fmesh
end

function sigma_dcount(lr::Logger)
end

function sigma_split(lr::Logger)
end

function sigma_gather(lr::Logger)
end
