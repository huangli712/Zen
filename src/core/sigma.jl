#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/23
#

"""
    sigma_reset(lr::Logger)
"""
function sigma_reset(lr::Logger)
    prompt(lr.log, "sigma::reset")
    prompt("Sigma : Reset")
    println("here")
end
