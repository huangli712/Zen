#
# project : pansy
# source  : plo.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/18
#

"""
    plo_ovlp()
"""
function plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1}, occupy::Array{F64,3})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    ovlp = zeros(F64, nproj, nproj, nspin)
end

"""
    plo_dm()
"""
function plo_dm()
end
