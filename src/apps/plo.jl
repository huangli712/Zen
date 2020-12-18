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
function plo_ovlp(chipsi::Array{C64,4}, weight::Array{F64,1})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    ovlp = zeros(F64, nproj, nproj, nspin)

    for s = 1:nspin
        for k = 1:nkpt
            wght = weight[k] / nkpt
            ovlp[:, :, s] = ovlp[:, :, s] + real(chipsi[:, :, k, s] * transpose(chipsi[:, :, k, s])) * wght
        end
    end
end

"""
    plo_dm()
"""
function plo_dm()
end
