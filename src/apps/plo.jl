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

Calculate the overlap of projectors
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
            B = conj(transpose(A))
            ovlp[:, :, s] = ovlp[:, :, s] + real(A * B) * wght
        end
    end

    @show ovlp[:, :, 1]
end

"""
    plo_dm()
"""
function plo_dm()
end
