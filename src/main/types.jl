mutable struct IterInfo
    total_iter
    dmft1_iter
    dmft2_iter
    dft_dmft_iter
end

function IterInfo(iter::Int64 = 0)
    IterInfo(iter, iter, iter, iter)
end

mutable struct Lattice
end

mutable struct DFTData
end
