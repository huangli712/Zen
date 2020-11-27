PCASE = Dict{String,Any}
        (
            "case" => ["case", "system's name"]
        )

PDFT  = Dict{String,Any}
        (
            "engine" => [],
            "smear"  => [],
            "kmesh"  => [],
        )

PDMFT = Dict{String,Any}
        (
        )

PCORE = Dict{String,Any}
        (
        )


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
