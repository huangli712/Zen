"""
    PCASE

Dictionary for parameters: case
"""
PCASE = Dict{String,Any}
        (
            "case" => ["unknown", "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}
        (
            "engine" => [],
            "smear"  => [],
            "kmesh"  => [],
            "magmom" => [],
            "lsymm"  => [],
            "lspins" => [],
            "lspinorb" => [],
            "window" => [],
            "loptim" => [],
            "lproj"  => [],
            "nproj"  => [],
            "sproj"  => []
        )

"""
    PDMFT

Dictionary for parameters: dynamical mean-field theory calculations
"""
PDMFT = Dict{String,Any}
        (
        )

"""
    PCORE

Dictionary for parameters: self-consistent calculations
"""
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
