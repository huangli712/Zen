"""
    make_trees()

Build the working directories at advance 
"""
function make_trees(d::Dict{String,Any})
    mkdir("dft")
    mkdir("dmft1")
    mkdir("dmft2")
    for i = 1:d["nimp"]
        mkdir("impurity.$i")
    end
end

function dft_driver(d::Dict{String,Any})
end
