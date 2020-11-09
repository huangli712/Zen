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

function dft_driver(it::IterInfo, d::Dict{String,Any})
    cd("dft")

    println(it.dft_iter)
    println(it.dmft_iter)
    println(it.dft_dmft_iter)

    it.dft_iter += 1

    cd("..")
end

function dmft_driver(d::Dict{String,Any})
end
