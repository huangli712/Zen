#!/usr/bin/env julia

include("zen.jl")
using .Zen

# check the version of julia runtime environment
require()

# print the welcome message 
welcome()

# parse the file case.toml to extract parameters
case = parse_toml(query_args(), "case", true) 
dft = parse_toml(query_args(), "dft", true)
dmft = parse_toml(query_args(), "dmft", true)
dft_dmft = parse_toml(query_args(), "dft_dmft", true)

# validate the parameters
# plug your codes here

# check the input files (which are essential for the calculation)
query_cars(dft)

# prepare the working directories
make_trees(dmft["impurity"])

# create a IterInfo object
it = IterInfo(0, 0, 0, 0)

if dft_dmft["mode"] == 1

    dft_init(it, case, dft)
    dft_run(it, dft)
    dft_save(it, dft)
    for iter in 1:dft_dmft["niter"]
        println("iter: $iter")
    end

else
    sorry()
end

goodbye()
