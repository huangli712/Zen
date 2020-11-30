#!/usr/bin/env julia

include("Zen.jl")
using .Zen

# check the version of julia runtime environment
require()

# print the welcome message 
welcome()

# parse the file case.toml to extract parameters
message("zen", "parse the configuration file")
cfg = parse_toml(query_args(), true)
@show cfg
exit(-1)
case = cfg["case"]
dft = cfg["dft"]
dmft = cfg["dmft"]
dft_dmft = cfg["dft_dmft"]

# validate the parameters
message("zen", "check the configuration parameters")
# plug your codes here

# write the parameters to stdout
message("zen", "display the configuration parameters")
param_case(case)
param_dft(dft)
param_dmft(dmft)
param_dft_dmft(dft_dmft)

# check the input files (which are essential for the calculation)
message("zen", "examine the essential input files")
query_cars(dft)

# prepare the working directories
message("zen", "create the working directories")
make_trees(dmft["impurity"])

# create a IterInfo object
message("zen", "make the instance of iterator")
it = IterInfo()

if dft_dmft["mode"] == 1

    message("zen", "enter one-shot mode")
    message("zen", "begin < dft block >")
    message("zen", "dft -> init")
    dft_init(it, case, dft)
    message("zen", "dft -> run")
    dft_run(it, dft)
    message("zen", "dft -> save")
    message("zen", "e_n_d < dft block >")
    dft_save(it, dft)
    for iter in 1:dft_dmft["niter"]
        message("zen", "dft_dmft_iter -> 0  dmft1_iter -> $iter dmft2_iter -> 0")
    end

else
    sorry()
end

goodbye()
