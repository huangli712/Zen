#!/usr/bin/env julia

include("Zen.jl")
using .Zen

# check the version of julia runtime environment
require()

# print the welcome message
welcome()

# parse the file case.toml to extract configuration parameters
cfg = parse_toml(query_args(), true)

# validate the configuration parameters
renew_config(cfg)
check_config()

# write the configuration parameters to stdout
message("Job Summary")
view_case()
view_dft()
view_dmft()
view_impurity()
view_solver()

# check the input files (which are essential for the calculation)
message("zen", "examine the essential input files")
query_cars()

# prepare the working directories
message("zen", "create the working directories")
make_trees()

# create a IterInfo object
message("zen", "make the instance of scf-consistent iterator")
it = IterInfo()

exit(-1)

if _m("mode") === 1

    message("zen", "enter one-shot mode")
    message("zen", "begin < dft block >")
    message("zen", "dft -> init")
    dft_init(it)
    message("zen", "dft -> run")
    dft_run(it)
    message("zen", "dft -> save")
    message("zen", "e_n_d < dft block >")
    dft_save(it)
    for iter = 1:_m("niter")
        message("zen", "dmft_cycle -> 0  dmft1_iter -> $iter dmft2_iter -> 0")
    end

else
    sorry()
end

goodbye()
