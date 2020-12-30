#!/usr/bin/env julia

#
# project : pansy
# source  : run.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/30
#

include("Zen.jl")
using .Zen

# S-1: check the version of julia runtime environment
require()

# S00: print the welcome message
welcome()

# S01: print the overview of zen
message("ZEN", "Overview")
overview()

# S02: parse the configuration file, get job's description
message("ZEN", "Parsing")
#
# S02.1: parse the case.toml file to extract configuration parameters
cfg = parse_toml(query_args(), true)
println("Extract configuration parameters from $(query_args())")
#
# S02.2: build the configuration dictionaries
renew_config(cfg)
println("Update configuration dictionaries")
#
# S02.3: validate the configuration dictionaries
check_config()
println("Verify configuration dictionaries\n")

# S03: write the configuration parameters to stdout
message("ZEN", "Listing")
#
# S03.1: show dict PCASE
view_case()
#
# S03.2: show dict PDFT
view_dft()
#
# S03.3: show dict PDMFT
view_dmft()
#
# S03.4: show dict PIMP
view_impurity()
#
# S03.5: show dict PSOLVER
view_solver()

# S04: initialize the job 
message("ZEN", "Initializing")
#
# S04.1: check the input files (which are essential for the calculation)
query_inps()
println("Check essential input files")
#
# S04.2: prepare the working directories
make_trees()
println("Create working directories")
#
# S04.3: create an IterInfo object
it = IterInfo()
println("Create self-consistent iterator")
#
# S04.4: check self-consistent mode
mode = _m("mode") === 1 ? "one-shot" : "fully self-consistent"
println("Check self-consistent mode: $mode\n")

message("ZEN", "Launching")

if _m("mode") === 1
    message("ZEN", "DFT")
    dft_init(it)
    dft_run(it)
    dft_save(it)
else
    sorry()
end

goodbye()
