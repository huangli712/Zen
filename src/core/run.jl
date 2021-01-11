#!/usr/bin/env julia

#
# project : pansy
# source  : run.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/04
#

include("Zen.jl")
using .Zen

# S-1: check the version of julia runtime environment
require()

# S00: print the welcome message
welcome()

# S01: print the overview for zen
message("ZEN", "OVERVIEW")
overview()

# S02: parse the configuration file, get job's description
message("ZEN", "PARSER")
#
# S02.1: parse the case.toml file to extract configuration parameters
println("Extract configuration parameters from $(query_args())")
cfg = parse_toml(query_args(), true)
#
# S02.2: build the configuration dictionaries
println("Build configuration dictionaries")
renew_config(cfg)
#
# S02.3: validate the configuration dictionaries
println("Verify configuration dictionaries\n")
check_config()

# S03: write the configuration parameters to stdout
message("ZEN", "VIEWER")
#
# S03.1: show dict PCASE
@timev list_case()
#
# S03.2: show dict PDFT
@timev list_dft()
#
# S03.3: show dict PDMFT
@timev list_dmft()
#
# S03.4: show dict PIMP
@timev list_impurity()
#
# S03.5: show dict PSOLVER
@timev list_solver()
exit(-1)
# S04: initialize the job
message("ZEN", "CREATOR")
#
# S04.1: check the input files (which are essential for the calculation)
println("Ensure essential input files")
query_inps()
#
# S04.2: prepare the working directories
println("Create working directories")
make_trees()
#
# S04.3: create an IterInfo object
println("Create self-consistent iterator")
it = IterInfo()
#
# S04.4: check self-consistent mode
mode = _m("mode") === 1 ? "one-shot" : "fully self-consistent"
println("Verify self-consistent mode: $mode\n")

message("ZEN", "START")

# well, we choose the one-shot mode
if _m("mode") === 1

#
# remarks:
#
# we would like to perform two successive dft runs if _d("loptim") is true.
# the purpose of the first dft run is to evaluate the fermi level. then a
# energy window is determined. we will use this window to generate optimal
# projectors in the second dft run. On the other hand, if _d("loptim") is
# false, only the first dft run is enough.
#

    # S05: perform dft calculation (for the first time).
    message("ZEN", "DFT")
    #
    # S05.1: prepare and check essential files for the dft engine
    println("Initialize everything needed by the dft engine")
    dft_init(it)
    #
    # S05.2: perform a self-consitent calculation at the dft level
    println("Launch the dft engine")
    dft_run(it)
    #
    # S05.3: backup the output files of the dft engine
    println("Save the output files\n")
    dft_save(it)

    # we want better optimal projectors
    # in the previous dft run, initial fermi level = 0 -> wrong energy
    # window -> wrong optimial projectors. but at this point, the fermi
    # level is updated, so we have to generate the optimal projectors
    # again within this window
    if _d("loptim")

        # S06: perform dft calculation (for the second time).
        message("ZEN", "DFT")
        #
        # S06.1: prepare and check essential files for the dft engine
        println("Initialize everything needed by the dft engine")
        dft_init(it)
        #
        # S06.2: perform a self-consitent calculation at the dft level
        println("Launch the dft engine")
        dft_run(it)
        #
        # S06.3: backup the output files of the dft engine
        println("Save the output files\n")
        dft_save(it)

    end

#
# remarks:
#
# the key Kohn-Sham data inclue lattice structures, k-mesh and its weights,
# tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. at
# first the adaptor will read in these data from the output files of dft
# engine. and then it will process the raw projectors (parse, label, group,
# filter, and rotate). finally, the adaptor will write down the processed
# data to some specified files within the IR format.
#

    # S07: To bridge the gap between dft engine and dmft engine by adaptor
    message("ZEN", "ADAPTOR")
    #
    # S07.1: prepare and check essential files for the adaptor
    println("Initialize everything needed by the adaptor")
    adaptor_init(it)
    #
    # S07.2: launch the adaptor. it will read the Kohn-Sham data from the
    # dft engine, postprocess them, and then write them to external files
    # with the IR format
    println("Launch the adaptor")
    adaptor_run(it)
    #
    # S07.3: backup the output files of the adaptor
    println("Save the output files\n")
    adaptor_save(it)

    for iter = 1:_m("niter")
        message("ZEN", "ITER : $iter")
        message("ZEN", "DMFT1")
        message("ZEN", "SOLVER")
    end

else
    sorry()
end

goodbye()
