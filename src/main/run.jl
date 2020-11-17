#!/usr/bin/env julia

using Printf

include("Zen.jl")
using .Zen

# check the version of julia runtime environment
require()

# print the welcome message 
welcome()

# parse the file case.toml to extract parameters
message("zen", "parse the configuration")
cfg = parse_toml(query_args(), true)
case = cfg["case"]
dft = cfg["dft"]
dmft = cfg["dmft"]
dft_dmft = cfg["dft_dmft"]

# validate the parameters
message("zen", "check the configuration parameters")
# plug your codes here

# write the parameters to stdout
println()
println("case -> ", case)
println("------"^8)

println()
println("dft parameters:")
println("------"^8)
println("dft  -> engine    -> ", dft["engine"])
println("dft  -> lspins    -> ", dft["lspins"])
println("dft  -> lspinorb  -> ", dft["lspinorb"])
println("dft  -> adaptor   -> lortho -> ", dft["adaptor"]["lortho"])
println("dft  -> projector -> lproj  -> ", dft["projector"]["lproj"])
println("dft  -> projector -> nproj  -> ", dft["projector"]["nproj"])
println("dft  -> projector -> sproj  -> ", dft["projector"]["sproj"])

println()
println("dmft parameters:")
println("------"^8)
println("dmft -> dcount    -> ", dmft["dcount"])
println("dmft -> nominal   -> ", dmft["nominal"])
println("dmft -> impurity  -> nimp   -> ", dmft["impurity"]["nimp"])
println("dmft -> impurity  -> atoms  -> ", dmft["impurity"]["atoms"])
println("dmft -> impurity  -> equiv  -> ", dmft["impurity"]["equiv"])
println("dmft -> impurity  -> shell  -> ", dmft["impurity"]["shell"])
println("dmft -> impurity  -> upara  -> ", dmft["impurity"]["upara"])
println("dmft -> impurity  -> jpara  -> ", dmft["impurity"]["jpara"])
println("dmft -> impurity  -> lpara  -> ", dmft["impurity"]["lpara"])
println("dmft -> solver    -> engine -> ", dmft["solver"]["engine"])

println()
println("dft_dmft parameters:")
println("------"^8)
println("dft_dmft -> mode    -> ", dft_dmft["mode"])
println("dft_dmft -> axis    -> ", dft_dmft["axis"])
println("dft_dmft -> beta    -> ", dft_dmft["beta"])
println("dft_dmft -> niter   -> ", dft_dmft["niter"])
println("dft_dmft -> mixer   -> ", dft_dmft["mixer"])
println("dft_dmft -> cc      -> ", dft_dmft["cc"])
println("dft_dmft -> ec      -> ", dft_dmft["ec"])
println("dft_dmft -> fc      -> ", dft_dmft["fc"])
println("dft_dmft -> lforce  -> ", dft_dmft["lforce"])
println("dft_dmft -> lcharge -> ", dft_dmft["lcharge"])
println("dft_dmft -> lenergy -> ", dft_dmft["lenergy"])
exit(-1)


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
