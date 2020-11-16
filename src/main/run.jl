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
case = parse_toml(query_args(), "case", true) 
dft = parse_toml(query_args(), "dft", true)
dmft = parse_toml(query_args(), "dmft", true)
dft_dmft = parse_toml(query_args(), "dft_dmft", true)

# validate the parameters
message("zen", "check the configuration parameters")
# plug your codes here

@printf("%s\n", "------"^8)
@printf("%s %s\n", "case ->", case)
@printf("%s\n", "------"^8)
@printf("%s %s\n", "dft  -> engine    ->", dft["engine"])
@printf("%s %s\n", "dft  -> lspins    ->", dft["lspins"])
@printf("%s %s\n", "dft  -> lspinorb  ->", dft["lspinorb"])
@printf("%s %s\n", "dft  -> adaptor   -> lortho ->", dft["adaptor"]["lortho"])
@printf("%s %s\n", "dft  -> projector -> lproj  ->", dft["projector"]["lproj"])
@printf("%s %s\n", "dft  -> projector -> nproj  ->", dft["projector"]["nproj"])
@printf("%s %s\n", "dft  -> projector -> sproj  ->", dft["projector"]["sproj"])
@printf("%s\n", "------"^8)
@printf("%s %s\n", "dmft -> dcount    ->", dmft["dcount"])
@printf("%s %s\n", "dmft -> nominal   ->", dmft["nominal"])
@printf("%s %s\n", "dmft -> impurity  -> nimp   ->", dmft["impurity"]["nimp"])
@printf("%s %s\n", "dmft -> impurity  -> atoms  ->", dmft["impurity"]["atoms"])
@printf("%s %s\n", "dmft -> impurity  -> equiv  ->", dmft["impurity"]["equiv"])
@printf("%s %s\n", "dmft -> impurity  -> shell  ->", dmft["impurity"]["shell"])
@printf("%s %s\n", "dmft -> impurity  -> upara  ->", dmft["impurity"]["upara"])
@printf("%s %s\n", "dmft -> impurity  -> jpara  ->", dmft["impurity"]["jpara"])
@printf("%s %s\n", "dmft -> impurity  -> lpara  ->", dmft["impurity"]["lpara"])
@printf("%s\n", "------"^8)
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
