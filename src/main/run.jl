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
@printf("%s %s\n", "dmft -> solver    -> engine ->", dmft["solver"]["engine"])
@printf("%s\n", "------"^8)
@printf("%s %s\n", "dft_dmft -> mode    ->", dft_dmft["mode"])
@printf("%s %s\n", "dft_dmft -> axis    ->", dft_dmft["axis"])
@printf("%s %s\n", "dft_dmft -> beta    ->", dft_dmft["beta"])
@printf("%s %s\n", "dft_dmft -> niter   ->", dft_dmft["niter"])
@printf("%s %s\n", "dft_dmft -> mixer   ->", dft_dmft["mixer"])
@printf("%s %s\n", "dft_dmft -> cc      ->", dft_dmft["cc"])
@printf("%s %s\n", "dft_dmft -> ec      ->", dft_dmft["ec"])
@printf("%s %s\n", "dft_dmft -> fc      ->", dft_dmft["fc"])
@printf("%s %s\n", "dft_dmft -> lforce  ->", dft_dmft["lforce"])
@printf("%s %s\n", "dft_dmft -> lcharge ->", dft_dmft["lcharge"])
@printf("%s %s\n", "dft_dmft -> lenergy ->", dft_dmft["lenergy"])
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
