#!/usr/bin/env julia

include("Zen.jl")
using .Zen

PT, PG, chipsi = vaspio_projs(joinpath(pwd(), "dft"))
kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
fermi = vaspio_fermi(joinpath(pwd(), "dft"))
cfg = parse_toml(query_args(), true)
renew_config(cfg)
plo_group(PG)
PGT, chipsi_r = plo_rotate(PG, chipsi)
enk = enk .- fermi
bmin, bmax, ib_window, chipsi_w = plo_window(enk, 2.0, -1.4, chipsi_r)
plo_orthog(ib_window, PGT, chipsi_w)
ovlp = plo_ovlp(PGT, chipsi_w, weight)
view_ovlp(PGT, ovlp)
dm = plo_dm(bmin, bmax, PGT, chipsi_w, weight, occupy)
view_dm(PGT, dm)
