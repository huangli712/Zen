#!/usr/bin/env julia

include("Zen.jl")
using .Zen


@show pkgdir(Zen), typeof(Zen)
exit(-1)


# parse the configuration parameter
cfg = parse_toml(query_args(), true)
renew_config(cfg)
check_config()
#@timev list_case()
@timev list_dft()
#@timev list_dmft()
#@timev list_impurity()
#@timev list_solver()
exit(-1)


# read the Kohn-Sham data
PT, PG, chipsi = vaspio_projs(joinpath(pwd(), "dft"))
kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
fermi = vaspio_fermi(joinpath(pwd(), "dft"))
volt, itet = vaspio_tetra(joinpath(pwd(), "dft"))

# process the Kohn-Sham data
plo_group(PG)
PGT, chipsi_r = plo_rotate(PG, chipsi)
enk = enk .- fermi
bmin, bmax, ib_window, chipsi_w = plo_window(enk, 2.0, -1.4, chipsi_r)
plo_orthog(ib_window, PGT, chipsi_w)

# qualify the Kohn-Sham data
plo_dos(bmin, bmax, PGT, chipsi_w, itet, enk)
#ovlp = plo_ovlp(PGT, chipsi_w, weight)
#dm = plo_dm(bmin, bmax, PGT, chipsi_w, weight, occupy)
#hamk = plo_hamk(bmin, bmax, PGT, chipsi_w, weight, enk)
#view_ovlp(PGT, ovlp)
#view_dm(PGT, dm)
#view_hamk(PGT, hamk)
