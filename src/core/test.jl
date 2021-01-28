#!/usr/bin/env julia

include("Zen.jl")
using .Zen

require()
welcome()
overview()
setup()
exhibit()
it = IterInfo()
adaptor_run(it)
exit(-1)

# read the Kohn-Sham data
#PT, PG, chipsi = vaspio_projs(joinpath(pwd(), "dft"))
#kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
#enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
#fermi = vaspio_fermi(joinpath(pwd(), "dft"))
#volt, itet = vaspio_tetra(joinpath(pwd(), "dft"))

# process the Kohn-Sham data
#plo_group(PG)
#PGT, chipsi_r = plo_rotate(PG, chipsi)
#enk = enk .- fermi
#bmin, bmax, ib_window, chipsi_w = plo_window(enk, 2.0, -1.4, chipsi_r)
#plo_orthog(ib_window, PGT, chipsi_w)

# qualify the Kohn-Sham data
#plo_dos(bmin, bmax, PGT, chipsi_w, itet, enk)
#ovlp = plo_ovlp(PGT, chipsi_w, weight)
#dm = plo_dm(bmin, bmax, PGT, chipsi_w, weight, occupy)
#hamk = plo_hamk(bmin, bmax, PGT, chipsi_w, weight, enk)
#view_ovlp(PGT, ovlp)
#view_dm(PGT, dm)
#view_hamk(PGT, hamk)
