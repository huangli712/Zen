#!/usr/bin/env julia

include("Zen.jl")
using .Zen

# test vaspio_kmesh() and irio_kmesh()
#--------------------------------------
#kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
#irio_kmesh(pwd(), kmesh, weight)

# test vaspio_tetra() and irio_tetra()
#--------------------------------------
#volt, itet = vaspio_tetra(joinpath(pwd(), "dft"))
#irio_tetra(pwd(), volt, itet)

# test vaspio_eigen() and irio_eigen()
#--------------------------------------
#enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
#irio_eigen(pwd(), enk, occupy)

# test vaspio_projs() and irio_projs()
#--------------------------------------
#chipsi1 = vaspio_projs(joinpath(pwd(), "dft"))
#chipsi2 = vaspio_projs(joinpath(pwd(), "dft"), false)
#irio_projs(pwd(), chipsi2)

# test vaspio_fermi() and irio_fermi()
#--------------------------------------
#fermi = vaspio_fermi(joinpath(pwd(), "dft"))
#irio_fermi(pwd(), fermi)

# test vaspio_lattice()
#--------------------------------------
#nsorts, natoms, symbols, atom_list, posi_list = vaspio_lattice(joinpath(pwd(), "dft"))
#@show nsorts, natoms, symbols, atom_list, posi_list

chipsi = vaspio_projs(joinpath(pwd(), "dft"), false)
kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
ovlp = plo_ovlp(chipsi, weight)
dm = plo_dm(chipsi, weight, occupy)
view_ovlp(ovlp)
view_dm(dm)
