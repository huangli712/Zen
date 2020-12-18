#!/usr/bin/env julia

include("Zen.jl")
using .Zen

#kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
#irio_kmesh(pwd(), kmesh, weight)

#volt, itet = vaspio_tetra(joinpath(pwd(), "dft"))
#irio_tetra(pwd(), volt, itet)

#enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
#irio_eigen(pwd(), enk, occupy)

#chipsi1 = vaspio_projs(joinpath(pwd(), "dft"))
#chipsi2 = vaspio_projs(joinpath(pwd(), "dft"), false)
#irio_projs(pwd(), chipsi2)

#fermi = vaspio_fermi(joinpath(pwd(), "dft"))
#irio_fermi(pwd(), fermi)

#nsorts, natoms, symbols, atom_list, posi_list = vaspio_poscar(pwd() * "/dft")
#@show nsorts, natoms, symbols, atom_list, posi_list
