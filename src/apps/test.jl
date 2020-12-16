#!/usr/bin/env julia

include("Zen.jl")
using .Zen

#kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
#irio_kmesh(pwd(), kmesh, weight)

#volt, itet = vaspio_tetra(joinpath(pwd(), "dft"))
#irio_tetra(pwd(), volt, itet)

enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
irio_eigen(pwd(), enk, occupy)

#nsorts, natoms, symbols, atom_list, posi_list = vaspio_poscar(pwd() * "/dft")
#@show nsorts, natoms, symbols, atom_list, posi_list

#chipsi = vaspio_projcar(pwd() * "/dft")
#@show chipsi[:,19,1724,1]

#chipsi = vaspio_locproj(pwd() * "/dft")
#@show chipsi[:,19,1724,1]

#kmesh, weight, ntet, volt, itet = vaspio_ibzkpt(pwd() * "/dft", true)
#irio_tetra(pwd(), ntet, volt, itet)

#enk, occupy = vaspio_eigenval(pwd() * "/dft")
#irio_eigen(pwd(), enk, occupy)
