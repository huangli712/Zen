#!/usr/bin/env julia

include("Zen.jl")
using .Zen

kmesh, weight = vaspio_ibzkpt(pwd() * "/dft", true)
for i in eachindex(weight)
    println(i, " ", kmesh[i,:], " ", weight[i])
end

#nsorts, natoms, symbols, atom_list, posi_list = vaspio_poscar(pwd() * "/dft")
#@show nsorts, natoms, symbols, atom_list, posi_list

#vaspio_eigenval(pwd() * "/dft")

#chipsi = vaspio_projcar(pwd() * "/dft")
#@show chipsi[:,19,1724,1]

#chipsi = vaspio_locproj(pwd() * "/dft")
#@show chipsi[:,19,1724,1]
