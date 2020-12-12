#!/usr/bin/env julia

include("Zen.jl")
using .Zen

#kmesh, weight, ntet, volt, itet = vaspio_ibzkpt(pwd() * "/dft", true)
#for i in eachindex(weight)
#    println(i, " ", kmesh[i,:], " ", weight[i])
#end
#@show ntet, volt
#for t in 1:ntet
#    @show t, itet[t,:]
#end

#nsorts, natoms, symbols, atom_list, posi_list = vaspio_poscar(pwd() * "/dft")
#@show nsorts, natoms, symbols, atom_list, posi_list

#enk, occupy = vaspio_eigenval(pwd() * "/dft")
#@show enk[60,20,1]

#chipsi = vaspio_projcar(pwd() * "/dft")
#@show chipsi[:,19,1724,1]

#chipsi = vaspio_locproj(pwd() * "/dft")
#@show chipsi[:,19,1724,1]

kmesh, weight, ntet, volt, itet = vaspio_ibzkpt(pwd() * "/dft", true)
irio_tetra(ntet, volt)
