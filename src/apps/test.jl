#!/usr/bin/env julia

include("Zen.jl")
using .Zen

#kmesh, weight = vaspio_ibzkpt(pwd() * "/dft")
#
#for i in eachindex(weight)
#    println(i, " ", kmesh[i,:], " ", weight[i])
#end

#vaspio_poscar(pwd() * "/dft")
#vaspio_eigenval(pwd() * "/dft")
#vaspio_locproj(pwd() * "/dft", true)
vaspio_projcar(pwd() * "/dft")
