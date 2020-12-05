#!/usr/bin/env julia

include("Zen.jl")
using .Zen

kmesh, weight = from_ibzkpt(pwd() * "/dft")

for i in eachindex(kmesh)
    #println(i, " ", kmesh[i,:], " ", weight[i])
    println(i)
end
