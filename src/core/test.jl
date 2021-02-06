#!/usr/bin/env julia

#
# Project : Pansy
# Source  : test.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/05
#

include("Zen.jl")
using .Zen

#setup()
#exhibit()
#it = IterInfo()
#adaptor_run(it)

z = 1.2
enk, occupy = vaspio_eigen("../03_SrVO3/dft")
fermi = vaspio_fermi("../03_SrVO3/dft")
volt, itet = vaspio_tetra("../03_SrVO3/dft")
@. enk = enk - fermi
wght = tetra_dos(z, itet, enk)
println(wght[:,280,1])
