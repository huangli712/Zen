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
chipsi1 = vaspio_projs(joinpath(pwd(), "dft"))
#chipsi2 = vaspio_projs(joinpath(pwd(), "dft"), false)
#irio_projs(pwd(), chipsi2)

# test vaspio_fermi() and irio_fermi()
#--------------------------------------
#fermi = vaspio_fermi(joinpath(pwd(), "dft"))
#irio_fermi(pwd(), fermi)

# test vaspio_lattice()
#--------------------------------------
#latt = vaspio_lattice(pwd())
#@show latt._case
#@show latt.scale
#@show latt.lvect
#@show latt.nsort
#@show latt.natom
#@show latt.sorts
#@show latt.atoms
#@show latt.coord

#chipsi = vaspio_projs(joinpath(pwd(), "dft"), false)
#kmesh, weight = vaspio_kmesh(joinpath(pwd(), "dft"))
#enk, occupy = vaspio_eigen(joinpath(pwd(), "dft"))
#ovlp = plo_ovlp(chipsi, weight)
#dm = plo_dm(chipsi, weight, occupy)
#view_ovlp(ovlp)
#view_dm(dm)

#orb_labels = ("s", 
#              "py", "pz", "px",
#              "dxy", "dyz", "dz2", "dxz", "dx2-y2",
#              "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)")
#
#for i in eachindex(orb_labels)
#    PrTrait(2, orb_labels[i])
#end
