#!/usr/bin/env julia

#
# Project : Begonia
# Source  : analyze.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/28
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_HOME"])

# Use the ZEN Framework
using Zen

print("Please specify the folder that contains the PROCAR file: ")
path = readline(stdin)
vaspio_procar(path)

exit()

    # Define orbital labels
    orb_labels = ["s",
                  "py", "pz", "px",
                  "dxy", "dyz", "dz2", "dxz", "dx2-y2"]

# The data are ready, then this function will interact with the users.
# Print essential information
println()
println("Number of spins: $nspin")
println("Number of k-points: $nkpt")
println("Number of bands: $nband")
println("Number of atoms: $natom")
println("Number of orbitals: $norbs")

# Enter a infinite loop until the users enter `q`
while true
    # Get atom index
    print("Please input atom index (integer, from 1 to $natom): ")
    atom_index = parse(I64, readline(stdin))

    # Get orbital index
    println("Orbitals: ", orb_labels)
    print("Please input orbital index (integer, from 1 to $norbs): ")
    orbital_index = parse(I64, readline(stdin))

    # Output the gathered information
    println("Selected atom index: $atom_index")
    println("Selected orbital index: $orbital_index")
    println("Selected orbital label: $(orb_labels[orbital_index])")

    # Sort, find out the most relevant orbitals
    v = sortperm(oab[orbital_index, atom_index, :], rev = true)

    # Output the band indices and weights
    print("Band index :")
    foreach(x -> @printf("%12i", x), v[1:5])
    println()
    print("Band weight:")
    foreach(x -> @printf("%12.7f", x), oab[orbital_index, atom_index, v[1:5]])
    println()

    # Prompt whether the users want to continue or quit
    println("If you want to continue, please enter `c` key, or else press `q` key")
    q = readline(stdin)

    # Quit the loop
    if q === "q"
        break
    end
end
