#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/04/01
#

#
# CT-HYB1 Quantum Impurity Solver
#

"""
    s_qmc1_init(it::IterInfo)

Check runtime environment of the CT-HYB1 quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_exec`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_init(it::IterInfo)
end

"""
    s_qmc1_exec(it::IterInfo)

Launch the CT-HYB1 quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB1")

    # Print the footer for a better visualization
    println()
end

"""
    s_qmc1_save(it::IterInfo)

Backup output files of the CT-HYB1 quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_exec`](@ref).
"""
function s_qmc1_save(it::IterInfo)
end

#
# CT-HYB2 Quantum Impurity Solver
#

"""
    s_qmc2_init(it::IterInfo)

Check runtime environment of the CT-HYB2 quantum impurity solver. Prepare
the necessary input files.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_exec`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_init(it::IterInfo)
    sorry()
end

"""
    s_qmc2_exec(it::IterInfo)

Launch the CT-HYB2 quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_save`](@ref).
"""
function s_qmc2_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB2")

    # Print the footer for a better visualization
    println()
end

"""
    s_qmc2_save(it::IterInfo)

Backup output files of the CT-HYB2 quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc2_init`](@ref), [`s_qmc2_exec`](@ref).
"""
function s_qmc2_save(it::IterInfo)
    sorry()
end

#
# HUB-I Quantum Impurity Solver
#

"""
    s_hub1_init(it::IterInfo)

Check runtime environment of the HIA quantum impurity solver. Prepare
the necessary input files.

See also: [`s_hub1_exec`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_init(it::IterInfo)
    sorry()
end

"""
    s_hub1_exec(it::IterInfo)

Launch the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_save`](@ref).
"""
function s_hub1_exec(it::IterInfo)
    # Print the header
    println("Engine : HIA")

    # Print the footer for a better visualization
    println()
end

"""
    s_hub1_save(it::IterInfo)

Backup output files of the HIA quantum impurity solver.

See also: [`s_hub1_init`](@ref), [`s_hub1_exec`](@ref).
"""
function s_hub1_save(it::IterInfo)
    sorry()
end

#
# NORG Quantum Impurity Solver
#

"""
    s_norg_init(it::IterInfo)

Check runtime environment of the NORG quantum impurity solver. Prepare
the necessary input files.

See also: [`s_norg_exec`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_init(it::IterInfo)
    sorry()
end

"""
    s_norg_exec(it::IterInfo)

Launch the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_save`](@ref).
"""
function s_norg_exec(it::IterInfo)
    # Print the header
    println("Engine : NORG")

    # Print the footer for a better visualization
    println()
end

"""
    s_norg_save(it::IterInfo)

Backup output files of the NORG quantum impurity solver.

See also: [`s_norg_init`](@ref), [`s_norg_exec`](@ref).
"""
function s_norg_save(it::IterInfo)
    sorry()
end
