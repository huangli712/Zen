#
# Project : Pansy
# Source  : solver.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/26
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
    ctqmc_setup()
    ctqmc_hyb_l()
    ctqmc_eimpx()
end

"""
    s_qmc1_exec(it::IterInfo)

Launch the CT-HYB1 quantum impurity solver.

This quantum impurity solver is from the `iQIST` software package.

See also: [`s_qmc1_init`](@ref), [`s_qmc1_save`](@ref).
"""
function s_qmc1_exec(it::IterInfo)
    # Print the header
    println("Engine : CT-HYB$(subscript(1))")

    # Get the home directory of quantum impurity solver
    solver_home = query_solver("ct_hyb1")

    # Determine mpi prefix (whether the solver is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "solver", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  Para : Using $numproc processors")

    # Select suitable solver program
    solver_exe = "$solver_home/ctqmc"
    @assert isfile(solver_exe)
    println("  Exec : $solver_exe")

    # Assemble command
    if isnothing(mpi_prefix)
        solver_cmd = solver_exe
    else
        solver_cmd = split("$mpi_prefix $solver_exe", " ")
    end

    # Launch it, the terminal output is redirected to solver.out
    run(pipeline(`$solver_cmd`, stdout = "solver.out"))

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
    println("Engine : CT-HYB$(subscript(2))")

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

function ctqmc_setup()
    open("solver.ctqmc.in", "w") do fout
        println(fout, "isscf = 1")
        println(fout, "isspn = 1")
        println(fout, "isort = 2")
        println(fout, "nband = 5")
        println(fout, "norbs = 10")
        println(fout, "ncfgs = 1024")

        println(fout, "Uc    = 4.0")
        println(fout, "Jz    = 0.7")
        println(fout, "mune  = 0.0")
        println(fout, "beta  = 40.0")
    end
end

function ctqmc_hyb_l()
    fmesh = []
    Delta = []
    cdim  = 0

    open("dmft.hyb_l", "r") do fin

        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, qdim, qdim, nmesh, nspin)

        for s = 1:nspin
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            cdim = parse(I64, line_to_array(strs)[7])
            @assert _t == 1 && _s == s
            for m = 1:nmesh
                fmesh[m] = parse(F64, line_to_array(fin)[3])
                # Parse hybridization functions
                for q = 1:cdim
                    for p = 1:cdim
                        _re, _im = parse.(F64, line_to_array(fin)[3:4])
                        Delta[p,q,m,s] = _re + _im * im
                    end
                end
            end
            # Skip two lines
            readline(fin)
            readline(fin)
        end
    end

    _, nband, nmesh, nspin = size(Delta)

    open("solver.hyb.in", "w") do fout
        for s = 1:nspin
            for p = 1:cdim
                orb = (s - 1) * cdim + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n", orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end
        end

        if nspin == 1
        for s = 1:nspin
            for p = 1:cdim
                orb = (2 - 1) * cdim + p
                for m = 1:nmesh
                    z = Delta[p,p,m,s]
                    @printf(fout, "%6i%16.8f%16.8f%16.8f%16.8f%16.8f\n", orb, fmesh[m], real(z), imag(z), 0.0, 0.0)
                end
                println(fout)
                println(fout)
            end
        end
        end
    end
end

function ctqmc_eimpx()
    Eimpx = []
    cdim  = 0

    open("dmft.eimpx", "r") do fin
    
        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        qdim = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for local impurity levels
        Eimpx = zeros(C64, qdim, qdim, nspin)

        for s = 1:nspin
            strs = readline(fin)
            _t = parse(I64, line_to_array(strs)[3])
            _s = parse(I64, line_to_array(strs)[5])
            cdim = parse(I64, line_to_array(strs)[7])
            @assert _t == 1 && _s == s

            # Parse local impurity levels
            for q = 1:cdim
                for p = 1:cdim
                    _re, _im = parse.(F64, line_to_array(fin)[3:4])
                    Eimpx[p,q,s] = _re + _im * im
                end
            end

            # Skip two lines
            readline(fin)
            readline(fin)
        end

    end

    _, nband, nspin = size(Eimpx)

    open("solver.eimp.in", "w") do fout

    end

end
