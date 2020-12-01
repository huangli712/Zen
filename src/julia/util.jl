"""
    require()

Check the version of julia runtime environment
"""
function require()
    if VERSION < v"1.5-"
        error("Please upgrade your julia to v1.5.0 or higher")
    end
end

"""
    query_args()

Check whether the configuration file (case.toml) is provided
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
       error("Please specify the configuration file at least")
    else
       ARGS[1]
    end
end

"""
    query_cars()

Check whether the essential input files exist
"""
function query_cars()
    #if PDFT["engine"][1] == "vasp"
    if Param(PDFT, "engine") === "vasp"
        if !isfile("POSCAR") || !isfile("POTCAR")
            error("Please provide both POSCAR and POTCAR files")
        end
    else
        sorry()
    end
end

"""
    query_zen()

Query the home directory for zen
"""
function query_zen()
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_dft(dft::Dict{String,Any})

Query the home directory of the dft engine
"""
function query_dft(dft::Dict{String,Any})
    if dft["engine"] == "vasp"
        if haskey(ENV, "VASP_HOME")
            ENV["VASP_HOME"]
        else
            error("VASP_HOME is undefined")
        end
    else
        sorry()
    end
end

"""
    view_case()

Print the parameters to stdout: for case
"""
function view_case()
    println()
    message("case summary")
    println("case -> ", Param(PCASE, "case"))
    println()
end

"""
    view_dft()

Print the parameters to stdout: for dft
"""
function view_dft()
    message("dft parameters")
    println("dft  -> engine   -> ", Param(PDFT, "engine"))
    println("dft  -> smear    -> ", Param(PDFT, "smear"))
    println("dft  -> kmesh    -> ", Param(PDFT, "kmesh"))
    println("dft  -> magmom   -> ", Param(PDFT, "magmom"))
    println("dft  -> lsymm    -> ", Param(PDFT, "lsymm"))
    println("dft  -> lspins   -> ", Param(PDFT, "lspins"))
    println("dft  -> lspinorb -> ", Param(PDFT, "lspinorb"))
    println("dft  -> window   -> ", Param(PDFT, "window"))
    println("dft  -> loptim   -> ", Param(PDFT, "loptim"))
    println("dft  -> lproj    -> ", Param(PDFT, "lproj"))
    println("dft  -> nproj    -> ", Param(PDFT, "nproj"))
    println("dft  -> sproj    -> ", Param(PDFT, "sproj"))
    println()
end

"""
    view_dmft()

Print the parameters to stdout: for dmft
"""
function view_dmft()
    message("dmft parameters")
    println("dmft -> mode     -> ", Param(PDMFT, "mode"))
    println("dmft -> axis     -> ", Param(PDMFT, "axis"))
    println("dmft -> beta     -> ", Param(PDMFT, "beta"))
    println("dmft -> niter    -> ", Param(PDMFT, "niter"))
    println("dmft -> mixer    -> ", Param(PDMFT, "mixer"))
    println("dmft -> dcount   -> ", Param(PDMFT, "dcount"))
    println("dmft -> cc       -> ", Param(PDMFT, "cc"))
    println("dmft -> ec       -> ", Param(PDMFT, "ec"))
    println("dmft -> fc       -> ", Param(PDMFT, "fc"))
    println("dmft -> lcharge  -> ", Param(PDMFT, "lcharge"))
    println("dmft -> lenergy  -> ", Param(PDMFT, "lenergy"))
    println("dmft -> lforce   -> ", Param(PDMFT, "lforce"))
    println()
end

"""
    view_impurity()

Print the parameters to stdout: for impurity
"""
function view_impurity()
    message("impurity parameters")
    println("impurity -> nsite  -> ", PIMP["nsite"][1])
    println("impurity -> atoms  -> ", PIMP["atoms"][1])
    println("impurity -> equiv  -> ", PIMP["equiv"][1])
    println("impurity -> shell  -> ", PIMP["shell"][1])
    println("impurity -> ising  -> ", PIMP["ising"][1])
    println("impurity -> occup  -> ", PIMP["occup"][1])
    println("impurity -> upara  -> ", PIMP["upara"][1])
    println("impurity -> jpara  -> ", PIMP["jpara"][1])
    println("impurity -> lpara  -> ", PIMP["lpara"][1])
    println()
end

"""
    view_solver()

Print the parameters to stdout: for solver
"""
function view_solver()
    message("solver parameters")
    println("solver   -> engine -> ", PSOLVER["engine"][1])
    println("solver   -> params -> ", PSOLVER["params"][1])
    println()
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                   |\n", color = :green)
    printstyled("========== ========== ===     ===  | A Modern DFT + DMFT Computation Framework\n", color = :green)
    printstyled("        // ||         || n     ||  |\n", color = :green)
    printstyled("       //  ||         ||  n    ||  |\n", color = :green)
    printstyled("  //zz//   ||eeeeeeee ||   n   ||  |\n", color = :green)
    printstyled(" //        ||         ||    n  ||  | Version: 0.0.3@d\n", color = :green)
    printstyled("//         ||         ||     n ||  | Release: 2020/12\n", color = :green)
    printstyled("========== ========== ===     ===  | Powered by the julia programming language\n", color = :green)
    printstyled("                                   |\n", color = :green)
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen
"""
function goodbye() 
    println("See you later")
end

function sorry()
    error("Sorry, this feature has not been implemented")
end

function message(from::String, msg::String)
    printstyled("["*from*"]: ", color = :green)
    println(msg)
end

function message(from::String)
    message(from,"")
end
