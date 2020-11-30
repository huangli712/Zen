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
    query_cars(dft::Dict{String,Any})

Check whether the essential input files exist
"""
function query_cars(dft::Dict{String,Any})
    if dft["engine"] == "vasp"
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

Print the parameters to stdout
"""
function view_case()
    println()
    message("case summary")
    println("case -> ", PCASE["case"][1])
    println()
end

"""
    view_dft()

Print the parameters to stdout
"""
function view_dft()
    message("dft parameters")
    println("dft  -> engine   -> ", PDFT["engine"][1])
    println("dft  -> smear    -> ", PDFT["smear"][1])
    println("dft  -> kmesh    -> ", PDFT["kmesh"][1])
    println("dft  -> magmom   -> ", PDFT["magmom"][1])
    println("dft  -> lsymm    -> ", PDFT["lsymm"][1])
    println("dft  -> lspins   -> ", PDFT["lspins"][1])
    println("dft  -> lspinorb -> ", PDFT["lspinorb"][1])
    println("dft  -> window   -> ", PDFT["window"][1])
    println("dft  -> loptim   -> ", PDFT["loptim"][1])
    println("dft  -> lproj    -> ", PDFT["lproj"][1])
    println("dft  -> nproj    -> ", PDFT["nproj"][1])
    println("dft  -> sproj    -> ", PDFT["sproj"][1])
    println()
end

"""
    view_dmft()

Print the parameters to stdout
"""
function view_dmft()
    message("dmft parameters")
    println("dmft -> mode     -> ", PDMFT["mode"][1])
    println("dmft -> axis     -> ", PDMFT["axis"][1])
    println("dmft -> beta     -> ", PDMFT["beta"][1])
    println("dmft -> niter    -> ", PDMFT["niter"][1])
    println("dmft -> mixer    -> ", PDMFT["mixer"][1])
    println("dmft -> dcount   -> ", PDMFT["dcount"][1])
    println("dmft -> cc       -> ", PDMFT["cc"][1])
    println("dmft -> ec       -> ", PDMFT["ec"][1])
    println("dmft -> fc       -> ", PDMFT["fc"][1])
    println("dmft -> lcharge  -> ", PDMFT["lcharge"][1])
    println("dmft -> lenergy  -> ", PDMFT["lenergy"][1])
    println("dmft -> lforce   -> ", PDMFT["lforce"][1])
    println()
end

"""
    view_impurity()

Print the parameters to stdout
"""
function view_impurity()
    message("impurity parameters")
    println("impurity -> nsite  -> ", PIMP["nsite"][1])
    println("impurity -> atoms  -> ", PIMP["atoms"][1])
    println("impurity -> equiv  -> ", PIMP["equiv"][1])
    println("impurity -> shell  -> ", PIMP["shell"][1])
    println("impurity -> occup  -> ", PIMP["occup"][1])
    println("impurity -> upara  -> ", PIMP["upara"][1])
    println("impurity -> jpara  -> ", PIMP["jpara"][1])
    println("impurity -> lpara  -> ", PIMP["lpara"][1])
    println("dmft -> solver    -> engine -> ", dmft["solver"]["engine"])

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
