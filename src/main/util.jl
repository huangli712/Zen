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

function param_case(case::String)
    println()
    println("case summary:")
    println("------"^8)
    println("case -> ", case)
end

function param_dft()
end

function param_dmft()
end

function param_dft_dmft()
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                   |\n", color = :green)
    printstyled("========== ========== ===     ===  | A Modern DFT + DMFT Simulation Framework\n", color = :green)
    printstyled("        // ||         || n     ||  |\n", color = :green)
    printstyled("       //  ||         ||  n    ||  |\n", color = :green)
    printstyled("  //zz//   ||eeeeeeee ||   n   ||  |\n", color = :green)
    printstyled(" //        ||         ||    n  ||  | Version: 0.0.2@d\n", color = :green)
    printstyled("//         ||         ||     n ||  | Release: 2020/11\n", color = :green)
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
