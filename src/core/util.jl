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

Check whether a case.toml file is provided
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
       error("Please provide a case.toml file to configure your calculation")
    else
       ARGS[1]
    end
end

"""
    query_zen()

Query the home directory for zen
"""
function query_zen()
    ENV["ZEN_HOME"]
end

"""
    check_dft(dft::Dict{String,Any})

Return the home directory of the dft engine
"""
function check_dft(dft::Dict{String,Any})
    if dft["engine"] == "vasp"
        ENV["VASP_HOME"]
    else
        sorry()
    end
end

"""
    check_inputs(dft::Dict{String,Any})

Check whether the essential input files exist
"""
function check_inputs(dft::Dict{String,Any})
    if dft["engine"] == "vasp"
        if !isfile("POSCAR") || !isfile("POTCAR")
            error("Please provide both POSCAR and POTCAR files")
        end
    else
        sorry()
    end
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
