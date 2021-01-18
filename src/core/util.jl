#
# project : pansy
# source  : util.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/19
#

"""
    @cswitch(constexpr, body)

Provides C-like switch statement with the ``falling through'' behavior. This
implement is borrowed from the following github repo.:
    https://github.com/Gnimuc/CSyntax.jl

# Examples
```julia
engine = get_d("engine")
@cswitch engine begin
    @case "vasp"
        vasp_init(it)
        break

    @default
        sorry()
        break
end
```
"""
macro cswitch(constexpr, body)
    case2label = Dict{Any,Symbol}()
    flow = Expr(:block)
    end_label = gensym("end")
    default_label = end_label

    for arg in body.args
        if Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@case")
            label = gensym("case")
            case2label[arg.args[3]] = label
            labelexpr = Expr(:symboliclabel, label)
            push!(flow.args, labelexpr)
        elseif Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@default")
            default_label = gensym("default")
            labelexpr = Expr(:symboliclabel, default_label)
            push!(flow.args, labelexpr)
        elseif arg == Expr(:break)
            labelexpr = Expr(:symbolicgoto, end_label)
            push!(flow.args, labelexpr)
        else
            push!(flow.args, arg)
        end
    end
    push!(flow.args, Expr(:symboliclabel, end_label))

    jumptable = Expr(:block)
    for (case, label) in case2label
        condition = Expr(:call, :(==), constexpr, case)
        push!(jumptable.args, Expr(:if, condition, Expr(:symbolicgoto, label)))
    end
    push!(jumptable.args[end].args, Expr(:symbolicgoto, default_label))

    return esc(Expr(:block, jumptable, flow))
end

"""
    require()

Check the version of julia runtime environment.
"""
function require()
    if VERSION < v"1.5-"
        error("Please upgrade your julia to v1.5.0 or higher")
    end
end

"""
    query_args()

Check whether the configuration file (case.toml) is provided.
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
        error("Please specify the configuration file")
    else
        ARGS[1]
    end
end

"""
    query_case()

Return case, in other words, the job's name.
"""
function query_case()
    basename( splitext(query_args())[1] )
end

"""
    query_inps()

Check whether the essential input files exist.
"""
function query_inps()
    # If the DFT engine is vasp, we have to ensure that the required input
    # files (POSCAR and POTCAR) are present.
    if get_d("engine") === "vasp"
        if !isfile("POSCAR") || !isfile("POTCAR")
            error("Please provide both POSCAR and POTCAR files")
        end
    else
        sorry()
    end
end

"""
    query_stop()

Query whether the case.stop file exists.
"""
function query_stop()
    isfile(query_case()*".stop")
end

"""
    query_zen()

Query the home directory of Zen.
"""
function query_zen()
    # we have to setup environment variable ZEN_HOME
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_dft()

Query the home directory of the dft engine
"""
function query_dft()
    # we have to setup environment variable VASP_HOME
    if get_d("engine") === "vasp"
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
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                        |\n", color = :green)
    printstyled("ZZZZZZZZZZZZ EEEEEEEEEEEE NNNNNNNNNNNN  | ", color = :green)
    printstyled("A Modern DFT + DMFT Computation Framework\n", color = :magenta)
    printstyled("          Z               N          N  |\n", color = :green)
    printstyled("         Z                N          N  |\n", color = :green)
    printstyled("   ZZZZZZ    EEEEEEEEEEEE N          N  |\n", color = :green)
    printstyled("  Z                       N          N  | ", color = :green)
    printstyled("Version: $__VERSION__\n", color = :magenta)
    printstyled(" Z                        N          N  | ", color = :green)
    printstyled("Release: $__RELEASE__\n", color = :magenta)
    printstyled("ZZZZZZZZZZZZ EEEEEEEEEEEE N          N  | ", color = :green)
    printstyled("Powered by the julia programming language\n", color = :magenta)
    printstyled("                                        |\n", color = :green)
    println()
end

"""
    overview()

Print out the overview of zen to the screen
"""
function overview()
    # build strings
    str1 = nprocs() === 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # write the information
    println("Starting time: ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Parallel execution: using ", nprocs(), str1, str2)
    println("Current working directory: ", pwd())
    println("Job description file: ", query_args())
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen
"""
function goodbye()
    println("See you later")
end

"""
    sorry()

Print an error message to the screen
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    message(from::String, msg::String)

Print an standard zen message to the screen
"""
function message(from::String, msg::String)
    printstyled(from * " > ", color = :green)
    printstyled(msg, color = :magenta)
    println()
end

"""
    line_to_array(io::IOStream)

Convert a line (reading from an iostream) to a string array
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end
