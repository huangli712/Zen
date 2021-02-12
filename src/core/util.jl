#
# project : pansy
# source  : util.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/02/12
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
    @ps1(str, c)

Wrapper for printstyled function.
"""
macro ps1(str, c)
    return :( printstyled($str, color = $c) )
end

"""
    @ps2(str1, c1, str2, c2)

Wrapper for printstyled function.
"""
macro ps2(str1, c1, str2, c2)
    ex = quote
        printstyled($str1, color = $c1)
        printstyled($str2, color = $c2)
    end
    return :( $(esc(ex)) )
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
    setup_args(x::Vararg{String})

Setup `ARGS` manually. This function is used only in REPL environment.
We can use this function to update `ARGS`, so that the `query_args()`
and the other related functions can work correctly.

# Examples
```julia-repl
julia > setup_args("SrVO3.toml")
1-element Array{String,1}:
 "SrVO3.toml"
```

See also: [`query_args`](@ref).
"""
function setup_args(x::Vararg{String})
    # Clean `ARGS`
    empty!(ARGS)

    # Convert the arguments to an array of strings
    X = collect(x)

    # Push them into `ARGS` one by one
    for i in eachindex(X)
        push!(ARGS, X[i])
    end

    # Return `ARGS`, only for debug
    ARGS
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

See also: [`query_stop`](@ref).
"""
function query_case()
    basename( splitext(query_args())[1] )
end

"""
    query_inps()

Check whether the essential input files exist.
"""
function query_inps()

#
# Remarks:
#
# For vasp, the essential input files are:
#    1. POSCAR
#    2. POTCAR
# As for the INCAR and KPOINTS, they will be generated automatically.
#

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

See also: [`query_case`](@ref).
"""
function query_stop()
    isfile(query_case() * ".stop")
end

"""
    query_zen()

Query the home directory of Zen.

See also: [`query_dft`](@ref).
"""
function query_zen()
    # We have to setup the environment variable ZEN_HOME
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_dft()

Query the home directory of the DFT engine.

See also: [`query_zen`](@ref).
"""
function query_dft()
    # We have to setup the environment variable VASP_HOME
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

Print out the welcome messages to the screen.
"""
function welcome()
    @ps1 "                                        |\n" :green
    @ps2 "ZZZZZZZZZZZZ EEEEEEEEEEEE NNNNNNNNNNNN  | "  :green "A Modern DFT + DMFT Computation Framework\n" :magenta
    @ps1 "          Z               N          N  |\n" :green
    @ps1 "         Z                N          N  |\n" :green
    @ps1 "   ZZZZZZ    EEEEEEEEEEEE N          N  |\n" :green
    @ps2 "  Z                       N          N  | "  :green "Version: $__VERSION__\n" :magenta
    @ps2 " Z                        N          N  | "  :green "Release: $__RELEASE__\n" :magenta
    @ps2 "ZZZZZZZZZZZZ EEEEEEEEEEEE N          N  | "  :green "Powered by the julia programming language\n" :magenta
    @ps1 "                                        |\n" :green
    println()
end

"""
    overview()

Print out the overview of Zen to the screen.
"""
function overview()
    # Build strings
    str1 = nprocs() === 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # Write the information
    prompt("ZEN", "Overview")
    println("Time: ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Para: Using ", nprocs(), str1, str2)
    println("Dirs: ", pwd())
    println("Task: ", query_args())
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen.
"""
function goodbye()
    println("See you later")
end

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    prompt(from::String, msg::String)

Print a format Zen message to the screen.
"""
function prompt(from::String, msg::String)
    @ps2 "$from > " :green msg :magenta
    println()
end

"""
    prompt(msg::String)

Print a format Zen message to the screen.
"""
function prompt(msg::String)
    @ps2 "Task -> " :blue msg :light_red
    println()
end

"""
    prompt(io::IOStream, msg::String)

Print a format Zen message to the given IOStream.
"""
function prompt(io::IOStream, msg::String)
    date = Dates.format(now(), "yyyy-mm-dd / HH:MM:SS")
    println(io, "$date  $msg")
end

"""
    line_to_array(io::IOStream)

Convert a line (reading from an iostream) to a string array.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end

"""
    line_to_cmplx(io::IOStream)

Convert a line (reading from an iostream) to a cmplx number. It is used
to parse the `LOCPROJ` file.

See also: [`vaspio_projs`](@ref).
"""
@inline function line_to_cmplx(io::IOStream)
    str = readline(io)
    _re = str[10:28]
    _im = str[30:end]
    return parse(F64, _re) + parse(F64, _im) * im
end
