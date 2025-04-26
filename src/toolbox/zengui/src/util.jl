#
# Project : Camellia
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/25
#

#=
### *Basic Macros*
=#

"""
    @cswitch(constexpr, body)

Provides a C-like switch statement with the *falling through* behavior.
This implementation was borrowed from the following github repository:

* https://github.com/Gnimuc/CSyntax.jl

### Examples
```julia
engine = get_d("engine")
@cswitch engine begin
    @case "vasp"
        just_do_it()
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

#=
### *Error Handler*
=#

"""
    sorry()

Print an error message to the screen.

### Arguments
N/A

### Returns
N/A
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

#=
### *Dictionary Utility*
=#

"""
    dict_to_toml(d::AbstractDict)

Convert an ordered dictionary to toml file (actually String).
"""
function dict_to_toml(d::AbstractDict)
    io = IOBuffer()
    TOML.print(io,d)
    return String(take!(io))
end

"""
    dict_to_ini(d::AbstractDict)

Convert an ordered dictionary to ini file (actually String).
"""
function dict_to_ini(d::AbstractDict)
    io = IOBuffer()
    for (key, value) in d
        println(io, "$key = $value")
    end
    return String(take!(io))
end

#=
### *Miscellaneous Utility*
=#

"""
    open_url(url::String)

Invoke the default web browser to open the given url. It only supports the
windows, macos, and linux systems.
"""
function open_url(url::String)
    if Sys.iswindows()
        run(`start $url`)
    elseif Sys.islinux()
        run(`xdg-open $url`)
    elseif Sys.isapple()
        run(`open $url`)
    else
        sorry()
    end
end
