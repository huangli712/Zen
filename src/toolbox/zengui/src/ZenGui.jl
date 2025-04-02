#
# Project : Camellia
# Source  : ZenGui.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/31
#

"""
    ZenGui

ZenGui is a general-purpose graphic user interface for ab initio dynamical
mean-field theory codes. It can be used to prepare necessary configuration
files for them. Now it supports the following codes:

* All-in-one DFT + DMFT package (**Zen**)
* Quantum impurity solver toolkit (**iQIST**)
* Analytic continuation tools (**ACFlow** and **ACTest**)
* Dynamical mean-field theory engines (**Dyson** and **DFermion**)

This code is under heavy development. **PLEASE USE IT AT YOUR OWN RISK**.

For more details about how to obtain, install and use the ZenGui code,
please visit the following website:

* `https://huangli712.github.io/projects/zengui/index.html`

Any suggestions, comments, and feedbacks are welcome. Enjoy it!
"""
module ZenGui

#=
### *Using Standard Libraries*
=#

using Dates
using Printf

#=
### *Using Third-Party Libraries*
=#

using CImGui
using CImGui.lib
using CImGui.CSyntax
using CImGui.CSyntax.CStatic

import GLFW
import ModernGL as GL

#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the ACFlow toolkit.

*Members* :

```text
I32, I64, API -> Numerical types (Integer).
F32, F64, APF -> Numerical types (Float).
C32, C64, APC -> Numerical types (Complex).
R32, R64, APR -> Numerical types (Union of Integer and Float).
N32, N64, APN -> Numerical types (Union of Integer, Float, and Complex).
#
__LIBNAME__   -> Name of this julia toolkit.
__VERSION__   -> Version of this julia toolkit.
__RELEASE__   -> Released date of this julia toolkit.
__AUTHORS__   -> Authors of this julia toolkit.
#
authors       -> Print the authors of ZenGui to screen.
```
=#

#
include("global.jl")
#
export I32, I64, API
export F32, F64, APF
export C32, C64, APC
export R32, R64, APR
export N32, N64, APN
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors

include("types.jl")

include("util.jl")

include("menu.jl")
export create_menu
export set_menu_file
export set_menu_edit
export set_menu_style
export set_menu_help

include("zen.jl")

include("dyson.jl")
include("dfermion.jl")

include("ctseg.jl")
include("cthyb.jl")
include("atomic.jl")

include("acflow.jl")
include("actest.jl")

include("about.jl")
export create_app_about

include("base.jl")
export zeng_run

end
