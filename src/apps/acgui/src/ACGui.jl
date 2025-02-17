#
# Project : Tulip
# Source  : ACGui.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/26
#

"""
    ACGui

ACGui is a simple web app for the `ACFlow` package. It depends on the Dash
web framework, and provide an useful ui to facilite analytic continuation
calculations. Now ACGui supports three analytic continuation solvers:

* Maximum Entropy Method (`MaxEnt` solver, `recommended`)
* Barycentric Rational Function Method (`BarRat` solver, `recommended`)
* Stochastic Pole eXpansion (`StochPX` solver, `recommended`)

The MaxEnt and BarRat solvers are extremely fast, so users can obtain the
calculated results quickly. However, the StochPX solver is quite slow. It
is not a good idea to start calculations with it through ACGui. In such
cases, users can download the relevant `ac.toml` files from the app, and
then submit tasks manually.
"""
module ACGui

using TOML
using Base64
using Dash
using ACFlow

include("layout.jl")
#
export acg_layout!
export layout_header_block
export layout_data_block
export layout_base_block
export layout_maxent_block
export layout_barrat_block
export layout_stochpx_block
export layout_calc_block
export layout_about_block

include("callback.jl")
#
export callbacks_in_data_tab
export callbacks_in_general_tab
export callbacks_in_solver_tab
export callbacks_in_run_tab
export callbacks_in_about_tab
export parse_parameters
export register_callback

include("base.jl")
#
export acg_clean
export acg_run

end
