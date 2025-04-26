#
# Project : Tulip
# Source  : ACGui.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/05
#

"""
    ACGui

ACGui is a simple web app for the `ACFlow` toolkit. It depends on the Dash
web framework, and provide an useful ui to facilite analytic continuation
calculations. Now ACGui supports seven analytic continuation solvers:

* Maximum Entropy Method (`MaxEnt` solver, `recommended`)
* Barycentric Rational Function Method (`BarRat` solver, `recommended`)
*  Nevanlinna Analytical Continuation (`NevanAC` solver)
* Stochastic Analytic Continuation (`StochAC` solver, Beach's algorithm)
* Stochastic Analytic Continuation (`StochSK` solver, Sandvik's algorithm)
* Stochastic Optimization Method (`StochOM` solver)
* Stochastic Pole eXpansion (`StochPX` solver, `recommended`)

The MaxEnt, BarRat, and NevanAC solvers are extremely fast, so users can
obtain the calculated results quickly. However, other stochastic solvers
are quite slow (they could spend a few hours solving analytic continuation
problems). It is not a good idea to perform calculations with them through
ACGui. In such cases, users can download the relevant `ac.toml` files from
this app, and then submit their tasks manually.
"""
module ACGui

using TOML
using Base64
using Dash
using ACFlow

# Setup frontend ui
include("layout.jl")
#
export acg_layout!
export layout_header_block
export layout_data_block
export layout_base_block
export layout_maxent_block
export layout_barrat_block
export layout_nevanac_block
export layout_stochac_block
export layout_stochsk_block
export layout_stochom_block
export layout_stochpx_block
export layout_calc_block
export layout_about_block

# Register important callbacks to respond users' input
include("callback.jl")
#
export callbacks_in_data_tab
export callbacks_in_general_tab
export callbacks_in_solver_tab
export callbacks_in_run_tab
export callbacks_in_about_tab
export parse_parameters
export register_callback

# Launch the server side
include("base.jl")
#
export acg_clean
export acg_run

end
