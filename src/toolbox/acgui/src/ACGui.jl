#
# Project : Tulip
# Source  : ACGui.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/24
#

"""
    ACGui

ACGui is a simple web app for the `ACFlow` package. It depends on the Dash
web framework, and provide an useful ui to facilite analytic continuation
calculations. Now ACGui supports six analytic continuation solvers:

* Maximum Entropy Method (`MaxEnt` solver, `recommended`)
* Barycentric Rational Function Method (`BarRat` solver, `recommended`)
* Stochastic Analytic Continuation (`StochAC` solver, Beach's algorithm)
* Stochastic Analytic Continuation (`StochSK` solver, Sandvik's algorithm)
* Stochastic Optimization Method (`StochOM` solver)
* Stochastic Pole eXpansion (`StochPX` solver, `recommended`)

The MaxEnt and BarRat solvers are extremely fast, so users can obtain the
calculated results quickly. However, the other stochastic solvers are
quite slow (they could spend several hours solving analytic continuation
problems). It is not a good idea to start calculations with them through
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
