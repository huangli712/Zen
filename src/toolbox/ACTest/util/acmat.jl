#!/usr/bin/env julia

#
# This script is used to build the spectral functions with negative
# weights, and the corresponding matrix-value Green's functions. This
# script will launch only 1 process.
#
# Usage:
#
#     $ ./acmat.jl act.toml
#

haskey(ENV,"ACTEST_HOME") && pushfirst!(LOAD_PATH, ENV["ACTEST_HOME"])

using ACTest

welcome()
overview()
read_param()
make_data_mat()
goodbye()
