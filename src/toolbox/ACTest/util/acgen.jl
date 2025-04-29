#!/usr/bin/env julia

#
# This script is used to build the spectral functions and the related
# correlation functions. It will launch only 1 process.
#
# Usage:
#
#     $ ./acgen.jl act.toml
#

haskey(ENV,"ACTEST_HOME") && pushfirst!(LOAD_PATH, ENV["ACTEST_HOME"])

using ACTest

welcome()
overview()
read_param()
make_data()
goodbye()
