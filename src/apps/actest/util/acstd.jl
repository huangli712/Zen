#!/usr/bin/env julia

#
# This script is used to create a standard dataset, which contains 100
# spectral functions and related correlation functions. This dataset
# is called ACT100. It should be used to benchmark the newly developed
# analytic continuation tools and methods. This script will launch only
# 1 process.
#
# Usage:
#
#     $ acstd.jl act.toml
#

using ACTest

welcome()
overview()
read_param()
make_data_std()
goodbye()
