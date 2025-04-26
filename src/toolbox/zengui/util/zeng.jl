#!/usr/bin/env julia

#
# This script is used to launch the ZenGui app. To make it work correctly,
# please make sure that the ZenGui is correctly installed via `Pkg.add()`,
# or the environment variable `ZEN_GUI` is correctly set. It should point
# to the `zengui/src` directory.
#
# Usage:
#
#     $ ./zeng.jl
#

haskey(ENV,"ZEN_GUI") && pushfirst!(LOAD_PATH, ENV["ZEN_GUI"])

using ZenGui

zeng_run()
