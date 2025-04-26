#!/usr/bin/env julia

#
# This julia script will start a web app. Then the users can visit this
# app with any modern web browser, such as Chrome, Edge, and Safari. The
# default URL is
#
#    http://127.0.0.1:8848
#
# This app provides a simple yet userful graphic user interface for the
# `ACFlow` toolkit. The users can use it to do analytic continuations
# and visualize the calculated results on the fly. Notice that this app
# suits short analytic continuation simulations.
#
# Usage:
#
#     $ ./acg.jl
#

haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])
haskey(ENV,"ACGUI_HOME") && pushfirst!(LOAD_PATH, ENV["ACGUI_HOME"])

using ACGui

acg_run()
