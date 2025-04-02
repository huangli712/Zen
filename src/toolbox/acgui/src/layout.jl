
# Project : Tulip
# Source  : layout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/26
#

"""
    acg_layout!(app::Dash.DashApp)

Global layout for the ACGui app.
"""
function acg_layout!(app::Dash.DashApp)
    app.layout = html_div() do
        layout_header_block(),
        html_br(),
        dcc_tabs(children = [
            dcc_tab(
                label = "Data",
                children = layout_data_block(),
                className = "custom-tab",
            ),
            dcc_tab(
                label = "General",
                children = layout_base_block(),
                className = "custom-tab",
            ),
            dcc_tab(
                label = "Solver",
                children = [
                    layout_maxent_block(),
                    layout_barrat_block(),
                    layout_stochpx_block(),
                ],
                className = "custom-tab",
            ),
            dcc_tab(
                label = "Run",
                children = layout_calc_block(),
                className = "custom-tab",
            ),
            dcc_tab(
                label = "About",
                children = layout_about_block(),
                className = "custom-tab",
            ),
        ]),
        html_br(),
        html_br(),
        html_br()
    end
end

"""
    layout_header_block()

Layout for the header part.
"""
function layout_header_block()
    html_center([
        html_h2("ACGui"),
        html_h3("A Graphic User Interface For ACFlow"),
    ])
end

"""
    layout_data_block()

Layout for the `data` tab. Users can upload imaginary-time or Matsubara
Green's functions via this tab.
"""
function layout_data_block()
    html_div([
        html_br(),
        # Table for basic information of the uploaded file.
        html_table([
            html_caption(
                html_b("Basic Information About the Uploaded File")
            ),
            html_thead(
                html_tr([
                    html_th("Filename"),
                    html_th("Type"),
                    html_th("Number of rows"),
                    html_th("Number of columns"),
                ])
            ),
            #
            html_tbody(
                html_tr([
                    html_td("N/A", id = "upload-file-name"),
                    html_td("N/A", id = "upload-file-type"),
                    html_td("N/A", id = "upload-file-nrow"),
                    html_td("N/A", id = "upload-file-ncol"),
                ])
            ),
        ]),
        html_br(),
        # The first 4 rows of the uploaded file
        html_div(id = "upload-file-head"),
        html_br(),
        # The last 4 rows of the uploaded file
        html_div(id = "upload-file-tail"),
        html_br(),
        # Upload control
        html_center(
            dcc_upload(
                children = html_div([
                    "Drag and Drop or ",
                    html_a("Select Files"),
                ]),
                id = "upload-data",
                multiple = false,
                className = "custom-upload",
            )
        ),
    ])
end

"""
    layout_base_block()

Layout for the `general` tab. Users should configure the basic parameters
for analytic continuation calculations in this tab.
"""
function layout_base_block()
    html_table([
        html_thead(
            html_tr(
                html_th(
                    html_label("Configuration parameters: general setup"),
                    colSpan = 3
                )
            )
        ),
        #
        html_tbody([
            html_tr([
                html_th(html_label("Filename for input data")),
                html_td(html_label("finput")),
                html_td(
                    dcc_input(
                        id = "base-finput",
                        type = "text",
                        value = "giw.data",
                        readOnly = true
                    )
                ),
            ]),
            # Now ACGui only supports three solvers.
            html_tr([
                html_th(html_label("Solver for the analytic continuation problem")),
                html_td(html_label("solver")),
                html_td(
                    dcc_dropdown(
                        id = "base-solver",
                        options = [
                            (label = "MaxEnt", value = "MaxEnt"),
                            (label = "BarRat", value = "BarRat"),
                            (label = "StochPX", value = "StochPX"),
                        ],
                        value = "MaxEnt",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Type of kernel function")),
                html_td(html_label("ktype")),
                html_td(
                    dcc_dropdown(
                        id = "base-ktype",
                        options = [
                            (label = "fermi", value = "fermi"),
                            (label = "boson", value = "boson"),
                            (label = "bsymm", value = "bsymm"),
                        ],
                        value = "fermi",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Type of default model function")),
                html_td(html_label("mtype")),
                html_td(
                    dcc_dropdown(
                        id = "base-mtype",
                        options = [
                            (label = "flat", value = "flat"),
                            (label = "gauss", value = "gauss"),
                            (label = "1gauss", value = "1gauss"),
                            (label = "2gauss", value = "2gauss"),
                            (label = "lorentz", value = "lorentz"),
                            (label = "1lorentz", value = "1lorentz"),
                            (label = "2lorentz", value = "2lorentz"),
                            (label = "risedecay", value = "risedecay"),
                        ],
                        value = "flat",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Grid for input data (imaginary axis)")),
                html_td(html_label("grid")),
                html_td(
                    dcc_dropdown(
                        id = "base-grid",
                        options = [
                            (label = "ftime", value = "ftime"),
                            (label = "fpart", value = "fpart"),
                            (label = "btime", value = "btime"),
                            (label = "bpart", value = "bpart"),
                            (label = "ffreq", value = "ffreq"),
                            (label = "ffrag", value = "ffrag"),
                            (label = "bfreq", value = "bfreq"),
                            (label = "bfrag", value = "bfrag"),
                        ],
                        value = "ftime",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Mesh for output data (real axis)")),
                html_td(html_label("mesh")),
                html_td(
                    dcc_dropdown(
                        id = "base-mesh",
                        options = [
                            (label = "linear", value = "linear"),
                            (label = "tangent", value = "tangent"),
                            (label = "lorentz", value = "lorentz"),
                            (label = "halflorentz", value = "halflorentz"),
                        ],
                        value = "linear",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of grid points")),
                html_td(html_label("ngrid")),
                html_td(
                    dcc_input(
                        id = "base-ngrid",
                        type = "text",
                        value = "10"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of mesh points")),
                html_td(html_label("nmesh")),
                html_td(
                    dcc_input(
                        id = "base-nmesh",
                        type = "text",
                        value = "501"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Right boundary (maximum value) of output mesh")),
                html_td(html_label("wmax")),
                html_td(
                    dcc_input(
                        id = "base-wmax",
                        type = "text",
                        value = "5.0"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Left boundary (minimum value) of output mesh")),
                html_td(html_label("wmin")),
                html_td(
                    dcc_input(
                        id = "base-wmin",
                        type = "text",
                        value = "-5.0"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Inverse temperature")),
                html_td(html_label("beta")),
                html_td(
                    dcc_input(
                        id = "base-beta",
                        type = "text",
                        value = "10.0"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Is it the offdiagonal part in matrix-valued function")),
                html_td(html_label("offdiag")),
                html_td(
                    dcc_radioitems(
                        id = "base-offdiag",
                        options = [
                            (label = "Yes", value = "true"),
                            (label = "No", value = "false"),
                        ],
                        value = "false",
                        labelStyle = Dict("display" => "inline-block")
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Are the analytic continuation results written into files")),
                html_td(html_label("fwrite")),
                html_td(
                    dcc_radioitems(
                        id = "base-fwrite",
                        options = [
                            (label = "Yes", value = "true"),
                            (label = "No", value = "false"),
                        ],
                        value = "false",
                        labelStyle = Dict("display" => "inline-block")
                    )
                ),
            ]),
        ]),
    ])
end

"""
    layout_maxent_block()

Layout for the `solver` tab. It is the panel for the `MaxEnt` solver. Note
that this panel can be hidden, if `solver` in `general` tab is not equal
to `MaxEnt`.
"""
function layout_maxent_block()
    html_table([
        html_thead(
            html_tr(
                html_th(html_label("[MaxEnt] block"), colSpan = 3)
            )
        ),
        #
        html_tbody([
            html_tr([
                html_th(html_label("How to determine the optimized α parameter")),
                html_td(html_label("method")),
                html_td(
                    dcc_dropdown(
                        id = "maxent-method",
                        options = [
                            (label = "historic", value = "historic"),
                            (label = "classic", value = "classic"),
                            (label = "bryan", value = "bryan"),
                            (label = "chi2kink", value = "chi2kink"),
                        ],
                        value = "chi2kink",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Type of the entropy term")),
                html_td(html_label("stype")),
                html_td(
                    dcc_dropdown(
                        id = "maxent-stype",
                        options = [
                            (label = "sj", value = "sj"),
                            (label = "br", value = "br"),
                        ],
                        value = "sj",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Total number of the chosen α parameters")),
                html_td(html_label("nalph")),
                html_td(
                    dcc_input(
                        id = "maxent-nalph",
                        type = "text",
                        value = "12"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Starting value for the α parameter")),
                html_td(html_label("alpha")),
                html_td(
                    dcc_input(
                        id = "maxent-alpha",
                        type = "text",
                        value = "1e9"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Scaling factor for the α parameter")),
                html_td(html_label("ratio")),
                html_td(
                    dcc_input(
                        id = "maxent-ratio",
                        type = "text",
                        value = "10.0"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Shall we preblur the kernel and spectrum")),
                html_td(html_label("blur")),
                html_td(
                    dcc_input(
                        id = "maxent-blur",
                        type = "text",
                        value = "-1.0"
                    )
                ),
            ]),
        ]),
    ], id = "maxent-block", hidden = true)
end

"""
    layout_barrat_block()

Layout for the `solver` tab. It is the panel for the `BarRat` solver. Note
that this panel can be hidden, if `solver` in `general` tab is not equal
to `BarRat`.
"""
function layout_barrat_block()
    html_table([
        html_thead(
            html_tr(
                html_th(html_label("[BarRat] block"), colSpan = 3)
            )
        ),
        #
        html_tbody([
            html_tr([
                html_th(html_label("Possible type of the spectrum")),
                html_td(html_label("atype")),
                html_td(
                    dcc_dropdown(
                        id = "barrat-atype",
                        options = [
                            (label = "cont", value = "cont"),
                            (label = "delta", value = "delta"),
                        ],
                        value = "cont",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("How to denoise the input data")),
                html_td(html_label("denoise")),
                html_td(
                    dcc_dropdown(
                        id = "barrat-denoise",
                        options = [
                            (label = "none", value = "none"),
                            (label = "prony_s", value = "prony_s"),
                            (label = "prony_o", value = "prony_o"),
                        ],
                        value = "none",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Threshold for the Prony approximation")),
                html_td(html_label("epsilon")),
                html_td(
                    dcc_input(
                        id = "barrat-epsilon",
                        type = "text",
                        value = "1e-10"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Cutoff for unphysical poles")),
                html_td(html_label("pcut")),
                html_td(
                    dcc_input(
                        id = "barrat-pcut",
                        type = "text",
                        value = "1e-3"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Tiny distance from the real axis")),
                html_td(html_label("eta")),
                html_td(
                    dcc_input(
                        id = "barrat-eta",
                        type = "text",
                        value = "1e-2"
                    )
                ),
            ]),
        ]),
    ], id = "barrat-block", hidden = true)
end

"""
    layout_stochpx_block()

Layout for the `solver` tab. It is the panel for the `StochPX` solver. Note
that this panel can be hidden, if `solver` in `general` tab is not equal
to `StochPX`.
"""
function layout_stochpx_block()
    html_table([
        html_thead(
            html_tr(
                html_th(html_label("[StochPX] block"), colSpan = 3)
            )
        ),
        #
        html_tbody([
            html_tr([
                html_th(html_label("How to evaluate the final spectral density")),
                html_td(html_label("method")),
                html_td(
                    dcc_dropdown(
                        id = "stochpx-method",
                        options = [
                            (label = "best", value = "best"),
                            (label = "mean", value = "mean"),
                        ],
                        value = "mean",
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of points of a very fine linear mesh")),
                html_td(html_label("nfine")),
                html_td(
                    dcc_input(
                        id = "stochpx-nfine",
                        type = "text",
                        value = "100000"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of poles")),
                html_td(html_label("npole")),
                html_td(
                    dcc_input(
                        id = "stochpx-npole",
                        type = "text",
                        value = "200"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of attempts (tries) to seek the solution")),
                html_td(html_label("ntry")),
                html_td(
                    dcc_input(
                        id = "stochpx-ntry",
                        type = "text",
                        value = "1000"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Number of Monte Carlo steps per attempt / try")),
                html_td(html_label("nstep")),
                html_td(
                    dcc_input(
                        id = "stochpx-nstep",
                        type = "text",
                        value = "1000000"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Artificial inverse temperature")),
                html_td(html_label("theta")),
                html_td(
                    dcc_input(
                        id = "stochpx-theta",
                        type = "text",
                        value = "1e+6"
                    )
                ),
            ]),
            html_tr([
                html_th(html_label("Tiny distance from the real axis")),
                html_td(html_label("eta")),
                html_td(
                    dcc_input(
                        id = "stochpx-eta",
                        type = "text",
                        value = "1e-4"
                    )
                ),
            ]),
        ]),
    ], id = "stochpx-block", hidden = true)
end

"""
    layout_calc_block()

Layout for the `run` tab. The users can start analytic continuation
simulations and visualize the calculated results via this tab.
"""
function layout_calc_block()
    html_div([
        html_br(),
        # The following four labels are always invisible. They are used to
        # collect the configuration parameters for ACFlow.
        html_label(
            children = "N/A",
            id = "dict-base",
            hidden = true,
        ),
        html_label(
            children = "N/A",
            id = "dict-maxent",
            hidden = true,
        ),
        html_label(
            children = "N/A",
            id = "dict-barrat",
            hidden = true,
        ),
        html_label(
            children = "N/A",
            id = "dict-stochpx",
            hidden = true,
        ),
        html_br(),
        html_center(
            html_button(
                children = "Start Analytic Continuation",
                id = "calc",
                n_clicks = 0
            )
        ),
        html_br(),
        html_center(
            html_button(
                children = "Get ac.toml only",
                id = "get-ac-toml",
                n_clicks = 0
            )
        ),
        html_div(id = "download-data"),
        html_br(),
        html_center(
            html_button(
                children = "Check err.out",
                id = "check-err-out",
                n_clicks = 0
            )
        ),
        html_br(),
        # To display err.out
        html_div(id = "err-out", hidden = true),
        html_br(),
        # To display the calculated results
        html_div(id = "canvas"),
    ])
end

"""
    layout_about_block()

Layout for the `about` tab. It is used to display the version number and
author of ACGui.
"""
function layout_about_block()
    html_div([
        html_h4("Version : v0.5.0-devel.241030"),
        html_h4("Release : 2024/10"),
        html_h4("Developed by Li Huang (hungli@caep.cn)"),
    ])
end
