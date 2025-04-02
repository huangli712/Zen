
# Project : Tulip
# Source  : callback.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/24
#

# The following global arrays are used to define the configure parameters
# for the ACFlow package.
#
# For the [BASE] block
const _PBASE = [
    "finput",
    "solver",
    "ktype",
    "mtype",
    "grid",
    "mesh",
    "ngrid",
    "nmesh",
    "wmax",
    "wmin",
    "beta",
    "offdiag",
    "fwrite"
]
#
# For the [MaxEnt] block
const _PMaxEnt = [
    "method",
    "stype",
    "nalph",
    "alpha",
    "ratio",
    "blur"
]
#
# For the [BarRat] block
const _PBarRat = [
    "atype",
    "denoise",
    "epsilon",
    "pcut",
    "eta"
]
#
# For the [StochAC] block
const _PStochAC = [
    "nfine",
    "ngamm",
    "nwarm",
    "nstep",
    "ndump",
    "nalph",
    "alpha",
    "ratio"
]
#
# For the [StochSK] block
const _PStochSK = [
    "method",
    "nfine",
    "ngamm",
    "nwarm",
    "nstep",
    "ndump",
    "retry",
    "theta",
    "ratio"
]
#
# For the [StochOM] block
const _PStochOM = [
    "ntry",
    "nstep",
    "nbox",
    "sbox",
    "wbox",
    "norm"
]
#
# For the [StochPX] block
const _PStochPX = [
    "method",
    "nfine",
    "npole",
    "ntry",
    "nstep",
    "theta",
    "eta"
]

"""
    callbacks_in_data_tab(app::Dash.DashApp)

Callbacks for the `data` tab. It only includes a callback, which is used
to upload files from client side to server side.
"""
function callbacks_in_data_tab(app::Dash.DashApp)
    callback!(
        app,
        Output("upload-file-name", "children"),
        Output("upload-file-type", "children"),
        Output("upload-file-nrow", "children"),
        Output("upload-file-ncol", "children"),
        Output("upload-file-head", "children"),
        Output("upload-file-tail", "children"),
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
    ) do contents, filename
        if !isnothing(filename)
            # At first, we should decode the contents, and convert them to
            # a readable string.
            content_type, content_string = split(contents, ',')
            decoded = base64decode(content_string)
            str = String(decoded)
            #
            # Write the upload data to file.
            open(filename, "w") do f
                write(f, str)
            end
            #
            # Get the first and last 4 rows of upload data
            str_vec = split(str, "\n", keepempty = false)
            str_head = [split(x, " ", keepempty = false) for x in str_vec[1:4]]
            str_tail = [split(x, " ", keepempty = false) for x in str_vec[end-3:end]]
            #
            # Get number of rows and columns of upload data
            nrow = length(str_vec)
            ncol = length(str_head[1])
            #
            # Build datatable for the first 4 rows of upload data
            dt_head = dash_datatable(
                columns = [
                    Dict(
                        "name" => "Column $i",
                        "id" => "column-$i"
                    ) for i in 1:ncol
                ],
                data = [
                    Dict(
                        "column-$i" => str_head[j][i] for i in 1:ncol
                    ) for j = 1:4
                ]
            )
            #
            # Build datatable for the last 4 rows of upload data
            dt_tail = dash_datatable(
                columns = [
                    Dict(
                        "name" => "Column $i",
                        "id" => "column-$i"
                    ) for i in 1:ncol
                ],
                data = [
                    Dict(
                        "column-$i" => str_tail[j][i] for i in 1:ncol
                    ) for j = 1:4
                ]
            )
            #
            return (filename, content_type, nrow, ncol, dt_head, dt_tail)
        else
            # Create an empty datatable
            dt = dash_datatable()
            return ("N/A", "N/A", "N/A", "N/A", dt, dt)
        end
    end
end

"""
    callbacks_in_general_tab(app::Dash.DashApp)

Callbacks for the `general` tab. It includes two callbacks. One is used to
control the `solver` tab. The other is used to gather parameters from this
tab, and then update `dict-base` in `run` tab.
"""
function callbacks_in_general_tab(app::Dash.DashApp)
    # Callback 1:
    #
    # Enable or disable specific solver panel in the `solver` tab.
    callback!(
        app,
        Output("maxent-block", "hidden"),
        Output("barrat-block", "hidden"),
        Output("stochac-block", "hidden"),
        Output("stochsk-block", "hidden"),
        Output("stochom-block", "hidden"),
        Output("stochpx-block", "hidden"),
        Input("base-solver", "value"),
    ) do solver
        # Enable `MaxEnt` solver
        if solver == "MaxEnt"
            return (false, true, true, true, true, true)
        end

        # Enable `BarRat` solver
        if solver == "BarRat"
            return (true, false, true, true, true, true)
        end

        # Enable `StochAC` solver
        if solver == "StochAC"
            return (true, true, false, true, true, true)
        end

        # Enable `StochSK` solver
        if solver == "StochSK"
            return (true, true, true, false, true, true)
        end

        # Enable `StochOM` solver
        if solver == "StochOM"
            return (true, true, true, true, false, true)
        end

        # Enable `StochPX` solver
        if solver == "StochPX"
            return (true, true, true, true, true, false)
        end
    end

    # Callback 2
    #
    # Collect parameters from all inputs in this tab. And then `dict-base`
    # in `run` tab is updated. Note that `dict-base` is hidden.
    callback!(
        app,
        Output("dict-base", "children"),
        [Input("base-$i", "value") for i in _PBASE],
    ) do vals...
        return join(vals, "|")
    end
end

"""
    callbacks_in_solver_tab(app::Dash.DashApp)

Callbacks for the `solver` tab. It includes six callbacks. All of them
are used to collect parameters that are relevant to analytic continuation
solvers.
"""
function callbacks_in_solver_tab(app::Dash.DashApp)
    # Callback 1
    #
    # Collect parameters from the `MaxEnt` panel. Then `dict-maxent` in
    # `run` tab will be updated. Note that `dict-maxent` is hidden.
    callback!(
        app,
        Output("dict-maxent", "children"),
        [Input("maxent-$i", "value") for i in _PMaxEnt],
    ) do vals...
        return join(vals, "|")
    end

    # Callback 2
    #
    # Collect parameters from the `BarRat` panel. Then `dict-barrat` in
    # `run` tab will be updated. Note that `dict-barrat` is hidden.
    callback!(
        app,
        Output("dict-barrat", "children"),
        [Input("barrat-$i", "value") for i in _PBarRat],
    ) do vals...
        return join(vals, "|")
    end

    # Callback 3
    #
    # Collect parameters from the `StochAC` panel. Then `dict-stochac` in
    # `run` tab will be updated. Note that `dict-stochac` is hidden.
    callback!(
        app,
        Output("dict-stochac", "children"),
        [Input("stochac-$i", "value") for i in _PStochAC],
    ) do vals...
        return join(vals, "|")
    end

    # Callback 4
    #
    # Collect parameters from the `StochSK` panel. Then `dict-stochsk` in
    # `run` tab will be updated. Note that `dict-stochsk` is hidden.
    callback!(
        app,
        Output("dict-stochsk", "children"),
        [Input("stochsk-$i", "value") for i in _PStochSK],
    ) do vals...
        return join(vals, "|")
    end

    # Callback 5
    #
    # Collect parameters from the `StochOM` panel. Then `dict-stochom` in
    # `run` tab will be updated. Note that `dict-stochom` is hidden.
    callback!(
        app,
        Output("dict-stochom", "children"),
        [Input("stochom-$i", "value") for i in _PStochOM],
    ) do vals...
        return join(vals, "|")
    end

    # Callback 6
    #
    # Collect parameters from the `StochPX` panel. Then `dict-stochpx` in
    # `run` tab will be updated. Note that `dict-stochpx` is hidden.
    callback!(
        app,
        Output("dict-stochpx", "children"),
        [Input("stochpx-$i", "value") for i in _PStochPX],
    ) do vals...
        return join(vals, "|")
    end
end

"""
    callbacks_in_run_tab(app::Dash.DashApp)

Callbacks for the `run` tab. It contains three callbacks. The first one
is for the `Start Analytic Continuation` button. The second one is for
the `Get ac.toml only` button. The third one is for the `Check err.out`
button.
"""
function callbacks_in_run_tab(app::Dash.DashApp)
    # Callback 1
    #
    # For the `Start Analytic Continuation` button. It will collect key
    # parameters, construct Dict structs, and launch the `ACFlow` package
    # to do the simulation. Finally, it will show the calculated results
    # in `canvas`.
    callback!(
        app,
        Output("canvas", "children"),
        Input("calc", "n_clicks"),
        State("dict-base", "children"),
        State("dict-maxent", "children"),
        State("dict-barrat", "children"),
        State("dict-stochac", "children"),
        State("dict-stochsk", "children"),
        State("dict-stochom", "children"),
        State("dict-stochpx", "children"),
    ) do btn, pbase, pmaxent, pbarrat, pstochac, pstochsk, pstochom, pstochpx
        if btn > 0
            # Convert parameters to dictionary
            B, S, solver = parse_parameters(
                pbase,
                pmaxent,
                pbarrat,
                pstochac,
                pstochsk,
                pstochom,
                pstochpx
            )

            # Print the resulting TOML file in terminal
            X = Dict("BASE"=>B, solver=>S)
            TOML.print(X)

            # Launch the `ACFlow` package to do analytic continuation.
            welcome()
            setup_param(B,S)
            mesh, Aout, _ = ACFlow.solve(ACFlow.read_data())

            # Visualize the calculated results
            fig = dcc_graph(
                figure = (
                    data = [(x = mesh, y = Aout),],
                )
            )

            return fig
        else
            return dcc_graph()
        end
    end

    # Callback 2
    #
    # Try to generate a TOML file from the collected parameters, and then
    # download it.
    callback!(
        app,
        Output("download-data", "children"),
        Input("get-ac-toml", "n_clicks"),
        State("dict-base", "children"),
        State("dict-maxent", "children"),
        State("dict-barrat", "children"),
        State("dict-stochac", "children"),
        State("dict-stochsk", "children"),
        State("dict-stochom", "children"),
        State("dict-stochpx", "children"),
    ) do btn, pbase, pmaxent, pbarrat, pstochac, pstochsk, pstochom, pstochpx
        if btn > 0
            # Convert parameters to dictionary
            B, S, solver = parse_parameters(
                pbase,
                pmaxent,
                pbarrat,
                pstochac,
                pstochsk,
                pstochom,
                pstochpx
            )

            # Print it to a TOML file
            X = Dict("BASE"=>B, solver=>S)
            io = IOBuffer()
            TOML.print(io,X)
            content = String(take!(io))

            return dcc_download(
                    data = Dict(
                        "content"=>content,
                        "filename" => "ac.toml"
                    )
                )
        else
            return nothing
        end
    end

    # Callback 3
    #
    # Show err.out in `err-out`.
    callback!(
        app,
        Output("err-out", "hidden"),
        Output("err-out", "children"),
        Input("check-err-out", "n_clicks"),
    ) do btn
        fn = "./err.out"
        if isfile(fn)
            err = read(fn, String)
        else
            err = "N/A"
        end
        #
        if iseven(btn)
            return(true, err)
        else
            return(false, err)
        end
    end
end

"""
    callbacks_in_about_tab(app::Dash.DashApp)

Callbacks for the `about` tab. Now it is empty.
"""
function callbacks_in_about_tab(app::Dash.DashApp)
end

"""
    parse_parameters(
        pbase::String,
        pmaxent::String,
        pbarrat::String,
        pstochac::String,
        pstochsk::String,
        pstochom::String,
        pstochpx::String
    )

Convert parameters to dictionary.
"""
function parse_parameters(
    pbase::String,
    pmaxent::String,
    pbarrat::String,
    pstochac::String,
    pstochsk::String,
    pstochom::String,
    pstochpx::String
)
    # For [BASE] block, it is necessary.
    array_base = split(pbase,"|")
    B = Dict{String,Any}(
        "finput" => string(array_base[1]),
        "solver" => string(array_base[2]),
        "ktype"  => string(array_base[3]),
        "mtype"  => string(array_base[4]),
        "grid"   => string(array_base[5]),
        "mesh"   => string(array_base[6]),
        "ngrid"  => parse(I64, array_base[7]),
        "nmesh"  => parse(I64, array_base[8]),
        "wmax"   => parse(F64, array_base[9]),
        "wmin"   => parse(F64, array_base[10]),
        "beta"   => parse(F64, array_base[11]),
        "offdiag" => parse(Bool, array_base[12]),
        "fwrite"  => parse(Bool, array_base[13]),
    )

    # For [MaxEnt] block, it is optional.
    if array_base[2] == "MaxEnt"
        array_maxent = split(pmaxent,"|")
        S = Dict{String,Any}(
            "method" => string(array_maxent[1]),
            "stype"  => string(array_maxent[2]),
            "nalph"  => parse(I64, array_maxent[3]),
            "alpha"  => parse(F64, array_maxent[4]),
            "ratio"  => parse(F64, array_maxent[5]),
            "blur"   => parse(F64, array_maxent[6]),
        )
    end

    # For [BarRat] block, it is optional.
    if array_base[2] == "BarRat"
        array_barrat = split(pbarrat,"|")
        S = Dict{String,Any}(
            "atype"   => string(array_barrat[1]),
            "denoise" => string(array_barrat[2]),
            "epsilon" => parse(F64, array_barrat[3]),
            "pcut"    => parse(F64, array_barrat[4]),
            "eta"     => parse(F64, array_barrat[5]),
        )
    end

    # For [StochAC] block, it is optional.
    if array_base[2] == "StochAC"
        array_stochac = split(pstochac,"|")
        S = Dict{String,Any}(
            "nfine"  => parse(I64, array_stochac[1]),
            "ngamm"  => parse(I64, array_stochac[2]),
            "nwarm"  => parse(I64, array_stochac[3]),
            "nstep"  => parse(I64, array_stochac[4]),
            "ndump"  => parse(I64, array_stochac[5]),
            "nalph"  => parse(I64, array_stochac[6]),
            "alpha"  => parse(F64, array_stochac[7]),
            "ratio"  => parse(F64, array_stochac[8]),
        )
    end

    # For [StochSK] block, it is optional.
    if array_base[2] == "StochSK"
        array_stochsk = split(pstochsk,"|")
        S = Dict{String,Any}(
            "method" => string(array_stochsk[1]),
            "nfine"  => parse(I64, array_stochsk[2]),
            "ngamm"  => parse(I64, array_stochsk[3]),
            "nwarm"  => parse(I64, array_stochsk[4]),
            "nstep"  => parse(I64, array_stochsk[5]),
            "ndump"  => parse(I64, array_stochsk[6]),
            "retry"  => parse(I64, array_stochsk[7]),
            "theta"  => parse(F64, array_stochsk[8]),
            "ratio"  => parse(F64, array_stochsk[9]),
        )
    end

    # For [StochOM] block, it is optional.
    if array_base[2] == "StochOM"
        array_stochom = split(pstochom,"|")
        S = Dict{String,Any}(
            "ntry"   => parse(I64, array_stochom[1]),
            "nstep"  => parse(I64, array_stochom[2]),
            "nbox"   => parse(I64, array_stochom[3]),
            "sbox"   => parse(F64, array_stochom[4]),
            "wbox"   => parse(F64, array_stochom[5]),
            "norm"   => parse(F64, array_stochom[6]),
        )
    end

    # For [StochPX] block, it is optional.
    if array_base[2] == "StochPX"
        array_stochpx = split(pstochpx,"|")
        S = Dict{String,Any}(
            "method" => string(array_stochpx[1]),
            "nfine"  => parse(I64, array_stochpx[2]),
            "npole"  => parse(I64, array_stochpx[3]),
            "ntry"   => parse(I64, array_stochpx[4]),
            "nstep"  => parse(I64, array_stochpx[5]),
            "theta"  => parse(F64, array_stochpx[6]),
            "eta"    => parse(F64, array_stochpx[7]),
        )
    end

    return B, S, array_base[2]
end

"""
    register_callback(app::Dash.DashApp)

Register all callbacks for the ACGui app.
"""
function register_callback(app::Dash.DashApp)
    callbacks_in_data_tab(app)
    callbacks_in_general_tab(app)
    callbacks_in_solver_tab(app)
    callbacks_in_run_tab(app)
    callbacks_in_about_tab(app)
end
