#
# Project : Camellia
# Source  : acflow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/26
#

"""
    create_app_acflow(p_open::Ref{Bool})

Create an UI window for the ACFlow toolkit, which provides some analytic
continuation tools.
"""
function create_app_acflow(p_open::Ref{Bool})
    # Create the ACFlow window, which can not be resized.
    CImGui.Begin(
        "ACFlow",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "ACFLOW"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _acflow_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the ac.toml file
    _acflow_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _acflow_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _acflow_top_block()

Setup widgets in the top of the window for the ACFlow toolkit.
"""
function _acflow_top_block()
    CImGui.Text("ACFlow: A Full-Fledged Analytic Continuation Toolkit In Julia")
end

"""
    _acflow_main_block()

Setup widgets associated with the parameters in the `ac.toml` file.
"""
function _acflow_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    if CImGui.BeginTabBar("ACFlowTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        if CImGui.BeginTabItem("general")
            _acflow_general_block()
            #
            CImGui.EndTabItem()
        end
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        if CImGui.BeginTabItem("solver")
            _acflow_solver_block()
            #
            CImGui.EndTabItem()
        end
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _acflow_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the ACFlow toolkit.
"""
function _acflow_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 80.0
    widget_button_height = 25.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View ac.toml")
    end
    #
    if CImGui.BeginPopupModal("View ac.toml", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_toml(build_acflow_dict())
            flags = CImGui.ImGuiInputTextFlags_ReadOnly
            flags = CImGui.ImGuiInputTextFlags_AllowTabInput | flags
            CImGui.InputTextMultiline("##source", text, 10000, ImVec2(400, 600), flags)
        end
        #
        if CImGui.Button("OK", ImVec2(widget_button_width, widget_button_height))
            CImGui.CloseCurrentPopup()
        end
        #
        CImGui.EndPopup()
    end
    #
    CImGui.SameLine()
    #
    if CImGui.Button("Close", ImVec2(widget_button_width, widget_button_height))
        p_open[] = false
    end
end

"""
    _acflow_general_block()

Setup widgets for the `general` tab.
"""
function _acflow_general_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    CImGui.Text("Configure [BASE] Block")

    # Input: finput
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic buf = "giw.data" * "\0"^60 begin
        CImGui.InputText(" Filename for input data", buf, length(buf))
        PBASE.finput = rstrip(buf,'\0')
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(finput)$(PBASE.finput)")
    #
    # Input: solver
    CImGui.SetNextItemWidth(widget_combo_width)
    solver_list = ["MaxEnt", "BarRat", "NevanAC", "StochAC", "StochSK", "StochOM", "StochPX"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Solver for the analytic continuation problem", &id, solver_list)
        PBASE.solver = solver_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(solver)$(PBASE.solver)")
    #
    # Input: ktype
    CImGui.SetNextItemWidth(widget_combo_width)
    ktype_list = ["fermi", "boson", "bsymm"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Type of kernel function", &id, ktype_list)
        PBASE.ktype = ktype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ktype)$(PBASE.ktype)")
    #
    # Input: mtype
    CImGui.SetNextItemWidth(widget_combo_width)
    mtype_list = ["flat", "gauss", "1gauss", "2gauss", "lorentz", "1lorentz", "2lorentz", "risedecay", "file"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Type of default model function", &id, mtype_list)
        PBASE.mtype = mtype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(mtype)$(PBASE.mtype)")
    #
    # Input: grid
    CImGui.SetNextItemWidth(widget_combo_width)
    grid_list = ["ftime", "fpart", "btime", "bpart", "ffreq", "ffrag", "bfreq", "bfrag"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Grid for input data (imaginary axis)", &id, grid_list)
        PBASE.grid = grid_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(grid)$(PBASE.grid)")
    #
    # Input: mesh
    CImGui.SetNextItemWidth(widget_combo_width)
    mesh_list = ["linear", "tangent", "lorentz", "halflorentz"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Mesh for output data (real axis)", &id, mesh_list)
        PBASE.mesh = mesh_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(mesh)$(PBASE.mesh)")
    #
    # Input: ngrid
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(10) begin
        @c CImGui.InputInt(" Number of grid points", &_i)
        PBASE.ngrid = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ngrid)$(PBASE.ngrid)")
    #
    # Input: nmesh
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(501) begin
        @c CImGui.InputInt(" Number of mesh points", &_i)
        PBASE.nmesh = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nmesh)$(PBASE.nmesh)")
    #
    # Input: wmax
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(5.0) begin
        @c CImGui.InputDouble(" Right boundary (maximum value) of output mesh", &_f)
        PBASE.wmax = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wmax)$(PBASE.wmax)")
    #
    # Input: wmin
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(-5.0) begin
        @c CImGui.InputDouble(" Left boundary (minimum value) of output mesh", &_f)
        PBASE.wmin = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wmin)$(PBASE.wmin)")
    #
    # Input: beta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(10.0) begin
        @c CImGui.InputDouble(" Inverse temperature", &_f)
        PBASE.beta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(beta)$(PBASE.beta)")
    #
    # Input: offdiag
    CImGui.SetNextItemWidth(widget_combo_width)
    offdiag_list = ["Yes", "No"]
    @cstatic id = Cint(1) begin
        @c CImGui.Combo(" Is it the offdiagonal correlation function", &id, offdiag_list)
        if id == 0
            PBASE.offdiag = true
        else
            PBASE.offdiag = false
        end
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(offdiag)$(PBASE.offdiag)")
    #
    # Input: fwrite
    CImGui.SetNextItemWidth(widget_combo_width)
    fwrite_list = ["Yes", "No"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Are the calculated results written into files", &id, fwrite_list)
        if id == 0
            PBASE.fwrite = true
        else
            PBASE.fwrite = false
        end
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(fwrite)$(PBASE.fwrite)")
end

"""
    _acflow_solver_block()

Setup widgets for the `solver` tab.
"""
function _acflow_solver_block()
    # It should change upon the selection of analytic continuation solver.
    CImGui.Text("Configure [$(PBASE.solver)] Block")

    @cswitch PBASE.solver begin

        @case "MaxEnt"
            _acflow_maxent_block()
            break

        @case "BarRat"
            _acflow_barrat_block()
            break

        @case "NevanAC"
            _acflow_nevanac_block()
            break

        @case "StochAC"
            _acflow_stochac_block()
            break

        @case "StochSK"
            _acflow_stochsk_block()
            break

        @case "StochOM"
            _acflow_stochom_block()
            break

        @case "StochPX"
            _acflow_stochpx_block()
            break

        @default
            sorry()
            break

    end
end

"""
    _acflow_maxent_block()

Setup widgets for the [MaxEnt] block in the ac.toml file.
"""
function _acflow_maxent_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: method
    CImGui.SetNextItemWidth(widget_combo_width)
    method_list = ["historic", "classic", "bryan", "chi2kink"]
    @cstatic id = Cint(3) begin
        @c CImGui.Combo(" How to determine the optimized α parameter", &id, method_list)
        PMaxEnt.method = method_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(method)$(PMaxEnt.method)")
    #
    # Input: stype
    CImGui.SetNextItemWidth(widget_combo_width)
    stype_list = ["sj", "br"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Type of the entropic term", &id, stype_list)
        PMaxEnt.stype = stype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(stype)$(PMaxEnt.stype)")
    #
    # Input: nalph
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(12) begin
        @c CImGui.InputInt(" Total number of the chosen α parameters", &_i)
        PMaxEnt.nalph = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nalph)$(PMaxEnt.nalph)")
    #
    # Input: alpha
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e9) begin
        @c CImGui.InputDouble(" Starting value for the α parameter", &_f)
        PMaxEnt.alpha = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(alpha)$(PMaxEnt.alpha)")
    #
    # Input: ratio
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(10.0) begin
        @c CImGui.InputDouble(" Scaling factor for the α parameter", &_f)
        PMaxEnt.ratio = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ratio)$(PMaxEnt.ratio)")
    #
    # Input: blur
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(-1.0) begin
        @c CImGui.InputDouble(" Shall we preblur the kernel and spectrum", &_f)
        PMaxEnt.blur = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(blur)$(PMaxEnt.blur)")
end

"""
    _acflow_barrat_block()

Setup widgets for the [BarRat] block in the ac.toml file.
"""
function _acflow_barrat_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: atype
    CImGui.SetNextItemWidth(widget_combo_width)
    atype_list = ["cont", "delta"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Possible type of the spectrum", &id, atype_list)
        PBarRat.atype = atype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(atype)$(PBarRat.atype)")
    #
    # Input: denoise
    CImGui.SetNextItemWidth(widget_combo_width)
    denoise_list = ["none", "prony_s", "prony_o"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" How to denoise the input data", &id, denoise_list)
        PBarRat.denoise = denoise_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(denoise)$(PBarRat.denoise)")
    #
    # Input: epsilon
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-10) begin
        @c CImGui.InputDouble(" Threshold for the Prony approximation", &_f)
        PBarRat.epsilon = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(epsilon)$(PBarRat.epsilon)")
    #
    # Input: pcut
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-3) begin
        @c CImGui.InputDouble(" Cutoff for unphysical poles", &_f)
        PBarRat.pcut = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(pcut)$(PBarRat.pcut)")
    #
    # Input: eta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-2) begin
        @c CImGui.InputDouble(" Tiny distance from the real axis", &_f)
        PBarRat.eta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(eta)$(PBarRat.eta)")
end

"""
    _acflow_nevanac_block()

Setup widgets for the [NevanAC] block in the ac.toml file.
"""
function _acflow_nevanac_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: pick
    CImGui.SetNextItemWidth(widget_combo_width)
    pick_list = ["Yes", "No"]
    @cstatic id = Cint(1) begin
        @c CImGui.Combo(" Check the Pick criterion or not", &id, pick_list)
        if id == 0
            PNevanAC.pick = true
        else
            PNevanAC.pick = false
        end
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(pick)$(PNevanAC.pick)")
    #
    # Input: pick
    CImGui.SetNextItemWidth(widget_combo_width)
    hardy_list = ["Yes", "No"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Perform Hardy basis optimization or not", &id, hardy_list)
        if id == 0
            PNevanAC.hardy = true
        else
            PNevanAC.hardy = false
        end
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(hardy)$(PNevanAC.hardy)")
    #
    # Input: hmax
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(50) begin
        @c CImGui.InputInt(" Upper cut off of Hardy order", &_i)
        PNevanAC.hmax = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(hmax)$(PNevanAC.hmax)")
    #
    # Input: alpha
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-4) begin
        @c CImGui.InputDouble(" Regulation parameter for smooth norm", &_f)
        PNevanAC.alpha = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(alpha)$(PNevanAC.alpha)")
    #
    # Input: eta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-2) begin
        @c CImGui.InputDouble(" Tiny distance from the real axis", &_f)
        PNevanAC.eta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(eta)$(PNevanAC.eta)")
end

"""
    _acflow_stochac_block()

Setup widgets for the [StochAC] block in the ac.toml file.
"""
function _acflow_stochac_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: nfine
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(10000) begin
        @c CImGui.InputInt(" Number of points of a very fine linear mesh", &_i)
        PStochAC.nfine = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nfine)$(PStochAC.nfine)")
    #
    # Input: ngamm
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(512) begin
        @c CImGui.InputInt(" Number of δ functions", &_i)
        PStochAC.ngamm = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ngamm)$(PStochAC.ngamm)")
    #
    # Input: nwarm
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(4000) begin
        @c CImGui.InputInt(" Number of Monte Carlo thermalization steps", &_i)
        PStochAC.nwarm = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nwarm)$(PStochAC.nwarm)")
    #
    # Input: nstep
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(4000000) begin
        @c CImGui.InputInt(" Number of Monte Carlo sweeping steps", &_i)
        PStochAC.nstep = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nstep)$(PStochAC.nstep)")
    #
    # Input: ndump
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(40000) begin
        @c CImGui.InputInt(" Intervals for monitoring Monte Carlo sweeps", &_i)
        PStochAC.ndump = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ndump)$(PStochAC.ndump)")
    #
    # Input: nalph
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(20) begin
        @c CImGui.InputInt(" Total number of the chosen α parameters", &_i)
        PStochAC.nalph = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nalph)$(PStochAC.nalph)")
    #
    # Input: alpha
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1.0) begin
        @c CImGui.InputDouble(" Starting value for the α parameter", &_f)
        PStochAC.alpha = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(alpha)$(PStochAC.alpha)")
    #
    # Input: ratio
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1.2) begin
        @c CImGui.InputDouble(" Scaling factor for the α parameter", &_f)
        PStochAC.ratio = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ratio)$(PStochAC.ratio)")
end

"""
    _acflow_stochsk_block()

Setup widgets for the [StochSK] block in the ac.toml file.
"""
function _acflow_stochsk_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: method
    CImGui.SetNextItemWidth(widget_combo_width)
    method_list = ["chi2min", "chi2kink"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" How to determine the optimized Θ parameter", &id, method_list)
        PStochSK.method = method_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(method)$(PStochSK.method)")
    #
    # Input: nfine
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(100000) begin
        @c CImGui.InputInt(" Number of points of a very fine linear mesh", &_i)
        PStochSK.nfine = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nfine)$(PStochSK.nfine)")
    #
    # Input: ngamm
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(1000) begin
        @c CImGui.InputInt(" Number of δ functions", &_i)
        PStochSK.ngamm = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ngamm)$(PStochSK.ngamm)")
    #
    # Input: nwarm
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(1000) begin
        @c CImGui.InputInt(" Number of Monte Carlo thermalization steps", &_i)
        PStochSK.nwarm = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nwarm)$(PStochSK.nwarm)")
    #
    # Input: nstep
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(20000) begin
        @c CImGui.InputInt(" Number of Monte Carlo sweeping steps", &_i)
        PStochSK.nstep = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nstep)$(PStochSK.nstep)")
    #
    # Input: ndump
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(200) begin
        @c CImGui.InputInt(" Intervals for monitoring Monte Carlo sweeps", &_i)
        PStochSK.ndump = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ndump)$(PStochSK.ndump)")
    #
    # Input: retry
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(10) begin
        @c CImGui.InputInt(" How often to recalculate the goodness function", &_i)
        PStochSK.retry = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(retry)$(PStochSK.retry)")
    #
    # Input: theta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e+6) begin
        @c CImGui.InputDouble(" Starting value for the Θ parameter", &_f)
        PStochSK.theta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(theta)$(PStochSK.theta)")
    #
    # Input: ratio
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(0.9) begin
        @c CImGui.InputDouble(" Scaling factor for the Θ parameter", &_f)
        PStochSK.ratio = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ratio)$(PStochSK.ratio)")
end

"""
    _acflow_stochom_block()

Setup widgets for the [StochOM] block in the ac.toml file.
"""
function _acflow_stochom_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: ntry
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(2000) begin
        @c CImGui.InputInt(" Number of attempts (tries) to seek the solution", &_i)
        PStochOM.ntry = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ntry)$(PStochOM.ntry)")
    #
    # Input: nstep
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(1000) begin
        @c CImGui.InputInt(" Number of Monte Carlo steps per attempt / try", &_i)
        PStochOM.nstep = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nstep)$(PStochOM.nstep)")
    #
    # Input: nbox
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(100) begin
        @c CImGui.InputInt(" Number of boxes to construct the spectrum", &_i)
        PStochOM.nbox = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nbox)$(PStochOM.nbox)")
    #
    # Input: sbox
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(0.005) begin
        @c CImGui.InputDouble(" Minimum area of the randomly generated boxes", &_f)
        PStochOM.sbox = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(sbox)$(PStochOM.sbox)")
    #
    # Input: wbox
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(0.02) begin
        @c CImGui.InputDouble(" Minimum width of the randomly generated boxes", &_f)
        PStochOM.wbox = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wbox)$(PStochOM.wbox)")
    #
    # Input: norm
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(-1.0) begin
        @c CImGui.InputDouble(" Is the norm calculated", &_f)
        PStochOM.norm = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(norm)$(PStochOM.norm)")
end

"""
    _acflow_stochpx_block()

Setup widgets for the [StochPX] block in the ac.toml file.
"""
function _acflow_stochpx_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: method
    CImGui.SetNextItemWidth(widget_combo_width)
    method_list = ["best", "mean"]
    @cstatic id = Cint(1) begin
        @c CImGui.Combo(" How to evaluate the final spectral density", &id, method_list)
        PStochPX.method = method_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(method)$(PStochPX.method)")
    #
    # Input: nfine
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(100000) begin
        @c CImGui.InputInt(" Number of grid points for a very fine mesh", &_i)
        PStochPX.nfine = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nfine)$(PStochPX.nfine)")
    #
    # Input: npole
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(200) begin
        @c CImGui.InputInt(" Number of poles", &_i)
        PStochPX.npole = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(npole)$(PStochPX.npole)")
    #
    # Input: ntry
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(1000) begin
        @c CImGui.InputInt(" Number of attempts (tries) to seek the solution", &_i)
        PStochPX.ntry = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ntry)$(PStochPX.ntry)")
    #
    # Input: nstep
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(1000000) begin
        @c CImGui.InputInt(" Number of Monte Carlo steps per attempt / try", &_i)
        PStochPX.nstep = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nstep)$(PStochPX.nstep)")
    #
    # Input: theta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e+6) begin
        @c CImGui.InputDouble(" Artificial inverse temperature", &_f)
        PStochPX.theta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(theta)$(PStochPX.theta)")
    #
    # Input: eta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1e-4) begin
        @c CImGui.InputDouble(" Tiny distance from the real axis", &_f)
        PStochPX.eta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(eta)$(PStochPX.eta)")
end
