#
# Project : Camellia
# Source  : actest.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/26
#

"""
    create_app_actest(p_open::Ref{Bool})

Create an UI window for the ACTest toolkit, which is used to benchmark the
analytic continuation tools as implemented in the ACFlow package.
"""
function create_app_actest(p_open::Ref{Bool})
    # Create the ACTest window, which can not be resized.
    CImGui.Begin(
        "ACTest",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "ACTEST"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _actest_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the act.toml file
    _actest_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _actest_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _actest_top_block()

Setup widgets in the top of the window for the ACTest toolkit.
"""
function _actest_top_block()
    CImGui.Text("ACTest: A Testing Toolkit For Analytic Continuation Methods And Codes")
end

"""
    _actest_main_block()

Setup widgets associated with the parameters in the `act.toml` file.
"""
function _actest_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    if CImGui.BeginTabBar("ACTestTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        if CImGui.BeginTabItem("general")
            _actest_general_block()
            #
            CImGui.EndTabItem()
        end
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        if CImGui.BeginTabItem("solver")
            _actest_solver_block()
            #
            CImGui.EndTabItem()
        end
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _actest_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the ACTest toolkit.
"""
function _actest_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 80.0
    widget_button_height = 25.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View act.toml")
    end
    #
    if CImGui.BeginPopupModal("View act.toml", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_toml(build_actest_dict())
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
    _actest_general_block()

Setup widgets for the `general` tab.
"""
function _actest_general_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    CImGui.Text("Configure [Test] Block")

    # Input: solver
    CImGui.SetNextItemWidth(widget_combo_width)
    solver_list = ["MaxEnt", "BarRat", "NevanAC", "StochAC", "StochSK", "StochOM", "StochPX"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Analytic continuation solver", &id, solver_list)
        PTEST.solver = solver_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(solver)$(PTEST.solver)")
    #
    # Input: ptype
    CImGui.SetNextItemWidth(widget_combo_width)
    ptype_list = ["gauss", "lorentz", "delta", "rectangle", "risedecay"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Type of peaks in the spectrum", &id, ptype_list)
        PTEST.ptype = ptype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ptype)$(PTEST.ptype)")
    #
    # Input: ktype
    CImGui.SetNextItemWidth(widget_combo_width)
    ktype_list = ["fermi", "boson", "bsymm"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Type of kernel function", &id, ktype_list)
        PTEST.ktype = ktype_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ktype)$(PTEST.ktype)")
    #
    # Input: grid
    CImGui.SetNextItemWidth(widget_combo_width)
    grid_list = ["ftime", "btime", "ffreq", "bfreq"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Grid for correlation function", &id, grid_list)
        PTEST.grid = grid_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(grid)$(PTEST.grid)")
    #
    # Input: mesh
    CImGui.SetNextItemWidth(widget_combo_width)
    mesh_list = ["linear", "tangent", "lorentz", "halflorentz"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Mesh for spectral function", &id, mesh_list)
        PTEST.mesh = mesh_list[id + 1]
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(mesh)$(PTEST.mesh)")
    #
    # Input: ngrid
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(10) begin
        @c CImGui.InputInt(" Number of grid points", &_i)
        PTEST.ngrid = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ngrid)$(PTEST.ngrid)")
    #
    # Input: nmesh
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(501) begin
        @c CImGui.InputInt(" Number of mesh points", &_i)
        PTEST.nmesh = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(nmesh)$(PTEST.nmesh)")
    #
    # Input: ntest
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _i = Cint(100) begin
        @c CImGui.InputInt(" Number of tests", &_i)
        PTEST.ntest = _i
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ntest)$(PTEST.ntest)")
    #
    # Input: wmax
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(5.0) begin
        @c CImGui.InputDouble(" Right boundary (maximum value) of real mesh", &_f)
        PTEST.wmax = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wmax)$(PTEST.wmax)")
    #
    # Input: wmin
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(-5.0) begin
        @c CImGui.InputDouble(" Left boundary (minimum value) of real mesh", &_f)
        PTEST.wmin = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wmin)$(PTEST.wmin)")
    #
    # Input: pmax
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(4.0) begin
        @c CImGui.InputDouble(" Right boundary (maximum value) for possible peaks", &_f)
        PTEST.pmax = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(pmax)$(PTEST.pmax)")
    #
    # Input: pmin
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(-4.0) begin
        @c CImGui.InputDouble(" Left boundary (minimum value) for possible peaks", &_f)
        PTEST.pmin = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(wmin)$(PTEST.pmin)")
    #
    # Input: beta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(10.0) begin
        @c CImGui.InputDouble(" Inverse temperature", &_f)
        PTEST.beta = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(beta)$(PTEST.beta)")
    #
    # Input: noise
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(1.0e-6) begin
        @c CImGui.InputDouble(" Noise level", &_f)
        PTEST.noise = _f
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(noise)$(PTEST.noise)")
    #
    # Input: offdiag
    CImGui.SetNextItemWidth(widget_combo_width)
    offdiag_list = ["Yes", "No"]
    @cstatic id = Cint(1) begin
        @c CImGui.Combo(" Is it the offdiagonal correlation function", &id, offdiag_list)
        if id == 0
            PTEST.offdiag = true
        else
            PTEST.offdiag = false
        end
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(offdiag)$(PTEST.offdiag)")
    #
    # Input: lpeak
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic buf = "1,2,3" * "\0"^60 begin
        CImGui.InputText(" Number of peaks in the spectrum", buf, length(buf))
        _buf = split(rstrip(buf,'\0'), ",")   # Get rid of blanks
        filter!(x -> length(x) > 0, _buf)     # Remove empty strings
        lpeak = map(x -> parse(I64, x), _buf) # Convert strings to integers
        unique!(lpeak) # Remove redundant numbers
        PTEST.lpeak = lpeak
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(lpeak)$(PTEST.lpeak)")
end

"""
    _actest_solver_block()

Setup widgets for the `solver` tab.
"""
function _actest_solver_block()
    # It should change upon the selection of analytic continuation solver.
    CImGui.Text("Configure [Solver] Block / Analytic Continuation Solver: $(PTEST.solver)")

    @cswitch PTEST.solver begin

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
