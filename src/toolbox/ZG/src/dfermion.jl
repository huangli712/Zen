#
# Project : Camellia
# Source  : dfermion.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/26
#

"""
    create_app_dfermion(p_open::Ref{Bool})

Create an UI window for the DFermion code, which is a dual fermion engine.
"""
function create_app_dfermion(p_open::Ref{Bool})
    # Create the DFermion window, which can not be resized.
    CImGui.Begin(
        "DFermion",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "DFERMION"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _dfermion_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the dfa.in file
    _dfermion_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _dfermion_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _dfermion_top_block()

Setup widgets in the top of the window for the DFermion code.
"""
function _dfermion_top_block()
    CImGui.Text("DFermion: Dual Fermion Application")
end

"""
    _dfermion_main_block()

Setup widgets associated with the parameters in the `dfa.in` file.
"""
function _dfermion_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    # There are four tabs
    if CImGui.BeginTabBar("dfermionTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        _dfermion_model_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        _dfermion_dimension_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_PURPLE)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_PURPLE)
        _dfermion_kmesh_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_SALMON)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_SALMON)
        _dfermion_cycle_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _dfermion_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the DFermion code.
"""
function _dfermion_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 80.0
    widget_button_height = 25.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View dfa.in")
    end
    #
    if CImGui.BeginPopupModal("View dfa.in", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_ini(build_dfermion_dict())
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
    _dfermion_model_block()

Setup widgets for the `model` tab.
"""
function _dfermion_model_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("model")
        CImGui.Text("Configure [model] Part")

        # Input: nband
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1) begin
            @c CImGui.SliderInt(" Number of correlated bands", &_i, 1, 7)
            PDFERMION.nband = _i
            _i != 1 && push!(_DFERMION, "nband")
            _i == 1 && delete!(_DFERMION, "nband")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nband)$(PDFERMION.nband)")
        #
        # Input: nspin
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of spin projections", &_i)
            _i = Cint(PDFERMION.nspin) # This parameter should not be changed.
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nspin)$(PDFERMION.nspin)")
        #
        # Input: norbs
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of correlated orbitals", &_i)
            _i = Cint(PDFERMION.nspin * PDFERMION.nband)
            PDFERMION.norbs = _i
            _i != 2 && push!(_DFERMION, "norbs")
            _i == 2 && delete!(_DFERMION, "norbs")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(norbs)$(PDFERMION.norbs)")
        #
        # Input: mune
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Chemical potential or fermi level", &_f)
            PDFERMION.mune = _f
            _f != 0.0 && push!(_DFERMION, "mune")
            _f == 0.0 && delete!(_DFERMION, "mune")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mune)$(PDFERMION.mune)")
        #
        # Input: beta
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0) begin
            @c CImGui.InputDouble(" Inversion of temperature", &_f)
            PDFERMION.beta = _f
            _f != 1.0 && push!(_DFERMION, "beta")
            _f == 1.0 && delete!(_DFERMION, "beta")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(beta)$(PDFERMION.beta)")
        #
        # Input: part
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0) begin
            @c CImGui.InputDouble(" Hopping parameter t for Hubbard model", &_f)
            PDFERMION.part = _f
            _f != 1.0 && push!(_DFERMION, "part")
            _f == 1.0 && delete!(_DFERMION, "part")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(part)$(PDFERMION.part)")

        CImGui.EndTabItem()
    end
end

"""
    _dfermion_dimension_block()

Setup widgets for the `dimension` tab.
"""
function _dfermion_dimension_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("dimension")
        CImGui.Text("Configure [dimension] Part")

        # Input: nffrq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(16) begin
            @c CImGui.SliderInt(" Number of fermionic frequencies for 2P GF", &_i, 8, 1024)
            PDFERMION.nffrq = _i
            _i != 16 && push!(_DFERMION, "nffrq")
            _i == 16 && delete!(_DFERMION, "nffrq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nffrq)$(PDFERMION.nffrq)")
        #
        # Input: nbfrq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(7) begin
            @c CImGui.SliderInt(" Number of bosonic frequncies for 2P GF", &_i, 4, 512)
            PDFERMION.nbfrq = _i
            _i != 7 && push!(_DFERMION, "nbfrq")
            _i == 7 && delete!(_DFERMION, "nbfrq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nbfrq)$(PDFERMION.nbfrq)")

        CImGui.EndTabItem()
    end
end

"""
    _dfermion_kmesh_block()

Setup widgets for the `k-mesh` tab.
"""
function _dfermion_kmesh_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("k-mesh")
        CImGui.Text("Configure [k-mesh] Part")

        # Input: nkpts
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(64) begin
            @c CImGui.InputInt(" Number of k-points (totally)", &_i)
            PDFERMION.nkpts = _i
            _i != 64 && push!(_DFERMION, "nkpts")
            _i == 64 && delete!(_DFERMION, "nkpts")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nkpts)$(PDFERMION.nkpts)")
        #
        # Input: nkp_x
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8) begin
            @c CImGui.InputInt(" Number of k-points (along x-axis)", &_i)
            PDFERMION.nkp_x = _i
            _i != 8 && push!(_DFERMION, "nkp_x")
            _i == 8 && delete!(_DFERMION, "nkp_x")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nkp_x)$(PDFERMION.nkp_x)")
        #
        # Input: nkp_y
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8) begin
            @c CImGui.InputInt(" Number of k-points (along y-axis)", &_i)
            PDFERMION.nkp_y = _i
            _i != 8 && push!(_DFERMION, "nkp_y")
            _i == 8 && delete!(_DFERMION, "nkp_y")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nkp_y)$(PDFERMION.nkp_y)")
        #
        # Input: nkp_z
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8) begin
            @c CImGui.InputInt(" Number of k-points (along z-axis)", &_i)
            PDFERMION.nkp_z = _i
            _i != 8 && push!(_DFERMION, "nkp_z")
            _i == 8 && delete!(_DFERMION, "nkp_z")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nkp_z)$(PDFERMION.nkp_z)")

        CImGui.EndTabItem()
    end
end

"""
    _dfermion_cycle_block()

Setup widgets for the `cycle` tab.
"""
function _dfermion_cycle_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("cycle")
        CImGui.Text("Configure [cycle] Part")

        # Input: isdia
        CImGui.SetNextItemWidth(widget_combo_width)
        isdia_list = ["second order", "ladder", "to be done"]
        @cstatic id = Cint(1) begin
            @c CImGui.Combo(" Running scheme of the code", &id, isdia_list)
            PDFERMION.isdia = id + 1
            id != 1 && push!(_DFERMION, "isdia")
            id == 1 && delete!(_DFERMION, "isdia")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isdia)$(PDFERMION.isdia)")
        #
        # Input: ndfit
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(10) begin
            @c CImGui.SliderInt(" Number of dual fermion iterations", &_i, 1, 100)
            PDFERMION.ndfit = _i
            _i != 10 && push!(_DFERMION, "ndfit")
            _i == 10 && delete!(_DFERMION, "ndfit")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ndfit)$(PDFERMION.ndfit)")
        #
        # Input: nbsit
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(10) begin
            @c CImGui.SliderInt(" Number of iterations for solving the BSE", &_i, 1, 100)
            PDFERMION.nbsit = _i
            _i != 10 && push!(_DFERMION, "nbsit")
            _i == 10 && delete!(_DFERMION, "nbsit")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nbsit)$(PDFERMION.nbsit)")
        #
        # Input: dfmix
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0) vmin = Cdouble(0.0) vmax = Cdouble(1.0) begin
            @c CImGui.SliderScalar(
                " Mixing parameter for dual fermion iteration",
                CImGui.ImGuiDataType_Double,
                &_f,
                &vmin, &vmax
            )
            PDFERMION.dfmix = _f
            _f != 1.0 && push!(_DFERMION, "dfmix")
            _f == 1.0 && delete!(_DFERMION, "dfmix")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(dfmix)$(PDFERMION.dfmix)")
        #
        # Input: bsmix
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.7) vmin = Cdouble(0.0) vmax = Cdouble(1.0) begin
            @c CImGui.SliderScalar(
                " Mixing parameter for solving the BSE",
                CImGui.ImGuiDataType_Double,
                &_f,
                &vmin, &vmax
            )
            PDFERMION.bsmix = _f
            _f != 0.7 && push!(_DFERMION, "bsmix")
            _f == 0.7 && delete!(_DFERMION, "bsmix")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(bsmix)$(PDFERMION.bsmix)")

        CImGui.EndTabItem()
    end
end
