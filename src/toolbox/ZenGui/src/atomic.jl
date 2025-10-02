#
# Project : Camellia
# Source  : atomic.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/10/02
#

"""
    create_app_atomic(p_open::Ref{Bool})

Create an UI window for the atomic code, which is an atomic eigenvalue
problem solver in the iQIST package.
"""
function create_app_atomic(p_open::Ref{Bool})
    # Create the atomic window, which can not be resized.
    CImGui.Begin(
        "iQIST | atomic",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "ATOMIC"
    end

    # Fix size of the window
    window_width = 1400.0
    window_height = 700.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _atomic_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the solver.atomic.in file
    _atomic_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _atomic_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _atomic_top_block()

Setup widgets in the top of the window for the iQIST/atomic code.
"""
function _atomic_top_block()
    CImGui.Text("atomic: An Atomic Eigenvalue Problem Solver")
end

"""
    _atomic_main_block()

Setup widgets associated with the parameters in the `solver.atomic.in` file.
"""
function _atomic_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    if CImGui.BeginTabBar("atomicTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        _atomic_model_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        _atomic_interaction_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_PURPLE)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_PURPLE)
        _atomic_natural_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_SALMON)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_SALMON)
        _atomic_algorithm_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _atomic_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the iQIST/atomic code.
"""
function _atomic_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 120.0
    widget_button_height = 60.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View solver.atomic.in")
    end
    #
    if CImGui.BeginPopupModal("View solver.atomic.in", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_ini(build_iqist_dict("atomic"))
            flags = CImGui.ImGuiInputTextFlags_ReadOnly
            flags = CImGui.ImGuiInputTextFlags_AllowTabInput | flags
            CImGui.InputTextMultiline("##source", text, 10000, ImVec2(600, 600), flags)
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
    _atomic_model_block()

Setup widgets for the `model` tab.
"""
function _atomic_model_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("model")
        CImGui.Text("Configure [model] Part")

        # Input: nband
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1) begin
            @c CImGui.SliderInt(" Number of correlated bands", &_i, 1, 7)
            PATOMIC.nband = _i
            _i != 1 && push!(_ATOMIC, "nband")
            _i == 1 && delete!(_ATOMIC, "nband")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nband)$(PATOMIC.nband)")
        #
        # Input: nspin
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of spin projections", &_i)
            _i = Cint(PATOMIC.nspin) # This parameter should not be changed.
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nspin)$(PATOMIC.nspin)")
        #
        # Input: norbs
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of correlated orbitals", &_i)
            _i = Cint(PATOMIC.nspin * PATOMIC.nband)
            PATOMIC.norbs = _i
            _i != 2 && push!(_ATOMIC, "norbs")
            _i == 2 && delete!(_ATOMIC, "norbs")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(norbs)$(PATOMIC.norbs)")
        #
        # Input: ncfgs
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(4) begin
            @c CImGui.InputInt(" Number of many-body configurations", &_i)
            _i = Cint(2 ^ PATOMIC.norbs)
            PATOMIC.ncfgs = _i
            _i != 4 && push!(_ATOMIC, "ncfgs")
            _i == 4 && delete!(_ATOMIC, "ncfgs")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ncfgs)$(PATOMIC.ncfgs)")

        CImGui.EndTabItem()
    end
end

"""
    _atomic_interaction_block()

Setup widgets for the `interaction` tab.
"""
function _atomic_interaction_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("interaction")
        CImGui.Text("Configure [interaction] Part")

        # Input: icu
        CImGui.SetNextItemWidth(widget_combo_width)
        icu_list = ["kanamori", "slater-cordon"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Type of Coulomb interaction matrix", &id, icu_list)
            PATOMIC.icu = id + 1
            id != 0 && push!(_ATOMIC, "icu")
            id == 0 && delete!(_ATOMIC, "icu")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(icu)$(PATOMIC.icu)")
        #
        # Input: Uc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(2.0) begin
            @c CImGui.InputDouble(" Intra-orbital Coulomb interaction", &_f)
            PATOMIC.Uc = _f
            _f != 2.0 && push!(_ATOMIC, "Uc")
            _f == 2.0 && delete!(_ATOMIC, "Uc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Uc)$(PATOMIC.Uc)")
        #
        # Input: Uv
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(2.0) begin
            @c CImGui.InputDouble(" Inter-orbital Coulomb interaction", &_f)
            PATOMIC.Uv = _f
            _f != 2.0 && push!(_ATOMIC, "Uv")
            _f == 2.0 && delete!(_ATOMIC, "Uv")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Uv)$(PATOMIC.Uv)")
        #
        # Input: Jz
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Hund's exchange interaction in z axis", &_f)
            PATOMIC.Jz = _f
            _f != 0.0 && push!(_ATOMIC, "Jz")
            _f == 0.0 && delete!(_ATOMIC, "Jz")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Jz)$(PATOMIC.Jz)")
        #
        # Input: Js
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Spin-flip interaction", &_f)
            PATOMIC.Js = _f
            _f != 0.0 && push!(_ATOMIC, "Js")
            _f == 0.0 && delete!(_ATOMIC, "Js")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Js)$(PATOMIC.Js)")
        #
        # Input: Jp
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Pair-hopping interaction", &_f)
            PATOMIC.Jp = _f
            _f != 0.0 && push!(_ATOMIC, "Jp")
            _f == 0.0 && delete!(_ATOMIC, "Jp")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Jp)$(PATOMIC.Jp)")
        #
        # Input: Ud
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(2.0) begin
            @c CImGui.InputDouble(" Coulomb interaction parameter (Slater type)", &_f)
            PATOMIC.Ud = _f
            _f != 2.0 && push!(_ATOMIC, "Ud")
            _f == 2.0 && delete!(_ATOMIC, "Ud")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Ud)$(PATOMIC.Ud)")
        #
        # Input: Jh
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Hund's exchange parameter (Slater type)", &_f)
            PATOMIC.Jh = _f
            _f != 0.0 && push!(_ATOMIC, "Jh")
            _f == 0.0 && delete!(_ATOMIC, "Jh")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Jh)$(PATOMIC.Jh)")

        CImGui.EndTabItem()
    end
end

"""
    _atomic_natural_block()

Setup widgets for the `natural eigenbasis` tab.
"""
function _atomic_natural_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("natural eigenbasis")
        CImGui.Text("Configure [natural eigenbasis] Part")

        # Input: ibasis
        CImGui.SetNextItemWidth(widget_combo_width)
        ibasis_list = ["builtin", "external"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" How to build the natural eigenbasis", &id, ibasis_list)
            PATOMIC.ibasis = id + 1
            id != 0 && push!(_ATOMIC, "ibasis")
            id == 0 && delete!(_ATOMIC, "ibasis")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ibasis)$(PATOMIC.ibasis)")
        #
        # Input: icf
        CImGui.SetNextItemWidth(widget_combo_width)
        icf_list = ["none", "diagonal", "off-diagonal"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Type of crystal field splitting", &id, icf_list)
            PATOMIC.icf = id
            id != 0 && push!(_ATOMIC, "icf")
            id == 0 && delete!(_ATOMIC, "icf")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(icf)$(PATOMIC.icf)")
        #
        # Input: isoc
        CImGui.SetNextItemWidth(widget_combo_width)
        isoc_list = ["none", "onsite"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Type of spin-orbit coupling", &id, isoc_list)
            PATOMIC.isoc = id
            id != 0 && push!(_ATOMIC, "isoc")
            id == 0 && delete!(_ATOMIC, "isoc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isoc)$(PATOMIC.isoc)")
        #
        # Input: mune
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Chemical potential or fermi level", &_f)
            PATOMIC.mune = _f
            _f != 0.0 && push!(_ATOMIC, "mune")
            _f == 0.0 && delete!(_ATOMIC, "mune")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mune)$(PATOMIC.mune)")
        #
        # Input: lambda
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Strength of spin-orbit coupling", &_f)
            PATOMIC.lambda = _f
            _f != 0.0 && push!(_ATOMIC, "lambda")
            _f == 0.0 && delete!(_ATOMIC, "lambda")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lambda)$(PATOMIC.lambda)")

        CImGui.EndTabItem()
    end
end

"""
    _atomic_algorithm_block()

Setup widgets for the `algorithm` tab.
"""
function _atomic_algorithm_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("algorithm")
        CImGui.Text("Configure [algorithm] Part")

        # Input: ictqmc
        CImGui.SetNextItemWidth(widget_combo_width)
        ictqmc_list = ["direct", "GQN (N)", "GQN (N+Sz)", "GQN (N+Sz+PS)", "GQN (N+Jz)", "auto"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" How to diagonalize the atomic Hamiltonian matrix", &id, ictqmc_list)
            PATOMIC.ictqmc = id + 1
            id != 0 && push!(_ATOMIC, "ictqmc")
            id == 0 && delete!(_ATOMIC, "ictqmc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ictqmc)$(PATOMIC.ictqmc)")
        #
        # Input: nmini
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(0) begin
            @c CImGui.SliderInt(" Lower boundary of occupancy N", &_i, 1, 14)
            PATOMIC.nmini = _i
            _i != 0 && push!(_ATOMIC, "nmini")
            _i == 0 && delete!(_ATOMIC, "nmini")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nmini)$(PATOMIC.nmini)")
        #
        # Input: nmaxi
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.SliderInt(" Upper boundary of occupancy N", &_i, 1, 14)
            PATOMIC.nmaxi = _i
            _i != 2 && push!(_ATOMIC, "nmaxi")
            _i == 2 && delete!(_ATOMIC, "nmaxi")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nmaxi)$(PATOMIC.nmaxi)")

        CImGui.EndTabItem()
    end
end
