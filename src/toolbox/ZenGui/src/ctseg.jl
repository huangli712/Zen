#
# Project : Camellia
# Source  : ctseg.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/10/02
#

"""
    create_app_ctseg(p_open::Ref{Bool})

Create an UI window for the ctseg code, which is a continuous-time quantum
impurity solver in the iQIST package.
"""
function create_app_ctseg(p_open::Ref{Bool})
    # Create the ctseg window, which can not be resized.
    CImGui.Begin(
        "iQIST | ctseg",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "CTSEG"
    end

    # Fix size of the window
    window_width = 1500.0
    window_height = 900.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _ctseg_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the solver.ctqmc.in file
    _ctseg_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _ctseg_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _ctseg_top_block()

Setup widgets in the top of the window for the iQIST/ctseg code.
"""
function _ctseg_top_block()
    CImGui.Text("ctseg: A Continuous-Time Quantum Impurity Solver")
end

"""
    _ctseg_main_block()

Setup widgets associated with the parameters in the `solver.ctqmc.in` file.
"""
function _ctseg_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    if CImGui.BeginTabBar("ctsegTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        _ctseg_model_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        _ctseg_dimension_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_PURPLE)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_PURPLE)
        _ctseg_symmetry_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_SALMON)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_SALMON)
        _ctseg_represent_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TAN)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TAN)
        _ctseg_monte_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TOMATO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TOMATO)
        _ctseg_measure_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_HOTPINK)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_HOTPINK)
        _ctseg_cycle_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _ctseg_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the iQIST/ctseg code.
"""
function _ctseg_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 120.0
    widget_button_height = 60.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View solver.ctqmc.in")
    end
    #
    if CImGui.BeginPopupModal("View solver.ctqmc.in", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_ini(build_iqist_dict("ctseg"))
            flags = CImGui.ImGuiInputTextFlags_ReadOnly
            flags = CImGui.ImGuiInputTextFlags_AllowTabInput | flags
            CImGui.InputTextMultiline("##source", text, 10000, ImVec2(600, 800), flags)
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
    _ctseg_model_block()

Setup widgets for the `model` tab.
"""
function _ctseg_model_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("model")
        CImGui.Text("Configure [model] Part")

        # Input: isscr
        CImGui.SetNextItemWidth(widget_combo_width)
        isscr_list = ["static", "plasmon pole", "ohmic", "realistic"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Whether the Coulomb interaction U is dynamic", &id, isscr_list)
            PCTSEG.isscr = id + 1
            id != 0 && push!(_CTSEG, "isscr")
            id == 0 && delete!(_CTSEG, "isscr")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isscr)$(PCTSEG.isscr)")
        #
        # Input: nband
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1) begin
            @c CImGui.SliderInt(" Number of correlated bands", &_i, 1, 7)
            PCTSEG.nband = _i
            _i != 1 && push!(_CTSEG, "nband")
            _i == 1 && delete!(_CTSEG, "nband")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nband)$(PCTSEG.nband)")
        #
        # Input: nspin
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of spin projections", &_i)
            _i = Cint(PCTSEG.nspin) # This parameter should not be changed.
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nspin)$(PCTSEG.nspin)")
        #
        # Input: norbs
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of correlated orbitals", &_i)
            _i = Cint(PCTSEG.nspin * PCTSEG.nband)
            PCTSEG.norbs = _i
            _i != 2 && push!(_CTSEG, "norbs")
            _i == 2 && delete!(_CTSEG, "norbs")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(norbs)$(PCTSEG.norbs)")
        #
        # Input: ncfgs
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(4) begin
            @c CImGui.InputInt(" Number of atomic eigenstates", &_i)
            _i = Cint(2 ^ PCTSEG.norbs)
            PCTSEG.ncfgs = _i
            _i != 4 && push!(_CTSEG, "ncfgs")
            _i == 4 && delete!(_CTSEG, "ncfgs")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ncfgs)$(PCTSEG.ncfgs)")
        #
        # Input: Uc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(4.0) begin
            @c CImGui.InputDouble(" Intra-orbital Coulomb interaction", &_f)
            PCTSEG.Uc = _f
            _f != 4.0 && push!(_CTSEG, "Uc")
            _f == 4.0 && delete!(_CTSEG, "Uc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Uc)$(PCTSEG.Uc)")
        #
        # Input: Jz
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0) begin
            @c CImGui.InputDouble(" Hund's exchange interaction in z axis", &_f)
            PCTSEG.Jz = _f
            _f != 0.0 && push!(_CTSEG, "Jz")
            _f == 0.0 && delete!(_CTSEG, "Jz")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(Jz)$(PCTSEG.Jz)")
        #
        # Input: lc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0) begin
            @c CImGui.InputDouble(" Strength of dynamical screening effect", &_f)
            PCTSEG.lc = _f
            _f != 1.0 && push!(_CTSEG, "lc")
            _f == 1.0 && delete!(_CTSEG, "lc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lc)$(PCTSEG.lc)")
        #
        # Input: wc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0) begin
            @c CImGui.InputDouble(" Screening frequency", &_f)
            PCTSEG.wc = _f
            _f != 1.0 && push!(_CTSEG, "wc")
            _f == 1.0 && delete!(_CTSEG, "wc")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(wc)$(PCTSEG.wc)")
        #
        # Input: mune
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(2.0) begin
            @c CImGui.InputDouble(" Chemical potential or fermi level", &_f)
            PCTSEG.mune = _f
            _f != 2.0 && push!(_CTSEG, "mune")
            _f == 2.0 && delete!(_CTSEG, "mune")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mune)$(PCTSEG.mune)")
        #
        # Input: beta
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(8.0) begin
            @c CImGui.InputDouble(" Inversion of temperature", &_f)
            PCTSEG.beta = _f
            _f != 8.0 && push!(_CTSEG, "beta")
            _f == 8.0 && delete!(_CTSEG, "beta")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(beta)$(PCTSEG.beta)")
        #
        # Input: part
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.5) begin
            @c CImGui.InputDouble(" Hopping parameter t for Hubbard model", &_f)
            PCTSEG.part = _f
            _f != 0.5 && push!(_CTSEG, "part")
            _f == 0.5 && delete!(_CTSEG, "part")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(part)$(PCTSEG.part)")

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_dimension_block()

Setup widgets for the `dimension` tab.
"""
function _ctseg_dimension_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("dimension")
        CImGui.Text("Configure [dimension] Part")

        # Input: mfreq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8193) begin
            @c CImGui.SliderInt(" Number of correlated bands", &_i, 2^10+1, 2^14+1)
            PCTSEG.mfreq = _i
            _i != 8193 && push!(_CTSEG, "mfreq")
            _i == 8193 && delete!(_CTSEG, "mfreq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mfreq)$(PCTSEG.mfreq)")
        #
        # Input: nffrq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(32) begin
            @c CImGui.SliderInt(" Number of fermionic frequencies for 2P function", &_i, 8, 1024)
            PCTSEG.nffrq = _i
            _i != 32 && push!(_CTSEG, "nffrq")
            _i == 32 && delete!(_CTSEG, "nffrq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nffrq)$(PCTSEG.nffrq)")
        #
        # Input: nbfrq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8) begin
            @c CImGui.SliderInt(" Number of bosonic frequncies for 2P function", &_i, 4, 512)
            PCTSEG.nbfrq = _i
            _i != 8 && push!(_CTSEG, "nbfrq")
            _i == 8 && delete!(_CTSEG, "nbfrq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nbfrq)$(PCTSEG.nbfrq)")
        #
        # Input: nfreq
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(128) begin
            @c CImGui.SliderInt(" Number of matsubara frequencies sampled by solver", &_i, 64, 1024)
            PCTSEG.nfreq = _i
            _i != 128 && push!(_CTSEG, "nfreq")
            _i == 128 && delete!(_CTSEG, "nfreq")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nfreq)$(PCTSEG.nfreq)")
        #
        # Input: ntime
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1024) begin
            @c CImGui.SliderInt(" Number of imaginary time slices sampled by solver", &_i, 256, 10240)
            PCTSEG.ntime = _i
            _i != 1024 && push!(_CTSEG, "ntime")
            _i == 1024 && delete!(_CTSEG, "ntime")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ntime)$(PCTSEG.ntime)")

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_symmetry_block()

Setup widgets for the `symmetry` tab.
"""
function _ctseg_symmetry_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("symmetry")
        CImGui.Text("Configure [symmetry] Part")

        # Input: isbnd
        CImGui.SetNextItemWidth(widget_combo_width)
        isbnd_list = ["no", "yes"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Symmetry of the impurity model (band part)", &id, isbnd_list)
            PCTSEG.isbnd = id + 1
            id != 0 && push!(_CTSEG, "isbnd")
            id == 0 && delete!(_CTSEG, "isbnd")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isbnd)$(PCTSEG.isbnd)")
        #
        # Input: isspn
        CImGui.SetNextItemWidth(widget_combo_width)
        isspn_list = ["no", "yes"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Symmetry of the impurity model (spin part)", &id, isspn_list)
            PCTSEG.isspn = id + 1
            id != 0 && push!(_CTSEG, "isspn")
            id == 0 && delete!(_CTSEG, "isspn")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isspn)$(PCTSEG.isspn)")

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_represent_block()

Setup widgets for the `representation` tab.
"""
function _ctseg_represent_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("representation")
        CImGui.Text("Configure [representation] Part")

        # Input: isort
        CImGui.SetNextItemWidth(widget_combo_width)
        isort_list = ["standard", "legendre", "intermediate"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Which basis will be used to do the measurement", &id, isort_list)
            PCTSEG.isort = id + 1
            id != 0 && push!(_CTSEG, "isort")
            id == 0 && delete!(_CTSEG, "isort")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isort)$(PCTSEG.isort)")
        #
        # Input: lemax
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(32) begin
            @c CImGui.InputInt(" Maximum expansion order for legendre polynomial", &_i)
            PCTSEG.lemax = _i
            _i != 32 && push!(_CTSEG, "lemax")
            _i == 32 && delete!(_CTSEG, "lemax")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lemax)$(PCTSEG.lemax)")
        #
        # Input: legrd
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(20001) begin
            @c CImGui.InputInt(" Number of mesh points for legendre polynomial", &_i)
            PCTSEG.legrd = _i
            _i != 20001 && push!(_CTSEG, "legrd")
            _i == 20001 && delete!(_CTSEG, "legrd")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(legrd)$(PCTSEG.legrd)")
        #
        # Input: svmax
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(32) begin
            @c CImGui.InputInt(" Maximum expansion order for svd polynomial", &_i)
            PCTSEG.svmax = _i
            _i != 32 && push!(_CTSEG, "svmax")
            _i == 32 && delete!(_CTSEG, "svmax")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(svmax)$(PCTSEG.svmax)")
        #
        # Input: svgrd
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2001) begin
            @c CImGui.InputInt(" Number of mesh points for svd polynomial", &_i)
            PCTSEG.svgrd = _i
            _i != 2001 && push!(_CTSEG, "svgrd")
            _i == 2001 && delete!(_CTSEG, "svgrd")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(svgrd)$(PCTSEG.svgrd)")

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_monte_block()

Setup widgets for the `monte carlo` tab.
"""
function _ctseg_monte_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("monte carlo")
        CImGui.Text("Configure [monte carlo] Part")

        # Input: mkink
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1024) begin
            @c CImGui.InputInt(" Maximum perturbation expansion order", &_i)
            PCTSEG.mkink = _i
            _i != 1024 && push!(_CTSEG, "mkink")
            _i == 1024 && delete!(_CTSEG, "mkink")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mkink)$(PCTSEG.mkink)")
        #
        # Input: nflip
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(20000) begin
            @c CImGui.InputInt(" Flip period for spin up and spin down states", &_i)
            PCTSEG.nflip = _i
            _i != 20000 && push!(_CTSEG, "nflip")
            _i == 20000 && delete!(_CTSEG, "nflip")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nflip)$(PCTSEG.nflip)")
        #
        # Input: ntherm
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(200000) begin
            @c CImGui.InputInt(" Number of thermalization steps", &_i)
            PCTSEG.ntherm = _i
            _i != 200000 && push!(_CTSEG, "ntherm")
            _i == 200000 && delete!(_CTSEG, "ntherm")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ntherm)$(PCTSEG.ntherm)")
        #
        # Input: nsweep
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(20000000) begin
            @c CImGui.InputInt(" Number of Monte Carlo sweeping steps", &_i)
            PCTSEG.nsweep = _i
            _i != 20000000 && push!(_CTSEG, "nsweep")
            _i == 20000000 && delete!(_CTSEG, "nsweep")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nsweep)$(PCTSEG.nsweep)")
        #
        # Input: nwrite
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2000000) begin
            @c CImGui.InputInt(" Output period for quantum impurity solver", &_i)
            PCTSEG.nwrite = _i
            _i != 2000000 && push!(_CTSEG, "nwrite")
            _i == 2000000 && delete!(_CTSEG, "nwrite")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nwrite)$(PCTSEG.nwrite)")
        #
        # Input: nclean
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(100000) begin
            @c CImGui.InputInt(" Clean update period for quantum impurity solver", &_i)
            PCTSEG.nclean = _i
            _i != 100000 && push!(_CTSEG, "nclean")
            _i == 100000 && delete!(_CTSEG, "nclean")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nclean)$(PCTSEG.nclean)")
        #
        # Input: nmonte
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(10) begin
            @c CImGui.InputInt(" How often to sample the physical observables 1", &_i)
            PCTSEG.nmonte = _i
            _i != 10 && push!(_CTSEG, "nmonte")
            _i == 10 && delete!(_CTSEG, "nmonte")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nmonte)$(PCTSEG.nmonte)")
        #
        # Input: ncarlo
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(10) begin
            @c CImGui.InputInt(" How often to sample the physical observables 2", &_i)
            PCTSEG.ncarlo = _i
            _i != 10 && push!(_CTSEG, "ncarlo")
            _i == 10 && delete!(_CTSEG, "ncarlo")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ncarlo)$(PCTSEG.ncarlo)")

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_measure_block()

Setup widgets for the `measurement` tab.
"""
function _ctseg_measure_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("measurement")
        CImGui.Text("Configure [measurement] Part")

        # Input: iswor
        CImGui.SetNextItemWidth(widget_combo_width)
        iswor_list = ["standard", "worm"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Which algorithm will be used to do the measurement", &id, iswor_list)
            PCTSEG.iswor = id + 1
            id != 0 && push!(_CTSEG, "iswor")
            id == 0 && delete!(_CTSEG, "iswor")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(iswor)$(PCTSEG.iswor)")
        #
        # Input: isobs
        if CImGui.CollapsingHeader("Calculate physical observables")
            @cstatic kin = false fid = false mag = false begin
                @c CImGui.Checkbox("kinetic energy fluctuation", &kin)
                @c CImGui.Checkbox("Fidelity susceptibility", &fid)
                @c CImGui.Checkbox("Powers of local magnetization", &mag)
                isobs = 1 + Int(kin) * 2^1 + Int(fid) * 2^2 + Int(mag) * 2^3
                PCTSEG.isobs = isobs
                isobs != 1 && push!(_CTSEG, "isobs")
                isobs == 1 && delete!(_CTSEG, "isobs")
            end
        end
        #
        # Input: issus
        if CImGui.CollapsingHeader("Calculate susceptibilities")
            @cstatic st = false ct = false sf = false cf = false begin
                @c CImGui.Checkbox("Spin-spin correlation function (time space)", &st)
                @c CImGui.Checkbox("Charge-charge correlation function (time space)", &ct)
                @c CImGui.Checkbox("Spin-spin correlation function (frequency space)", &sf)
                @c CImGui.Checkbox("Charge-charge correlation function (frequency space)", &cf)
                issus = 1 + Int(st) * 2^1 + Int(ct) * 2^2 + Int(sf) * 2^3 + Int(cf) * 2^4
                PCTSEG.issus = issus
                issus != 1 && push!(_CTSEG, "issus")
                issus == 1 && delete!(_CTSEG, "issus")
            end
        end
        #
        # Input: isvrt
        if CImGui.CollapsingHeader("Calculate two-particle Green's functions")
            @cstatic p1 = Cint(0) p2 = Cint(0) begin
                @c CImGui.RadioButton("Block: AABB / Channel: particle-hole", &p1, 1)
                @c CImGui.RadioButton("Block: ABBA / Channel: particle-hole", &p1, 2)
                @c CImGui.RadioButton("Block: AABB / Channel: particle-particle", &p2, 1)
                @c CImGui.RadioButton("Block: ABBA / Channel: particle-particle", &p2, 2)
                isvrt = 1
                p1 == 1 && (isvrt = isvrt + 2^1)
                p1 == 2 && (isvrt = isvrt + 2^2)
                p2 == 1 && (isvrt = isvrt + 2^3)
                p2 == 2 && (isvrt = isvrt + 2^4)
                PCTSEG.isvrt = isvrt
                isvrt != 1 && push!(_CTSEG, "isvrt")
                isvrt == 1 && delete!(_CTSEG, "isvrt")
            end
        end

        CImGui.EndTabItem()
    end
end

"""
    _ctseg_cycle_block()

Setup widgets for the `cycle` tab.
"""
function _ctseg_cycle_block()
    # Define default size for widgets
    widget_input_width = 180
    widget_combo_width = 180

    if CImGui.BeginTabItem("cycle")
        CImGui.Text("Configure [cycle] Part")

        # Input: isscf
        CImGui.SetNextItemWidth(widget_combo_width)
        isscf_list = ["one-shot", "self-consistent"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Scheme of dynamical mean-field theory calculations", &id, isscf_list)
            PCTSEG.isscf = id + 1
            id != 0 && push!(_CTSEG, "isscf")
            id == 0 && delete!(_CTSEG, "isscf")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isscf)$(PCTSEG.isscf)")
        #
        # Input: niter
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(20) begin
            @c CImGui.SliderInt(" Number of self-consistent iterations", &_i, 1, 100)
            PCTSEG.niter = _i
            _i != 20 && push!(_CTSEG, "niter")
            _i == 20 && delete!(_CTSEG, "niter")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(niter)$(PCTSEG.niter)")
        #
        # Input: alpha
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.7) vmin = Cdouble(0.0) vmax = Cdouble(1.0) begin
            @c CImGui.SliderScalar(
                " Mixing factor for self-consistent engine",
                CImGui.ImGuiDataType_Double,
                &_f,
                &vmin, &vmax
            )
            PCTSEG.alpha = _f
            _f != 0.7 && push!(_CTSEG, "alpha")
            _f == 0.7 && delete!(_CTSEG, "alpha")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(alpha)$(PCTSEG.alpha)")

        CImGui.EndTabItem()
    end
end
