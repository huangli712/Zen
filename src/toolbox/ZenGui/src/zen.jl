#
# Project : Camellia
# Source  : zen.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/4/26
#

"""
    @widgets_generator_dft

Macro for generating codes for the `dft` tab.

See also: [`_zen_dft_block`](@ref).
"""
macro widgets_generator_dft(x)
    ex = quote
        i = $x

        # Input: sproj
        @cstatic buf = "1 : d : Pr" * "\0"^60 begin
            CImGui.SetNextItemWidth(widget_input_width)
            CImGui.InputText(" Specifications for generating projector $i", buf, length(buf))
            PDFT.sproj[i] = rstrip(buf,'\0')
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(sproj_$i)$(PDFT.sproj[i])")
        end
        #
        # Input: window
        @cstatic vec = Cdouble[-1.4,6.0] begin
            CImGui.SetNextItemWidth(widget_input_width * 2)
            CImGui.InputScalarN(
                " Band window for normalizing projector $i",
                CImGui.ImGuiDataType_Double,
                vec,
                2
            )
            PDFT.window[2*i-1] = vec[1]
            PDFT.window[2*i] = vec[2]
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(window_$i)$(PDFT.window[2*i-1:2*i])")
        end
    end

    return :( $(esc(ex)) )
end

"""
    @widgets_generator_impurity

Macro for generating codes for the `impurity` tab.

See also: [`_zen_impurity_block`](@ref).
"""
macro widgets_generator_impurity(x)
    ex = quote
        i = $x

        # Input: atoms
        @cstatic buf = "V : 2" * "\0"^60 begin
            CImGui.SetNextItemWidth(widget_input_width)
            CImGui.InputText(" Chemical symbols of impurity atom $i", buf, length(buf))
            PIMPURITY.atoms[i] = rstrip(buf,'\0')
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(atoms_$i)$(PIMPURITY.atoms[i])")
        end
        #
        # Input: equiv
        @cstatic _i = Cint(1) begin
            CImGui.SetNextItemWidth(widget_input_width)
            @c CImGui.InputInt(" Equivalency of quantum impurity atom $i", &_i)
            PIMPURITY.equiv[i] = _i
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(equiv_$i)$(PIMPURITY.equiv[i])")
        end
        #
        # Input: shell
        @cstatic id = Cint(0) begin
            CImGui.SetNextItemWidth(widget_combo_width)
            shell_list = ["s", "p", "d", "f", "d_t2g", "d_eg"]
            @c CImGui.Combo(" Angular momenta of correlated orbital $i", &id, shell_list)
            PIMPURITY.shell[i] = shell_list[id + 1]
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(shell_$i)$(PIMPURITY.shell[i])")
        end
        #
        # Input: ising
        @cstatic id = Cint(0) begin
            CImGui.SetNextItemWidth(widget_combo_width)
            ising_list = ["ising", "full"]
            @c CImGui.Combo(" Interaction types of correlated orbital $i", &id, ising_list)
            PIMPURITY.ising[i] = ising_list[id + 1]
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(ising_$i)$(PIMPURITY.ising[i])")
        end
        #
        # Input: occup
        @cstatic _f = Cdouble(1.0) begin
            CImGui.SetNextItemWidth(widget_input_width)
            @c CImGui.InputDouble(" Nominal impurity occupancy $i", &_f)
            PIMPURITY.occup[i] = _f
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(occup_$i)$(PIMPURITY.occup[i])")
        end
        #
        # Input: upara
        @cstatic _f = Cdouble(4.0) begin
            CImGui.SetNextItemWidth(widget_input_width)
            @c CImGui.InputDouble(" Coulomb interaction parameter $i", &_f)
            PIMPURITY.upara[i] = _f
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(upara_$i)$(PIMPURITY.upara[i])")
        end
        #
        # Input: jpara
        @cstatic _f = Cdouble(0.7) begin
            CImGui.SetNextItemWidth(widget_input_width)
            @c CImGui.InputDouble(" Hund's coupling parameter $i", &_f)
            PIMPURITY.jpara[i] = _f
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(jpara_$i)$(PIMPURITY.jpara[i])")
        end
        #
        # Input: lpara
        @cstatic _f = Cdouble(0.0) begin
            CImGui.SetNextItemWidth(widget_input_width)
            @c CImGui.InputDouble(" Spin-orbit coupling parameter $i", &_f)
            PIMPURITY.lpara[i] = _f
            CImGui.SameLine()
            CImGui.TextColored(COL_MAGENTA, "(lpara_$i)$(PIMPURITY.lpara[i])")
        end
    end

    return :( $(esc(ex)) )
end

"""
    create_app_zen(p_open::Ref{Bool})

Create an UI window for the Zen toolkit, which is an integrated package
for ab initio dynamical mean-field theory calculations.
"""
function create_app_zen(p_open::Ref{Bool})
    # Create the Zen window, which can not be resized.
    CImGui.Begin(
        "Zen",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "ZEN"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _zen_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the case.toml file
    _zen_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _zen_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _zen_top_block()

Setup widgets in the top of the window for the Zen package.
"""
function _zen_top_block()
    CImGui.Text("Zen: A modern DFT + DMFT computation framework")
end

"""
    _zen_main_block()

Setup widgets associated with the parameters in the `case.toml` file.
"""
function _zen_main_block()
    tab_bar_flags = CImGui.ImGuiTabBarFlags_None
    #
    if CImGui.BeginTabBar("ZenTabBar", tab_bar_flags)
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TEAL)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TEAL)
        _zen_case_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_INDIGO)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_INDIGO)
        _zen_dft_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_PURPLE)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_PURPLE)
        _zen_dmft_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_SALMON)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_SALMON)
        _zen_impurity_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.PushStyleColor(CImGui.ImGuiCol_Tab, COL_TAN)
        CImGui.PushStyleColor(CImGui.ImGuiCol_TabSelected, COL_TAN)
        _zen_solver_block()
        CImGui.PopStyleColor(2)
        #
        CImGui.EndTabBar()
    end
end

"""
    _zen_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the Zen package.
"""
function _zen_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 80.0
    widget_button_height = 25.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View case.toml")
    end
    #
    if CImGui.BeginPopupModal("View case.toml", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_toml(build_zen_dict())
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
    _zen_case_block()

Setup widgets for the `case` tab.
"""
function _zen_case_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("case")
        CImGui.Text("Configure [case] Block")

        # Input: case
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic buf = "SrVO3" * "\0"^60 begin
            CImGui.InputText(" System's name or seedname", buf, length(buf))
            PCASE.case = rstrip(buf,'\0')
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(case)$(PCASE.case)")

        CImGui.EndTabItem()
    end
end

"""
    _zen_dft_block()

Setup widgets for the `dft` tab.
"""
function _zen_dft_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("dft")
        CImGui.Text("Configure [dft] Block")

        # Input: engine
        CImGui.SetNextItemWidth(widget_combo_width)
        engine_list = ["vasp", "qe"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Name of density functional theory code", &id, engine_list)
            PDFT.engine = engine_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(engine)$(PDFT.engine)")
        #
        # Input: projtype
        CImGui.SetNextItemWidth(widget_combo_width)
        projtype_list = ["plo", "wannier"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Types of projectors", &id, projtype_list)
            PDFT.projtype = projtype_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(projtype)$(PDFT.projtype)")
        #
        # Input: smear
        CImGui.SetNextItemWidth(widget_combo_width)
        smear_list = ["mp2", "mp1", "gauss", "tetra"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Scheme for smearing", &id, smear_list)
            PDFT.smear = smear_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(smear)$(PDFT.smear)")
        #
        # Input: kmesh
        CImGui.SetNextItemWidth(widget_combo_width)
        kmesh_list = ["accurate", "medium", "coarse", "file"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" K-mesh for brillouin zone integration", &id, kmesh_list)
            PDFT.kmesh = kmesh_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(kmesh)$(PDFT.kmesh)")
        #
        # Input: magmom
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic buf = "0.0" * "\0"^60 begin
            CImGui.InputText(" Initial magnetic moments", buf, length(buf))
            PDFT.magmom = rstrip(buf,'\0')
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(magmom)$(PDFT.magmom)")
        #
        # Input: ncycle
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8) begin
            @c CImGui.InputInt(" Number of DFT iterations per DFT + DMFT cycle", &_i)
            PDFT.ncycle = _i
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ncycle)$(PDFT.ncycle)")
        #
        # Input: lsymm
        CImGui.SetNextItemWidth(widget_combo_width)
        lsymm_list = ["yes", "no"]
        @cstatic id = Cint(1) begin
            @c CImGui.Combo(" Is the symmetry turned on or off", &id, lsymm_list)
            if id == 0
                PDFT.lsymm = true
            else
                PDFT.lsymm = false
            end
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lsymm)$(PDFT.lsymm)")
        #
        # Input: lspins
        CImGui.SetNextItemWidth(widget_combo_width)
        lspins_list = ["yes", "no"]
        @cstatic id = Cint(1) begin
            @c CImGui.Combo(" Are the spin orientations polarized or not", &id, lspins_list)
            if id == 0
                PDFT.lspins = true
            else
                PDFT.lspins = false
            end
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lspins)$(PDFT.lspins)")
        #
        # Input: lspinorb
        CImGui.SetNextItemWidth(widget_combo_width)
        lspinorb_list = ["yes", "no"]
        @cstatic id = Cint(1) begin
            @c CImGui.Combo(" Is the spin-orbit coupling considered or not", &id, lspinorb_list)
            if id == 0
                PDFT.lspinorb = true
            else
                PDFT.lspinorb = false
            end
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lspinorb)$(PDFT.lspinorb)")
        #
        # Input: lproj
        CImGui.SetNextItemWidth(widget_combo_width)
        lproj_list = ["yes", "no"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Are the projectors generated or not", &id, lproj_list)
            if id == 0
                PDFT.lproj = true
            else
                PDFT.lproj = false
            end
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lproj)$(PDFT.lproj)")

        # For the separator
        CImGui.Spacing()
        CImGui.Separator()
        CImGui.Spacing()

        # Input: nsite
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(1) begin
            @c CImGui.SliderInt(" Number of (correlated) impurity sites", &_i, 1, 9)
            PIMPURITY.nsite = _i
            #
            resize!(PDFT.sproj, PIMPURITY.nsite)
            fill!(PDFT.sproj, PDFT.sproj[1])
            resize!(PDFT.window, PIMPURITY.nsite * 2)
            for i = 2:PIMPURITY.nsite
                PDFT.window[2*i-1] = PDFT.window[1]
                PDFT.window[2*i] = PDFT.window[2]
            end
            #
            resize!(PIMPURITY.atoms, PIMPURITY.nsite)
            fill!(PIMPURITY.atoms, PIMPURITY.atoms[1])
            resize!(PIMPURITY.equiv, PIMPURITY.nsite)
            fill!(PIMPURITY.equiv, PIMPURITY.equiv[1])
            resize!(PIMPURITY.shell, PIMPURITY.nsite)
            fill!(PIMPURITY.shell, PIMPURITY.shell[1])
            resize!(PIMPURITY.ising, PIMPURITY.nsite)
            fill!(PIMPURITY.ising, PIMPURITY.ising[1])
            resize!(PIMPURITY.occup, PIMPURITY.nsite)
            fill!(PIMPURITY.occup, PIMPURITY.occup[1])
            resize!(PIMPURITY.upara, PIMPURITY.nsite)
            fill!(PIMPURITY.upara, PIMPURITY.upara[1])
            resize!(PIMPURITY.jpara, PIMPURITY.nsite)
            fill!(PIMPURITY.jpara, PIMPURITY.jpara[1])
            resize!(PIMPURITY.lpara, PIMPURITY.nsite)
            fill!(PIMPURITY.lpara, PIMPURITY.lpara[1])
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nsite)$(PIMPURITY.nsite)")

        # For the separator
        CImGui.Spacing()
        CImGui.Separator()
        CImGui.Spacing()

        # Input: sproj and window
        for i = 1:PIMPURITY.nsite
            if CImGui.CollapsingHeader("impurity $i")
                i == 1 && @widgets_generator_dft 1
                i == 2 && @widgets_generator_dft 2
                i == 3 && @widgets_generator_dft 3
                i == 4 && @widgets_generator_dft 4
                i == 5 && @widgets_generator_dft 5
                i == 6 && @widgets_generator_dft 6
                i == 7 && @widgets_generator_dft 7
                i == 8 && @widgets_generator_dft 8
                i == 9 && @widgets_generator_dft 9
            end
        end

        CImGui.EndTabItem()
    end
end

"""
    _zen_dmft_block()

Setup widgets for the `dmft` tab.
"""
function _zen_dmft_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("dmft")
        CImGui.Text("Configure [dmft] Block")

        # Input: mode
        CImGui.SetNextItemWidth(widget_combo_width)
        mode_list = ["one-shot", "self-consistent"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Scheme of dynamical mean-field theory calculations", &id, mode_list)
            PDMFT.mode = id + 1
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mode)$(PDMFT.mode)")
        #
        # Input: axis
        CImGui.SetNextItemWidth(widget_combo_width)
        axis_list = ["imaginary", "real"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Imaginary-time axis or real-frequency axis", &id, axis_list)
            PDMFT.axis = id + 1
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(axis)$(PDMFT.axis)")
        #
        # Input: niter
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(60) begin
            @c CImGui.SliderInt(" Maximum allowed number of DFT + DMFT iterations", &_i, 1, 100)
            PDMFT.niter = _i
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(niter)$(PDMFT.niter)")
        #
        # Input: nmesh
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(8193) begin
            @c CImGui.SliderInt(" Number of frequency points", &_i, 2^7+1, 2^14+1)
            PDMFT.nmesh = _i
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(nmesh)$(PDMFT.nmesh)")
        #
        # Input: dcount
        CImGui.SetNextItemWidth(widget_combo_width)
        dcount_list = ["fll1", "fll2", "amf", "held", "exact"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Scheme of double counting term", &id, dcount_list)
            PDMFT.dcount = dcount_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(dcount)$(PDMFT.dcount)")
        #
        # Input: beta
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(40.0) begin
            @c CImGui.InputDouble(" Inverse system temperature", &_f)
            PDMFT.beta = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(beta)$(PDMFT.beta)")
        #
        # Input: mixer
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.1) vmin = Cdouble(0.0) vmax = Cdouble(1.0) begin
            @c CImGui.SliderScalar(
                " Mixing factor",
                CImGui.ImGuiDataType_Double,
                &_f,
                &vmin, &vmax
            )
            PDMFT.mixer = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mixer)$(PDMFT.mixer)")
        #
        # Input: mc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0001) begin
            @c CImGui.InputDouble(" Convergence criterion of chemical potential", &_f)
            PDMFT.mc = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(mc)$(PDMFT.mc)")
        #
        # Input: cc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(1.0e-6) begin
            @c CImGui.InputDouble(" Convergence criterion of charge", &_f)
            PDMFT.cc = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(cc)$(PDMFT.cc)")
        #
        # Input: ec
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0001) begin
            @c CImGui.InputDouble(" Convergence criterion of total energy", &_f)
            PDMFT.ec = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ec)$(PDMFT.ec)")
        #
        # Input: sc
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _f = Cdouble(0.0001) begin
            @c CImGui.InputDouble(" Convergence criterion of self-energy function", &_f)
            PDMFT.sc = _f
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(sc)$(PDMFT.sc)")
        #
        # Input: lfermi
        CImGui.SetNextItemWidth(widget_combo_width)
        lfermi_list = ["yes", "no"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Whether chemical potential should be updated", &id, lfermi_list)
            if id == 0
                PDMFT.lfermi = true
            else
                PDMFT.lfermi = false
            end
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(lfermi)$(PDMFT.lfermi)")

        CImGui.EndTabItem()
    end
end

"""
    _zen_impurity_block()

Setup widgets for the `impurity` tab.
"""
function _zen_impurity_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("impurity")
        CImGui.Text("Configure [impurity] Block")
        #
        for i = 1:PIMPURITY.nsite
            if CImGui.CollapsingHeader("impurity $i")
                i == 1 && @widgets_generator_impurity 1
                i == 2 && @widgets_generator_impurity 2
                i == 3 && @widgets_generator_impurity 3
                i == 4 && @widgets_generator_impurity 4
                i == 5 && @widgets_generator_impurity 5
                i == 6 && @widgets_generator_impurity 6
                i == 7 && @widgets_generator_impurity 7
                i == 8 && @widgets_generator_impurity 8
                i == 9 && @widgets_generator_impurity 9
            end
        end
        #
        CImGui.EndTabItem()
    end
end

"""
    _zen_solver_block()

Setup widgets for the `solver` tab.
"""
function _zen_solver_block()
    # Widgets for the ctseg quantum impurity solver
    function _layout_ctseg()
        empty!(PSOLVER.params)

        # Input: isscr
        CImGui.SetNextItemWidth(widget_combo_width)
        isscr_list = ["static", "plasmon pole", "ohmic", "realistic"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" whether the Coulomb interaction U is dynamic", &id, isscr_list)
            push!(PSOLVER.params, "isscr = $(id + 1)")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isscr)$(last(PSOLVER.params))")
        #
        # Input: isort
        CImGui.SetNextItemWidth(widget_combo_width)
        isort_list = ["standard", "legendre", "intermediate"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Which basis will be used to do the measurement", &id, isort_list)
            push!(PSOLVER.params, "isort = $(id + 1)")
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(isort)$(last(PSOLVER.params))")
    end

    # Widgets for the cthyb quantum impurity solver
    function _layout_cthyb() end

    # Widgets for the hia quantum impurity solver
    function _layout_hia() end

    # Widgets for the norg quantum impurity solver
    function _layout_norg() end

    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    if CImGui.BeginTabItem("solver")
        CImGui.Text("Configure [solver] Block")

        # Input: engine
        CImGui.SetNextItemWidth(widget_combo_width)
        engine_list = ["ctseg", "cthyb", "hia", "norg"]
        @cstatic id = Cint(0) begin
            @c CImGui.Combo(" Name of quantum impurity solver", &id, engine_list)
            PSOLVER.engine = engine_list[id + 1]
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(engine)$(PSOLVER.engine)")
        #
        # Input: ncycle
        CImGui.SetNextItemWidth(widget_input_width)
        @cstatic _i = Cint(2) begin
            @c CImGui.InputInt(" Number of solver iterations per DFT + DMFT cycle", &_i)
            PSOLVER.ncycle = _i
        end
        CImGui.SameLine()
        CImGui.TextColored(COL_MAGENTA, "(ncycle)$(PSOLVER.ncycle)")

        # For the separator
        CImGui.Spacing()
        CImGui.Separator()
        CImGui.Spacing()

        # It should change upon the selection of quantum impurity solver.
        CImGui.Text("Quantum Impurity Solver: $(PSOLVER.engine)")

        @cswitch PSOLVER.engine begin

            @case "ctseg"
                _layout_ctseg()
                break

            @case "cthyb"
                _layout_cthyb()
                break

            @case "hia"
                _layout_hia()
                break

            @case "norg"
                _layout_norg()
                break

            @default
                break

        end

        CImGui.EndTabItem()
    end
end
