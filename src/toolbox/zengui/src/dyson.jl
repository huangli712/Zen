#
# Project : Camellia
# Source  : dyson.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/26
#

"""
    create_app_dyson(p_open::Ref{Bool})

Create an UI window for the Dyson code, which is a dynamical mean-field
theory engine.
"""
function create_app_dyson(p_open::Ref{Bool})
    # Create the Dyson window, which can not be resized.
    CImGui.Begin(
        "Dyson",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "DYSON"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # For the widgets in the top of this window
    _dyson_top_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For all the blocks in the dmft.in file
    _dyson_main_block()

    # For the separator
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Spacing()

    # For the widgets in the bottom of this window
    _dyson_bottom_block(p_open)

    # End of this window
    CImGui.End()
end

"""
    _dyson_top_block()

Setup widgets in the top of the window for the Dyson code.
"""
function _dyson_top_block()
    CImGui.Text("Dyson: Dyson's Equation Solver For The Zen Computation Framework")
end

"""
    _dyson_main_block()

Setup widgets associated with the parameters in the `dmft.in` file.
"""
function _dyson_main_block()
    # Define default size for widgets
    widget_input_width = 100
    widget_combo_width = 100

    # Input: task
    CImGui.SetNextItemWidth(widget_combo_width)
    task_list = ["calc. hybridization",
                 "calc. density matrix",
                 "calc. fermi level",
                 "calc. impurity level",
                 "calc. eigenvalues",
                 "calc. spectrum",
                 "calc. density of states",
                "to be done"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Running mode of the code", &id, task_list)
        PDYSON.task = id + 1
        id != 0 && push!(_DYSON, "task")
        id == 0 && delete!(_DYSON, "task")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(task)$(PDYSON.task)")
    #
    # Input: axis
    CImGui.SetNextItemWidth(widget_combo_width)
    axis_list = ["imaginary axis", "real axis"]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Axis for Brillouin zone integration", &id, axis_list)
        PDYSON.axis = id + 1
        id != 0 && push!(_DYSON, "axis")
        id == 0 && delete!(_DYSON, "axis")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(axis)$(PDYSON.axis)")
    #
    # Input: beta
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(8.0) begin
        @c CImGui.InputDouble(" Inverse temperature", &_f)
        PDYSON.beta = _f
        _f != 8.0 && push!(_DYSON, "beta")
        _f == 8.0 && delete!(_DYSON, "beta")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(beta)$(PDYSON.beta)")
    #
    # Input: mc
    CImGui.SetNextItemWidth(widget_input_width)
    @cstatic _f = Cdouble(0.0001) begin
        @c CImGui.InputDouble(" Convergence criterion for fermi level search", &_f)
        PDYSON.mc = _f
        _f != 0.0001 && push!(_DYSON, "mc")
        _f == 0.0001 && delete!(_DYSON, "mc")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(mc)$(PDYSON.mc)")
    #
    # Input: lfermi
    CImGui.SetNextItemWidth(widget_combo_width)
    lfermi_list = [".true.", ".false."]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Whether the fermi level should be updated", &id, lfermi_list)
        PDYSON.lfermi = lfermi_list[id+1]
        id != 0 && push!(_DYSON, "lfermi")
        id == 0 && delete!(_DYSON, "lfermi")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(lfermi)$(PDYSON.lfermi)")
    #
    # Input: ltetra
    CImGui.SetNextItemWidth(widget_combo_width)
    ltetra_list = [".true.", ".false."]
    @cstatic id = Cint(0) begin
        @c CImGui.Combo(" Whether the analytical tetrahedron method is used", &id, ltetra_list)
        PDYSON.ltetra = ltetra_list[id+1]
        id != 0 && push!(_DYSON, "ltetra")
        id == 0 && delete!(_DYSON, "ltetra")
    end
    CImGui.SameLine()
    CImGui.TextColored(COL_MAGENTA, "(ltetra)$(PDYSON.ltetra)")
end

"""
    _dyson_bottom_block(p_open::Ref{Bool})

Setup widgets in the bottom of the window for the Dyson code.
"""
function _dyson_bottom_block(p_open::Ref{Bool})
    # Define default size for widgets
    widget_button_width = 80.0
    widget_button_height = 25.0

    # For the buttons
    if CImGui.Button("View", ImVec2(widget_button_width, widget_button_height))
        CImGui.OpenPopup("View dmft.in")
    end
    #
    if CImGui.BeginPopupModal("View dmft.in", C_NULL, CImGui.ImGuiWindowFlags_AlwaysAutoResize)
        @cstatic text="Hello World!" begin
            text = dict_to_ini(build_dyson_dict())
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
