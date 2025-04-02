#
# Project : Camellia
# Source  : menu.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/01
#

"""
    create_menu()

Generate menu items in the main window.
"""
function create_menu()
    if CImGui.BeginMainMenuBar()
        set_menu_file()
        set_menu_edit()
        set_menu_style()
        set_menu_help()
        #
        CImGui.EndMainMenuBar()
    end
end

"""
    set_menu_file()

Setup items in menu ``File''.
"""
function set_menu_file()
    if CImGui.BeginMenu("File")
        @c CImGui.MenuItem("Open...", C_NULL, &FMENU._OPEN)
        @c CImGui.MenuItem("Save As", C_NULL, &FMENU._SAVE)
        #
        CImGui.Separator()
        #
        @c CImGui.MenuItem("Exit", C_NULL, &FMENU._EXIT)
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_edit()

Setup items in menu ``Edit''.
"""
function set_menu_edit()
    if CImGui.BeginMenu("Edit")
        if CImGui.BeginMenu("Integrated Package")
            @c CImGui.MenuItem("Zen", C_NULL, &FMENU.ZEN)
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Dynamical Mean-Field Theory")
            @c CImGui.MenuItem("Dyson", C_NULL, &FMENU.DYSON)
            @c CImGui.MenuItem("DFermion", C_NULL, &FMENU.DFERMION)
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Quantum Impurity Solvers")
            @c CImGui.MenuItem("iQIST | CTSEG", C_NULL, &FMENU.CTSEG)
            @c CImGui.MenuItem("iQIST | CTHYB", C_NULL, &FMENU.CTHYB)
            @c CImGui.MenuItem("iQIST | ATOMIC", C_NULL, &FMENU.ATOMIC)
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Analytic Continuation")
            @c CImGui.MenuItem("ACFlow", C_NULL, &FMENU.ACFLOW)
            @c CImGui.MenuItem("ACTest", C_NULL, &FMENU.ACTEST)
            #
            CImGui.EndMenu()
        end
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_style()

Setup items in menu ``Style''.
"""
function set_menu_style()
    if CImGui.BeginMenu("Style")
        @c CImGui.MenuItem("Classic", C_NULL, &FMENU._CLASSIC)
        @c CImGui.MenuItem("Dark", C_NULL, &FMENU._DARK)
        @c CImGui.MenuItem("Light", C_NULL, &FMENU._LIGHT)
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_help()

Setup items in menu ``Help''.
"""
function set_menu_help()
    if CImGui.BeginMenu("Help")
        if CImGui.BeginMenu("Documentation")
            @c CImGui.MenuItem("Zen", C_NULL, &FMENU._ZEN)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("Dyson", C_NULL, &FMENU._DYSON)
            @c CImGui.MenuItem("DFermion", C_NULL, &FMENU._DFERMION)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("iQIST", C_NULL, &FMENU._IQIST)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("ACFlow", C_NULL, &FMENU._ACFLOW)
            @c CImGui.MenuItem("ACTest", C_NULL, &FMENU._ACTEST)
            #
            CImGui.EndMenu()
        end
        #
        CImGui.Separator()
        #
        @c CImGui.MenuItem("User's manual", C_NULL, &FMENU._ZENGUI)
        @c CImGui.MenuItem("About ZenGui", C_NULL, &FMENU._ABOUT)
        #
        CImGui.EndMenu()
    end
end
