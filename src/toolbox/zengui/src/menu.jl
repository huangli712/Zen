#
# Project : Camellia
# Source  : menu.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/25
#

"""
    create_menu()

Generate all menu items in the main window. Note that the `FMENU` struct
should be modified here.

See also: [`MenuFlags`](@ref) and [`FMENU`](@ref).
"""
function create_menu()
    if CImGui.BeginMainMenuBar()
        set_menu_file()  # For the ``File'' menu
        set_menu_edit()  # For the ``Edit'' menu
        set_menu_style() # For the ``Style'' menu
        set_menu_help()  # For the ``Help'' menu
        #
        CImGui.EndMainMenuBar()
    end
end

"""
    set_menu_file()

Setup menu items in ``File''. There are only two items: Save and Exit.
"""
function set_menu_file()
    if CImGui.BeginMenu("File")
        @c CImGui.MenuItem("Save", C_NULL, &FMENU.F_SAVE)
        #
        CImGui.Separator()
        #
        @c CImGui.MenuItem("Exit", C_NULL, &FMENU.F_EXIT)
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_edit()

Setup menu items in ``Edit''. They are related to the apps developed by
myself. Now only the Zen, Dyson, DFermion, iQIST, ACFlow, and ACTest codes
are supported.
"""
function set_menu_edit()
    if CImGui.BeginMenu("Edit")
        if CImGui.BeginMenu("Integrated Package")
            @c CImGui.MenuItem("Zen", C_NULL, &FMENU.E_ZEN)
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Quantum Many-Body Theory Engines")
            @c CImGui.MenuItem("Dyson", C_NULL, &FMENU.E_DYSON)
            @c CImGui.MenuItem("DFermion", C_NULL, &FMENU.E_DFERMION)
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Quantum Impurity Solvers")
            if CImGui.BeginMenu("iQIST")
                @c CImGui.MenuItem("ctseg", C_NULL, &FMENU.E_CTSEG)
                @c CImGui.MenuItem("cthyb", C_NULL, &FMENU.E_CTHYB)
                @c CImGui.MenuItem("atomic", C_NULL, &FMENU.E_ATOMIC)
                #
                CImGui.EndMenu()
            end
            #
            CImGui.EndMenu()
        end
        #
        if CImGui.BeginMenu("Analytic Continuation Tools")
            @c CImGui.MenuItem("ACFlow", C_NULL, &FMENU.E_ACFLOW)
            @c CImGui.MenuItem("ACTest", C_NULL, &FMENU.E_ACTEST)
            #
            CImGui.EndMenu()
        end
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_style()

Setup menu items in ``Style''. They are used to modify the appearance,
including background image and color styles, of this window-based app.
"""
function set_menu_style()
    if CImGui.BeginMenu("Style")
        @c CImGui.MenuItem("Change Background", C_NULL, &FMENU.S_BGIMAGE)
        #
        CImGui.Separator()
        #
        @c CImGui.MenuItem("Classic", C_NULL, &FMENU.S_CLASSIC)
        @c CImGui.MenuItem("Dark", C_NULL, &FMENU.S_DARK)
        @c CImGui.MenuItem("Light", C_NULL, &FMENU.S_LIGHT)
        #
        CImGui.EndMenu()
    end
end

"""
   set_menu_help()

Setup menu items in ``Help''. They are related to the documentation and
user guides of all the apps.
"""
function set_menu_help()
    if CImGui.BeginMenu("Help")
        if CImGui.BeginMenu("Documentation")
            @c CImGui.MenuItem("Zen", C_NULL, &FMENU.H_ZEN)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("Dyson", C_NULL, &FMENU.H_DYSON)
            @c CImGui.MenuItem("DFermion", C_NULL, &FMENU.H_DFERMION)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("iQIST", C_NULL, &FMENU.H_IQIST)
            #
            CImGui.Separator()
            #
            @c CImGui.MenuItem("ACFlow", C_NULL, &FMENU.H_ACFLOW)
            @c CImGui.MenuItem("ACTest", C_NULL, &FMENU.H_ACTEST)
            #
            CImGui.EndMenu()
        end
        #
        CImGui.Separator()
        #
        @c CImGui.MenuItem("User's Manual", C_NULL, &FMENU.H_ZENGUI)
        @c CImGui.MenuItem("About ZenGui", C_NULL, &FMENU.H_ABOUT)
        #
        CImGui.EndMenu()
    end
end
