#
# Project : Camellia
# Source  : save.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/26
#

"""
    save_zen(p_open::Ref{Bool})

Save configuration file (`case.toml`) for the `Zen` package.
"""
function save_zen(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save Zen",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "case.toml")
    CImGui.Text("The configurtion file for Zen will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the case.toml file will be stored
    # in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_zen_dict()
        open("case.toml", "w") do fout
            TOML.print(fout, D)
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_dyson(p_open::Ref{Bool})

Save configuration file (`dmft.in`) for the `Dyson` code.
"""
function save_dyson(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save Dyson",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "dmft.in")
    CImGui.Text("The configurtion file for Dyson will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the dmft.in file will be stored
    # in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_dyson_dict()
        open("dmft.in", "w") do fout
            for (key, value) in D
                println(fout, "$key = $value")
            end
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_dfermion(p_open::Ref{Bool})

Save configuration file (`dfa.in`) for the `DFermion` code.
"""
function save_dfermion(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save DFermion",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "dfa.in")
    CImGui.Text("The configurtion file for DFermion will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the dfa.in file will be stored
    # in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_dfermion_dict()
        open("dfa.in", "w") do fout
            for (key, value) in D
                println(fout, "$key = $value")
            end
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_ctseg(p_open::Ref{Bool})

Save configuration file (`solver.ctqmc.in`) for the `iQIST/ctseg` code.
"""
function save_ctseg(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save iQIST | ctseg",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "solver.ctqmc.in")
    CImGui.Text("The configurtion file for iQIST/ctseg will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the solver.ctqmc.in file will be
    # stored in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_iqist_dict("ctseg")
        open("solver.ctqmc.in", "w") do fout
            for (key, value) in D
                println(fout, "$key = $value")
            end
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_cthyb(p_open::Ref{Bool})

Save configuration file (`solver.ctqmc.in`) for the `iQIST/cthyb` code.
"""
function save_cthyb(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save iQIST | cthyb",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "solver.ctqmc.in")
    CImGui.Text("The configurtion file for iQIST/cthyb will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the solver.ctqmc.in file will be
    # stored in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_iqist_dict("cthyb")
        open("solver.ctqmc.in", "w") do fout
            for (key, value) in D
                println(fout, "$key = $value")
            end
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_atomic(p_open::Ref{Bool})

Save configuration file (`solver.atomic.in`) for the `iQIST/atomic` code.
"""
function save_atomic(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save iQIST | atomic",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "solver.atomic.in")
    CImGui.Text("The configurtion file for iQIST/atomic will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the solver.atomic.in file will be
    # stored in the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_iqist_dict("atomic")
        open("solver.atomic.in", "w") do fout
            for (key, value) in D
                println(fout, "$key = $value")
            end
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_acflow(p_open::Ref{Bool})

Save configuration file (`ac.toml`) for the `ACFlow` toolkit.
"""
function save_acflow(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save ACFlow",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "ac.toml")
    CImGui.Text("The configurtion file for ACFlow will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the ac.toml file will be stored in
    # the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_acflow_dict()
        open("ac.toml", "w") do fout
            TOML.print(fout, D)
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_actest(p_open::Ref{Bool})

Save configuration file (`act.toml`) for the `ACTest` toolkit.
"""
function save_actest(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save ACTest",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    file = joinpath(pwd(), "act.toml")
    CImGui.Text("The configurtion file for ACTest will be saved at:")
    CImGui.TextColored(COL_MAGENTA, "  $file")

    # If the button is pressed, then the act.toml file will be stored in
    # the current directory.
    if CImGui.Button("Save It")
        p_open[] = false
        #
        D = build_actest_dict()
        open("act.toml", "w") do fout
            TOML.print(fout, D)
        end
    end

    # Close the popup window
    CImGui.End()
end

"""
    save_nothing(p_open::Ref{Bool})

This function will create an empty window and do nothing.
"""
function save_nothing(p_open::Ref{Bool})
    # Create a popup window
    CImGui.Begin(
        "Save Nothing",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    CImGui.TextColored(COL_MAGENTA, "Nothing to be saved!")

    # If the button is pressed, then close this window.
    if CImGui.Button("Close")
        p_open[] = false
    end

    # Close the popup window
    CImGui.End()
end
