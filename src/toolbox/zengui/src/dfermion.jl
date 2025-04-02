#
# Project : Camellia
# Source  : dfermion.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/31
#

function create_app_dfermion(p_open::Ref{Bool})
    # Create the DFermion window, which is modal and can not be resized.
    CImGui.Begin(
        "DFermion",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    # Fix size of the window
    window_width = 400.0
    window_height = 300.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    CImGui.Text("DFermion")

    # End of this window
    CImGui.End()
end
