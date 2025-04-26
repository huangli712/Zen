#
# Project : Camellia
# Source  : cthyb.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/21
#

"""
    create_app_cthyb(p_open::Ref{Bool})

Create an UI window for the cthyb code, which is a continuous-time quantum
impurity solver in the iQIST package.
"""
function create_app_cthyb(p_open::Ref{Bool})
    # Create the cthyb window, which can not be resized.
    CImGui.Begin(
        "iQIST | CTHYB",
        p_open,
        CImGui.ImGuiWindowFlags_NoResize
    )

    # Setup the flag for active window
    if CImGui.IsWindowFocused()
        CWIN.name = "CTHYB"
    end

    # Fix size of the window
    window_width = 600.0
    window_height = 600.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    CImGui.Text("CTHYB")

    # End of this window
    CImGui.End()
end
