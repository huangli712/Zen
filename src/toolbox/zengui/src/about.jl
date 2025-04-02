#
# Project : Camellia
# Source  : about.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/31
#

"""
    create_app_about(p_open::Ref{Bool})

Display the `About` window, which is used to show some userful information
for users.
"""
function create_app_about(p_open::Ref{Bool})
    # Create the about window, which is modal and can not be resized.
    CImGui.Begin(
        "About ZenGui",
        p_open,
        CImGui.ImGuiWindowFlags_Modal | CImGui.ImGuiWindowFlags_NoResize
    )

    # Fix size of the window
    window_width = 400.0
    window_height = 300.0
    CImGui.SetWindowSize(ImVec2(window_width, window_height))

    # Below are some widges

    # We want to make sure `ZenGui` is shown in the middle of the window.
    txt_width = CImGui.CalcTextSize("ZenGui").x
    offset = (window_width - txt_width) / 2.0
    CImGui.SameLine(offset)
    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "ZenGui")
    #
    CImGui.Spacing()
    CImGui.TextWrapped("A general-purposed graphic user interface for " *
        "ab initio dynamical mean-field theory codes")
    #
    CImGui.Spacing()
    CImGui.Separator()
    CImGui.Dummy(ImVec2(0.0,10.0))

    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "Author :")
    CImGui.SameLine()
    CImGui.Text("Li Huang")
    #
    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "Contact:")
    CImGui.SameLine()
    CImGui.Text("huangli at caep.cn")
    #
    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "Version:")
    CImGui.SameLine()
    CImGui.Text("v0.2.0-devel.250331")
    #
    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "License:")
    CImGui.SameLine()
    CImGui.Text("GNU General Public License Version 3")
    #
    CImGui.TextColored(ImVec4(1.0,0.0,1.0,1.0), "Github :")
    CImGui.SameLine()
    CImGui.Text("https://github.com/huangli712/ZenGui")

    CImGui.Spacing()
    CImGui.TextWrapped("Powered by the Julia language (v$VERSION) " *
        "and the Dear ImGui library (v$(CImGui.IMGUI_VERSION)).")

    # Create a `OK` button. It will reset `p_open`.
    CImGui.Dummy(ImVec2(0.0,10.0))
    #
    # Change the default color for the button
    CImGui.PushStyleColor(CImGui.ImGuiCol_Button, ImVec4(1.0,0.0,1.0,1.0))
    #
    button_width = 80.0
    button_height = 25.0
    if CImGui.Button("OK", ImVec2(button_width, button_height))
        p_open[] = false
    end
    #
    # Reset the color
    CImGui.PopStyleColor(1)

    # End of this window
    CImGui.End()
end
