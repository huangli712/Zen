#
# Project : Camellia
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/01
#

function setup_config_flags()
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable
    io.IniFilename = C_NULL
end

function tweak_window()
    io = CImGui.GetIO()
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    style.AntiAliasedLines = true
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(style.Colors, CImGui.ImGuiCol_WindowBg, ImVec4(col.x, col.y, col.z, 1.0f0))
    end
end

function setup_fonts()
    fonts_dir = "/Users/lihuang/Library/Fonts"
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    CImGui.AddFontFromFileTTF(fonts, joinpath(fonts_dir, "FiraCode-Regular.ttf"), 16)
end
