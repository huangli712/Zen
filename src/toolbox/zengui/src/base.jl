#
# Project : Camellia
# Source  : base.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/25
#

# Basic algebraic operation for ImVec2.
Base.:+(v1::ImVec2, v2::ImVec2) = ImVec2(v1.x + v2.x, v1.y + v2.y)

#=
### *Main Loop*
=#

"""
    zeng_run()

Main function. It launchs the graphic user interface and respond to user
inputs unitl the main window is closed.
"""
function zeng_run()
    # Setup backend for Dear ImGui. Now CImGui.jl only supports opengl.
    CImGui.set_backend(:GlfwOpenGL3)

    # Setup context for Dear ImGui
    ctx = CImGui.CreateContext()

    # Setup color's style for Dear ImGui
    CImGui.StyleColorsDark()

    # Setup flags for Dear ImGui, enabling docking and multi-viewport.
    setup_flags()

    # Setup window's style
    #
    # When viewports are enabled, we tweak WindowRounding and WindowBg so
    # platform windows can look identical to regular ones.
    setup_window()

    # Load special fonts if available
    setup_fonts()

    # Global id for texture.
    # When texture_id is nothing, it means that the texture has not been loaded.
    texture_id = nothing
    CImGui.render(ctx; window_title = "ZenGui") do
        # If texture_id is nothing, we should try to load the image and
        # setup the texture's id.
        if isnothing(texture_id)
            texture_id = load_texture()
        end

        # Setup global menu in the main window
        create_menu()

        # Setup the background images for the main viewport
        setup_background(texture_id)

        # Respond to menu events
        #
        # For File menu
        FMENU.F_SAVE     && @c handle_menu_save(&FMENU.F_SAVE)
        FMENU.F_EXIT     && return :imgui_exit_loop
        #
        # For Edit menu
        FMENU.E_ZEN      && @c create_app_zen(&FMENU.E_ZEN)
        FMENU.E_DYSON    && @c create_app_dyson(&FMENU.E_DYSON)
        FMENU.E_DFERMION && @c create_app_dfermion(&FMENU.E_DFERMION)
        FMENU.E_CTSEG    && @c create_app_ctseg(&FMENU.E_CTSEG)
        FMENU.E_CTHYB    && @c create_app_cthyb(&FMENU.E_CTHYB)
        FMENU.E_ATOMIC   && @c create_app_atomic(&FMENU.E_ATOMIC)
        FMENU.E_ACFLOW   && @c create_app_acflow(&FMENU.E_ACFLOW)
        FMENU.E_ACTEST   && @c create_app_actest(&FMENU.E_ACTEST)
        #
        # For Style menu
        # Once the menu `Change Background` is clicked, texture_id must be
        # reset to nothing. Then the load_texture() function is called
        # again to load a new image and reassign the texture_id.
        FMENU.S_BGIMAGE  && (texture_id = nothing)
        FMENU.S_BGIMAGE  && @c handle_menu_background(&FMENU.S_BGIMAGE)
        FMENU.S_CLASSIC  && @c handle_menu_classic(&FMENU.S_CLASSIC)
        FMENU.S_DARK     && @c handle_menu_dark(&FMENU.S_DARK)
        FMENU.S_LIGHT    && @c handle_menu_light(&FMENU.S_LIGHT)
        #
        # For Help menu
        FMENU.H_ZEN      && @c handle_menu_zen(&FMENU.H_ZEN)
        FMENU.H_DYSON    && @c handle_menu_dyson(&FMENU.H_DYSON)
        FMENU.H_DFERMION && @c handle_menu_dfermion(&FMENU.H_DFERMION)
        FMENU.H_IQIST    && @c handle_menu_iqist(&FMENU.H_IQIST)
        FMENU.H_ACFLOW   && @c handle_menu_acflow(&FMENU.H_ACFLOW)
        FMENU.H_ACTEST   && @c handle_menu_actest(&FMENU.H_ZENGUI)
        FMENU.H_ZENGUI   && @c handle_menu_zengui(&FMENU.H_ZENGUI)
        FMENU.H_ABOUT    && @c create_app_about(&FMENU.H_ABOUT)
    end
end

#=
### *Configure Application*
=#

"""
    load_texture()

Load images from the ZenGui/src/.images directory. Note that there are
12 images now. This function will pick one image randomly and load it.
Finally, it will return an `ImTextureID` object which is associted with
the selected image.

See also: [`setup_background`](@ref).
"""
function load_texture()
    # Prepare images
    #
    # Setup image list
    img_list = ["bg1.jpg", "bg2.jpg", "bg3.jpg",
                "bg4.jpg", "bg5.png", "bg6.jpg",
                "bg7.jpg", "bg8.jpg", "bg9.jpg",
                "bg10.jpg",
                "bg11.jpg",
                "bg12.jpg"]
    #
    # Select one image randomly.
    img_indx = rand(MersenneTwister(), 1:length(img_list))
    #
    # Setup directory and path for the selected image
    img_dir = joinpath(ENV["ZEN_GUI"], ".images")
    img_path = joinpath(img_dir, img_list[img_indx])

    # Load the selected image by the FileIO and Images packages.
    img = RGBA.(rotr90(FileIO.load(img_path)))
    width, height = size(img)

    # Setup texture and the related properties by opengl
    texture_id = GLuint(0)
    @c glGenTextures(1, &texture_id)
    glBindTexture(GL_TEXTURE_2D, texture_id)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, reinterpret(UInt8, img))

    # Return an ImTextureID object
    return CImGui.ImTextureID(texture_id)
end

"""
    setup_flags()

Setup configuration flags for the Dear ImGui library.
"""
function setup_flags()
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_ViewportsEnable
    io.IniFilename = C_NULL # We never generate the ini file for this app.
end

"""
    setup_fonts()

Setup fonts for this graphic user interface. Note that the files for fonts
should be saved in the ZenGui/src/.fonts directory.
"""
function setup_fonts()
    fonts_dir = joinpath(ENV["ZEN_GUI"], ".fonts")
    fonts = unsafe_load(CImGui.GetIO().Fonts)
    CImGui.AddFontFromFileTTF(
        fonts,
        joinpath(fonts_dir, "FiraCode-Regular.ttf"),
        16,
        C_NULL,
        CImGui.GetGlyphRangesGreek(fonts) # To display the Greek letters
    )
end

"""
    setup_window()

Tweak the window's style in this graphic user interface.
"""
function setup_window()
    style = Ptr{ImGuiStyle}(CImGui.GetStyle())
    style.AntiAliasedLines = true
    #
    # Setup background color for window
    io = CImGui.GetIO()
    if unsafe_load(io.ConfigFlags) & ImGuiConfigFlags_ViewportsEnable == ImGuiConfigFlags_ViewportsEnable
        style.WindowRounding = 5.0f0
        col = CImGui.c_get(style.Colors, CImGui.ImGuiCol_WindowBg)
        CImGui.c_set!(
            style.Colors,
            CImGui.ImGuiCol_WindowBg,
            ImVec4(col.x, col.y, col.z, 1.0f0)
        )
    end
    #
    # Setup default colors for buttons
    CImGui.c_set!(style.Colors, CImGui.ImGuiCol_Button, COL_ORANGE)
    CImGui.c_set!(style.Colors, CImGui.ImGuiCol_ButtonHovered, COL_LIGHTGREEN)
    #
    # Setup default colors for sliders
    CImGui.c_set!(style.Colors, CImGui.ImGuiCol_SliderGrab, COL_PINK)
    CImGui.c_set!(style.Colors, CImGui.ImGuiCol_SliderGrabActive, COL_PINK)
end

"""
    setup_background(texture_id)

Setup the background image for this app. `texture_id` is an `ImTextureID`
object. which is initialized by `load_texture()`.

See also: [`load_texture`](@ref).
"""
function setup_background(texture_id)
    # Get default viewpoint, determine its size and position.
    viewport = unsafe_load(CImGui.GetMainViewport())
    pos = viewport.Pos
    size = viewport.Size

    # Get draw list for background
    #
    # The selected image will be drawn on it directly without creating a
    # new window.
    drawlist = CImGui.GetBackgroundDrawList()

    # Draw the image
    CImGui.AddImage(drawlist, texture_id, pos, pos + size, (0.0, 1.0), (1.0, 0.0))
end

#=
### *Menu Handler*
=#

"""
    handle_menu_save(p_open::Ref{Bool})

Respond the menu event: save. Try to save configurtion files for various
tools or codes.
"""
function handle_menu_save(p_open::Ref{Bool})
    @cswitch CWIN.name begin

        @case "ZEN"
            if FMENU.E_ZEN
                save_zen(p_open)
            else
                p_open[] = false
            end
            break

        @case "DYSON"
            if FMENU.E_DYSON
                save_dyson(p_open)
            else
                p_open[] = false
            end
            break

        @case "DFERMION"
            if FMENU.E_DFERMION
                save_dfermion(p_open)
            else
                p_open[] = false
            end
            break

        @case "CTSEG"
            if FMENU.E_CTSEG
                save_ctseg(p_open)
            else
                p_open[] = false
            end
            break

        @case "CTHYB"
            if FMENU.E_CTHYB
                save_cthyb(p_open)
            else
                p_open[] = false
            end
            break

        @case "ATOMIC"
            if FMENU.E_ATOMIC
                save_atomic(p_open)
            else
                p_open[] = false
            end
            break

        @case "ACFLOW"
            if FMENU.E_ACFLOW
                save_acflow(p_open)
            else
                p_open[] = false
            end
            break

        @case "ACTEST"
            if FMENU.E_ACTEST
                save_actest(p_open)
            else
                p_open[] = false
            end
            break

        @default
            save_nothing(p_open)
            break

    end
end

"""
    handle_menu_background(p_open::Ref{Bool})

Respond the menu event: change background. It will change the background
image randomly.

See also: [`setup_background`](@ref).
"""
function handle_menu_background(p_open::Ref{Bool})
    p_open[] = false
end

"""
    handle_menu_classic(p_open::Ref{Bool})

Respond the menu event: classic. Change the appearance of graphic user
interface to classic style.

See also: [`setup_window`](@ref).
"""
function handle_menu_classic(p_open::Ref{Bool})
    CImGui.StyleColorsClassic()
    setup_window() # We should reset color style for the buttons
    p_open[] = false
end

"""
    handle_menu_dark(p_open::Ref{Bool})

Respond the menu event: dark. Change the appearance of graphic user
interface to dark style. Note that the defalt style is dark.

See also: [`setup_window`](@ref).
"""
function handle_menu_dark(p_open::Ref{Bool})
    CImGui.StyleColorsDark()
    setup_window() # We should reset color style for the buttons
    p_open[] = false
end

"""
    handle_menu_light(p_open::Ref{Bool})

Respond the menu event: light. Change the appearance of graphic user
interface to light style.

See also: [`setup_window`](@ref).
"""
function handle_menu_light(p_open::Ref{Bool})
    CImGui.StyleColorsLight()
    setup_window() # We should reset color style for the buttons
    p_open[] = false
end

"""
    handle_menu_zen(p_open::Ref{Bool})

Respond the menu event: zen. Try to open documentation for the Zen package.
"""
function handle_menu_zen(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/zen/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_dyson(p_open::Ref{Bool})

Respond the menu event: dyson. Try to open documentation for the Dyson
code.
"""
function handle_menu_dyson(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/dyson/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_dfermion(p_open::Ref{Bool})

Respond the menu event: dfermion. Try to open documentation for the
DFermion code.
"""
function handle_menu_dfermion(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/dfermion/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_iqist(p_open::Ref{Bool})

Respond the menu event: iqist. Try to open documentation for the iQIST
package.
"""
function handle_menu_iqist(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/iqist_new/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_acflow(p_open::Ref{Bool})

Respond the menu event: acflow. Try to open documentation for the ACFlow
toolkit.
"""
function handle_menu_acflow(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/acflow/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_actest(p_open::Ref{Bool})

Respond the menu event: actest. Try to open documentation for the ACTest
toolkit.
"""
function handle_menu_actest(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/actest/index.html"
    open_url(url)
    p_open[] = false
end

"""
    handle_menu_zengui(p_open::Ref{Bool})

Respond the menu event: zengui. Try to open documentation for the ZenGui
application.
"""
function handle_menu_zengui(p_open::Ref{Bool})
    url = "https://huangli712.github.io/projects/zengui/index.html"
    open_url(url)
    p_open[] = false
end
