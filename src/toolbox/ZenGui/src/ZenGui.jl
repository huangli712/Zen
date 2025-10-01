#
# Project : Camellia
# Source  : ZenGui.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/06/27
#

"""
    ZenGui

ZenGui is a general-purpose graphic user interface for ab initio dynamical
mean-field theory codes, which have been developed at the China Academy of
Engineeing Physics. It should be used to generate necessary configuration
files for them. Now it supports the following codes:

* All-in-one DFT + DMFT package (**Zen**)
* Dynamical mean-field theory engines (**Dyson** and **DFermion**)
* Quantum impurity solver toolkit (**iQIST**)
* Analytic continuation tools (**ACFlow** and **ACTest**)

More codes will be supported in the future. Now this code is under heavy
development. **PLEASE USE IT AT YOUR OWN RISK**.

For more details about how to obtain, install and use this ZenGui app,
please visit the following website:

* `https://huangli712.github.io/projects/zengui/index.html`

Any suggestions, comments, and feedbacks are welcome. Enjoy it!
"""
module ZenGui

#=
### *Using Standard Libraries*
=#

using Random
using TOML
using Base.Sys

#=
### *Using Third-Party Libraries*
=#

#=
The **OrderedCollections** library provides support to the OrderedDict struct.
=#

using OrderedCollections

#=
The **FileIO** and **Images** packages are used to load and parse images.
=#

using FileIO
using Images

#=
The ZenGui application relies on the **Dear ImGui library**, which is a
C++ immediate mode graphic user interface library. Note that **CImGui**
is a C-API wrapper for **Dear ImGui**, and **CImGui.jl** provides a Julia
interface to **CImGui**.
=#

using CImGui
using CImGui.lib
using CImGui.CSyntax
using CImGui.CSyntax.CStatic

#=
**GLFW** is an open source multi-platform library for OpenGL, OpenGL ES
and Vulkan development on the desktop. It just provides a simple API for
creating windows, contexts and surfaces, receiving user's input and events.
**ModernGL** is a OpenGL binding for Julia.

Here, **GLFW** and **ModernGL** are backends for **Dear ImGui**.
=#

using GLFW
using ModernGL

#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the ZenGui app.

*Members* :

```text
I32, I64, API -> Numerical types (Integer).
F32, F64, APF -> Numerical types (Float).
C32, C64, APC -> Numerical types (Complex).
R32, R64, APR -> Numerical types (Union of Integer and Float).
N32, N64, APN -> Numerical types (Union of Integer, Float, and Complex).
#
__LIBNAME__   -> Name of this julia toolkit.
__VERSION__   -> Version of this julia toolkit.
__RELEASE__   -> Released date of this julia toolkit.
__AUTHORS__   -> Authors of this julia toolkit.
#
authors       -> Print the authors of ZenGui to screen.
```
=#

#
include("global.jl")
#
export I32, I64, API
export F32, F64, APF
export C32, C64, APC
export R32, R64, APR
export N32, N64, APN
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors

#=
### *Includes And Exports* : *util.jl*
=#

#=
*Summary* :

To provide some useful utility macros and functions.

*Members* :

```text
@cswitch     -> C-style switch.
#
sorry        -> Say sorry.
#
dict_to_toml -> Convert dictionary to toml file.
dict_to_ini  -> Convert dictionary to ini file.
#
open_url     -> Invoke the default web browser to open the given url.
```
=#

#
include("util.jl")
#
export @cswitch
#
export sorry
#
export dict_to_toml
export dict_to_ini
#
export open_url

#=
### *Includes And Exports* : *types.jl*
=#

#=
*Summary* :

Define some dicts and structs, which are used to store the configuration
parameters or represent some essential data structures.

*Members* :

```text
CURRENT_WINDOW  -> A struct used to keep the name of the current window.
CWIN            -> An instance for the `CURRENT_WINDOW` struct.
#
MenuFlags       -> A struct used to track the status of all the menu items.
FMENU           -> An instance for the `MenuFlags` struct.
#
ZEN_PCASE       -> It represents the `[case]` block in the `case.toml` file.
ZEN_PDFT        -> It represents the `[dft]` block in the `case.toml` file.
ZEN_PDMFT       -> It represents the `[dmft]` block in the `case.toml` file.
ZEN_PIMPURITY   -> It represents the `[impurity]` block in the `case.toml` file.
ZEN_PSOLVER     -> It represents the `[solver]` block in the `case.toml` file.
PCASE           -> An instance for the `ZEN_PCASE` struct.
PDFT            -> An instance for the `ZEN_PDFT` struct.
PDMFT           -> An instance for the `ZEN_PDMFT` struct.
PIMPURITY       -> An instance for the `ZEN_PIMPURITY` struct.
PSOLVER         -> An instance for the `ZEN_PSOLVER` struct.
#
DYSON_PDYSON    -> It contains the parameters in the `dmft.in` file.
_DYSON          -> It records the names of modified parameters for the Dyson code.
PDYSON          -> An instance for the `DYSON_PDYSON` struct.
#
DFERMION_PDFERMION -> It contains the parameters in the `dfa.in` file.
_DFERMION       -> It records the names of modified parameters for the DFermion code.
PDFERMION       -> An instance for the `DFERMION_PDFERMION` struct.
#
IQIST_PCTSEG    -> It contains the parameters in the `solver.ctqmc.in` file.
IQIST_PCTHYB    -> It contains the parameters in the `solver.ctqmc.in` file.
IQIST_PATOMIC   -> It contains the parameters in the `solver.atomic.in` file.
_CTSEG          -> It records the names of modified parameters for the ctseg code.
_CTHYB          -> It records the names of modified parameters for the cthyb code.
_ATOMIC         -> It records the names of modified parameters for the atomic code.
PCTSEG          -> An instance for the `IQIST_PCTSEG` struct.
PCTHYB          -> An instance for the `IQIST_PCTHYB` struct.
PATOMIC         -> An instance for the `IQIST_PATOMIC` struct.
#
ACFLOW_PBASE    -> It represents the `[BASE]` block in the `ac.toml` file.
ACFLOW_PMaxEnt  -> It represents the `[MaxEnt]` block in the `ac.toml` file.
ACFLOW_PBarRat  -> It represents the `[BarRat]` block in the `ac.toml` file.
ACFLOW_PNevanAC -> It represents the `[NevanAC]` block in the `ac.toml` file.
ACFLOW_PStochAC -> It represents the `[StochAC]` block in the `ac.toml` file.
ACFLOW_PStochSK -> It represents the `[StochSK]` block in the `ac.toml` file.
ACFLOW_PStochOM -> It represents the `[StochOM]` block in the `ac.toml` file.
ACFLOW_PStochPX -> It represents the `[StochPX]` block in the `ac.toml` file.
PBASE           -> An instance for the `ACFLOW_PBASE` struct.
PMaxEnt         -> An instance for the `ACFLOW_PMaxEnt` struct.
PBarRat         -> An instance for the `ACFLOW_PBarRat` struct.
PNevanAC        -> An instance for the `ACFLOW_PNevanAC` struct.
PStochAC        -> An instance for the `ACFLOW_PStochAC` struct.
PStochSK        -> An instance for the `ACFLOW_PStochSK` struct.
PStochOM        -> An instance for the `ACFLOW_PStochOM` struct.
PStochPX        -> An instance for the `ACFLOW_PStochPX` struct.
#
ACTEST_PTEST    -> It represents the `[Test]` block in the `act.toml` file.
PTEST           -> An instance for the `ACTEST_PTEST` struct.
#
struct_to_dict      -> Convert a struct to an ordered dictionary.
#
build_zen_dict      -> Assemble dictionary for the Zen package.
build_dyson_dict    -> Assemble dictionary for the Dyson code.
build_dfermion_dict -> Assemble dictionary for the DFermion code.
build_iqist_dict    -> Assemble dictionary for the iQIST package.
build_acflow_dict   -> Assemble dictionary for the ACFlow toolkit.
build_actest_dict   -> Assemble dictionary for the ACTest toolkit.
```
=#

#
include("types.jl")
#
export CURRENT_WINDOW
export CWIN
#
export MenuFlags
export FMENU
#
export ZEN_PCASE
export ZEN_PDFT
export ZEN_PDMFT
export ZEN_PIMPURITY
export ZEN_PSOLVER
export PCASE
export PDFT
export PDMFT
export PIMPURITY
export PSOLVER
#
export DYSON_PDYSON
export _DYSON
export PDYSON
#
export DFERMION_PDFERMION
export _DFERMION
export PDFERMION
#
export IQIST_PCTSEG
export IQIST_PCTHYB
export IQIST_PATOMIC
export _CTSEG
export _CTHYB
export _ATOMIC
export PCTSEG
export PCTHYB
export PATOMIC
#
export ACFLOW_PBASE
export ACFLOW_PMaxEnt
export ACFLOW_PBarRat
export ACFLOW_PNevanAC
export ACFLOW_PStochAC
export ACFLOW_PStochSK
export ACFLOW_PStochOM
export ACFLOW_PStochPX
export PBASE
export PMaxEnt
export PBarRat
export PNevanAC
export PStochAC
export PStochSK
export PStochOM
export PStochPX
#
export ACTEST_PTEST
export PTEST
#
export struct_to_dict
#
export build_zen_dict
export build_dyson_dict
export build_dfermion_dict
export build_iqist_dict
export build_acflow_dict
export build_actest_dict

#=
### *Includes And Exports* : *save.jl*
=#

#=
*Summary* :

Save configuration files for various quantum many-body codes.

*Members* :

```text
save_zen      -> Write case.toml for the Zen package.
save_dyson    -> Write dmft.in for the Dyson code.
save_dfermion -> Write dfa.in for the DFermion code.
save_ctseg    -> Write solver.ctqmc.in for the iQIST/ctseg code.
save_cthyb    -> Write solver.ctqmc.in for the iQIST/cthyb code.
save_atomic   -> Write solver.ctqmc.in for the iQIST/atomic code.
save_acflow   -> Write ac.toml for the ACFlow toolkit.
save_actest   -> Write act.toml for the ACTest toolkit.
save_nothing  -> Write nothing.
```
=#

#
include("save.jl")
#
export save_zen
export save_dyson
export save_dfermion
export save_ctseg
export save_cthyb
export save_atomic
export save_acflow
export save_actest
export save_nothing

#=
### *Includes And Exports* : *menu.jl*
=#

#=
*Summary* :

Setup global menu for the ZenGui app.

*Members* :

```text
create_menu    -> Create all the menu.
#
set_menu_file  -> Create the menu items in ``File''.
set_menu_edit  -> Create the menu items in ``Edit''.
set_menu_style -> Create the menu items in ``Style''.
set_menu_help  -> Create the menu items in ``Help''.
```
=#

#
include("menu.jl")
#
export create_menu
#
export set_menu_file
export set_menu_edit
export set_menu_style
export set_menu_help

#=
### *Includes And Exports* : *zen.jl*
=#

#=
*Summary* :

Create and display a window for configuring the Zen package.

*Members* :

```text
@widgets_generator_dft -> Macro for generating codes for the `dft` tab.
@widgets_generator_impurity -> Macro for generating codes for the `impurity` tab.
#
create_app_zen      -> Create and display the `Zen` window.
#
_zen_top_block      -> Setup widgets in the top of the window.
_zen_main_block     -> Setup widgets associated with the parameters.
_zen_bottom_block   -> Setup widgets in the bottom of the window.
#
_zen_case_block     -> Setup widgets for the `case` tab.
_zen_dft_block      -> Setup widgets for the `dft` tab.
_zen_dmft_block     -> Setup widgets for the `dmft` tab.
_zen_impurity_block -> Setup widgets for the `impurity` tab.
_zen_solver_block   -> Setup widgets for the `solver` tab.
```
=#

#
include("zen.jl")
#
export @widgets_generator_dft
export @widgets_generator_impurity
#
export create_app_zen
#
export _zen_top_block
export _zen_main_block
export _zen_bottom_block
#
export _zen_case_block
export _zen_dft_block
export _zen_dmft_block
export _zen_impurity_block
export _zen_solver_block

#=
### *Includes And Exports* : *dyson.jl*
=#

#=
*Summary* :

Create and display a window for configuring the Dyson code.

*Members* :

```text
create_app_dyson    -> Create and display the `Dyson` window.
#
_dyson_top_block    -> Setup widgets in the top of the window.
_dyson_main_block   -> Setup widgets associated with the parameters.
_dyson_bottom_block -> Setup widgets in the bottom of the window.
```
=#

#
include("dyson.jl")
#
export create_app_dyson
#
export _dyson_top_block
export _dyson_main_block
export _dyson_bottom_block

#=
### *Includes And Exports* : *dfermion.jl*
=#

#=
*Summary* :

Create and display a window for configuring the DFermion code.

*Members* :

```text
create_app_dfermion    -> Create and display the `DFermion` window.
#
_dfermion_top_block    -> Setup widgets in the top of the window.
_dfermion_main_block   -> Setup widgets associated with the parameters.
_dfermion_bottom_block -> Setup widgets in the bottom of the window.
#
_dfermion_model_block  -> Setup widgets for the `model` tab.
_dfermion_dimension_block -> Setup widgets for the `dimension` tab.
_dfermion_kmesh_block  -> Setup widgets for the `k-mesh` tab.
_dfermion_cycle_block  -> Setup widgets for the `cycle` tab.
```
=#

#
include("dfermion.jl")
#
export create_app_dfermion
#
export _dfermion_top_block
export _dfermion_main_block
export _dfermion_bottom_block
#
export _dfermion_model_block
export _dfermion_dimension_block
export _dfermion_kmesh_block
export _dfermion_cycle_block

#=
### *Includes And Exports* : *ctseg.jl*
=#

#=
*Summary* :

Create and display a window for configuring the iQIST/ctseg code.

*Members* :

```text
create_app_ctseg       -> Create and display the `ctseg` window.
#
_ctseg_top_block       -> Setup widgets in the top of the window.
_ctseg_main_block      -> Setup widgets associated with the parameters.
_ctseg_bottom_block    -> Setup widgets in the bottom of the window.
#
_ctseg_model_block     -> Setup widgets for the `model` tab.
_ctseg_dimension_block -> Setup widgets for the `dimension` tab.
_ctseg_symmetry_block  -> Setup widgets for the `symmetry` tab.
_ctseg_represent_block -> Setup widgets for the `representation` tab.
_ctseg_monte_block     -> Setup widgets for the `monte carlo` tab.
_ctseg_measure_block   -> Setup widgets for the `measurement` tab.
_ctseg_cycle_block     -> Setup widgets for the `cycle` tab.
```
=#

#
include("ctseg.jl")
#
export create_app_ctseg
#
export _ctseg_top_block
export _ctseg_main_block
export _ctseg_bottom_block
#
export _ctseg_model_block
export _ctseg_dimension_block
export _ctseg_symmetry_block
export _ctseg_represent_block
export _ctseg_monte_block
export _ctseg_measure_block
export _ctseg_cycle_block

#=
### *Includes And Exports* : *cthyb.jl*
=#

#=
*Summary* :

Create and display a window for configuring the iQIST/cthyb code.

*Members* :

```text
create_app_cthyb -> Create and display the `cthyb` window.
```
=#

#
include("cthyb.jl")
#
export create_app_cthyb

#=
### *Includes And Exports* : *atomic.jl*
=#

#=
*Summary* :

Create and display a window for configuring the iQIST/atomic code.

*Members* :

```text
create_app_atomic         -> Create and display the `atomic` window.
#
_atomic_top_block         -> Setup widgets in the top of the window.
_atomic_main_block        -> Setup widgets associated with the parameters.
_atomic_bottom_block      -> Setup widgets in the bottom of the window.
#
_atomic_model_block       -> Setup widgets for the `model` tab.
_atomic_interaction_block -> Setup widgets for the `interaction` tab.
_atomic_natural_block     -> Setup widgets for the `natural eigenbasis` tab.
_atomic_algorithm_block   -> Setup widgets for the `algorithm` tab.
```
=#

#
include("atomic.jl")
#
export create_app_atomic
#
export _atomic_top_block
export _atomic_main_block
export _atomic_bottom_block
#
export _atomic_model_block
export _atomic_interaction_block
export _atomic_natural_block
export _atomic_algorithm_block

#=
### *Includes And Exports* : *acflow.jl*
=#

#=
*Summary* :

Create and display a window for configuring the ACFlow toolkit.

*Members* :

```text
create_app_acflow     -> Create and display the `ACFlow` window.
#
_acflow_top_block     -> Setup widgets in the top of the window.
_acflow_main_block    -> Setup widgets associated with the parameters.
_acflow_bottom_block  -> Setup widgets in the bottom of the window.
#
_acflow_general_block -> Setup widgets for the `general` tab.
_acflow_solver_block  -> Setup widgets for the `solver` tab.
#
_acflow_maxent_block  -> Setup widgets for the [MaxEnt] block in the ac.toml file.
_acflow_barrat_block  -> Setup widgets for the [BarRat] block in the ac.toml file.
_acflow_nevanac_block -> Setup widgets for the [NevanAC] block in the ac.toml file.
_acflow_stochac_block -> Setup widgets for the [StochAC] block in the ac.toml file.
_acflow_stochsk_block -> Setup widgets for the [StochSK] block in the ac.toml file.
_acflow_stochom_block -> Setup widgets for the [StochOM] block in the ac.toml file.
_acflow_stochpx_block -> Setup widgets for the [StochPX] block in the ac.toml file.
```
=#

#
include("acflow.jl")
#
export create_app_acflow
#
export _acflow_top_block
export _acflow_main_block
export _acflow_bottom_block
#
export _acflow_general_block
export _acflow_solver_block
#
export _acflow_maxent_block
export _acflow_barrat_block
export _acflow_nevanac_block
export _acflow_stochac_block
export _acflow_stochsk_block
export _acflow_stochom_block
export _acflow_stochpx_block

#=
### *Includes And Exports* : *actest.jl*
=#

#=
*Summary* :

Create and display a window for configuring the ACTest toolkit.

*Members* :

```text
create_app_actest     -> Create and display the `ACTest` window.
#
_actest_top_block     -> Setup widgets in the top of the window.
_actest_main_block    -> Setup widgets associated with the parameters.
_actest_bottom_block  -> Setup widgets in the bottom of the window.
#
_actest_general_block -> Setup widgets for the `general` tab.
_actest_solver_block  -> Setup widgets for the `solver` tab.
```
=#

#
include("actest.jl")
#
export create_app_actest
#
export _actest_top_block
export _actest_main_block
export _actest_bottom_block
#
export _actest_general_block
export _actest_solver_block

#=
### *Includes And Exports* : *about.jl*
=#

#=
*Summary* :

Create and display the `About` window for the ZenGui app.

*Members* :

```text
create_app_about -> Create and display the `About` window.
```
=#

#
include("about.jl")
#
export create_app_about

#=
### *Includes And Exports* : *base.jl*
=#

#=
*Summary* :

To create main window and respond to users inputs for the ZenGui app.

*Members* :

```text
zeng_run         -> Main function for the ZenGui app.
#
load_texture     -> Load figures from the ZenGui/src/.images directory.
load_logo        -> Load logo image from the ZenGui/src/.images directory.
#
setup_flags      -> Setup configuration flags for the Dear ImGui library.
setup_fonts      -> Setup fonts for this graphic user interface.
setup_window     -> Tweak the window's style in this graphic user interface.
setup_background -> Setup the background figure for this app.
#
handle_menu_save       -> Respond the menu event: save.
handle_menu_background -> Respond the menu event: change background.
handle_menu_classic    -> Respond the menu event: classic.
handle_menu_dark       -> Respond the menu event: dark.
handle_menu_light      -> Respond the menu event: light.
handle_menu_zen        -> Respond the menu event: zen.
handle_menu_dyson      -> Respond the menu event: dyson.
handle_menu_dfermion   -> Respond the menu event: dfermion.
handle_menu_iqist      -> Respond the menu event: iqist.
handle_menu_acflow     -> Respond the menu event: acflow.
handle_menu_actest     -> Respond the menu event: actest.
handle_menu_zengui     -> Respond the menu event: zengui.
```
=#

#
include("base.jl")
#
export zeng_run
#
export load_texture
export load_logo
#
export setup_flags
export setup_fonts
export setup_window
export setup_background
#
export handle_menu_save
export handle_menu_background
export handle_menu_classic
export handle_menu_dark
export handle_menu_light
export handle_menu_zen
export handle_menu_dyson
export handle_menu_dfermion
export handle_menu_iqist
export handle_menu_acflow
export handle_menu_actest
export handle_menu_zengui

end
