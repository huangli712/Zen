## Index

```@index
Pages = ["library.md"]
```

## Modules

```@docs
ZenGui
```

## Numerical Types

```@docs
I32
I64
API
F32
F64
APF
C32
C64
APC
R32
R64
APR
N32
N64
APN
```

## String Constants

```@docs
__LIBNAME__
__VERSION__
__RELEASE__
__AUTHORS__
authors
```

## General Macros

```@docs
@cswitch
```

## General Utilities

```@docs
sorry
open_url
```

## Structs: Menu

```@docs
CURRENT_WINDOW
CWIN
MenuFlags
FMENU
```

## Widgets: Menu

```@docs
create_menu
set_menu_file
set_menu_edit
set_menu_style
set_menu_help
handle_menu_save
handle_menu_background
handle_menu_dark
handle_menu_classic
handle_menu_light
handle_menu_zen
handle_menu_dyson
handle_menu_dfermion
handle_menu_iqist
handle_menu_acflow
handle_menu_actest
handle_menu_zengui
```

## Structs: Zen

```@docs
ZEN_PCASE
ZEN_PDFT
ZEN_PDMFT
ZEN_PIMPURITY
ZEN_PSOLVER
PCASE
PDFT
PDMFT
PIMPURITY
PSOLVER
```

## Widgets: Zen

```@docs
@widgets_generator_dft
@widgets_generator_impurity
create_app_zen
_zen_top_block
_zen_main_block
_zen_bottom_block
_zen_case_block
_zen_dft_block
_zen_dmft_block
_zen_impurity_block
_zen_solver_block
build_zen_dict
```

## Structs: Dyson

```@docs
DYSON_PDYSON
_DYSON
PDYSON
```

## Widgets: Dyson

```@docs
create_app_dyson
_dyson_top_block
_dyson_main_block
_dyson_bottom_block
build_dyson_dict
```

## Structs: DFermion

```@docs
DFERMION_PDFERMION
_DFERMION
PDFERMION
```

## Widgets: DFermion

```@docs
create_app_dfermion
_dfermion_top_block
_dfermion_main_block
_dfermion_bottom_block
_dfermion_model_block
_dfermion_dimension_block
_dfermion_kmesh_block
_dfermion_cycle_block
build_dfermion_dict
```

## Structs: iQIST

```@docs
IQIST_PCTSEG
IQIST_PCTHYB
IQIST_PATOMIC
_CTSEG
_CTHYB
_ATOMIC
PCTSEG
PCTHYB
PATOMIC
```

## Widgets: iQIST

```@docs
create_app_ctseg
create_app_cthyb
create_app_atomic
_ctseg_top_block
_ctseg_main_block
_ctseg_bottom_block
_ctseg_model_block
_ctseg_dimension_block
_ctseg_symmetry_block
_ctseg_represent_block
_ctseg_measure_block
_ctseg_monte_block
_ctseg_cycle_block
_atomic_top_block
_atomic_main_block
_atomic_bottom_block
_atomic_model_block
_atomic_interaction_block
_atomic_natural_block
_atomic_algorithm_block
build_iqist_dict
```

## Structs: ACFlow

```@docs
ACFLOW_PBASE
ACFLOW_PMaxEnt
ACFLOW_PBarRat
ACFLOW_PNevanAC
ACFLOW_PStochAC
ACFLOW_PStochSK
ACFLOW_PStochOM
ACFLOW_PStochPX
PBASE
PMaxEnt
PBarRat
PNevanAC
PStochAC
PStochSK
PStochOM
PStochPX
```

## Widgets: ACFlow

```@docs
create_app_acflow
_acflow_top_block
_acflow_main_block
_acflow_bottom_block
_acflow_general_block
_acflow_solver_block
_acflow_maxent_block
_acflow_barrat_block
_acflow_nevanac_block
_acflow_stochac_block
_acflow_stochsk_block
_acflow_stochom_block
_acflow_stochpx_block
build_acflow_dict
```

## Structs: ACTest

```@docs
ACTEST_PTEST
PTEST
```

## Widgets: ACTest

```@docs
create_app_actest
_actest_top_block
_actest_main_block
_actest_bottom_block
_actest_general_block
_actest_solver_block
build_actest_dict
```

## Widgets: About

```@docs
create_app_about
```

## Functions: Utilities

```@docs
dict_to_toml
dict_to_ini
struct_to_dict
```

## Functions: Setup

```@docs
load_texture
load_logo
setup_flags
setup_fonts
setup_window
setup_background
```

## Functions: Save

```@docs
save_zen
save_dyson
save_dfermion
save_ctseg
save_cthyb
save_atomic
save_acflow
save_actest
save_nothing
```

## Functions: Main

```@docs
zeng_run
```
