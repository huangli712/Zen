# Base

*To provide the core functions to control the DFT engine, DMFT engine, quantum impurity solvers, Kohn-Sham adaptor, self-energy engine, and mixer engine. The DFT + DMFT iteration (one-shot mode or charge fully self-consistent mode) is also implemented in this file. This file also includes some functions to watch and manipulate the IterInfo struct.*

*Source: base.jl*

## Contents

```@contents
Pages = ["base.md"]
```

## Index

```@index
Pages = ["base.md"]
```

## Functions

```@docs
ready
go
final
cycle1
cycle2
try_dft
try_dmft
try_solver
try_adaptor
try_sigma
try_mixer
monitor
suspend
suicide
dft_core
dmft_core
solver_core
adaptor_core
sigma_core
mixer_core
energy_core
build_trees
clear_trees
incr_it
zero_it
prev_it
cntr_it
show_it
conv_it
```
