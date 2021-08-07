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
cycle3
cycle4
cycle5
cycle6
cycle7
cycle8
monitor
suspend
suicide
dft_run
dmft_run
solver_run
adaptor_run
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
