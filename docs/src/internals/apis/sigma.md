# Self-energy functions

*Tools for treating the self-energy functions ``\Sigma``, double counting terms ``\Sigma_{\text{dc}}``. Note that the function `sigma_split()` is designed for the hybridization functions ``\Delta`` and local impurity levels ``\epsilon_i``, instead of the self-energy functions.*

*Source: sigma.jl*

## Contents

```@contents
Pages = ["sigma.md"]
```

## Index

```@index
Pages = ["sigma.md"]
```

## Functions

```@docs
sigma_reset
sigma_dcount
sigma_split
sigma_gather
cal_dc_fll
cal_dc_amf
cal_dc_held
cal_dc_exact
read_sigma
read_sigdc
write_sigma
write_sigdc
```
