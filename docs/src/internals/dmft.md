# Dynamical Mean-Field Theory

*Wrapper for dynamical mean-field theory engine. It also provides some essential tools to deal with (read and write) the fermi level ``\mu``, hybridization functions ``\Delta``, local impurity levels ``\epsilon_i``, and correlation-induced correction for density matrix ``\Gamma``.*

*Source: dmft.jl*

## Contents

```@contents
Pages = ["dmft.md"]
```

## Index

```@index
Pages = ["dmft.md"]
```

## Functions

```@docs
dmft_init
dmft_exec
dmft_save
read_fermi
read_delta
read_eimpx
read_gcorr
write_delta
write_eimpx
write_gcorr
```
