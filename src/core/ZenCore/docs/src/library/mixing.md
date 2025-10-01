*Tools for mixing the self-energy functions ``\Sigma``, hybridization functions ``\Delta``, and local impurity levels ``\epsilon_i``. They adopted the linear mixing algorithm. We also implement the so-called Kerker algorithm to mix the correlation-induced correction for density matrix ``\Gamma``.*

*Source: mixer.jl*

```@index
Pages = ["mixer.md"]
```

## Functions

```@docs
mixer_call
mixer_sigma
mixer_delta
mixer_eimpx
mixer_gcorr
amix
distance
```
