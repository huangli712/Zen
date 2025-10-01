*Wrapper for various quantum impurity solvers. Now only the CT-HYB₁, CT-HYB₂, HIA, and NORG quantum impurity solvers are supported.*

*Source: solver.jl*

```@index
Pages = ["solver.md"]
```

## Functions

```@docs
solver_call
solver_copy
solver_sigma
solver_nimpx
solver_edmft
s_qmc1_init
s_qmc1_exec
s_qmc1_save
s_qmc1_copy
s_qmc2_init
s_qmc2_exec
s_qmc2_save
s_qmc2_copy
s_hub1_init
s_hub1_exec
s_hub1_save
s_hub1_copy
s_norg_init
s_norg_exec
s_norg_save
s_norg_copy
ctqmc_setup
ctqmc_atomx
ctqmc_delta
ctqmc_eimpx
ctqmc_sigma
ctqmc_nimpx
ctqmc_edmft
norg_setup
norg_delta
norg_eimpx
norg_sigma
norg_nimpx
norg_edmft
GetSigma
GetNimpx
GetEdmft
GetSymmetry
GetImpurity
CatImpurity
FixImpurity
```
