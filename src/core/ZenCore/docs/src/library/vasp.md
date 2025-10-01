*Tools for the vasp software package (adaptor). It provide a lot of functions to deal with the vasp-related files.*

*Source: vasp.jl*

```@index
Pages = ["vasp.md"]
```

## Functions

```@docs
dft_call(::VASPEngine, ::IterInfo)
dft_stop(::VASPEngine)
dft_resume(::VASPEngine)
adaptor_call(::VASPEngine, ::Dict{Symbol,Any})
vasp_adaptor
vasp_init
vasp_exec
vasp_save
vasp_back
vasp_stop
vaspc_incar
vaspc_kpoints
vaspc_gcorr
vaspc_stopcar
vaspc_lock
vaspq_stopcar
vaspq_lock
vaspq_files
vaspio_nband
vaspio_valence
vaspio_energy
vaspio_procar
vaspio_lattice
vaspio_kmesh
vaspio_tetra
vaspio_eigen
vaspio_projs
vaspio_fermi
vaspio_charge
```
