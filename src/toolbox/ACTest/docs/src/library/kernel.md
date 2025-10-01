*Build kernel functions.*

The ACTest toolkit supports twelve types of kernel functions. They are:

* Fermionic imaginary time kernel (`ktype = "fermi", grid = "ftime"`)
* Fermionic Matsubara kernel (`ktype = "fermi", grid = "ffreq"`)
* Bosonic imaginary time kernel (`ktype = "boson", grid = "btime"`)
* Bosonc Matsubara kernel (`ktype = "boson", grid = "bfreq"`)
* Symmetric bosonic imaginary time kernel (`ktype = "bsymm", grid = "btime"`)
* Symmetric bosonic Matsubara kernel (`ktype = "bsymm", grid = "bfreq"`)

```@index
Pages = ["kernel.md"]
```

## Making Kernels

```@docs
build_kernel
build_kernel_symm
```
