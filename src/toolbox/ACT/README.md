# ACTest

The `ACTest` toolkit implements a spectral function / correlation function generator in Juila. It can generate multiple datasets for various spectral functions automatically. In addition, the corresponding imaginary time or Matsubara Green's functions (or the other correlation functions) with artifical Gaussian noises are also synthetized. It provides some useful scripts to perform extensive analytic continuation calculations and analyze the benchmark results. Now it is interfaced with the `ACFlow` and the `MiniPole` toolkits. But interfaced with the other analytic continuation methods or tools are also straightforward.

This toolkit is currently under developement. **PLEASE USE IT AT YOUR OWN RISK!**

## Version

v1.1.1-devel.250121

## License

GNU General Public License Version 3

## Installation

Please type the following commands in Julia's REPL to install the ACTest toolkit:

```julia-repl
julia> using Pkg
julia> Pkg.add("https://github.com/huangli712/ACTest")
```

To update the ACTest toolkit, please use the following commands:

```julia-repl
julia> using Pkg
julia> Pkg.update("ACTest")
```

## Documentation

See `actest/docs` or visit https://huangli712.github.io/projects/actest/index.html for the latest manual.
