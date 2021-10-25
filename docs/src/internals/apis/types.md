# Types

*Define some dicts and structs, which are used to store the config parameters or represent some essential data structures.*

*Source: types.jl*

## Contents

```@contents
Pages = ["types.md"]
```

## Index

```@index
Pages = ["types.md"]
```

## Customized Types

```@docs
DType
```

## Global Dicts

```@docs
PCASE
PDFT
PDMFT
PIMP
PSOLVER
```

## Customized Structs: DFT Engine

```@docs
AbstractEngine
NULLEngine
VASPEngine
QEEngine
WANNIEREngine
_engine_
nameof(::NULLEngine)
```

## Customized Structs: Quantum Impurity Solver

```@docs
AbstractSolver
NULLSolver
CTHYB₁Solver
CTHYB₂Solver
HIASolver
NORGSolver
ATOMSolver
_solver_
nameof(::NULLSolver)
```

## Customized Structs: Adaptor

```@docs
AbstractAdaptor
NULLAdaptor
PLOAdaptor
WANNIERAdaptor
_adaptor_
nameof(::NULLAdaptor)
```

## Customized Structs: Sigma Engine

```@docs
AbstractMode
NULLMode
RESETMode
DCOUNTMode
SPLITMode
GATHERMode
_mode_
nameof(::NULLMode)
```

## Customized Structs: Mixer Engine

```@docs
AbstractMixer
NULLMixer
ΣMixer
ΔMixer
EMixer
ΓMixer
_mixer_
nameof(::NULLMixer)
```

## Structs

```@docs
Logger
Energy
IterInfo
Lattice
Mapping
Impurity
PrTrait
PrGroup
PrWindow
```

## Constructors

```@docs
Logger()
Energy()
IterInfo()
Lattice(::String, ::F64, ::I64, ::I64)
Mapping(::I64, ::I64, ::I64)
Impurity(::I64, ::String, ::I64, ::I64, ::String, ::String, ::F64, ::F64, ::F64, ::F64, ::F64)
PrTrait(::I64, ::String)
PrGroup(::I64, ::I64)
PrWindow(::Array{I64,3}, ::Tuple{R64,R64})
```

## Operators

```@docs
==
```

## Traits

```@docs
show(io::IO, logger::Logger)
show(io::IO, ene::Energy)
show(io::IO, it::IterInfo)
show(io::IO, latt::Lattice)
show(io::IO, map::Mapping)
show(io::IO, imp::Impurity)
show(io::IO, PT::PrTrait)
show(io::IO, PG::PrGroup)
show(io::IO, PW::PrWindow)
getproperty(et::Energy, sym::Symbol)
```
