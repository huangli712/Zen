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
_engine_
```

## Customized Structs: Quantum Impurity Solver

```@docs
_solver_
```

## Customized Structs: Adaptor

```@docs
_adaptor_
```

## Customized Structs: Sigma Engine

```@docs
_mode_
```

## Customized Structs: Mixer Engine

```@docs
_mixer_
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
