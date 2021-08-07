# Types

## Contents

```@contents
Pages = ["types.md"]
```

## Index

```@index
Pages = ["types.md"]
```

## Global Dicts

```@docs
PCASE
PDFT
PDMFT
PIMP
PSOLVER
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
