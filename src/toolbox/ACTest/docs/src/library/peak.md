*Declare various peaks, which are used to build the spectral functions.*

Now the ACTest toolkit supports the following five types of peaks:

* Gaussian peak (`ptype = "gauss"`)
* Lorentzian peak (`ptype = "lorentz"`)
* ``\delta``-like peak (`ptype = "delta"`)
* Rectangle peak (`ptype = "rectangle"`)
* Rise-And-Decay peak (`ptype = "risedecay"`)

```@index
Pages = ["peak.md"]
```

## Types

```@docs
AbstractPeak
GaussianPeak
LorentzianPeak
DeltaPeak
RectanglePeak
RiseDecayPeak
```

## Base.* Functions

```@docs
Base.show(io::IO, 𝑝::GaussianPeak)
Base.show(io::IO, 𝑝::LorentzianPeak)
Base.show(io::IO, 𝑝::DeltaPeak)
Base.show(io::IO, 𝑝::RectanglePeak)
Base.show(io::IO, 𝑝::RiseDecayPeak)
```
