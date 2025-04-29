#
# Project : Lily
# Source  : peak.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/28
#

#=
*Remarks* :

Now the `ACTest` toolkit supports the following types of peaks.

**Gaussian peak**

```math
\begin{equation}
p(\omega) = A\exp{\left[-\frac{(\omega - ϵ)^2}{2\Gamma^2}\right]}
\end{equation}
```

**Lorentzian peak**

```math
\begin{equation}
p(\omega) = \frac{A}{\pi} \frac{\Gamma}{(\omega - \epsilon)^2 + \Gamma^2}
\end{equation}
```

**Delta-like peak**

```math
\begin{equation}
p(\omega) = A\exp{\left[-\frac{(\omega - ϵ)^2}{2\gamma^2}\right]},~
\text{where}~\gamma = 0.01
\end{equation}
```

**Rectangle peak**

```math
\begin{equation}
p(\omega) =
\begin{cases}
h, \quad \text{if}~\omega \in [c-w/2,c+w/2], \\
0, \quad \text{else}. \\
\end{cases}
\end{equation}
```

**Rise-And-Decay peak**

```math
\begin{equation}
p(\omega) = h \exp{(-|\omega - c|^{\gamma})}
\end{equation}
```

The spectral functions can be considered as a superposition of the above
peaks. That is to say: ``A(\omega) = \sum_i p_i(\omega)``.
=#

"""
    (𝑝::GaussianPeak)(ω::F64)

Evaluate the gaussian peak at ω.

### Arguments
* ω -> ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`GaussianPeak`](@ref).
"""
function (𝑝::GaussianPeak)(ω::F64)
    return 𝑝.A * exp( -(ω - 𝑝.ϵ) ^ 2.0 / (2.0 * 𝑝.Γ ^ 2.0) )
end

"""
    (𝑝::GaussianPeak)(ω::Vector{F64})

Evaluate the gaussian peak at real mesh.

### Arguments
* ω -> Real mesh, ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`GaussianPeak`](@ref).
"""
function (𝑝::GaussianPeak)(ω::Vector{F64})
    return @. 𝑝.A * exp( -(ω - 𝑝.ϵ) ^ 2.0 / (2.0 * 𝑝.Γ ^ 2.0) )
end

"""
    (𝑝::LorentzianPeak)(ω::F64)

Evaluate the lorentzian peak at ω.

### Arguments
* ω -> ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`LorentzianPeak`](@ref).
"""
function (𝑝::LorentzianPeak)(ω::F64)
    return 𝑝.A / π * 𝑝.Γ / ((ω - 𝑝.ϵ) ^ 2.0 + 𝑝.Γ ^ 2.0)
end

"""
    (𝑝::LorentzianPeak)(ω::Vector{F64})

Evaluate the lorentzian peak at real mesh.

### Arguments
* ω -> Real mesh, ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`LorentzianPeak`](@ref).
"""
function (𝑝::LorentzianPeak)(ω::Vector{F64})
    return @. 𝑝.A / π * 𝑝.Γ / ((ω - 𝑝.ϵ) ^ 2.0 + 𝑝.Γ ^ 2.0)
end

"""
    (𝑝::DeltaPeak)(ω::F64)

Evaluate the δ-like peak at ω.

### Arguments
* ω -> ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`DeltaPeak`](@ref).
"""
function (𝑝::DeltaPeak)(ω::F64)
    return 𝑝.A * exp( -(ω - 𝑝.ϵ) ^ 2.0 / (2.0 * 𝑝.Γ ^ 2.0) )
end

"""
    (𝑝::DeltaPeak)(ω::Vector{F64})

Evaluate the δ-like peak at real mesh.

### Arguments
* ω -> Real mesh, ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`DeltaPeak`](@ref).
"""
function (𝑝::DeltaPeak)(ω::Vector{F64})
    return @. 𝑝.A * exp( -(ω - 𝑝.ϵ) ^ 2.0 / (2.0 * 𝑝.Γ ^ 2.0) )
end

"""
    (𝑝::RectanglePeak)(ω::F64)

Evaluate the rectangle peak at ω.

### Arguments
* ω -> ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`RectanglePeak`](@ref).
"""
function (𝑝::RectanglePeak)(ω::F64)
    function f(x)
        𝑝.c - 𝑝.w / 2.0 ≤ x ≤ 𝑝.c + 𝑝.w / 2.0
    end
    return f(ω) ? 𝑝.h : zero(ω)
end

"""
    (𝑝::RectanglePeak)(ω::Vector{F64})

Evaluate the rectangle peak at real mesh.

### Arguments
* ω -> Real mesh, ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`RectanglePeak`](@ref).
"""
function (𝑝::RectanglePeak)(ω::Vector{F64})
    function f(x)
        𝑝.c - 𝑝.w / 2.0 ≤ x ≤ 𝑝.c + 𝑝.w / 2.0
    end
    return map(x -> f(x) ? 𝑝.h : zero(eltype(ω)), ω)
end

"""
    (𝑝::RiseDecayPeak)(ω::F64)

Evaluate the rise-and-decay peak at ω.

### Arguments
* ω -> ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`RiseDecayPeak`](@ref).
"""
function (𝑝::RiseDecayPeak)(ω::F64)
    return 𝑝.h * exp( - ( abs(ω - 𝑝.c) ) ^ 𝑝.γ )
end

"""
    (𝑝::RiseDecayPeak)(ω::Vector{F64})

Evaluate the rise-and-decay peak at real mesh.

### Arguments
* ω -> Real mesh, ω ∈ ℝ.

### Returns
* val -> 𝑝(ω).

See also: [`RiseDecayPeak`](@ref).
"""
function (𝑝::RiseDecayPeak)(ω::Vector{F64})
    return @. 𝑝.h * exp( - ( abs(ω - 𝑝.c) ) ^ 𝑝.γ )
end

# A general function to evaluate the peak at given mesh.
(𝑝::AbstractPeak)(ω::AbstractMesh) = 𝑝(ω.mesh)
