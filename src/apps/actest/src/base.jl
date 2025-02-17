#
# Project : Lily
# Source  : base.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/11/14
#

"""
    reprod(am::AbstractMesh, kernel::Matrix{F64}, A::Vector{F64})

Try to reproduce the input data, which can be compared with the raw data
to see whether the analytic continuation is reasonable.

### Arguments
* am -> Real frequency mesh.
* kernel -> The kernel function.
* A -> The calculated spectral function, A(ω).

### Returns
* G -> Reconstructed correlators, G(τ) or G(iωₙ), Vector{F64}.

See also: [`AbstractMesh`](@ref).
"""
function reprod(am::AbstractMesh, kernel::Matrix{F64}, A::Vector{F64})
    ndim, nmesh = size(kernel)
    @assert nmesh == length(am) == length(A)

    @einsum KA[i,j] := kernel[i,j] * A[j]

    G = zeros(F64, ndim)
    for i = 1:ndim
        G[i] = trapz(am, view(KA, i, :))
    end

    return G
end

"""
    setup_param(C::Dict{String,Any}, reset::Bool = true)

Setup the configuration dictionaries via function call. Here `C` contains
parameters for general setup. If `reset` is true, then the configuration
dictionaries will be reset to their default values at first. Later, `C`
will be used to customized the dictionaries further.

### Arguments
See above explanations.

### Returns
N/A

See also: [`read_param`](@ref).
"""
function setup_param(C::Dict{String,Any}, reset::Bool = true)
    # _PTEST contains the default parameters.
    # If reset is true, it will be used to update the PTEST dictionary.
    reset && rev_dict(_PTEST)
    rev_dict(C)
end

"""
    read_param()

Setup the configuration dictionaries via an external file. The valid
format of a configuration file is `toml`.

### Arguments
N/A

### Returns
N/A

See also: [`setup_param`](@ref).
"""
function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
    chk_dict()
    see_dict()
end

"""
    make_data_std()

Try to build a standard dataset (`ACT100`), which contains 100 typical
spectral functions and the corresponding correlation functions. These
spectral functions are defined in `src/dataset.jl`. Please use the script
`util/acstd.jl` to do this job.

### Arguments
N/A

### Returns
N/A
"""
function make_data_std()
    # Get dicts for the standard test (ACT100)
    ACT100 = union(STD_FG, STD_FD, STD_FRD, STD_BG, STD_BD, STD_BRD)

    # We have to make sure ntest == 100
    ntest = get_t("ntest")
    @assert ntest == length(ACT100) == 100

    # Initialize the random number generator
    # It is used to generate random noise only.
    seed = rand(1:10000) * myid() + 1981
    rng = MersenneTwister(seed)
    println("Random number seed: ", seed)

    # Start the loop
    println()
    for i = 1:ntest
        @printf("Test -> %4i / %4i\n", i, ntest)

        # Get dict for current test
        dict = ACT100[i]

        # Prepare grid for input data
        grid = make_grid(dict["grid"])
        println("Build grid for input data: ", length(grid), " points")

        # Prepare mesh for output spectrum
        mesh = make_mesh(dict["ktype"], dict["mesh"])
        println("Build mesh for spectrum: ", length(mesh), " points")

        # Prepare kernel function
        kernel = make_kernel(mesh, grid)
        println("Build default kernel: ", dict["ktype"])

        # Generate spectral functions
        sf = make_spectrum(mesh, dict["peaks"], dict["signs"])
        #
        # Generate Green's functions
        green = make_green(rng, sf, kernel, grid)
        #
        # Write generated data
        write_spectrum(i, sf)
        write_backward(i, green)
        #
        println()
    end
end

"""
    make_data()

Try to generate spectral functions and the corresponding Green's functions.

### Arguments
N/A

### Returns
N/A
"""
function make_data()
    # Get number of tests
    ntest = get_t("ntest")

    # Initialize the random number generator
    seed = 2002 #rand(1:10000) * myid() + 1981
    rng = MersenneTwister(seed)
    println("Random number seed: ", seed)

    # Prepare grid for input data
    grid = make_grid()
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    # Prepare kernel function
    kernel = make_kernel(mesh, grid)
    println("Build default kernel: ", get_t("ktype"))

    # Start the loop
    println()
    for i = 1:ntest
        @printf("Test -> %4i / %4i\n", i, ntest)
        #
        # Generate spectral functions
        sf = make_spectrum(rng, mesh)
        #
        # Generate Green's functions
        green = make_green(rng, sf, kernel, grid)
        #
        # Write generated data
        write_spectrum(i, sf)
        write_backward(i, green)
        #
        println()
    end
end

"""
    make_peak(rng::AbstractRNG)

Generate peak to build the final spectral function.

### Arguments
* rng -> Random number generator.

### Returns
* 𝑝 -> A Peak struct (subtype of AbstractPeak).

See also: [`AbstractPeak`](@ref).
"""
function make_peak(rng::AbstractRNG)
    ptype = get_t("ptype")
    pmax  = get_t("pmax")
    pmin  = get_t("pmin")

    @cswitch ptype begin
        @case "gauss"
            A = rand(rng)
            Γ = rand(rng)
            ϵ = rand(rng) * (pmax - pmin) + pmin
            𝑝 = GaussianPeak(A, Γ, ϵ)
            break

        @case "lorentz"
            A = rand(rng)
            Γ = rand(rng)
            ϵ = rand(rng) * (pmax - pmin) + pmin
            𝑝 = LorentzianPeak(A, Γ, ϵ)
            break

        @case "delta"
            A = rand(rng)
            Γ = 0.02 # Very sharp gaussian peak
            ϵ = rand(rng) * (pmax - pmin) + pmin
            𝑝 = DeltaPeak(A, Γ, ϵ)
            break

        @case "rectangle"
            c = rand(rng) * (pmax - pmin) + pmin
            w = rand(rng) * min(c - pmin, pmax - c) * 2.0
            h = rand(rng)
            @assert pmin ≤ c - w/2.0 ≤ c + w/2.0 ≤ pmax
            𝑝 = RectanglePeak(c, w, h)
            break

        @case "risedecay"
            c = rand(rng) * (pmax - pmin) + pmin
            γ = rand(rng) * 2.0 # γ ∈ (0.0,2.0)
            h = rand(rng)
            𝑝 = RiseDecayPeak(c, γ, h)
            break

        @default
            sorry()
            break
    end
    println(𝑝)

    return 𝑝
end

"""
    make_spectrum(rng::AbstractRNG, mesh::AbstractMesh)

Generate a spectral function randomly at given mesh.

### Arguments
* rng -> Random number generator.
* mesh -> Real frequency mesh, ω.

### Returns
* sf -> A SpectralFunction struct.

See also: [`make_data`](@ref).
"""
function make_spectrum(rng::AbstractRNG, mesh::AbstractMesh)
    # Extract essential parameters
    ktype = get_t("ktype")
    offdiag = get_t("offdiag")
    lpeak = get_t("lpeak")

    # Get number of peaks
    npeak, = rand(rng, lpeak, 1)
    @printf("number of peaks : %2i\n", npeak)

    # Determine signs for all peaks
    signs = ones(F64, npeak)
    if offdiag
        # How many negative signs are there?
        # We have to make sure that at least one sign is negative.
        nsign = rand(rng, 1:npeak)
        selected = rand(rng, 1:npeak, nsign)
        signs[selected] .= -1.0
        @assert count(x -> x < 0.0, signs) ≥ 1
    end

    image = zeros(F64, length(mesh))
    #
    for sign in signs
        # Get sign for the current peak
        @printf("sign : %4.2f\n", sign)
        #
        # Generate peak
        𝑝 = make_peak(rng)
        #
        # Add up to the spectrum
        image = image + sign * 𝑝(mesh)
    end
    #
    # Normalize the spectrum
    if !offdiag
        image = image ./ trapz(mesh,image)
    else
        # We have to make sure that A(ω) should exhibit negative weights
        # at some points.
        #
        # For fermionic systems
        if ktype == "fermi"
            # A(ω) > 0
            if count(x -> x < 0.0, image) == 0
                @. image = image * (-1.0)
            end
            #
            @assert any(x -> x < 0.0, image)
        # For bosonic systems
        else
            # A(ω) > 0
            if count(x -> x < 0.0, image .* mesh.mesh) == 0
                @. image = image * (-1.0)
            end
            #
            @assert any(x -> x < 0.0, image .* mesh.mesh)
        end
    end

    return SpectralFunction(mesh, image)
end

"""
    make_spectrum(
        mesh::AbstractMesh,
        pv::Vector{<:AbstractPeak},
        sv::Vector{F64}
    )

Generate a spectral function with given peaks and signs at given mesh.

### Arguments
* mesh -> Real frequency mesh, ω.
* pv -> List of peaks.
* sv -> List of signs.

### Returns
* sf -> A SpectralFunction struct.

See also: [`make_data_std`](@ref).
"""
function make_spectrum(
    mesh::AbstractMesh,
    pv::Vector{<:AbstractPeak},
    sv::Vector{F64}
    )
    # Get number of peaks
    npeak = length(pv)
    @assert length(pv) == length(sv)
    @printf("number of peaks : %2i\n", npeak)

    image = zeros(F64, length(mesh))
    #
    for i = 1:npeak
        # Determine sign of the current peak
        @printf("sign : %4.2f\n", sv[i])
        #
        # Generate peak
        𝑝 = pv[i]
        println(𝑝)
        #
        # Add up to the spectrum
        image = image + sv[i] * 𝑝(mesh)
    end
    #
    # Normalize the spectrum
    if count(x -> x > 0.0, sv) == npeak
        # All signs are positive. It is not for off-diagonal element.
        # We have to normalize the spectrum.
        image = image ./ trapz(mesh,image)
    end

    return SpectralFunction(mesh, image)
end

"""
    make_green(
        rng::AbstractRNG,
        sf::SpectralFunction,
        kernel::Matrix{F64},
        grid::AbstractGrid
    )

For given spectral function A and kernel matrix K, try to generate the
corresponding correlation function G (note that G ≡ KA).

### Arguments
* rng -> Random number generator.
* sf -> A SpectralFunction struct, A(ω).
* kernel -> Kernel matrix.
* grid -> Grid for correlation function.

### Returns
* gf -> A GreenFunction struct.
"""
function make_green(
    rng::AbstractRNG,
    sf::SpectralFunction,
    kernel::Matrix{F64},
    grid::AbstractGrid
    )
    # Get the noise level
    # If δ < 0, it means noise-free.
    δ = get_t("noise")

    # Calculate Green's function
    green = reprod(sf.mesh, kernel, sf.image)

    # Setup standard deviation
    ngrid = length(green)
    if δ < 0.0
        δ = 0.0
        error = fill(1.0e-4, ngrid)
    else
        error = fill(δ, ngrid)
    end

    # Setup random noise
    noise = randn(rng, F64, ngrid) * δ

    return GreenFunction(grid, green .+ noise, error)
end

"""
    make_grid(_grid = nothing)

To generate imaginary time grid or Masubara grid for many-body correlator.
It will return a sub-type of the AbstractGrid struct.

### Arguments
* _grid -> If it is present, we should use it to setup grid.

### Returns
* grid -> Imaginary time or imaginary frequency grid.

See also: [`AbstractGrid`](@ref).
"""
function make_grid(_grid = nothing)
    # Extract key parameters
    if isnothing(_grid)
        grid = get_t("grid")
    else
        grid = _grid
    end
    #
    ngrid = get_t("ngrid")
    β = get_t("beta")

    _grid = nothing
    @cswitch grid begin
        @case "ftime"
            _grid = FermionicImaginaryTimeGrid(ngrid, β)
            break

        @case "btime"
            _grid = BosonicImaginaryTimeGrid(ngrid, β)
            break

        @case "ffreq"
            _grid = FermionicMatsubaraGrid(ngrid, β)
            break

        @case "bfreq"
            _grid = BosonicMatsubaraGrid(ngrid, β)
            break

        @default
            sorry()
            break
    end

    return _grid
end

"""
    make_mesh(_ktype = nothing, _mesh = nothing)

Try to generate an uniform (linear) or non-uniform (non-linear) mesh for
the spectral function in real axis.

### Arguments
* _ktype -> If it is present, we should use it to setup ktype.
* _mesh -> If it is present, we should use it to setup mesh.

### Returns
* mesh -> Real frequency mesh. It should be a subtype of AbstractMesh.

See also: [`LinearMesh`](@ref), [`TangentMesh`](@ref), [`LorentzMesh`](@ref).
"""
function make_mesh(_ktype = nothing, _mesh = nothing)
    # Predefined parameters for mesh generation
    #
    # Note that the parameters `f1` and `cut` are only for the generation
    # of the non-uniform mesh.
    #
    f1::F64 = 2.1
    cut::F64 = 0.01

    # Setup parameters according to act.toml
    pmesh = get_t("pmesh")
    if !isa(pmesh, Missing)
        (length(pmesh) == 1) && begin
            Γ, = pmesh
            f1 = Γ
            cut = Γ
        end
    end

    # Get essential parameters
    if isnothing(_ktype)
        ktype = get_t("ktype")
    else
        ktype = _ktype
    end
    #
    if isnothing(_mesh)
        mesh = get_t("mesh")
    else
        mesh = _mesh
    end
    #
    nmesh = get_t("nmesh")
    wmax::F64 = get_t("wmax")
    wmin::F64 = get_t("wmin")
    #
    # For bosonic correlators of Hermitian operators, the spectral
    # function is defined in (0, ∞) only.
    if ktype == "bsymm"
        @assert wmin ≥ 0.0
        @assert wmax ≥ 0.0
        @assert wmax > wmin
    end

    # Try to generate the required mesh
    @cswitch mesh begin
        @case "linear"
            return LinearMesh(nmesh, wmin, wmax)
            break

        @case "tangent"
            return TangentMesh(nmesh, wmin, wmax, f1)
            break

        @case "lorentz"
            return LorentzMesh(nmesh, wmin, wmax, cut)
            break

        @case "halflorentz"
            return HalfLorentzMesh(nmesh, wmax, cut)
            break

        @default
            sorry()
            break
    end
end

"""
    make_kernel(am::AbstractMesh, ag::AbstractGrid)

Try to generate various kernel functions.

### Arguments
* am -> Real frequency mesh.
* ag -> Imaginary axis grid.

### Returns
* kernel -> Kernel function, a 2D array, (ntime,nmesh) or (nfreq,nmesh).

See also: [`AbstractMesh`](@ref), [`AbstractGrid`](@ref).
"""
function make_kernel(am::AbstractMesh, ag::AbstractGrid)
    ktype = get_t("ktype")
    grid = get_t("grid")

    @cswitch ktype begin
        @case "fermi"
            @assert grid in ("ftime", "ffreq")
            return build_kernel(am, ag)
            break

        @case "boson"
            @assert grid in ("btime", "bfreq")
            return build_kernel(am, ag)
            break

        @case "bsymm"
            @assert grid in ("btime", "bfreq")
            return build_kernel_symm(am, ag)
            break

        @default
            sorry()
            break
    end
end
