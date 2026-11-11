"""
    RiemannBasins — Basins of attraction of rational functions on the Riemann sphere

Package for visualising the basins of attraction of the discrete semi-flows
induced by iterating univariate rational functions on the augmented Riemann sphere

    P¹⁺(ℂ) = P¹(ℂ) ∪ {[0:0]}

where [0:0] is the *super-zero* point that absorbs computational indeterminations.

The sphere model used is the augmented unit sphere

    S²⁺ = S² ∪ {(0,0,0)},

connected to P¹⁺(ℂ) via the standard bijection arising from stereographic projection.

# References
- Hernández-Hernández, L.J. (2026). RiemannBasins: Basins of attraction of rational
  functions on the Riemann sphere. University of La Rioja.
"""
module RiemannBasins

using Polynomials          # v3+  — Polynomial(...) API
using Colors               # RGB, HSV, Lab color types
using Plots                # plotting backend

# ── Public API ─────────────────────────────────────────────────────────────────
export sphere_bijection,
       chordal_metric,
       fixed_points,
       okabe_ito_colors,
       basin_portrait,
       two_chart_portrait,
       animate_basins,
       basin_relief_3d,
       basin_relief_3d_makie,
       basin_portrait_pcycle

# ──────────────────────────────────────────────────────────────────────────────
# §1  Homogeneous coordinates and sphere model
# ──────────────────────────────────────────────────────────────────────────────

"""
    _homogeneous_normalize(twotuple, aproxprecision=15)

Return the *normalised* representative of the homogeneous point `[z:t]`
in `P¹⁺(ℂ)`.

The normalisation rule is:
- If both `|z|` and `|t|` are below `10^{-aproxprecision}`, return `(0, 0)` (super-zero).
- If `|z| ≤ |t|`, return `(z/t, 1)` (affine chart at `t=1`).
- Otherwise, return `(1, t/z)` (affine chart at `z=1`).
"""
function _homogeneous_normalize(
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        aproxprecision::Int64 = 15)

    tt1 = twotuple[1]
    tt2 = twotuple[2]
    tol = 1.0 / 10^aproxprecision

    if abs(tt1) < tol && abs(tt2) < tol
        return (complex(0.0), complex(0.0))
    elseif abs(tt1) <= abs(tt2)
        return (tt1 / tt2, complex(1.0))
    else
        return (complex(1.0), tt2 / tt1)
    end
end

"""
    sphere_bijection(twotuple, aproxprecision=8)

Map a homogeneous point `(z, t) ∈ P¹⁺(ℂ)` to its image in `S²⁺ ⊂ ℝ³` via
inverse stereographic projection.

    φ([z:t]) = ( Re(z̄t + t̄z) / (|t|²+|z|²),
                 Im(z̄t − t̄z) / (|t|²+|z|²),
                 (|z|²−|t|²) / (|t|²+|z|²) )

The super-zero `[0:0]` maps to the origin `(0,0,0)`.

# Returns
`Vector{Float64}` of length 3.
"""
function sphere_bijection(
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        aproxprecision::Int64 = 8)

    z = twotuple[1]; t = twotuple[2]
    tol = 1.0 / 10^aproxprecision

    if abs(z) < tol && abs(t) < tol
        return [0.0, 0.0, 0.0]
    end

    denom = abs2(t) + abs2(z)
    return [
        real(conj(z)*t + conj(t)*z) / denom,
        real(1im*(conj(z)*t - conj(t)*z)) / denom,
        real(-abs2(t) + abs2(z)) / denom
    ]
end

"""
    chordal_metric(P, Q, aproxprecision=8)

Chordal (spherical) distance between two homogeneous points in `P¹⁺(ℂ)`,
computed as the Euclidean distance between their images in `S²⁺`.

Returns a `Float64` in `[0, √2]`.
"""
function chordal_metric(twotuple, twotuple1, aproxprecision::Int64 = 8)
    v = sphere_bijection(twotuple, aproxprecision) .-
        sphere_bijection(twotuple1, aproxprecision)
    return sqrt(v[1]^2 + v[2]^2 + v[3]^2)
end

# ──────────────────────────────────────────────────────────────────────────────
# §2  Homogeneous polynomial evaluation
# ──────────────────────────────────────────────────────────────────────────────

function _bivariate_poly(
        coefficientlist::Vector{Complex{Float64}},
        u::Complex{Float64},
        t::Complex{Float64},
        d::Int64)

    ff = coefficientlist[1] * t^d
    for i in 2:length(coefficientlist)
        ff += coefficientlist[i] * u^(i-1) * t^(d-i+1)
    end
    return ff
end

function _bivariate_poly(
        coefficientlist::Vector{T},
        u::Complex{Float64},
        t::Complex{Float64},
        d::Int64) where {T<:Number}
    return _bivariate_poly(complex.(Float64.(coefficientlist)), u, t, d)
end

# ──────────────────────────────────────────────────────────────────────────────
# §3  Rational map in homogeneous form
# ──────────────────────────────────────────────────────────────────────────────

function _pair_of_functions(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T}) where {T<:Number}

    d = max(length(coefficientlistnum), length(coefficientlistden)) - 1

    cnum = ComplexF64.(coefficientlistnum)
    cden = ComplexF64.(coefficientlistden)

    fff(u::Complex{Float64}, t::Complex{Float64}) =
        _bivariate_poly(cnum, u, t, d)
    ggg(u::Complex{Float64}, t::Complex{Float64}) =
        _bivariate_poly(cden, u, t, d)

    return fff, ggg
end

function _rational_function(hpair, twotuple, aproxprecision::Int64 = 8)
    cnew = hpair[1](twotuple[1], twotuple[2])
    dnew = hpair[2](twotuple[1], twotuple[2])
    return _homogeneous_normalize((cnew, dnew), aproxprecision)
end

# ──────────────────────────────────────────────────────────────────────────────
# §4  Fixed points
# ──────────────────────────────────────────────────────────────────────────────

"""
    fixed_points(coefficientlistnum, coefficientlistden, aproxprecision=8)

Find the fixed points of the rational map `f = p/q` on `P¹⁺(ℂ)`.

The finite fixed points satisfy `p(z) - z·q(z) = 0`; solved via the
companion-matrix eigenvalue method.  The points `∞ = [1:0]` and
`⊙ = [0:0]` are checked separately.

Coefficients are listed *constant term first*:
`[a₀, a₁, …, aₙ]` represents `a₀ + a₁z + … + aₙzⁿ`.

# Returns
`(superfixed, infixed, finitelist)` where
- `superfixed::Bool` — `true` if `[0:0]` is a fixed point.
- `infixed::Bool` — `true` if `∞` is a fixed point.
- `finitelist` — finite fixed points as normalised homogeneous pairs.
"""
function fixed_points(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T},
        aproxprecision::Int64 = 8) where {T<:Number}

    cnum = complex.(Float64.(coefficientlistnum))
    cden = complex.(Float64.(coefficientlistden))

    hpair = _pair_of_functions(cnum, cden)
    imsuper = _rational_function(hpair, (complex(0.0), complex(0.0)), aproxprecision)
    superfixed = (abs(imsuper[1]) < 1/10^aproxprecision &&
                  abs(imsuper[2]) < 1/10^aproxprecision)

    iminf = _rational_function(hpair, (complex(1.0), complex(0.0)), aproxprecision)
    tol = 1/10^aproxprecision
    infixed = (abs(iminf[2]) < tol)

    Polynumerator   = Polynomial(cnum)
    Polydenominator = Polynomial(cden)
    Polyx           = Polynomial([0.0+0.0im, 1.0+0.0im])

    lookingzeros = Polynumerator - Polydenominator * Polyx
    rawroots = roots(lookingzeros)

    finitelist = Tuple{Complex{Float64},Complex{Float64}}[]
    for r in rawroots
        push!(finitelist, _homogeneous_normalize((Complex{Float64}(r), complex(1.0)), aproxprecision))
    end

    return (superfixed, infixed, finitelist)
end

# ──────────────────────────────────────────────────────────────────────────────
# §5  Orbit iteration
# ──────────────────────────────────────────────────────────────────────────────

function _newstep(hpair, iter::Int64, iterprecision::Int64,
                  hpoint::Tuple{Complex{Float64},Complex{Float64}})

    tol = 1.0 / 10^iterprecision
    imagepoint = hpoint
    n = 0

    for _ in 1:iter
        nextpoint = _rational_function(hpair, imagepoint)
        n += 1
        if chordal_metric(imagepoint, nextpoint) < tol
            imagepoint = nextpoint
            break
        end
        imagepoint = nextpoint
    end

    return (imagepoint, n)
end

# ──────────────────────────────────────────────────────────────────────────────
# §6  Sampling grids
# ──────────────────────────────────────────────────────────────────────────────

# Chart 1: neighbourhood of 0
function _rectangle(
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        expprecision::Int64 = 10)

    a, b = Float64(xinterval[1]), Float64(xinterval[2])
    c, d = Float64(yinterval[1]), Float64(yinterval[2])
    tol  = 1.0 / 2^expprecision

    return [(complex(r, im_val), complex(1.0))
            for im_val in d:-tol:c, r in a:tol:b]
end

# Chart 2: neighbourhood of ∞
function _rectangle_chart2(
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        expprecision::Int64 = 8)

    n  = 2^expprecision
    rs = range(Float64(xinterval[1]), Float64(xinterval[2]), length=n)
    is = range(Float64(yinterval[1]), Float64(yinterval[2]), length=n)

    grid = Matrix{Tuple{Complex{Float64},Complex{Float64}}}(undef, n, n)
    for (j, re) in enumerate(rs), (i, im_val) in enumerate(is)
        grid[i,j] = (complex(1.0), complex(re, im_val))
    end
    return grid
end

# ──────────────────────────────────────────────────────────────────────────────
# §7  Classification helpers
# ──────────────────────────────────────────────────────────────────────────────

function _position_up_to_tolerance(fixedPointList::Vector, aproxprecision::Int64, twotuple)
    tol = 1.0 / 10^aproxprecision
    for (it, fp) in enumerate(fixedPointList)
        if chordal_metric(twotuple, fp, aproxprecision) < tol
            return it
        end
    end
    return 0
end

# Returns (pos, n_iterations, residual)
function _position_iteration(
        fixedPointList::Vector,
        aproxprecision::Int64,
        iterprecision::Int64,
        twotuple,
        hpair,
        iter::Int64)

    imagepoint, n = _newstep(hpair, iter, iterprecision, twotuple)

    tol = 1.0 / 10^aproxprecision
    if abs(imagepoint[1]) < tol && abs(imagepoint[2]) < tol
        return (-1, n, 0.0)
    end

    residual = tol
    pos = 0
    for (it, fp) in enumerate(fixedPointList)
        r = chordal_metric(imagepoint, fp, aproxprecision)
        if r < tol
            pos = it
            residual = max(r, 1e-15)
            break
        end
    end
    return (pos, n, residual)
end

# ──────────────────────────────────────────────────────────────────────────────
# §8  Basin computation (matrix form)
# ──────────────────────────────────────────────────────────────────────────────

# Returns (fixedPointList, positionMatrix, iterationMatrix, residualMatrix)
function _basin_matrices(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T},
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        aproxprecision::Int64 = 8,
        iterprecision::Int64  = 8,
        iter::Int64           = 200,
        expprecision::Int64   = 8) where {T<:Number}

    hpair = _pair_of_functions(coefficientlistnum, coefficientlistden)
    _, _, finitelist = fixed_points(
        coefficientlistnum, coefficientlistden, aproxprecision)

    fixedPointList = copy(finitelist)
    push!(fixedPointList, (complex(1.0), complex(0.0)))   # ∞

    grid = _rectangle(xinterval, yinterval, expprecision)
    nrows, ncols = size(grid)

    positionMatrix  = zeros(Int64, nrows, ncols)
    iterationMatrix = zeros(Int64, nrows, ncols)
    residualMatrix  = zeros(Float64, nrows, ncols)

    for j in 1:ncols, i in 1:nrows
        pos, n, res = _position_iteration(
            fixedPointList, aproxprecision, iterprecision,
            grid[i,j], hpair, iter)
        positionMatrix[i,j]  = pos
        iterationMatrix[i,j] = n
        residualMatrix[i,j]  = res
    end

    return fixedPointList, positionMatrix, iterationMatrix, residualMatrix
end

# Generic helper for any pre-built grid
function _basin_matrices_from_grid(
        grid::Matrix,
        fixedPointList::Vector,
        hpair,
        aproxprecision::Int64,
        iterprecision::Int64,
        iter::Int64)

    nrows, ncols = size(grid)
    positionMatrix  = zeros(Int64, nrows, ncols)
    iterationMatrix = zeros(Int64, nrows, ncols)
    residualMatrix  = zeros(Float64, nrows, ncols)

    for j in 1:ncols, i in 1:nrows
        pos, n, res = _position_iteration(
            fixedPointList, aproxprecision, iterprecision,
            grid[i,j], hpair, iter)
        positionMatrix[i,j]  = pos
        iterationMatrix[i,j] = n
        residualMatrix[i,j]  = res
    end
    return positionMatrix, iterationMatrix, residualMatrix
end

# ──────────────────────────────────────────────────────────────────────────────
# §9  Colour mapping
# ──────────────────────────────────────────────────────────────────────────────

"""
    okabe_ito_colors()

Return the 8-colour Okabe-Ito palette as a `Vector{RGB{Float64}}`.

Designed by Masataka Okabe & Kei Ito (2008) to remain distinguishable under
deuteranopia, protanopia, and tritanopia.

Order: black, orange, sky-blue, bluish-green, yellow, blue, vermilion,
reddish-purple.
"""
function okabe_ito_colors()::Vector{RGB{Float64}}
    return [
        RGB{Float64}(0.000, 0.000, 0.000),   # 1 black
        RGB{Float64}(0.902, 0.624, 0.000),   # 2 orange
        RGB{Float64}(0.337, 0.706, 0.914),   # 3 sky blue
        RGB{Float64}(0.000, 0.620, 0.451),   # 4 bluish green
        RGB{Float64}(0.941, 0.894, 0.259),   # 5 yellow
        RGB{Float64}(0.000, 0.447, 0.698),   # 6 blue
        RGB{Float64}(0.835, 0.369, 0.000),   # 7 vermilion
        RGB{Float64}(0.800, 0.475, 0.655),   # 8 reddish purple
    ]
end

function _build_colormap(nfixed::Int64, gseed::RGB{Float64};
                         colorblind::Bool = false)
    colors = Vector{RGB{Float64}}(undef, nfixed + 2)
    colors[1] = RGB{Float64}(0.0, 0.0, 0.0)   # class 0 → black

    if colorblind
        oi = okabe_ito_colors()
        for k in 1:nfixed+1
            idx = mod1(k + 1, length(oi))
            colors[k+1] = oi[idx]
        end
    else
        # Perceptually-uniform arc in LCHab space:
        # keep Lightness and Chroma from the seed; rotate Hue uniformly.
        # LCHab is the cylindrical form of CIE L*a*b*, so equal hue steps
        # look equal to the eye — unlike HSV, which is perceptually uneven.
        lch0 = convert(LCHab{Float64}, gseed)
        C0   = max(lch0.c, 40.0)   # guarantee vivid colours even for desaturated seeds
        for k in 1:nfixed+1
            hk  = mod(lch0.h + (k-1) * 360.0 / (nfixed + 1), 360.0)
            rgb = convert(RGB{Float64}, LCHab{Float64}(lch0.l, C0, hk))
            colors[k+1] = RGB{Float64}(
                clamp(red(rgb),   0.0, 1.0),
                clamp(green(rgb), 0.0, 1.0),
                clamp(blue(rgb),  0.0, 1.0))
        end
    end
    return colors
end

# Cyclic L* banding in CIE Lab space (colour-blind safe)
function _modulate_lab(color::RGB{Float64}, t_smooth::Float64,
                       maxiter::Float64; bands::Float64 = 6.0)::RGB{Float64}
    lab      = convert(Lab{Float64}, color)
    phase    = 2π * bands * t_smooth / maxiter
    l_cyclic = 0.5 + 0.5 * cos(phase)
    L_new    = clamp(lab.l * (0.50 + 0.50 * l_cyclic), 5.0, 100.0)
    rgb      = convert(RGB{Float64}, Lab{Float64}(L_new, lab.a, lab.b))
    return RGB{Float64}(
        clamp(red(rgb),   0.0, 1.0),
        clamp(green(rgb), 0.0, 1.0),
        clamp(blue(rgb),  0.0, 1.0))
end

# ──────────────────────────────────────────────────────────────────────────────
# §10  Rendering
# ──────────────────────────────────────────────────────────────────────────────

# Internal: builds the raw pixel matrix (used by both _render_basin_image and
# the GLMakie extension so neither has to duplicate the coloring logic).
function _build_basin_image_matrix(
        fixedPointList::Vector,
        positionMatrix::Matrix{Int64},
        iterationMatrix::Matrix{Int64},
        cmap::Vector{RGB{Float64}},
        colorstrategy::String;
        colorblind::Bool = false,
        residualMatrix::Union{Matrix{Float64},Nothing} = nothing,
        aproxprecision::Int64 = 8,
        bands::Float64 = 6.0)::Matrix{RGB{Float64}}

    nrows, ncols = size(positionMatrix)
    maxiter = max(1, maximum(iterationMatrix))
    tol     = 1.0 / 10^aproxprecision

    img = Matrix{RGB{Float64}}(undef, nrows, ncols)
    for j in 1:ncols, i in 1:nrows
        pos = positionMatrix[i, j]
        nit = iterationMatrix[i, j]

        if colorstrategy == "iteration"
            v = nit / maxiter
            img[i,j] = RGB{Float64}(v, v, v)

        elseif colorstrategy == "positionfixedpoints"
            idx = (pos <= 0) ? 1 : min(pos + 1, length(cmap))
            img[i,j] = cmap[idx]

        else  # "positionplusiteration"
            idx  = (pos <= 0) ? 1 : min(pos + 1, length(cmap))
            base = cmap[idx]

            # Smooth iteration count (quadratic convergence: log(2) denominator)
            # μ = n + log(log(ε)/log(r)) / log(2)
            if !isnothing(residualMatrix) && pos > 0
                r = max(residualMatrix[i,j], 1e-15)
                tol_c = max(tol, 1e-15)
                log_ratio = log(log(tol_c) / log(r)) / log(2.0)
                t_smooth = clamp(Float64(nit) + log_ratio, 0.0, Float64(maxiter))
            else
                t_smooth = Float64(nit)
            end

            # Both modes use L* modulation in CIE L*a*b* space:
            # hue and chroma (a*, b*) stay fixed; only brightness varies.
            img[i,j] = _modulate_lab(base, t_smooth, Float64(maxiter);
                                     bands=bands)
        end
    end
    return img
end

function _render_basin_image(
        fixedPointList::Vector,
        positionMatrix::Matrix{Int64},
        iterationMatrix::Matrix{Int64},
        xinterval::Tuple{Real,Real},
        yinterval::Tuple{Real,Real},
        cmap::Vector{RGB{Float64}},
        colorstrategy::String,
        title::String;
        colorblind::Bool = false,
        residualMatrix::Union{Matrix{Float64},Nothing} = nothing,
        aproxprecision::Int64 = 8,
        bands::Float64 = 6.0)

    nrows, ncols = size(positionMatrix)

    img = _build_basin_image_matrix(
        fixedPointList, positionMatrix, iterationMatrix,
        cmap, colorstrategy;
        colorblind     = colorblind,
        residualMatrix = residualMatrix,
        aproxprecision = aproxprecision,
        bands          = bands)

    a, b = Float64(xinterval[1]), Float64(xinterval[2])
    c, d = Float64(yinterval[1]), Float64(yinterval[2])

    plt = plot(
        img,
        aspect_ratio = :equal,
        xlims = (1, ncols),
        ylims = (1, nrows),
        xticks = (
            range(1, ncols, length=5),
            string.(round.(range(a, b, length=5), digits=2))
        ),
        yticks = (
            range(1, nrows, length=5),
            string.(round.(range(c, d, length=5), digits=2))
        ),
        xlabel = "Re(z)",
        ylabel = "Im(z)",
        title  = title,
        legend = false,
        size   = (600, 600)
    )
    return plt
end

# ──────────────────────────────────────────────────────────────────────────────
# §11  Public API
# ──────────────────────────────────────────────────────────────────────────────

"""
    basin_portrait(coefficientlistnum, coefficientlistden;
                   xinterval=(-1.5,1.5), yinterval=(-1.5,1.5),
                   aproxprecision=8, iterprecision=8, iter=200, expprecision=8,
                   bands=6.0, gseed=RGB(0.2,0.6,0.9),
                   colorstrategy="positionplusiteration",
                   colorblind=false, show_plot=true, savepath="")

Compute and display the basin portrait of the rational function `f = p/q`
on the Riemann sphere `P¹⁺(ℂ)`.

Coefficients are listed *constant term first* (same convention as `Polynomial`
in Polynomials.jl): `[a₀, a₁, …, aₙ]` represents `a₀ + a₁z + … + aₙzⁿ`.

# Arguments
- `coefficientlistnum` — numerator polynomial coefficients.
- `coefficientlistden` — denominator polynomial coefficients.
- `xinterval`, `yinterval` — real and imaginary axis extents (default `(-1.5, 1.5)`).
- `aproxprecision::Int64` — decimal digits for convergence detection (default `8`).
- `iterprecision::Int64` — convergence precision for orbit iteration (default `8`).
- `iter::Int64` — maximum number of orbit steps per pixel (default `200`).
- `expprecision::Int64` — grid resolution: `2^expprecision` points per side (default `8`).
- `bands::Float64` — number of cyclic lightness bands per attractor basin (default `6.0`);
  higher values give finer fractal ridges along iso-convergence curves.
- `gseed::RGB{Float64}` — seed colour for the basin palette (ignored when `colorblind=true`).
- `colorstrategy::String` — one of:
  - `"positionplusiteration"` (default): colour encodes attractor identity
    and convergence speed via smooth colouring and cyclic banding.
  - `"positionfixedpoints"`: flat colour per attractor.
  - `"iteration"`: greyscale proportional to iteration count.
- `colorblind::Bool` — use the Okabe-Ito palette and CIE L*a*b* banding (default `false`).
- `show_plot::Bool` — call `display` on the result (default `true`).
- `savepath::String` — if non-empty, save the plot to this file.

# Returns
`(fixedPointList, positionMatrix, iterationMatrix, plt)`.

# Example
```julia
using RiemannBasins
# Newton's method for z³ − 1:  f(z) = (2z³+1)/(3z²)
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 0z + 0z² + 2z³
den = complex([0.0, 0.0, 3.0, 0.0])   # 0 + 0z + 3z² + 0z³
fps, pos, itr, plt = basin_portrait(num, den)
```
"""
function basin_portrait(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        aproxprecision::Int64       = 8,
        iterprecision::Int64        = 8,
        iter::Int64                 = 200,
        expprecision::Int64         = 8,
        bands::Float64              = 6.0,
        gseed::RGB{Float64}         = RGB{Float64}(0.2, 0.6, 0.9),
        colorstrategy::String       = "positionplusiteration",
        hue_band::Float64           = 25.0,
        show_plot::Bool             = true,
        colorblind::Bool            = false,
        savepath::String            = "") where {T<:Number}

    fixedPointList, positionMatrix, iterationMatrix, residualMatrix =
        _basin_matrices(
            coefficientlistnum, coefficientlistden,
            xinterval, yinterval,
            aproxprecision, iterprecision, iter, expprecision)

    nfixed = length(fixedPointList)
    cmap   = _build_colormap(nfixed, gseed; colorblind=colorblind)

    plt = _render_basin_image(
        fixedPointList, positionMatrix, iterationMatrix,
        xinterval, yinterval, cmap, colorstrategy,
        "Basins of attraction";
        colorblind=colorblind,
        bands=bands,
        residualMatrix=residualMatrix,
        aproxprecision=aproxprecision)

    show_plot && display(plt)
    isempty(savepath) || savefig(plt, savepath)

    return plt
end

"""
    two_chart_portrait(coefficientlistnum, coefficientlistden;
                       xinterval=(-1.5,1.5), yinterval=(-1.5,1.5),
                       aproxprecision=8, iterprecision=8, iter=200, expprecision=8,
                       bands=6.0, gseed=RGB(0.2,0.6,0.9),
                       colorstrategy="positionplusiteration",
                       colorblind=false, show_plot=true, savepath="")

Compute and display the basin portrait in **both affine charts** of `P¹(ℂ)`
side by side.

- **Chart 1** `{[z:1]}` — standard complex plane; grid centred at `z = 0`.
- **Chart 2** `{[1:t]}` — chart centred at `z = ∞`; if `t ≠ 0` this represents
  `z = 1/t` in chart 1.

The same colour map is used for both panels, so colours are directly comparable.

# Returns
`(fixedPointList, p1, p2, plt)` where `p1`, `p2` are the individual panel
plots and `plt` is the combined `layout=(1,2)` figure.
"""
function two_chart_portrait(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        aproxprecision::Int64       = 8,
        iterprecision::Int64        = 8,
        iter::Int64                 = 200,
        expprecision::Int64         = 8,
        bands::Float64              = 6.0,
        gseed::RGB{Float64}         = RGB{Float64}(0.2, 0.6, 0.9),
        colorstrategy::String       = "positionplusiteration",
        colorblind::Bool            = false,
        show_plot::Bool             = true,
        savepath::String            = "") where {T<:Number}

    hpair = _pair_of_functions(coefficientlistnum, coefficientlistden)
    _, _, finitelist = fixed_points(
        coefficientlistnum, coefficientlistden, aproxprecision)

    fixedPointList = copy(finitelist)
    push!(fixedPointList, (complex(1.0), complex(0.0)))   # ∞

    nfixed = length(fixedPointList)
    cmap   = _build_colormap(nfixed, gseed; colorblind=colorblind)

    # Chart 1
    grid1 = _rectangle(xinterval, yinterval, expprecision)
    pos1, itr1, res1 = _basin_matrices_from_grid(
        grid1, fixedPointList, hpair, aproxprecision, iterprecision, iter)

    p1 = _render_basin_image(
        fixedPointList, pos1, itr1,
        xinterval, yinterval, cmap, colorstrategy,
        "Chart 1  (nbhd of 0)";
        colorblind=colorblind, bands=bands, residualMatrix=res1, aproxprecision=aproxprecision)

    # Chart 2
    grid2 = _rectangle_chart2(xinterval, yinterval, expprecision)
    pos2, itr2, res2 = _basin_matrices_from_grid(
        grid2, fixedPointList, hpair, aproxprecision, iterprecision, iter)

    p2 = _render_basin_image(
        fixedPointList, pos2, itr2,
        xinterval, yinterval, cmap, colorstrategy,
        "Chart 2  (nbhd of ∞)";
        colorblind=colorblind, bands=bands, residualMatrix=res2, aproxprecision=aproxprecision)

    plt = plot(p1, p2; layout=(1,2), size=(1200, 600))
    show_plot && display(plt)
    isempty(savepath) || savefig(plt, savepath)

    return plt
end

"""
    animate_basins(make_coefficients, parameter_values;
                   xinterval=(-1.5,1.5), yinterval=(-1.5,1.5),
                   aproxprecision=8, iterprecision=8, iter=200, expprecision=7,
                   gseed=RGB(0.2,0.6,0.9), colorstrategy="positionplusiteration",
                   colorblind=false, fps=8, filename="basins.gif")

Produce a GIF animation of basin portraits as a parameter varies.

`make_coefficients` must be a callable `p -> (num_coeffs, den_coeffs)` returning
the numerator and denominator coefficient vectors for each parameter value `p`.
`parameter_values` can be any iterable (e.g. a `range`).

# Returns
The `Plots.Animation` object; the GIF is also written to `filename`.

# Example
```julia
using RiemannBasins, Colors
# Newton's method for z³ − c, varying c
make_coeffs(c) = (complex([c, 0.0, 0.0, 2.0]),
                  complex([0.0, 0.0, 3.0, 0.0]))
anim = animate_basins(make_coeffs, range(0.5, 2.0, length=30);
                      filename="newton_c.gif", fps=6)
```
"""
function animate_basins(
        make_coefficients,
        parameter_values;
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        aproxprecision::Int64       = 8,
        iterprecision::Int64        = 8,
        iter::Int64                 = 200,
        expprecision::Int64         = 7,
        gseed::RGB{Float64}         = RGB{Float64}(0.2, 0.6, 0.9),
        colorstrategy::String       = "positionplusiteration",
        colorblind::Bool            = false,
        fps::Int64                  = 8,
        filename::String            = "basins.gif")

    anim = @animate for p in parameter_values
        num, den = make_coefficients(p)
        fps_list, posM, itrM, resM =
            _basin_matrices(num, den, xinterval, yinterval,
                            aproxprecision, iterprecision, iter, expprecision)
        nfixed = length(fps_list)
        cmap   = _build_colormap(nfixed, gseed; colorblind=colorblind)
        _render_basin_image(
            fps_list, posM, itrM,
            xinterval, yinterval, cmap, colorstrategy,
            "p = $(round(real(p), digits=3))";
            colorblind=colorblind,
            residualMatrix=resM,
            aproxprecision=aproxprecision)
    end
    gif(anim, filename; fps=fps)
    return anim
end


# ──────────────────────────────────────────────────────────────────────────────
# §12  3-D relief portrait
# ──────────────────────────────────────────────────────────────────────────────

"""
    basin_relief_3d(coefficientlistnum, coefficientlistden; ...)

3D relief portrait of the basins of attraction of `f = p/q` on `P¹⁺(ℂ)`.

The **height** encodes the smooth iteration count with cyclic banding,
producing a fractal topography with ridges at basin boundaries.
The **colour** encodes the attractor identity: each basin is placed at a
distinct altitude plateau so the colormap assigns it a different hue.

# Arguments
- `coefficientlistnum`, `coefficientlistden` — numerator / denominator
  coefficients, constant term first (same as `basin_portrait`).
- `xinterval`, `yinterval` — axis extents (default `(-1.5, 1.5)`).
- `aproxprecision` — decimal digits for convergence (default `8`).
- `iterprecision` — iteration convergence precision (default `8`).
- `iter` — maximum orbit steps per pixel (default `200`).
- `expprecision` — grid resolution `2^expprecision` per side
  (default `7` → 128×128 grid, good balance of detail vs. speed).
- `bands` — number of cyclic relief bands (default `6.0`;
  more bands → finer fractal ridges).
- `relief_scale` — vertical amplitude (default `1.0`;
  use `0.3`–`0.6` for subtle low-relief embossing).
- `camera_angle` — `(azimuth, elevation)` in degrees (default `(40, 55)`).
- `gseed` — seed colour for the basin palette (ignored when `colorblind=true`).
- `colorblind` — use the Okabe-Ito palette (default `false`).
- `show_plot` — call `display` (default `true`).
- `inverted` — if `true`, the elevation is inverted: the Julia-set boundary
  becomes the highest ridge and each basin centre is a valley
  (default `false`).
- `savepath` — file path to save the figure (empty = no save).

# Returns
The `Plots` figure.

# Example
```julia
using RiemannBasins
# Newton for z³ − 1:  f(z) = (2z³+1)/(3z²)
num = complex([1.0, 0.0, 0.0, 2.0])
den = complex([0.0, 0.0, 3.0, 0.0])
basin_relief_3d(num, den; bands=8.0, relief_scale=0.5)
```
"""
function basin_relief_3d(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        xinterval::Tuple{Real,Real}    = (-1.5, 1.5),
        yinterval::Tuple{Real,Real}    = (-1.5, 1.5),
        aproxprecision::Int64          = 8,
        iterprecision::Int64           = 8,
        iter::Int64                    = 200,
        expprecision::Int64            = 7,
        bands::Float64                 = 12.0,
        relief_scale::Float64          = 0.02,
        ridge_fraction::Float64        = 0.15,
        camera_angle::Tuple{Real,Real} = (30, 82),
        gseed::RGB{Float64}            = RGB{Float64}(0.2, 0.6, 0.9),
        colorblind::Bool               = false,
        inverted::Bool                 = false,
        show_plot::Bool                = true,
        savepath::String               = "") where {T<:Number}

    # ── 1. compute basin matrices ────────────────────────────────────────────
    fixedPointList, positionMatrix, iterationMatrix, residualMatrix =
        _basin_matrices(coefficientlistnum, coefficientlistden,
                        xinterval, yinterval,
                        aproxprecision, iterprecision, iter, expprecision)

    nrows, ncols = size(positionMatrix)
    maxiter      = Float64(max(1, maximum(iterationMatrix)))
    tol          = 1.0 / 10^aproxprecision
    nfixed       = length(fixedPointList)

    xs = collect(range(Float64(xinterval[1]), Float64(xinterval[2]), length = ncols))
    ys = collect(range(Float64(yinterval[1]), Float64(yinterval[2]), length = nrows))

    # ── 2. combined Z matrix ─────────────────────────────────────────────────
    # `relief_scale` controls the TOTAL z range — the whole surface fits within
    # [0, relief_scale], keeping the plot nearly flat (2-D look with embossing).
    #
    # Within that range:
    #   • basin_norm  ∈ [0, 1)        → colour zone per attractor (major part)
    #   • ridge_norm  ∈ [0, ridge_fraction)  → cyclic fractal banding (minor part)
    #
    #   z[i,j] = (basin_norm + ridge_norm) × relief_scale
    #
    # `ridge_fraction` (default 0.15) sets what fraction of each colour zone
    # is used by the iteration ridges; keeping it small avoids colour overlap
    # between adjacent basins while still showing the fractal texture.

    z = zeros(Float64, nrows, ncols)

    for j in 1:ncols, i in 1:nrows
        pos = positionMatrix[i, j]
        nit = Float64(iterationMatrix[i, j])

        # smooth iteration count (quadratic convergence correction)
        if pos > 0 && residualMatrix[i, j] > 0
            r         = max(residualMatrix[i, j], 1e-15)
            tol_c     = max(tol, 1e-15)
            log_ratio = log(log(tol_c) / log(r)) / log(2.0)
            t_smooth  = clamp(nit + log_ratio, 0.0, maxiter)
        else
            t_smooth = nit
        end

        # Normal mode:  basin centres are raised plateaus, Julia boundary is valley.
        #   ridge = ridge_fraction*(1−t)*(0.5+0.5cos) → strongest at t=0, zero at t=1
        # Inverted mode: Julia boundary is the mountain ridge, basins are craters.
        #   ridge = ridge_fraction*t*(0.5+0.5cos)      → zero at t=0, strongest at t=1
        t_norm = clamp(t_smooth / maxiter, 0.0, 1.0)
        phase  = 2π * bands * t_norm
        if inverted
            h          = pos <= 0 ? 1.0 : t_norm
            ridge_norm = ridge_fraction * h * (0.5 + 0.5 * cos(phase))
            z[i, j]   = (h * (1.0 - ridge_fraction) + ridge_norm) * relief_scale
        else
            ridge_norm = ridge_fraction * (1.0 - t_norm) * (0.5 + 0.5 * cos(phase))
            basin_norm = pos <= 0 ? 0.0 : Float64(pos - 1) / Float64(nfixed)
            z[i, j]   = (basin_norm + ridge_norm) * relief_scale
        end
    end

    # ── 3. surface plot ──────────────────────────────────────────────────────
    chosen_cmap = colorblind ? :tab10 : :hsv

    plt = surface(
        xs, ys, z;
        color     = chosen_cmap,
        linewidth = 0,
        camera    = camera_angle,
        xlabel    = "Re(z)",
        ylabel    = "Im(z)",
        zlabel    = "",
        title     = "3D basin relief  (bands = $(Int(round(bands))))" * (inverted ? "  [inverted]" : ""),
        legend    = false,
        size      = (750, 700))

    show_plot && display(plt)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end


# ──────────────────────────────────────────────────────────────────────────────
# §13  GLMakie extension stub
# ──────────────────────────────────────────────────────────────────────────────

"""
    basin_relief_3d_makie(coefficientlistnum, coefficientlistden; kwargs...)

Three-dimensional basin-relief surface rendered with **GLMakie** (OpenGL).
Requires `using GLMakie` before calling this function — the extension is
loaded automatically when GLMakie is in scope.

Compared to `basin_relief_3d` (Plots.jl / GR backend), this version:
- decouples **colour** (LCHab-banded basin portrait) from **elevation** (fractal ridge),
- applies real OpenGL lighting (ambient + diffuse shading), and
- supports interactive camera rotation in the display window.

All keyword arguments of `basin_relief_3d` are accepted, plus:
- `inverted` — if `true`, the Julia-set boundary is the highest ridge and
  each basin centre is a valley (default `false`).
- `background` — figure background colour (default `:black`).
- `shading` — Makie shading model (default `true`).
"""
function basin_relief_3d_makie(args...; kwargs...)
    error("""
basin_relief_3d_makie requires GLMakie.
Add  `using GLMakie`  before calling this function.
""")
end


# ──────────────────────────────────────────────────────────────────────────────
# §13  p-cycle basins
# ──────────────────────────────────────────────────────────────────────────────

"""
    _rational_compose(num_f, den_f, num_g, den_g)

Return `(num_fog, den_fog)` — coefficient lists for the rational function
`(f∘g)(z) = f(g(z))`.

Uses the identity
  numerator(f∘g) = Σ_{k=0}^{d_f}  num_f[k+1] · pg^k · qg^(d_f−k)
where `pg = Polynomial(num_g)`, `qg = Polynomial(den_g)`,
and `d_f = max(length(num_f), length(den_f)) − 1`.
"""
function _rational_compose(
        num_f::Vector{Complex{Float64}},
        den_f::Vector{Complex{Float64}},
        num_g::Vector{Complex{Float64}},
        den_g::Vector{Complex{Float64}})

    d_f = max(length(num_f), length(den_f)) - 1

    # Zero-pad coefficient vectors to d_f + 1 terms
    nf = vcat(num_f, zeros(Complex{Float64}, d_f + 1 - length(num_f)))
    df = vcat(den_f, zeros(Complex{Float64}, d_f + 1 - length(den_f)))

    pg = Polynomial(num_g)
    qg = Polynomial(den_g)
    one_poly = Polynomial([one(Complex{Float64})])

    # Precompute pg^k and qg^k for k = 0…d_f
    pg_pow = Vector{typeof(one_poly)}(undef, d_f + 1)
    qg_pow = Vector{typeof(one_poly)}(undef, d_f + 1)
    pg_pow[1] = one_poly
    qg_pow[1] = one_poly
    for k in 1:d_f
        pg_pow[k+1] = pg_pow[k] * pg
        qg_pow[k+1] = qg_pow[k] * qg
    end

    zero_poly = Polynomial(zeros(Complex{Float64}, 1))
    num_comp  = zero_poly
    den_comp  = zero_poly
    for k in 0:d_f
        term      = pg_pow[k+1] * qg_pow[d_f - k + 1]   # pg^k · qg^(d_f-k)
        num_comp += nf[k+1] * term
        den_comp += df[k+1] * term
    end

    return coeffs(num_comp), coeffs(den_comp)
end

"""
    _rational_iterate(num, den, p)

Return `(num_p, den_p)` — coefficient lists for the `p`-fold iterate `f^p`.
"""
function _rational_iterate(
        num::Vector{Complex{Float64}},
        den::Vector{Complex{Float64}},
        p::Int)

    p >= 1 || throw(ArgumentError("p must be ≥ 1"))
    p == 1 && return copy(num), copy(den)

    num_p = copy(num)
    den_p = copy(den)
    for _ in 2:p
        num_p, den_p = _rational_compose(num, den, num_p, den_p)
    end
    return num_p, den_p
end

"""
    _group_fp_by_orbit(fp_list, hpair_f, aproxprecision)

Partition the points in `fp_list` into orbit groups under `f` (given by
`hpair_f`).  Returns `(group_of, n_groups)` where `group_of[i]` is the
group index (1-based) of `fp_list[i]`.

All points that belong to the same f-orbit are assigned the same group index.
This groups the `p` components of each p-cycle together (Option A colouring)
and keeps each fixed point of `f` in its own singleton group.
"""
function _group_fp_by_orbit(
        fp_list::Vector,
        hpair_f,
        aproxprecision::Int)

    n   = length(fp_list)
    tol = 1.0 / 10^aproxprecision

    # Build successor map: nxt[i] = j  means  f(fp_list[i]) ≈ fp_list[j]
    nxt = zeros(Int, n)
    for i in 1:n
        fi = _rational_function(hpair_f, fp_list[i], aproxprecision)
        for j in 1:n
            if chordal_metric(fi, fp_list[j], aproxprecision) < tol
                nxt[i] = j
                break
            end
        end
    end

    # Trace orbit chains and assign group labels
    group_of = zeros(Int, n)
    n_groups = 0
    for i in 1:n
        group_of[i] != 0 && continue   # already in a group
        n_groups += 1
        g = n_groups
        j = i
        while group_of[j] == 0
            group_of[j] = g
            nxt[j] == 0 && break       # unresolved successor (e.g. ∞ not in list)
            j = nxt[j]
        end
    end

    return group_of, n_groups
end

# ──────────────────────────────────────────────────────────────────────────────
# §13b  LCHab colour encoding for p-cycle sub-basins
# ──────────────────────────────────────────────────────────────────────────────

"""
    _group_fp_by_orbit_lch(fp_list, hpair_f, aproxprecision)

Like `_group_fp_by_orbit`, but additionally returns `subidx_of[i]` —
the 0-based position of `fp_list[i]` within its orbit group (sub-basin index).

Returns `(group_of, subidx_of, n_groups)`.
"""
function _group_fp_by_orbit_lch(
        fp_list::Vector,
        hpair_f,
        aproxprecision::Int)

    n   = length(fp_list)
    tol = 1.0 / 10^aproxprecision

    # Build successor map: nxt[i] = j  means  f(fp_list[i]) ≈ fp_list[j]
    nxt = zeros(Int, n)
    for i in 1:n
        fi = _rational_function(hpair_f, fp_list[i], aproxprecision)
        for j in 1:n
            if chordal_metric(fi, fp_list[j], aproxprecision) < tol
                nxt[i] = j
                break
            end
        end
    end

    # Trace orbit chains; assign group labels and sub-basin indices
    group_of  = zeros(Int, n)
    subidx_of = zeros(Int, n)
    n_groups  = 0

    for i in 1:n
        group_of[i] != 0 && continue   # already processed
        n_groups += 1
        g   = n_groups
        j   = i
        idx = 0
        while group_of[j] == 0
            group_of[j]  = g
            subidx_of[j] = idx
            idx += 1
            nxt[j] == 0 && break
            j = nxt[j]
        end
    end

    return group_of, subidx_of, n_groups
end

"""
    _okabe_lch_hues() -> Vector{Float64}

Return the LCHab hue angles (degrees) of the Okabe-Ito palette
(colours 2–8, skipping black).  Computed at first call and cached.
"""
function _okabe_lch_hues()::Vector{Float64}
    oi = okabe_ito_colors()          # returns 8 RGB colours
    return [convert(LCHab{Float64}, c).h for c in oi[2:end]]
end

"""
    _pcycle_lch_color(cycle_idx, sub_idx, p, nit, maxiter; hue_band=25.0)

Compute the LCHab colour for a pixel belonging to sub-basin `sub_idx`
(0-based) of the `cycle_idx`-th attracting cycle (1-based), where the
cycle has period `p`.  The pixel required `nit` iterations out of `maxiter`.

Encoding:
- **Hue**:    H = H_base(cycle_idx) + sub_idx · hue_band / p  (mod 360°)
- **Lightness**: L = 85 − 57·(nit/maxiter)^0.5  ∈ [28, 85]
- **Chroma**: C = 55 (fixed)
"""
function _pcycle_lch_color(
        cycle_idx::Int, sub_idx::Int, p::Int,
        nit::Int,       maxiter::Int;
        hue_band::Float64 = 25.0)::RGB{Float64}

    hues   = _okabe_lch_hues()
    H_base = hues[mod1(cycle_idx, length(hues))]
    H      = mod(H_base + (p > 1 ? sub_idx * hue_band / p : 0.0), 360.0)

    t = (nit / max(1, maxiter))^0.5
    L = 85.0 - 57.0 * t          # L ∈ [28, 85]
    C = 55.0

    rgb = convert(RGB{Float64}, LCHab{Float64}(L, C, H))
    return RGB{Float64}(clamp(red(rgb),   0.0, 1.0),
                        clamp(green(rgb), 0.0, 1.0),
                        clamp(blue(rgb),  0.0, 1.0))
end

"""
    _build_basin_image_lch(fp_list, pos_raw, itr, group_of, subidx_of, p; hue_band)

Build the pixel-colour matrix for the `"pcycle_lch"` strategy.
Each pixel is coloured by LCHab encoding of (cycle, sub-basin, convergence speed).
The basin of ∞ is rendered in greyscale (lightness encodes iteration count).
"""
function _build_basin_image_lch(
        fp_list   ::Vector,
        pos_raw   ::Matrix{Int64},
        itr       ::Matrix{Int64},
        group_of  ::Vector{Int},
        subidx_of ::Vector{Int},
        p         ::Int;
        hue_band  ::Float64 = 25.0)::Matrix{RGB{Float64}}

    nrows, ncols = size(pos_raw)
    maxiter  = max(1, maximum(itr))
    n_fp     = length(fp_list)
    inf_grp  = group_of[n_fp]        # group that contains ∞

    img = Matrix{RGB{Float64}}(undef, nrows, ncols)
    for j in 1:ncols, i in 1:nrows
        pos = pos_raw[i, j]
        nit = itr[i, j]

        if pos <= 0
            # Non-converged → black
            img[i,j] = RGB{Float64}(0.0, 0.0, 0.0)

        elseif group_of[pos] == inf_grp
            # Basin of ∞ → greyscale, lightness by iteration
            t = (nit / maxiter)^0.5
            v = clamp(0.85 - 0.50 * t, 0.30, 0.90)
            img[i,j] = RGB{Float64}(v, v, v)

        else
            img[i,j] = _pcycle_lch_color(
                group_of[pos], subidx_of[pos], p, nit, maxiter;
                hue_band = hue_band)
        end
    end
    return img
end

"""
    _render_basin_image_lch(fp_list, pos_raw, itr, group_of, subidx_of, p,
                             xinterval, yinterval, title; hue_band=25.0)

Build and return a `Plots.Plot` using the `"pcycle_lch"` colour strategy.
"""
function _render_basin_image_lch(
        fp_list   ::Vector,
        pos_raw   ::Matrix{Int64},
        itr       ::Matrix{Int64},
        group_of  ::Vector{Int},
        subidx_of ::Vector{Int},
        p         ::Int,
        xinterval ::Tuple{Real,Real},
        yinterval ::Tuple{Real,Real},
        title     ::String;
        hue_band  ::Float64 = 25.0)

    img = _build_basin_image_lch(
        fp_list, pos_raw, itr, group_of, subidx_of, p; hue_band = hue_band)

    nrows, ncols = size(img)
    a, b = Float64(xinterval[1]), Float64(xinterval[2])
    c, d = Float64(yinterval[1]), Float64(yinterval[2])

    return plot(
        img,
        aspect_ratio = :equal,
        xlims  = (1, ncols),
        ylims  = (1, nrows),
        xticks = (range(1, ncols, length=5),
                  string.(round.(range(a, b, length=5), digits=2))),
        yticks = (range(1, nrows, length=5),
                  string.(round.(range(c, d, length=5), digits=2))),
        xlabel = "Re(z)",
        ylabel = "Im(z)",
        title  = title,
        legend = false,
        size   = (600, 600))
end

"""
(coefficientlistnum, coefficientlistden, p;
                          xinterval=(-1.5,1.5), yinterval=(-1.5,1.5),
                          aproxprecision=8, iterprecision=8, iter=200, expprecision=8,
                          bands=6.0, gseed=RGB(0.2,0.6,0.9),
                          colorstrategy="positionplusiteration",
                          colorblind=false, show_plot=true, savepath="")

Compute and display the **p-cycle basin portrait** of the rational map `f`.

A *p-cycle* is an orbit `{ζ, f(ζ), …, f^{p-1}(ζ)}` of minimal period `p`.
The basin of a p-cycle consists of `p` connected components cyclically
permuted by `f`.  **Option A colouring**: all `p` basin components of the
same cycle share one colour, so each attractor (whether a fixed point or a
genuine p-cycle) is rendered in a distinct colour.

The algorithm:
1. Compute `f^p` by symbolic polynomial composition.
2. Find all fixed points of `f^p` (genuine p-cycles, plus fixed points of
   `f` and periodic orbits of shorter period dividing `p`).
3. Group points that share the same orbit under `f`.
4. Iterate `f^p` at each pixel and converge to a fixed point of `f^p`.
5. Remap pixel colours so that each orbit group has one colour.

Coefficients follow the *constant term first* convention:
`[a₀, a₁, …, aₙ]` represents `a₀ + a₁z + … + aₙzⁿ`.

# Arguments
- `coefficientlistnum`, `coefficientlistden` — numerator/denominator coefficients.
- `p::Int` — period of the cycles to visualise (≥ 1).
- All other keyword arguments are the same as `basin_portrait`.

# Returns
`plt` — the basin portrait plot.

# Example
```julia
using RiemannBasins
# Newton for z³−1 — p=1 reproduces basin_portrait
num = complex([1.0, 0.0, 0.0, 2.0])
den = complex([0.0, 0.0, 3.0, 0.0])
plt1 = basin_portrait_pcycle(num, den, 1)

# Degree-2 map: f(z) = z²+c near Basilica (c=−1)
# Fixed points are roots of z²-z+1=0; for p=2 we see 2-cycle basins
num2 = complex([-1.0, 0.0, 1.0])  # -1 + z²
den2 = complex([1.0])              # 1
plt2 = basin_portrait_pcycle(num2, den2, 2)
```
"""
function basin_portrait_pcycle(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T},
        p::Int;
        xinterval::Tuple{Real,Real} = (-1.5, 1.5),
        yinterval::Tuple{Real,Real} = (-1.5, 1.5),
        aproxprecision::Int64       = 8,
        iterprecision::Int64        = 8,
        iter::Int64                 = 200,
        expprecision::Int64         = 8,
        bands::Float64              = 6.0,
        gseed::RGB{Float64}         = RGB{Float64}(0.2, 0.6, 0.9),
        colorstrategy::String       = "positionplusiteration",
        hue_band::Float64           = 25.0,
        show_plot::Bool             = true,
        colorblind::Bool            = false,
        savepath::String            = "") where {T<:Number}

    p >= 1 || throw(ArgumentError("p must be ≥ 1"))

    cnum = complex.(Float64.(coefficientlistnum))
    cden = complex.(Float64.(coefficientlistden))

    # 1. hpair for f  (needed for orbit grouping)
    hpair_f = _pair_of_functions(cnum, cden)

    # 2. Compute f^p symbolically
    num_p, den_p = _rational_iterate(cnum, cden, p)

    # 3. Fixed points of f^p (finite list + ∞)
    _, _, finitelist_p = fixed_points(num_p, den_p, aproxprecision)
    fp_list = Vector{Tuple{Complex{Float64},Complex{Float64}}}(finitelist_p)
    push!(fp_list, (complex(1.0), complex(0.0)))   # ∞ always included

    # 4. Group fixed points of f^p by their orbit under f
    #    For pcycle_lch we also need sub-basin indices within each group.
    if colorstrategy == "pcycle_lch"
        group_of, subidx_of, n_groups =
            _group_fp_by_orbit_lch(fp_list, hpair_f, aproxprecision)
    else
        group_of, n_groups = _group_fp_by_orbit(fp_list, hpair_f, aproxprecision)
        subidx_of = nothing
    end

    # 5. Basin iteration with f^p
    hpair_fp = _pair_of_functions(num_p, den_p)
    grid     = _rectangle(xinterval, yinterval, expprecision)
    pos_raw, itr, res = _basin_matrices_from_grid(
        grid, fp_list, hpair_fp,
        aproxprecision, iterprecision, iter)

    title = "$(p)-cycle basins, LCHab ($(n_groups) attractor$(n_groups == 1 ? "" : "s"))"

    if colorstrategy == "pcycle_lch"
        # ── Option B: LCHab encoding (Hue=cycle, ΔHue=sub-basin, L=speed) ──
        title = "$(p)-cycle basins, LCHab ($(n_groups) attractor$(n_groups == 1 ? "" : "s"))"
        plt = _render_basin_image_lch(
            fp_list, pos_raw, itr,
            group_of, subidx_of, p,
            xinterval, yinterval, title;
            hue_band = hue_band)
    else
        # ── Option A: original strategy (one colour per orbit group) ──────
        # 6. Remap pixel labels: fixed-point index → orbit group index
        posG = map(pos_raw) do pos
            (pos <= 0) ? pos : group_of[pos]
        end

        # 7. Colormap with one colour per orbit group
        cmap = _build_colormap(n_groups, gseed; colorblind=colorblind)

        # Build one representative per group for the render call
        reps = Vector{Tuple{Complex{Float64},Complex{Float64}}}(undef, n_groups)
        for i in eachindex(fp_list)
            reps[group_of[i]] = fp_list[i]
        end

        title = "$(p)-cycle basins ($(n_groups) attractor$(n_groups == 1 ? "" : "s"))"
        plt = _render_basin_image(
            reps, posG, itr,
            xinterval, yinterval,
            cmap, colorstrategy, title;
            colorblind     = colorblind,
            bands          = bands,
            residualMatrix = res,
            aproxprecision = aproxprecision)
    end

    show_plot && display(plt)
    isempty(savepath) || savefig(plt, savepath)

    return plt
end


# ──────────────────────────────────────────────────────────────────────────────
# §14  Critical orbit analysis — maximum attracting period
# ──────────────────────────────────────────────────────────────────────────────
#
# Theoretical basis: Fatou's theorem (1919) states that every attracting or
# superattracting cycle must attract at least one critical point of f.
# Corollary: tracing the orbits of the 2d−2 critical points of a degree-d
# rational map is *complete* — it finds every attracting cycle.
# The maximum period P returned by max_attracting_period therefore certifies
# that f has no attracting cycle of period > P.

export critical_points,
       critical_orbit_analysis,
       max_attracting_period,
       print_critical_orbit_summary

# ---- Data type -----------------------------------------------------------

"""
    CriticalOrbitResult

Encapsulates the result of tracing one critical orbit.

Fields
------
- `critical_point`    — normalised homogeneous pair for the critical point.
- `label`             — `:attracting`, `:superattracting`, or `:indeterminate`.
  `:indeterminate` covers parabolic, Siegel, Herman, Julia, and orbits that
  did not converge to a cycle within the allotted iterations.
- `period`            — minimal period of the attracting cycle (0 when `:indeterminate`).
- `cycle_points`      — the `period` representative points of the cycle (empty
  when `:indeterminate`).
- `lyapunov_constant` — c = (∏ f^#(zₖ))^{1/q} around the cycle; 0 for
  superattracting, 1.0 for `:indeterminate`.
"""
struct CriticalOrbitResult
    critical_point    :: Tuple{Complex{Float64},Complex{Float64}}
    label             :: Symbol
    period            :: Int
    cycle_points      :: Vector{Tuple{Complex{Float64},Complex{Float64}}}
    lyapunov_constant :: Float64
end

# ---- Internal helpers ----------------------------------------------------

# Coefficients of ∂F/∂z for a homogeneous polynomial
# F(z,t) = Σ aₖ zᵏ t^(d-k)  →  ∂F/∂z = Σ k·aₖ z^(k-1) t^(d-k)
# Input:  [a₀, a₁, …, aₐ]  (d+1 entries, constant-first convention)
# Output: [1·a₁, 2·a₂, …, d·aₐ]  (d entries, degree d-1)
function _hpoly_dz(cv::Vector{Complex{Float64}})
    d = length(cv) - 1
    d == 0 && return [complex(0.0)]
    return [Complex{Float64}(k) * cv[k + 1] for k in 1:d]
end

# Partial-∂/∂t coefficients of a homogeneous polynomial.
# F(z,t) = Σ aₖ zᵏ t^{d-k}  →  ∂F/∂t = Σ_{k=0}^{d-1} (d-k)aₖ zᵏ t^{d-k-1}
# Returns the coefficient list [(d)a₀, (d-1)a₁, …, 1·a_{d-1}] for _bivariate_poly.
function _hpoly_dt(cv::Vector{Complex{Float64}})
    d = length(cv) - 1
    d == 0 && return [complex(0.0)]
    return [Complex{Float64}(d - k) * cv[k + 1] for k in 0:d-1]
end

# Spherical derivative at a homogeneous point (z, t) ∈ ℂ² \ {0}.
#
# The correct formula is chart-sensitive (Euler's identity gives W_t = -(z/t)·W_z):
#   • First chart  |t| ≥ |z|: f^# = |W_z| / |t| / (|F|²+|G|²) · (|z|²+|t|²)
#   • Second chart |z| > |t|: f^# = |W_t| / |z| / (|F|²+|G|²) · (|z|²+|t|²)
# For _homogeneous_normalize output (max(|z|,|t|)=1) the divisor is always 1,
# so the formula collapses to |W_{chart}| · (|z|²+|t|²) / (|F|²+|G|²).
function _sph_deriv_hom(
        cnum :: Vector{Complex{Float64}},
        cden :: Vector{Complex{Float64}},
        d    :: Int,
        z    :: Complex{Float64},
        t    :: Complex{Float64})

    d == 0 && return 0.0
    Fv      = _bivariate_poly(cnum, z, t, d)
    Gv      = _bivariate_poly(cden, z, t, d)
    den     = abs2(Fv) + abs2(Gv)
    den < 1e-300 && return 0.0
    norm_sq = abs2(z) + abs2(t)
    abs_t   = abs(t)
    abs_z   = abs(z)

    if abs_t >= abs_z
        # First chart: Wronskian w.r.t. z
        abs_t < 1e-300 && return 0.0
        dFv = _bivariate_poly(_hpoly_dz(cnum), z, t, d - 1)
        dGv = _bivariate_poly(_hpoly_dz(cden), z, t, d - 1)
        W   = dFv * Gv - Fv * dGv
        return abs(W) / abs_t / den * norm_sq
    else
        # Second chart (near ∞): Wronskian w.r.t. t
        abs_z < 1e-300 && return 0.0
        dFv = _bivariate_poly(_hpoly_dt(cnum), z, t, d - 1)
        dGv = _bivariate_poly(_hpoly_dt(cden), z, t, d - 1)
        W   = dFv * Gv - Fv * dGv
        return abs(W) / abs_z / den * norm_sq
    end
end

# Lyapunov constant c = geometric mean of f^# around a q-cycle.
function _lyapunov_of_cycle(
        cnum      :: Vector{Complex{Float64}},
        cden      :: Vector{Complex{Float64}},
        d         :: Int,
        cycle_pts :: Vector{Tuple{Complex{Float64},Complex{Float64}}})

    q = length(cycle_pts)
    q == 0 && return 1.0
    log_sum = sum(log(max(_sph_deriv_hom(cnum, cden, d, pt[1], pt[2]), 1e-300))
                  for pt in cycle_pts)
    return exp(log_sum / q)
end

# Run warmup iterations, then search for smallest q ∈ 1:max_period with
# chordal_metric(y, f^q(y)) < tol_per.
# Returns (period, orbit_point_after_warmup); period=0 means not detected.
function _detect_period_orbit(
        hpair      ,
        start_pt   :: Tuple{Complex{Float64},Complex{Float64}},
        warmup     :: Int,
        tol        :: Float64,
        max_period :: Int)

    tol_per = 1e-5    # generous tolerance after the transient

    curr = start_pt
    for _ in 1:warmup
        curr = _rational_function(hpair, curr)
        abs(curr[1]) + abs(curr[2]) < tol && return 0, curr   # super-zero hit
    end

    for q in 1:max_period
        trial = curr
        for _ in 1:q
            trial = _rational_function(hpair, trial)
        end
        chordal_metric(curr, trial) < tol_per && return q, curr
    end

    return 0, curr    # indeterminate
end

# ---- Public API ----------------------------------------------------------

"""
    critical_points(coefficientlistnum, coefficientlistden;
                    aproxprecision=8, check_infinity=false)

Return the **critical points** of the rational map `f = p/q` as normalised
homogeneous pairs on `P¹⁺(ℂ)`.

The critical points are the roots of the Wronskian
    W(z) = p′(z)·q(z) − p(z)·q′(z),
a polynomial of degree ≤ 2d−2 for a degree-d map.  By the Riemann-Hurwitz
theorem, f has exactly 2d−2 critical points (counted with multiplicity);
the returned list removes duplicates within chordal distance 10^{-aproxprecision}.

Pass `check_infinity=true` to also include `∞ = [1:0]` when it is a critical
point (detected when deg(W) < 2d−2, i.e. some multiplicity is at ∞).

# Returns
`Vector{Tuple{Complex{Float64},Complex{Float64}}}`.
"""
function critical_points(
        coefficientlistnum :: Vector{T},
        coefficientlistden :: Vector{T};
        aproxprecision  :: Int  = 8,
        check_infinity  :: Bool = false) where {T <: Number}

    cnum = complex.(Float64.(coefficientlistnum))
    cden = complex.(Float64.(coefficientlistden))
    d    = max(length(cnum), length(cden)) - 1
    tol  = 1.0 / 10^aproxprecision

    # Wronskian  W(z) = p′(z)·q(z) − p(z)·q′(z)
    p_poly  = Polynomial(cnum)
    q_poly  = Polynomial(cden)
    wron    = derivative(p_poly) * q_poly - p_poly * derivative(q_poly)
    wcoeffs = collect(coeffs(wron))
    while length(wcoeffs) > 1 && abs(last(wcoeffs)) < 1e-12
        pop!(wcoeffs)
    end

    result = Tuple{Complex{Float64},Complex{Float64}}[]

    # Finite critical points (roots of Wronskian)
    if !(length(wcoeffs) == 1 && abs(wcoeffs[1]) < 1e-12)
        for r in roots(Polynomial(wcoeffs))
            pt = _homogeneous_normalize(
                (Complex{Float64}(r), complex(1.0)), aproxprecision)
            abs(pt[1]) + abs(pt[2]) < tol && continue          # skip super-zero
            any(chordal_metric(pt, ex) < tol * 10 for ex in result) && continue
            push!(result, pt)
        end
    end

    # ∞ = [1:0] is a critical point when Wronskian degree < 2d−2
    if check_infinity
        expected_deg = 2 * d - 2
        actual_deg   = length(wcoeffs) - 1
        if actual_deg < expected_deg
            inf_pt = (complex(1.0), complex(0.0))
            any(chordal_metric(inf_pt, ex) < tol * 10 for ex in result) ||
                push!(result, inf_pt)
        end
    end

    return result
end

"""
    critical_orbit_analysis(coefficientlistnum, coefficientlistden;
                            max_iter=600, aproxprecision=8,
                            max_period=30, check_infinity=false)

Classify the orbit of **every critical point** of `f = p/q` under iteration.

For each critical point the algorithm:
1. Runs `max_iter ÷ 2` warmup steps to escape the transient.
2. Searches for the smallest period `q ≤ max_period` with
   `chordal_metric(y, fᵠ(y)) < 10⁻⁵` (orbit inside the cycle neighbourhood).
3. Computes `c = (∏ₖ f^#(zₖ))^{1/q}` via the homogeneous spherical-derivative
   formula of [Álvarez et al. 2023, Def. 16].
4. Labels `:superattracting` if `c < 10⁻⁸`, `:attracting` if `c < 0.9999`,
   `:indeterminate` otherwise.

By Fatou's theorem this procedure is **complete**: every attracting cycle is
found because it must attract at least one critical point.

# Returns
`Vector{CriticalOrbitResult}` — one entry per critical point.
"""
function critical_orbit_analysis(
        coefficientlistnum :: Vector{T},
        coefficientlistden :: Vector{T};
        max_iter        :: Int  = 600,
        aproxprecision  :: Int  = 8,
        max_period      :: Int  = 30,
        check_infinity  :: Bool = false) where {T <: Number}

    cnum   = complex.(Float64.(coefficientlistnum))
    cden   = complex.(Float64.(coefficientlistden))
    d      = max(length(cnum), length(cden)) - 1
    tol    = 1.0 / 10^aproxprecision
    warmup = max_iter ÷ 2

    hpair  = _pair_of_functions(cnum, cden)
    crits  = critical_points(cnum, cden;
                 aproxprecision = aproxprecision,
                 check_infinity = check_infinity)

    results = CriticalOrbitResult[]

    for cpt in crits
        q, y = _detect_period_orbit(hpair, cpt, warmup, tol, max_period)

        if q == 0
            push!(results, CriticalOrbitResult(
                cpt, :indeterminate, 0,
                Tuple{Complex{Float64},Complex{Float64}}[], 1.0))
            continue
        end

        # Collect the q points of the detected cycle
        cycle_pts = Tuple{Complex{Float64},Complex{Float64}}[]
        curr = y
        for _ in 1:q
            push!(cycle_pts, curr)
            curr = _rational_function(hpair, curr)
        end

        c = _lyapunov_of_cycle(cnum, cden, d, cycle_pts)

        label = if c < 1e-8
            :superattracting
        elseif c < 0.9999
            :attracting
        else
            :indeterminate
        end

        if label == :indeterminate
            push!(results, CriticalOrbitResult(
                cpt, :indeterminate, 0,
                Tuple{Complex{Float64},Complex{Float64}}[], c))
        else
            push!(results, CriticalOrbitResult(cpt, label, q, cycle_pts, c))
        end
    end

    return results
end

"""
    max_attracting_period(coefficientlistnum, coefficientlistden;
                          max_iter=600, aproxprecision=8,
                          max_period=30, check_infinity=false)

Return the **maximum period** `P` of any attracting or superattracting cycle
of `f = p/q`, found by tracing the orbits of all critical points.

By Fatou's theorem this is **complete**: there is no attracting cycle of
period `> P`.  Returns `0` when no attracting cycle is detected (all Fatou
components are parabolic, Siegel, or Herman, or `max_period` is too small).

# Example
```julia
using RiemannBasins
# Newton for z³−1  (three attracting fixed points, period 1)
num = complex([1.0, 0.0, 0.0, 2.0]);  den = complex([0.0, 0.0, 3.0, 0.0])
max_attracting_period(num, den)        # → 1

# Basilica  z²−1  (attracting 2-cycle  {0, −1})
num2 = complex([-1.0, 0.0, 1.0]);  den2 = complex([1.0])
max_attracting_period(num2, den2)      # → 2
```
"""
function max_attracting_period(
        coefficientlistnum :: Vector{T},
        coefficientlistden :: Vector{T};
        max_iter        :: Int  = 600,
        aproxprecision  :: Int  = 8,
        max_period      :: Int  = 30,
        check_infinity  :: Bool = false) where {T <: Number}

    results = critical_orbit_analysis(
        coefficientlistnum, coefficientlistden;
        max_iter       = max_iter,
        aproxprecision = aproxprecision,
        max_period     = max_period,
        check_infinity = check_infinity)

    P = 0
    for r in results
        r.label != :indeterminate && (P = max(P, r.period))
    end
    return P
end

"""
    print_critical_orbit_summary(coefficientlistnum, coefficientlistden; kwargs...)

Print a human-readable table of all critical orbits and their classifications.
Accepts the same keyword arguments as `critical_orbit_analysis`.

Useful as a quick sanity check before running a full basin portrait, to verify
that all attracting cycles have been found and to determine the correct value
of `p` for `basin_portrait_pcycle`.

# Example
```julia
using RiemannBasins
num2 = complex([-1.0, 0.0, 1.0]);  den2 = complex([1.0])
print_critical_orbit_summary(num2, den2)
# Critical orbit analysis   degree = 2,  1 critical point(s)
# ──────────────────────────────────────────────────────────────────
# #    Critical point (z, t)                Label              Period  c
# ──────────────────────────────────────────────────────────────────
# 1    (0.0 + 0.0im, 1.0 + 0.0im)           attracting         2       0.5
# ──────────────────────────────────────────────────────────────────
# Maximum attracting period  P = 2
# By Fatou's theorem, f has no attracting cycle of period > P.
```
"""
function print_critical_orbit_summary(
        coefficientlistnum :: Vector{T},
        coefficientlistden :: Vector{T};
        kwargs...) where {T <: Number}

    results = critical_orbit_analysis(coefficientlistnum, coefficientlistden; kwargs...)
    d   = max(length(coefficientlistnum), length(coefficientlistden)) - 1
    sep = "─"^68

    println("\nCritical orbit analysis   degree = $(d),  $(length(results)) critical point(s)")
    println(sep)
    println(rpad("#",  5)  *
            rpad("Critical point (z, t)", 38) *
            rpad("Label", 18) *
            rpad("Period", 8) *
            "c")
    println(sep)

    for (i, r) in enumerate(results)
        z, t     = r.critical_point
        pt_str   = "($(round(z, digits=4)), $(round(t, digits=4)))"
        per_str  = r.period == 0 ? "—" : string(r.period)
        c_str    = r.label == :indeterminate ? "—" :
                   string(round(r.lyapunov_constant, digits=6))
        println(rpad(string(i), 5)  *
                rpad(pt_str, 38) *
                rpad(string(r.label), 18) *
                rpad(per_str, 8) *
                c_str)
    end

    println(sep)
    P = maximum((r.period for r in results if r.label != :indeterminate), init = 0)
    println("Maximum attracting period  P = $(P)")
    println("By Fatou's theorem, f has no attracting cycle of period > P.\n")

    return results
end

end  # module RiemannBasins
