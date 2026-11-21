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
using Base64              # base64 encoding for HTML export

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
       basin_portrait_pcycle,
       mandelbrot_portrait,
       basin_sphere_data,
       basin_measure,
       basin_sphere_portrait_3d,
       basin_sphere_portrait_html

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
    return _bivariate_poly(ComplexF64.(coefficientlist), u, t, d)
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

    cnum = ComplexF64.(coefficientlistnum)
    cden = ComplexF64.(coefficientlistden)

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
            string.(round.(range(d, c, length=5), digits=2))
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
# ──────────────────────────────────────────────────────────────────────────────
# §  Mandelbrot set  –  parameter portrait for the family z² + c
# ──────────────────────────────────────────────────────────────────────────────

"""
    mandelbrot_portrait(; re_range, im_range, n_re, n_im, maxiter, dpi,
                          marked_cs, title, show_plot, savepath)

Compute and display the Mandelbrot set for the quadratic family  z ↦ z² + c.

For each pixel c the orbit of the critical point 0 is iterated up to `maxiter`
steps.  Points for which the orbit does not escape (|z| ≤ 2) are coloured
near-black (inside M).  Exterior points are coloured using a perceptually
uniform LCH ramp based on the smooth escape-time (normalized iteration count),
consistent with the LCHab style used in `basin_portrait_pcycle`.

# Arguments
- `re_range`   – real-axis extent, default `(-2.5, 1.0)`
- `im_range`   – imaginary-axis extent, default `(-1.25, 1.25)`
- `n_re`       – horizontal pixel count, default `800`
- `n_im`       – vertical pixel count, default `600`
- `maxiter`    – maximum iterations, default `256`
- `dpi`        – plot resolution, default `150`
- `marked_cs`  – vector of `(c::ComplexF64, label::String)` pairs to overlay
                  as labelled yellow dots (the four canonical examples by default)
- `title`      – plot title
- `show_plot`  – display the plot interactively (default `true`)
- `savepath`   – if non-empty, save the figure to this path

# Examples
```julia
# Full default view with the four canonical c values marked
mandelbrot_portrait()

# Zoom into the main cardioid region, no marked points
mandelbrot_portrait(re_range=(-0.8, 0.4), im_range=(-0.65, 0.65),
                    marked_cs=[], n_re=1200, n_im=1200)

# High-resolution for a poster
mandelbrot_portrait(n_re=3000, n_im=2250, dpi=300,
                    savepath="poster_mandelbrot.png")
```
"""
function mandelbrot_portrait(;
        re_range   ::Tuple{Real,Real} = (-2.5,  1.0),
        im_range   ::Tuple{Real,Real} = (-1.25, 1.25),
        n_re       ::Int    = 800,
        n_im       ::Int    = 600,
        maxiter    ::Int    = 256,
        dpi        ::Int    = 150,
        marked_cs          = [
            (-1.0 + 0.0im,            "Basílica"),
            (-1.7549 + 0.0im,         "Airplane"),
            (-0.1226 + 0.7449im,      "Rabbit"),
            (-0.50434 + 0.56277im,    "Period-5"),
        ],
        title      ::String = "Mandelbrot set  –  z² + c",
        show_plot  ::Bool   = true,
        savepath   ::String = "")

    re_vals = LinRange(Float64(re_range[1]), Float64(re_range[2]), n_re)
    im_vals = LinRange(Float64(im_range[1]), Float64(im_range[2]), n_im)

    # ── 1. Iteration ─────────────────────────────────────────────────────────
    # nit[i,j] = smooth escape count for pixel (i row from top, j col)
    #          = NaN if inside M
    nit = fill(NaN, n_im, n_re)

    for j in 1:n_re
        c = complex(re_vals[j], 0.0)   # will add im part per row
        for i in 1:n_im
            c_pix = complex(re_vals[j], im_vals[n_im + 1 - i])
            z     = zero(ComplexF64)
            for n in 1:maxiter
                z = z * z + c_pix
                if abs2(z) > 4.0
                    # Smooth (normalized) iteration count
                    nit[i, j] = n - log2(max(1.0, log2(abs(z))))
                    break
                end
            end
        end
    end

    # ── 2. Colour mapping ────────────────────────────────────────────────────
    # Interior  → near-black
    # Exterior  → LCHab  H=220° (blue), C=55, L ∈ [20,85] by escape speed
    inside_col = RGB{Float64}(0.04, 0.04, 0.10)

    ext_vals = filter(!isnan, vec(nit))
    mu_min   = isempty(ext_vals) ? 0.0 : minimum(ext_vals)
    mu_max   = isempty(ext_vals) ? 1.0 : maximum(ext_vals)
    mu_span  = max(mu_max - mu_min, 1e-9)

    img = Matrix{RGB{Float64}}(undef, n_im, n_re)
    for i in 1:n_im, j in 1:n_re
        if isnan(nit[i, j])
            img[i, j] = inside_col
        else
            t = clamp((nit[i, j] - mu_min) / mu_span, 0.0, 1.0)
            L = 85.0 - 65.0 * t^0.4   # 85 (fast escape) → 20 (slow)
            img[i, j] = convert(RGB{Float64}, LCHab(L, 55.0, 220.0))
        end
    end

    # ── 3. Base plot ─────────────────────────────────────────────────────────
    plt = plot(
        img,
        aspect_ratio = :equal,
        xlims = (1, n_re),
        ylims = (1, n_im),
        title  = title,
        xlabel = "Re(c)",
        ylabel = "Im(c)",
        dpi    = dpi,
        size   = (n_re, n_im);
        # Custom tick labels: convert pixel indices to coordinate values
        xticks = (
            range(1, n_re; length = 5),
            [string(round(v, digits=2)) for v in LinRange(Float64(re_range[1]), Float64(re_range[2]), 5)]
        ),
        yticks = (
            range(1, n_im; length = 5),
            [string(round(v, digits=2)) for v in LinRange(Float64(im_range[1]), Float64(im_range[2]), 5)]
        ),
    )

    # ── 4. Marked parameter values ───────────────────────────────────────────
    if !isempty(marked_cs)
        dot_color = RGB{Float64}(1.0, 0.85, 0.0)   # golden yellow

        px(re) = 1 + (Float64(re) - Float64(re_range[1])) /
                     (Float64(re_range[2]) - Float64(re_range[1])) * (n_re - 1)
        py(im) = n_im - (Float64(im) - Float64(im_range[1])) /
                        (Float64(im_range[2]) - Float64(im_range[1])) * (n_im - 1)

        xs = [px(real(c)) for (c, _) in marked_cs]
        ys = [py(imag(c)) for (c, _) in marked_cs]

        scatter!(plt, xs, ys;
            markersize  = 7,
            markercolor = dot_color,
            markerstrokecolor = :black,
            markerstrokewidth = 1.5,
            label       = "",
        )

        for ((c, lbl), x, y) in zip(marked_cs, xs, ys)
            # Place label above the dot
            annotate!(plt, x, y - 18,
                text(lbl, :white, :center, 9))
        end
    end

    show_plot && display(plt)
    isempty(savepath) || savefig(plt, savepath)
    return plt
end


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
    _cyclic_palette_inset!(img, p, cycle_idx; hue_band, palette_size, margin)

Embed a circular colour-wheel legend into the top-right corner of `img` **in place**.

The wheel has `p` equal sectors, arranged counter-clockwise starting from the
right (+x axis): sector `k` (0-indexed) covers [k·2π/p, (k+1)·2π/p).
Going counter-clockwise the colours appear in the cyclic order
B₀ → B₁ → … → B_{p-1} → B₀, mirroring the orbit permutation f: Bₖ → Bₖ₊₁.

A dark border ring and a one-pixel white outer halo separate the wheel
from the basin image.
"""
function _cyclic_palette_inset!(
        img          ::Matrix{RGB{Float64}},
        p            ::Int,
        cycle_idx    ::Int  = 1;
        hue_band     ::Float64 = 25.0,
        palette_size ::Int  = 90,     # diameter in pixels
        margin       ::Int  = 6)      # gap from image border

    nrows, ncols = size(img)
    r  = palette_size ÷ 2            # radius (integer)
    r2 = Float64(r)                  # float radius

    # Top-right corner anchor
    ci = margin + r + 1              # row of centre  (small = near top)
    cj = ncols - margin - r          # col of centre  (near right edge)

    white  = RGB{Float64}(1.0, 1.0, 1.0)
    border = RGB{Float64}(0.15, 0.15, 0.15)  # dark grey rim

    for di in -(r+2):(r+2), dj in -(r+2):(r+2)
        ii = ci + di
        jj = cj + dj
        (1 <= ii <= nrows && 1 <= jj <= ncols) || continue

        d = sqrt(Float64(di)^2 + Float64(dj)^2)

        if d <= r2 - 1.5
            # Interior: determine sector from CCW angle
            # Image rows increase downward, so negate di for math y-axis
            angle_ccw = mod(atan(Float64(-di), Float64(dj)), 2π)
            sub_k = clamp(floor(Int, angle_ccw / (2π / p)), 0, p - 1)
            img[ii, jj] = _pcycle_lch_color(cycle_idx, sub_k, p, 0, 1;
                                            hue_band = hue_band)
        elseif d <= r2 + 0.5
            img[ii, jj] = border          # dark border ring
        elseif d <= r2 + 1.5
            img[ii, jj] = white           # white separation halo
        end
    end
    return img
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
        fp_list            ::Vector,
        pos_raw            ::Matrix{Int64},
        itr                ::Matrix{Int64},
        group_of           ::Vector{Int},
        subidx_of          ::Vector{Int},
        p                  ::Int,
        xinterval          ::Tuple{Real,Real},
        yinterval          ::Tuple{Real,Real},
        title              ::String;
        hue_band           ::Float64 = 25.0,
        show_cyclic_palette::Bool    = false)

    img = _build_basin_image_lch(
        fp_list, pos_raw, itr, group_of, subidx_of, p; hue_band = hue_band)

    # Optionally embed circular palette legend.
    # We detect the attracting finite cycle by finding the non-∞ group with the
    # most pixels (the ∞ group is at group_of[end] since ∞ is always last in fp_list).
    if show_cyclic_palette
        inf_grp = group_of[end]
        grp_counts = Dict{Int,Int}()
        for pos in pos_raw
            pos > 0 || continue
            g = group_of[pos]
            g == inf_grp && continue
            grp_counts[g] = get(grp_counts, g, 0) + 1
        end
        attr_grp = isempty(grp_counts) ? 1 : argmax(grp_counts)
        _cyclic_palette_inset!(img, p, attr_grp; hue_band = hue_band)
    end

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
                  string.(round.(range(d, c, length=5), digits=2))),
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
        show_cyclic_palette::Bool   = false,
        show_plot::Bool             = true,
        colorblind::Bool            = false,
        savepath::String            = "") where {T<:Number}

    p >= 1 || throw(ArgumentError("p must be ≥ 1"))

    cnum = ComplexF64.(coefficientlistnum)
    cden = ComplexF64.(coefficientlistden)

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
            hue_band            = hue_band,
            show_cyclic_palette = show_cyclic_palette)
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
- `label`             — `:attracting`, `:superattracting`, `:parabolic`,
  `:siegel`, `:herman`, or `:indeterminate`.
  `:indeterminate` covers Julia-set orbits and orbits that did not converge
  within the allotted iterations.
- `period`            — minimal period of the attracting or neutral cycle (0
  when not applicable).
- `cycle_points`      — the `period` representative points of the cycle (empty
  when not applicable).
- `lyapunov_constant` — c = (∏ f^#(zₖ))^{1/q} around the cycle; 0 for
  superattracting, 1.0 otherwise.
- `rotation_number`   — estimated Brjuno rotation number ρ ∈ (0,1) for
  `:siegel` and `:herman` labels; `NaN` otherwise.
"""
struct CriticalOrbitResult
    critical_point    :: Tuple{Complex{Float64},Complex{Float64}}
    label             :: Symbol
    period            :: Int
    cycle_points      :: Vector{Tuple{Complex{Float64},Complex{Float64}}}
    lyapunov_constant :: Float64
    rotation_number   :: Float64

    # 5-argument backward-compatible constructor (rotation_number = NaN)
    CriticalOrbitResult(cp, lbl, per, cyc, lyap) =
        new(cp, lbl, per, cyc, lyap, NaN)
    # 6-argument constructor (explicit rotation number)
    CriticalOrbitResult(cp, lbl, per, cyc, lyap, rho) =
        new(cp, lbl, per, cyc, lyap, rho)
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

# ---- Birkhoff average tools for parabolic detection ----------------------

"""
    _birkhoff_partial_averages(hpair, cnum, cden, d, start_pt, checkpoints)

Compute partial Birkhoff averages of `log(f^#)` along the orbit of `start_pt`.

Returns a vector `B` where
    `B[i] = (1/checkpoints[i]) ∑_{k=0}^{checkpoints[i]-1} log(f^#(orbit[k]))`
in the same order as the sorted `checkpoints`.

**Sign convention**: `B[i] < 0` when the orbit is attracted toward a cycle
(where `f^# < 1`).  For a **parabolic** basin `B[i] → 0⁻` at the slow rate
`∼ log(N)/N` (Birkhoff's ergodic theorem applied to the neutral measure on the
parabolic cycle).
"""
function _birkhoff_partial_averages(
        hpair       ,
        cnum        :: Vector{Complex{Float64}},
        cden        :: Vector{Complex{Float64}},
        d           :: Int,
        start_pt    :: Tuple{Complex{Float64},Complex{Float64}},
        checkpoints :: Vector{Int})

    isempty(checkpoints) && return Float64[]
    checkpoints = sort(checkpoints)
    N_max    = last(checkpoints)
    averages = Float64[]

    curr    = start_pt
    log_sum = 0.0
    ci      = 1          # checkpoint index

    for n in 1:N_max
        s        = _sph_deriv_hom(cnum, cden, d, curr[1], curr[2])
        log_sum += log(max(s, 1e-300))
        curr     = _rational_function(hpair, curr)

        if ci <= length(checkpoints) && n == checkpoints[ci]
            push!(averages, log_sum / n)
            ci += 1
        end
    end

    return averages
end

"""
    _aitken_accelerate(seq)

Apply one pass of Aitken's Δ² process to `seq`, returning an accelerated
sequence of length `length(seq) − 2`.

For a sequence converging as `sₙ = L + a·ρⁿ` with `|ρ| < 1`, the transform
recovers `L` exactly from three consecutive terms:
```
sₙ* = sₙ − (sₙ₊₁ − sₙ)² / (sₙ₊₂ − 2sₙ₊₁ + sₙ)
```
Applied to Birkhoff partial averages it removes the leading-order error term
and dramatically accelerates convergence in both the geometric (attracting)
and the logarithmic (parabolic) regimes.
"""
function _aitken_accelerate(seq::Vector{Float64})
    n = length(seq)
    n < 3 && return copy(seq)
    result = Vector{Float64}(undef, n - 2)
    for i in 1:(n-2)
        Δ   = seq[i+1] - seq[i]
        Δ2  = seq[i+2] - 2*seq[i+1] + seq[i]
        result[i] = abs(Δ2) < 1e-30 * (abs(seq[i]) + 1e-100) ? seq[i] :
                    seq[i] - Δ^2 / Δ2
    end
    return result
end

"""
    _classify_marginal(hpair, cnum, cden, d, start_pt, q, cycle_pts, c;
                       N_birkhoff=800)

Attempt to separate `:parabolic` from generic `:indeterminate` for a critical
orbit whose Lyapunov constant `c` lies above the attracting threshold `0.9999`.

**Algorithm**:
1. Compute four partial Birkhoff averages of `log(f^#)` at
   `N = N_birkhoff÷8, N÷4, N÷2, N`.
2. Apply Aitken's Δ² process to accelerate convergence.
3. Use three criteria jointly:
   - **Sign**: all four averages are strictly negative (orbit is contracting).
   - **Monotone growth**: averages increase toward zero (slow convergence).
   - **Rate**: the ratio `|B_{N÷4}| / |B_N|` exceeds 1.3, indicating
     sub-geometric (logarithmic) convergence characteristic of parabolic basins.
4. If all three hold and the Aitken-accelerated estimate is within `1e-5` of
   zero, label the component `:parabolic`.  If the accelerated estimate is
   instead clearly negative (`< −1e-5`), re-label as `:attracting` with the
   corrected Lyapunov constant `exp(accelerated_estimate)`.
   Otherwise return `:indeterminate`.

# Returns
`(label, c_refined)` — `label ∈ {:parabolic, :attracting, :indeterminate}`
and `c_refined` is the best available estimate of the Lyapunov constant.
"""
function _classify_marginal(
        hpair      ,
        cnum       :: Vector{Complex{Float64}},
        cden       :: Vector{Complex{Float64}},
        d          :: Int,
        start_pt   :: Tuple{Complex{Float64},Complex{Float64}},
        q          :: Int,
        cycle_pts  :: Vector{Tuple{Complex{Float64},Complex{Float64}}},
        c          :: Float64;
        N_birkhoff :: Int = 800)

    # Four logarithmically-spaced checkpoints
    cps = [N_birkhoff ÷ 8, N_birkhoff ÷ 4, N_birkhoff ÷ 2, N_birkhoff]
    B   = _birkhoff_partial_averages(hpair, cnum, cden, d, start_pt, cps)

    length(B) < 4 && return :indeterminate, c

    # Criterion 1: all Birkhoff averages negative (orbit is attracted, f^# < 1 on avg)
    all_negative = all(b < 0.0 for b in B)
    # Criterion 2: averages increase monotonically toward zero (not oscillating)
    monotone_up  = all(B[i] < B[i+1] for i in 1:length(B)-1)
    # Criterion 3: sub-geometric convergence rate (parabolic gives log(N)/N decay)
    #   For B_N ≈ A·log(N)/N:  |B_{N/4}| / |B_N| ≈ 4·log(N/4)/log(N) ≈ 2–4
    #   For geometric:          ratio ≈ 1
    ratio        = abs(B[end]) > 1e-15 ? abs(B[2]) / abs(B[end]) : 0.0

    # Aitken-accelerated estimate of the limit
    B_accel      = _aitken_accelerate(B)
    accel_limit  = isempty(B_accel) ? B[end] : B_accel[end]

    # Parabolic: three structural criteria plus Aitken limit not clearly < 0.
    # Note: Aitken applied to a log(N)/N sequence overshoots slightly positive,
    # so we only require accel_limit > -0.01 (not clearly attracting), not ≈ 0.
    if all_negative && monotone_up && ratio > 1.3 && accel_limit > -0.01
        return :parabolic, 1.0
    end

    # If Aitken reveals a clearly attracting cycle (c was barely above threshold)
    if accel_limit < -1e-5
        return :attracting, exp(accel_limit)
    end

    return :indeterminate, c
end

# ---- Siegel disk and Herman ring detection --------------------------------

"""
    _hom_to_affine(pt)

Convert a normalised homogeneous pair `(z, t)` to the affine complex number
`z/t`.  Returns `Inf + 0im` when the point is at infinity (`|t| ≪ |z|`).
"""
function _hom_to_affine(pt::Tuple{Complex{Float64},Complex{Float64}})
    z, t = pt
    abs(t) < 1e-13 * max(abs(z), 1.0) && return complex(Inf)
    return z / t
end

"""
    _check_brjuno(rho, k_max = 14)

Return `true` if the rotation number `rho ∈ (0,1)` appears to satisfy the
Brjuno condition  `∑_{n≥1} log(q_{n+1})/q_n < ∞`, estimated via the first
`k_max` partial quotients of its continued-fraction expansion.

The golden ratio `(√5−1)/2` and all numbers of bounded-type satisfy the
condition; Liouville-type numbers do not.
"""
function _check_brjuno(rho::Float64, k_max::Int = 14)
    isnan(rho)                       && return false
    rho < 1e-7 || rho > 1.0 - 1e-7  && return false   # essentially rational

    x       = rho
    q_prev  = Int64(0)   # q_{-1}
    q_curr  = Int64(1)   # q_0
    brjuno  = 0.0

    for _ in 1:k_max
        x < 1e-13 && return false         # orbit closed up → rational
        a      = floor(Int64, 1.0 / x)
        a == 0 && return false
        x      = 1.0 / x - a
        q_next = a * q_curr + q_prev

        # Brjuno sum term: log(q_{n+1}) / q_n
        brjuno += log(Float64(q_next)) / Float64(q_curr)
        brjuno > 50.0 && return false     # sum diverging → not Brjuno

        q_prev  = q_curr
        q_curr  = q_next
        q_curr > 10^14 && break           # overflow guard; sum converged
    end

    return brjuno < 30.0
end

"""
    _winding_number_discrete(orbit_affine, center)

Approximate the winding number of the closed discrete curve `orbit_affine`
around `center` by summing argument increments (divided by 2π).

The curve is closed by connecting `orbit_affine[end]` back to `orbit_affine[1]`.
"""
function _winding_number_discrete(
        orbit_affine :: Vector{Complex{Float64}},
        center       :: Complex{Float64})

    w = 0.0
    N = length(orbit_affine)
    for k in 1:N
        z_k   = orbit_affine[k]
        z_kp1 = orbit_affine[k < N ? k + 1 : 1]
        d_k   = z_k   - center
        d_kp1 = z_kp1 - center
        (abs(d_k) < 1e-10 || abs(d_kp1) < 1e-10) && continue
        w += angle(d_kp1 / d_k) / (2π)
    end
    return w
end

"""
    _classify_quasi_periodic(hpair, cnum, cden, d, cpt_after_warmup,
                              all_crits; N_orbit, eps_var)

Attempt to classify a critical orbit that did not converge to a periodic cycle
as either `:siegel` (Siegel disk), `:herman` (Herman ring), or `:indeterminate`.

Also returns the estimated rotation number `rho`.

**Algorithm**
1. Iterate the map `N_orbit` more steps from `cpt_after_warmup`.
2. Discard if orbit escapes to infinity.
3. Compute variance of the affine orbit; if too small the orbit is still
   converging and is not quasi-periodic.
4. Estimate rotation number by Birkhoff average of argument increments
   around the orbit centroid.
5. Verify the Brjuno condition via the continued-fraction expansion.
6. Distinguish Siegel/Herman by the winding-number signature of the orbit
   around the other finite critical points.
   - **Siegel**: all other critical points on the same side of the orbit
     (all winding numbers equal, hence all inside or all outside).
   - **Herman**: at least one critical point on each side (mixed winding
     numbers, indicating the orbit separates the plane nontrivially).
"""
function _classify_quasi_periodic(
        hpair            ,
        cnum             :: Vector{Complex{Float64}},
        cden             :: Vector{Complex{Float64}},
        d                :: Int,
        cpt_after_warmup :: Tuple{Complex{Float64},Complex{Float64}},
        all_crits        :: Vector{Tuple{Complex{Float64},Complex{Float64}}};
        N_orbit          :: Int     = 1200,
        eps_var          :: Float64 = 1e-6)

    # 1. Run N_orbit additional iterates
    curr    = cpt_after_warmup
    orbit_h = Vector{Tuple{Complex{Float64},Complex{Float64}}}(undef, N_orbit)
    for k in 1:N_orbit
        curr       = _rational_function(hpair, curr)
        orbit_h[k] = curr
    end

    orbit_a = [_hom_to_affine(pt) for pt in orbit_h]

    # 2. Discard if orbit escaped to ∞
    any(isinf(real(z)) for z in orbit_a) && return :indeterminate, NaN

    # 3. Variance check: orbit must be spatially spread out
    n_a   = length(orbit_a)
    mu_z  = sum(orbit_a) / n_a
    var_z = sum(abs2(z - mu_z) for z in orbit_a) / n_a
    var_z < eps_var && return :indeterminate, NaN

    # 4. Estimate rotation number via Birkhoff average of argument increments
    psi_sum = 0.0
    psi_sq  = 0.0
    n_valid = 0
    for k in 1:N_orbit - 1
        d_k   = orbit_a[k]     - mu_z
        d_kp1 = orbit_a[k + 1] - mu_z
        (abs(d_k) < 1e-10 || abs(d_kp1) < 1e-10) && continue
        psi_k    = angle(d_kp1 / d_k) / (2π)   # ∈ (−0.5, 0.5]
        psi_sum += psi_k
        psi_sq  += psi_k^2
        n_valid += 1
    end

    n_valid < 20 && return :indeterminate, NaN

    psi_mean = psi_sum / n_valid
    psi_var  = psi_sq  / n_valid - psi_mean^2

    # Low variance of increments confirms quasi-periodicity (near-constant rotation)
    psi_var > 0.05 && return :indeterminate, NaN

    rho = mod(psi_mean, 1.0)
    abs(rho) < 1e-5 || abs(rho - 1.0) < 1e-5 && return :indeterminate, NaN

    # 5. Brjuno condition
    _check_brjuno(rho) || return :indeterminate, rho

    # 6. Siegel vs Herman: winding-number signature around other finite critical points
    cpt_a         = _hom_to_affine(cpt_after_warmup)
    other_crits_a = Complex{Float64}[]
    for c in all_crits
        ca = _hom_to_affine(c)
        isinf(real(ca))         && continue   # skip ∞
        abs(ca - cpt_a) < 1e-5 && continue   # skip self
        push!(other_crits_a, ca)
    end

    # Only one finite critical point (e.g. quadratic polynomial) → Siegel disk
    isempty(other_crits_a) && return :siegel, rho

    # Compute winding numbers of the orbit around each other critical point
    winds     = [_winding_number_discrete(orbit_a, c) for c in other_crits_a]
    wind_ints = round.(Int, winds)

    # Siegel: all winding numbers equal (all crits on the same side of the orbit)
    # Herman: mixed winding numbers (the orbit separates critical points)
    if all(w == wind_ints[1] for w in wind_ints)
        return :siegel, rho
    else
        return :herman, rho
    end
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

    cnum = ComplexF64.(coefficientlistnum)
    cden = ComplexF64.(coefficientlistden)
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
5. *(When `detect_parabolic=true`)* For marginal cycles (`c ≥ 0.9999`), applies
   Birkhoff average acceleration (`_classify_marginal`) to distinguish
   `:parabolic` components from quasi-periodic ones, then attempts to
   classify the latter as `:siegel` (Siegel disk) or `:herman` (Herman ring)
   via rotation-number estimation and winding-number tests.

By Fatou's theorem this procedure is **complete**: every attracting cycle is
found because it must attract at least one critical point.

# Returns
`Vector{CriticalOrbitResult}` — one entry per critical point.
Possible labels: `:superattracting`, `:attracting`, `:parabolic`,
`:siegel`, `:herman`, `:indeterminate`.
"""
function critical_orbit_analysis(
        coefficientlistnum :: Vector{T},
        coefficientlistden :: Vector{T};
        max_iter          :: Int  = 600,
        aproxprecision    :: Int  = 8,
        max_period        :: Int  = 30,
        check_infinity    :: Bool = false,
        detect_parabolic  :: Bool = true,
        N_birkhoff        :: Int  = 800) where {T <: Number}

    cnum   = ComplexF64.(coefficientlistnum)
    cden   = ComplexF64.(coefficientlistden)
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
            # Orbit did not converge to any periodic cycle within max_period.
            # This is the typical case for Siegel disks and Herman rings (the
            # critical orbit stays quasi-periodic rather than landing on a cycle).
            if detect_parabolic
                qp_lbl, rho = _classify_quasi_periodic(
                    hpair, cnum, cden, d, y, crits;
                    N_orbit = N_birkhoff)
                push!(results, CriticalOrbitResult(
                    cpt, qp_lbl, 0,
                    Tuple{Complex{Float64},Complex{Float64}}[], 1.0, rho))
            else
                push!(results, CriticalOrbitResult(
                    cpt, :indeterminate, 0,
                    Tuple{Complex{Float64},Complex{Float64}}[], 1.0))
            end
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
            # Attempt to separate parabolic from Siegel/Herman via Birkhoff averages
            if detect_parabolic
                label2, c2 = _classify_marginal(hpair, cnum, cden, d, cpt,
                                                q, cycle_pts, c;
                                                N_birkhoff = N_birkhoff)
                if label2 == :parabolic
                    push!(results, CriticalOrbitResult(cpt, :parabolic, q, cycle_pts, 1.0))
                elseif label2 == :attracting
                    push!(results, CriticalOrbitResult(cpt, :attracting, q, cycle_pts, c2))
                else
                    # _classify_marginal returned :indeterminate but a cycle WAS found.
                    # Could be a Siegel disk with a near-resonant period; try QP test.
                    qp_lbl, rho = _classify_quasi_periodic(
                        hpair, cnum, cden, d, y, crits;
                        N_orbit = N_birkhoff)
                    if qp_lbl !== :indeterminate
                        push!(results, CriticalOrbitResult(
                            cpt, qp_lbl, q, cycle_pts, c, rho))
                    else
                        push!(results, CriticalOrbitResult(
                            cpt, :indeterminate, 0,
                            Tuple{Complex{Float64},Complex{Float64}}[], c))
                    end
                end
            else
                push!(results, CriticalOrbitResult(
                    cpt, :indeterminate, 0,
                    Tuple{Complex{Float64},Complex{Float64}}[], c))
            end
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


# ──────────────────────────────────────────────────────────────────────────────
# §15  Spherical subdivision  —  basin_sphere_data, basin_measure
# ──────────────────────────────────────────────────────────────────────────────
#
# Reference: Gutiérrez, Hernández, Magreñán, Rivas (2016).
#   "Measures of the basins of attracting n-cycles for the relaxed Newton's method."
#
# The Riemann sphere S² is covered by  6 × 4^k  near-equal spherical
# quadrilaterals obtained by projecting the uniformly-subdivided faces of the
# circumscribed unit cube.  For each cell the orbit of its spherical centroid
# is iterated to identify the attractor; cell areas are approximated by
# 4π / (6 × 4^k), which converges to the true area as k → ∞.

"""
    SphereBasinData

Result of [`basin_sphere_data`](@ref).

Fields
------
- `points`           — cell-centre positions on S² as `NTuple{3,Float64}`.
- `basin_indices`    — attractor index for each cell (0 = not converged / Julia set,
                        k ≥ 1 = converged to `fixed_point_list[k]`).
- `fixed_point_list` — attractor list (homogeneous pairs); last entry is always
                        ∞ = `(1+0im, 0+0im)`.
- `cell_area`        — approximate spherical area per cell: `4π / total_cells`.
- `subdivision`      — subdivision level `k` used (`6 × 4^k` cells total).
- `total_cells`      — `6 × 4^subdivision`.
"""
struct SphereBasinData
    points           :: Vector{NTuple{3,Float64}}
    basin_indices    :: Vector{Int}
    fixed_point_list :: Vector{Tuple{Complex{Float64},Complex{Float64}}}
    cell_area        :: Float64
    subdivision      :: Int
    total_cells      :: Int
end

"""
    _cubic_sphere_points(subdivision) -> Vector{NTuple{3,Float64}}

Generate cell-centre positions for the cubic subdivision of S² at level
`subdivision` (`6 × 4^subdivision` cells).

The unit cube faces are subdivided into `n × n` grids with `n = 2^subdivision`;
each grid-cell centre is projected radially onto S².

The 6 faces (and their parameter directions):
- Face 1 `(u, v, +1)` — top   (+z)
- Face 2 `(u, v, −1)` — bottom (−z)
- Face 3 `(+1, u, v)` — right  (+x)
- Face 4 `(−1, u, v)` — left   (−x)
- Face 5 `(u, +1, v)` — front  (+y)
- Face 6 `(u, −1, v)` — back   (−y)
"""
function _cubic_sphere_points(subdivision::Int)::Vector{NTuple{3,Float64}}
    n      = 2^subdivision          # cells per edge per face
    step   = 2.0 / n                # cell width in parameter space [-1,1]
    c0     = -1.0 + step / 2.0     # centre of first cell

    # 6 face parametrisations: (u,v) ∈ [-1,1]² → cube point
    face_maps = (
        (u, v) -> ( u,  v,  1.0),  # +z top
        (u, v) -> ( u,  v, -1.0),  # -z bottom
        (u, v) -> ( 1.0,  u,  v),  # +x right
        (u, v) -> (-1.0,  u,  v),  # -x left
        (u, v) -> ( u,  1.0,  v),  # +y front
        (u, v) -> ( u, -1.0,  v),  # -y back
    )

    total  = 6 * n * n
    points = Vector{NTuple{3,Float64}}(undef, total)
    idx    = 0

    for face_map in face_maps
        for j in 0:(n-1), i in 0:(n-1)
            u = c0 + i * step
            v = c0 + j * step
            x, y, z = face_map(u, v)
            r = sqrt(x*x + y*y + z*z)
            idx += 1
            points[idx] = (x/r, y/r, z/r)
        end
    end
    return points
end

"""
    _sphere_to_homogeneous(pt) -> Tuple{ComplexF64, ComplexF64}

Convert a unit-sphere point `pt = (x₁, x₂, x₃)` to homogeneous coordinates
`[z : t]` via the inverse stereographic projection from the north pole.

Convention (consistent with `sphere_bijection`):
- North pole (0, 0, +1) → `∞` = `[1 : 0]`
- South pole (0, 0, −1) → `0` = `[0 : 1]`
- Generic point         → `z = (x₁ + i x₂)/(1 − x₃)`, `t = 1`
"""
function _sphere_to_homogeneous(
        pt::NTuple{3,Float64})::Tuple{Complex{Float64},Complex{Float64}}

    x1, x2, x3 = pt
    if x3 >= 1.0 - 1e-10          # north pole / near ∞
        return (complex(1.0), complex(0.0))
    else
        z = complex(x1, x2) / (1.0 - x3)
        return (z, complex(1.0))
    end
end

"""
    basin_sphere_data(coefficientlistnum, coefficientlistden;
                      subdivision=5, aproxprecision=8, iterprecision=8, iter=200)
                      -> SphereBasinData

Compute basin membership for the cubic subdivision of the Riemann sphere S².

For each of the `6 × 4^subdivision` cells the orbit of the cell-centre
representative is iterated until it converges to an attractor or the iteration
limit `iter` is reached (algorithm of Gutiérrez, Hernández, Magreñán, Rivas 2016).

**Subdivision levels and cell counts:**
| `subdivision` | cells    | approx. cell area (°) |
|:---:|----------|:---:|
|  4  | 1 536    | 8.2 msr  |
|  5  | 6 144    | 2.0 msr  |
|  6  | 24 576   | 0.51 msr |
|  7  | 98 304   | 0.13 msr |

# Arguments
- `coefficientlistnum`, `coefficientlistden` — numerator / denominator
  coefficients, **constant term first** (Polynomials.jl convention).
- `subdivision::Int` — subdivision level (default `5`; gives 6 144 cells).
- `aproxprecision::Int` — decimal digits for attractor identification
  (tolerance `c₂ = 10^{-aproxprecision}`, default `8`).
- `iterprecision::Int` — decimal digits for Cauchy convergence test
  (tolerance `c₁ = 10^{-iterprecision}`, default `8`).
- `iter::Int` — maximum orbit steps per cell (default `200`).

# Returns
A [`SphereBasinData`](@ref) struct.

# Example
```julia
using RiemannBasins
# Newton's method for z³ − 1:  f(z) = (2z³+1)/(3z²)
num = complex([1.0, 0.0, 0.0, 2.0])
den = complex([0.0, 0.0, 3.0, 0.0])
data = basin_sphere_data(num, den; subdivision=5)
areas = basin_measure(data)        # Dict: attractor index → area (steradians)
```
"""
function basin_sphere_data(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        subdivision   ::Int = 5,
        aproxprecision::Int = 8,
        iterprecision ::Int = 8,
        iter          ::Int = 200) where {T<:Number}

    cnum = ComplexF64.(coefficientlistnum)
    cden = ComplexF64.(coefficientlistden)

    # 1. Build rational function and find attractors
    hpair = _pair_of_functions(cnum, cden)
    _, _, finitelist = fixed_points(cnum, cden, aproxprecision)
    fp_list = Vector{Tuple{Complex{Float64},Complex{Float64}}}(finitelist)
    push!(fp_list, (complex(1.0), complex(0.0)))   # ∞ always last

    # 2. Generate cubic subdivision of S²
    pts   = _cubic_sphere_points(subdivision)
    total = length(pts)

    # 3. Classify each cell
    basin_idx = Vector{Int}(undef, total)
    for k in eachindex(pts)
        pt_hom = _sphere_to_homogeneous(pts[k])
        pos, _, _ = _position_iteration(
            fp_list, aproxprecision, iterprecision, pt_hom, hpair, iter)
        basin_idx[k] = max(0, pos)   # -1 (super-zero) → 0
    end

    cell_area = 4π / total
    return SphereBasinData(pts, basin_idx, fp_list, cell_area, subdivision, total)
end

"""
    basin_measure(data::SphereBasinData) -> Dict{Int,Float64}

Estimate the spherical area of each basin from a [`SphereBasinData`](@ref) result.

Returns a `Dict` mapping attractor index → estimated area (steradians).
Key `0` collects non-converged / Julia-set cells.

The estimate uses `cell_area × count` for each basin, where
`cell_area = 4π / total_cells` is uniform (the cubic projection is not
perfectly equal-area; the estimate converges to the true value as
`subdivision → ∞`).

**Initial efficiency** (Gutiérrez et al. 2016, Definition 2.14):
the probability that a uniformly-random point on S² converges to a
*finite* attracting fixed point equals `(sum of finite-basin areas) / 4π`.
"""
function basin_measure(data::SphereBasinData)::Dict{Int,Float64}
    counts = Dict{Int,Int}()
    for idx in data.basin_indices
        counts[idx] = get(counts, idx, 0) + 1
    end
    areas = Dict{Int,Float64}()
    for (k, v) in counts
        areas[k] = v * data.cell_area
    end
    return areas
end

"""
    basin_measure(coefficientlistnum, coefficientlistden; kwargs...) -> Dict{Int,Float64}

Convenience wrapper: compute [`basin_sphere_data`](@ref) and return basin areas.

All keyword arguments are forwarded to `basin_sphere_data`
(`subdivision`, `aproxprecision`, `iterprecision`, `iter`).

# Example
```julia
num = complex([1.0, 0.0, 0.0, 2.0])   # Newton for z³−1
den = complex([0.0, 0.0, 3.0, 0.0])
areas = basin_measure(num, den; subdivision=6)
# Basin of ∞ has index == length(fp_list)
```
"""
function basin_measure(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        kwargs...) where {T<:Number}
    data = basin_sphere_data(coefficientlistnum, coefficientlistden; kwargs...)
    return basin_measure(data)
end

# ──────────────────────────────────────────────────────────────────────────────
# §16  GLMakie sphere portrait stub
# ──────────────────────────────────────────────────────────────────────────────

"""
    basin_sphere_portrait_3d(coefficientlistnum, coefficientlistden; kwargs...)

Render the basin portrait on the actual Riemann sphere S² with **GLMakie**
(OpenGL, interactive rotation).

Requires `using GLMakie` before calling — the extension is loaded automatically
when GLMakie is in scope.

Keyword arguments (forwarded to the extension):
- `subdivision::Int` — cubic-subdivision level (default `5`; gives 6 144 cells).
- `aproxprecision`, `iterprecision`, `iter` — same meaning as `basin_sphere_data`.
- `gseed::RGB{Float64}` — seed colour for the LCHab basin palette.
- `colorblind::Bool` — use the Okabe-Ito palette (default `false`).
- `show_plot::Bool` — display window (default `true`).
- `savepath::String` — save figure to this path (empty = no save).
"""
function basin_sphere_portrait_3d(args...; kwargs...)
    error("""
basin_sphere_portrait_3d requires GLMakie.
Add  `using GLMakie`  before calling this function.
""")
end


# ──────────────────────────────────────────────────────────────────────────────
# §17  HTML sphere portrait (pure Julia – no GLMakie required)
# ──────────────────────────────────────────────────────────────────────────────

"""
    basin_sphere_portrait_html(coefficientlistnum, coefficientlistden; kwargs...)

Export the basin portrait on S² as a **self-contained HTML file** that uses
Three.js and OrbitControls so collaborators can rotate the sphere in any
modern browser — no Julia or GLMakie required.

## Keyword arguments
| Keyword          | Default            | Description                                  |
|------------------|--------------------|----------------------------------------------|
| `ntheta`         | `150`              | Latitude grid points                         |
| `nphi`           | `300`              | Longitude grid points (nphi+1 columns used)  |
| `aproxprecision` | `8`                | Convergence tolerance (10^-p)                |
| `iterprecision`  | `8`                | Same as `basin_sphere_data`                  |
| `iter`           | `200`              | Maximum iteration count                      |
| `bands`          | `6.0`              | Banding frequency (0.0 = no banding)         |
| `gseed`          | `RGB(0.2,0.6,0.9)` | Seed colour for LCHab basin palette          |
| `colorblind`     | `false`            | Use Okabe-Ito palette instead                |
| `htmlpath`       | `"basin_sphere.html"` | Output path for the HTML file             |
| `title`          | `"Basin portrait on S²"` | Page title shown in the browser        |

Returns the `htmlpath` string.
"""
function basin_sphere_portrait_html(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        ntheta        ::Int          = 150,
        nphi          ::Int          = 300,
        aproxprecision::Int          = 8,
        iterprecision ::Int          = 8,
        iter          ::Int          = 200,
        bands         ::Float64      = 6.0,
        gseed         ::RGB{Float64} = RGB{Float64}(0.2, 0.6, 0.9),
        colorblind    ::Bool         = false,
        htmlpath      ::String       = "basin_sphere.html",
        title         ::String       = "Basin portrait on S²") where {T<:Number}

    # ── 1. Build rational-map internals ──────────────────────────────────────
    cnum  = ComplexF64.(coefficientlistnum)
    cden  = ComplexF64.(coefficientlistden)
    hpair = _pair_of_functions(cnum, cden)
    _, _, finitelist = fixed_points(cnum, cden, aproxprecision)
    fp_list = Vector{Tuple{Complex{Float64},Complex{Float64}}}(finitelist)
    push!(fp_list, (complex(1.0), complex(0.0)))
    nfixed = length(fp_list)
    cmap   = _build_colormap(nfixed, gseed; colorblind = colorblind)

    # ── 2. Grid ──────────────────────────────────────────────────────────────
    thetas = range(1e-4, π - 1e-4, length = ntheta)
    ncols  = nphi + 1

    X = Float64[sin(θ) * cos(φ)
                for θ in thetas, φ in range(0.0, 2π, length = ncols)]
    Y = Float64[sin(θ) * sin(φ)
                for θ in thetas, φ in range(0.0, 2π, length = ncols)]
    Z = Float64[cos(θ)
                for θ in thetas, φ in range(0.0, 2π, length = ncols)]

    # ── 3. Classify ──────────────────────────────────────────────────────────
    basin_idx = Matrix{Int}(undef, ntheta, ncols)
    iter_mat  = Matrix{Int}(undef, ntheta, ncols)
    resid_mat = Matrix{Float64}(undef, ntheta, ncols)
    tol_c     = 1.0 / 10^aproxprecision

    for j in 1:nphi, i in 1:ntheta
        x3 = Z[i, j]
        pt_hom = x3 >= 1.0 - 1e-10 ?
            (complex(1.0), complex(0.0)) :
            (complex(X[i, j], Y[i, j]) / (1.0 - x3), complex(1.0))
        pos, n, residual = _position_iteration(
            fp_list, aproxprecision, iterprecision, pt_hom, hpair, iter)
        basin_idx[i, j] = max(0, pos)
        iter_mat[i, j]  = n
        resid_mat[i, j] = residual
    end
    # close the φ-seam (column ncols repeats column 1)
    basin_idx[:, ncols] = basin_idx[:, 1]
    iter_mat[:, ncols]  = iter_mat[:, 1]
    resid_mat[:, ncols] = resid_mat[:, 1]

    # ── 4. Per-vertex RGB bytes ───────────────────────────────────────────────
    nv        = ntheta * ncols
    rgb_bytes = Vector{UInt8}(undef, 3 * nv)
    for i in 1:ntheta, j in 1:ncols
        v = (i - 1) * ncols + (j - 1)   # 0-indexed vertex in JS
        k = basin_idx[i, j]
        c = cmap[k <= 0 ? 1 : min(k + 1, length(cmap))]
        if k > 0 && bands > 0.0
            r   = max(resid_mat[i, j], 1e-15)
            tv  = max(tol_c, 1e-15)
            lr  = log(log(tv) / log(r)) / log(2.0)
            ts  = clamp(Float64(iter_mat[i, j]) + lr, 0.0, Float64(iter))
            c   = _modulate_lab(c, ts, Float64(iter); bands = bands)
        end
        rgb_bytes[3v + 1] = round(UInt8, clamp(red(c),   0.0, 1.0) * 255)
        rgb_bytes[3v + 2] = round(UInt8, clamp(green(c), 0.0, 1.0) * 255)
        rgb_bytes[3v + 3] = round(UInt8, clamp(blue(c),  0.0, 1.0) * 255)
    end

    b64data = Base64.base64encode(rgb_bytes)

    # ── 5. Compose HTML ───────────────────────────────────────────────────────
    html = string(
        "<!DOCTYPE html>\n",
        "<html lang=\"en\">\n",
        "<head>\n",
        "<meta charset=\"utf-8\">\n",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
        "<title>", title, "</title>\n",
        "<style>\n",
        "  * { margin: 0; padding: 0; box-sizing: border-box; }\n",
        "  body { background: #000; overflow: hidden; }\n",
        "  canvas { display: block; }\n",
        "  #info {\n",
        "    position: absolute; top: 12px; left: 12px;\n",
        "    color: #ccc; font: 13px/1.5 sans-serif;\n",
        "    pointer-events: none; user-select: none;\n",
        "  }\n",
        "  #info b { color: #fff; font-size: 15px; }\n",
        "</style>\n",
        "</head>\n",
        "<body>\n",
        "<div id=\"info\"><b>", title, "</b><br>",
        "Drag · scroll to zoom · right-drag to pan</div>\n",
        "<script src=\"https://unpkg.com/three@0.128.0/build/three.min.js\"></script>\n",
        "<script src=\"https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js\"></script>\n",
        "<script>\n",
        "// ── data embedded by RiemannBasins.jl ────────────────────────────────────────────────────────────────────────────────\n",
        "const NTHETA = ", string(ntheta), ";\n",
        "const NPHI   = ", string(nphi), ";\n",
        "const B64    = \"", b64data, "\";\n",
        "\n",
        "// ── decode Base64 → Uint8Array ──────────────────────────────────────────────────────────────────────────────\n",
        "function decodeB64(b64) {\n",
        "  const bin = atob(b64);\n",
        "  const arr = new Uint8Array(bin.length);\n",
        "  for (let k = 0; k < bin.length; k++) arr[k] = bin.charCodeAt(k);\n",
        "  return arr;\n",
        "}\n",
        "const colorBytes = decodeB64(B64);\n",
        "\n",
        "// ── Three.js scene ──────────────────────────────────────────────────────────────────────────────────────\n",
        "const scene    = new THREE.Scene();\n",
        "const W = window.innerWidth, H = window.innerHeight;\n",
        "const camera   = new THREE.PerspectiveCamera(45, W / H, 0.01, 100);\n",
        "camera.position.set(0, 0, 2.6);\n",
        "\n",
        "const renderer = new THREE.WebGLRenderer({ antialias: true });\n",
        "renderer.setPixelRatio(window.devicePixelRatio);\n",
        "renderer.setSize(W, H);\n",
        "document.body.appendChild(renderer.domElement);\n",
        "\n",
        "const controls = new THREE.OrbitControls(camera, renderer.domElement);\n",
        "controls.enableDamping = true;\n",
        "controls.dampingFactor = 0.08;\n",
        "controls.minDistance   = 1.1;\n",
        "controls.maxDistance   = 8.0;\n",
        "\n",
        "// ── build BufferGeometry ───────────────────────────────────────────────────────────────────────────────────\n",
        "const NCOLS = NPHI + 1;\n",
        "const NV    = NTHETA * NCOLS;\n",
        "\n",
        "const positions = new Float32Array(NV * 3);\n",
        "const colorsArr = new Float32Array(NV * 3);\n",
        "\n",
        "for (let i = 0; i < NTHETA; i++) {\n",
        "  const theta = 1e-4 + (Math.PI - 2e-4) * i / (NTHETA - 1);\n",
        "  for (let j = 0; j < NCOLS; j++) {\n",
        "    const phi = 2 * Math.PI * j / NPHI;\n",
        "    const v   = i * NCOLS + j;\n",
        "    positions[v*3]     = Math.sin(theta) * Math.cos(phi);\n",
        "    positions[v*3 + 1] = Math.sin(theta) * Math.sin(phi);\n",
        "    positions[v*3 + 2] = Math.cos(theta);\n",
        "    colorsArr[v*3]     = colorBytes[v*3]     / 255;\n",
        "    colorsArr[v*3 + 1] = colorBytes[v*3 + 1] / 255;\n",
        "    colorsArr[v*3 + 2] = colorBytes[v*3 + 2] / 255;\n",
        "  }\n",
        "}\n",
        "\n",
        "// triangulate the grid (two triangles per cell)\n",
        "const indices = [];\n",
        "for (let i = 0; i < NTHETA - 1; i++) {\n",
        "  for (let j = 0; j < NPHI; j++) {\n",
        "    const a = i       * NCOLS + j;\n",
        "    const b = a + 1;\n",
        "    const c = (i + 1) * NCOLS + j;\n",
        "    const d = c + 1;\n",
        "    indices.push(a, c, b,  b, c, d);\n",
        "  }\n",
        "}\n",
        "\n",
        "const geo = new THREE.BufferGeometry();\n",
        "geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));\n",
        "geo.setAttribute('color',    new THREE.BufferAttribute(colorsArr, 3));\n",
        "geo.setIndex(indices);\n",
        "\n",
        "const mat  = new THREE.MeshBasicMaterial({\n",
        "  vertexColors: true,\n",
        "  side: THREE.DoubleSide\n",
        "});\n",
        "scene.add(new THREE.Mesh(geo, mat));\n",
        "\n",
        "// ── resize handler ──────────────────────────────────────────────────────────────────────────────────────────\n",
        "window.addEventListener('resize', () => {\n",
        "  camera.aspect = window.innerWidth / window.innerHeight;\n",
        "  camera.updateProjectionMatrix();\n",
        "  renderer.setSize(window.innerWidth, window.innerHeight);\n",
        "});\n",
        "\n",
        "// ── render loop ─────────────────────────────────────────────────────────────────────────────────────────────\n",
        "(function animate() {\n",
        "  requestAnimationFrame(animate);\n",
        "  controls.update();\n",
        "  renderer.render(scene, camera);\n",
        "})();\n",
        "</script>\n",
        "</body>\n",
        "</html>\n"
    )

    open(htmlpath, "w") do io
        write(io, html)
    end
    return htmlpath
end


end  # module RiemannBasins
