module RiemannBasinsGLMakieExt

using RiemannBasins
using GLMakie
using Colors

# ─────────────────────────────────────────────────────────────────────────────
# Override the stub defined in the parent module.
# Julia loads this extension automatically when the user does `using GLMakie`.
# ─────────────────────────────────────────────────────────────────────────────

function RiemannBasins.basin_relief_3d_makie(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        xinterval::Tuple{Real,Real}    = (-1.5, 1.5),
        yinterval::Tuple{Real,Real}    = (-1.5, 1.5),
        aproxprecision::Int64          = 8,
        iterprecision::Int64           = 8,
        iter::Int64                    = 200,
        expprecision::Int64            = 7,
        bands::Float64                 = 50.0,
        relief_scale::Float64          = 0.01,
        ridge_fraction::Float64        = 0.15,
        camera_angle::Tuple{Real,Real} = (30, 82),
        gseed::RGB{Float64}            = RGB{Float64}(0.2, 0.6, 0.9),
        colorblind::Bool               = false,
        colorstrategy::String          = "positionplusiteration",
        inverted::Bool                 = false,
        background                     = :black,
        shading::Bool                  = true,
        show_plot::Bool                = true,
        savepath::String               = "") where {T<:Number}

    # ── 1. Compute basin matrices (reuse existing internals) ─────────────────
    fixedPointList, positionMatrix, iterationMatrix, residualMatrix =
        RiemannBasins._basin_matrices(coefficientlistnum, coefficientlistden,
                                      xinterval, yinterval,
                                      aproxprecision, iterprecision, iter, expprecision)

    nrows, ncols = size(positionMatrix)
    maxiter      = Float64(max(1, maximum(iterationMatrix)))
    tol          = 1.0 / 10^aproxprecision
    nfixed       = length(fixedPointList)

    xs = collect(range(Float64(xinterval[1]), Float64(xinterval[2]), length = ncols))
    ys = collect(range(Float64(yinterval[1]), Float64(yinterval[2]), length = nrows))

    # ── 2. Colour matrix — full LCHab-banded basin portrait ─────────────────
    # Build the colour map (LCHab arc or Okabe-Ito)
    cmap = RiemannBasins._build_colormap(nfixed, gseed; colorblind = colorblind)

    # Build raw pixel matrix (returns Matrix{RGB{Float64}}, not a Plots object)
    img = RiemannBasins._build_basin_image_matrix(
        fixedPointList, positionMatrix, iterationMatrix,
        cmap, colorstrategy;
        colorblind     = colorblind,
        residualMatrix = residualMatrix,
        aproxprecision = aproxprecision,
        bands          = bands)

    # ── 3. Elevation matrix — decoupled fractal ridge ────────────────────────
    # basin_norm  → altitude plateau per attractor  (major steps)
    # ridge_norm  → damped cosine ripple            (fractal texture)
    # Inverted: basin centres are raised mounds, boundary is the valley.
    z = zeros(Float64, nrows, ncols)
    for j in 1:ncols, i in 1:nrows
        pos = positionMatrix[i, j]
        nit = Float64(iterationMatrix[i, j])

        if pos > 0 && residualMatrix[i, j] > 0
            r         = max(residualMatrix[i, j], 1e-15)
            tol_c     = max(tol, 1e-15)
            log_ratio = log(log(tol_c) / log(r)) / log(2.0)
            t_smooth  = clamp(nit + log_ratio, 0.0, maxiter)
        else
            t_smooth = nit
        end

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

    # ── 4. Convert to Makie-compatible types ─────────────────────────────────
    # Makie surface(xs, ys, zs): zs[i,j] ↔ (xs[i], ys[j])  →  need transpose
    z_mk   = collect(z')                           # (ncols × nrows)
    # Convert RGB{Float64} → RGBAf for Makie
    img_mk = [RGBAf(red(img[r,c]), green(img[r,c]), blue(img[r,c]), 1.0f0)
              for c in 1:ncols, r in 1:nrows]     # (ncols × nrows) after transpose

    # ── 5. GLMakie figure ────────────────────────────────────────────────────
    az = deg2rad(Float64(camera_angle[1]))
    el = deg2rad(Float64(camera_angle[2]))

    fig = Figure(backgroundcolor = background, size = (820, 780))
    ax  = Axis3(fig[1, 1];
        xlabel          = "Re(z)",
        ylabel          = "Im(z)",
        zlabel          = "",
        title           = "3D basin relief  (bands = $(Int(round(bands))))" * (inverted ? "  [inverted]" : ""),
        azimuth         = az,
        elevation       = el,
        titlecolor      = :white,
        xlabelcolor     = :white,
        ylabelcolor     = :white,
        xticklabelcolor = :white,
        yticklabelcolor = :white,
        zticklabelcolor = :white,
        xgridcolor      = (:white, 0.15),
        ygridcolor      = (:white, 0.15),
        zgridcolor      = (:white, 0.15),
        backgroundcolor = (:black, 0.0))

    surface!(ax, xs, ys, z_mk;
        color       = img_mk,
        shading     = shading,
        interpolate = true)

    show_plot  && display(GLMakie.Screen(), fig)
    isempty(savepath) || save(savepath, fig)
    return fig
end

end  # module RiemannBasinsGLMakieExt
