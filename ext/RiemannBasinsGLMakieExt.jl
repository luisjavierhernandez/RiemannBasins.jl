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
        bands::Float64                 = 10.0,
        relief_scale::Float64          = 0.01,
        ridge_fraction::Float64        = 0.15,
        smooth_z::Float64              = 1.5,
        camera_angle::Tuple{Real,Real} = (30, 82),
        gseed::RGB{Float64}            = RGB{Float64}(0.2, 0.6, 0.9),
        colorblind::Bool               = false,
        colorstrategy::String          = "positionplusiteration",
        inverted::Bool                 = false,
        background                     = :black,
        px_per_unit::Int               = 1,
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
    cmap = RiemannBasins._build_colormap(nfixed, gseed; colorblind = colorblind)

    img = RiemannBasins._build_basin_image_matrix(
        fixedPointList, positionMatrix, iterationMatrix,
        cmap, colorstrategy;
        colorblind     = colorblind,
        residualMatrix = residualMatrix,
        aproxprecision = aproxprecision,
        bands          = bands)

    # ── 3. Elevation matrix — decoupled fractal ridge ────────────────────────
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

    # ── 3b. Gaussian smoothing — eliminates sharp spikes ──────────────────────
    if smooth_z > 0.0
        r   = max(1, ceil(Int, 3smooth_z))
        ker = [exp(-i^2 / (2smooth_z^2)) for i in -r:r]
        ker ./= sum(ker)
        nr, nc = size(z)
        tmp = copy(z)
        for i in 1:nr, j in 1:nc    # horizontal pass
            s = 0.0
            for di in -r:r
                s += ker[di+r+1] * z[i, clamp(j+di, 1, nc)]
            end
            tmp[i, j] = s
        end
        for i in 1:nr, j in 1:nc    # vertical pass
            s = 0.0
            for di in -r:r
                s += ker[di+r+1] * tmp[clamp(i+di, 1, nr), j]
            end
            z[i, j] = s
        end
    end

    # ── 4. Convert to Makie-compatible types ─────────────────────────────────
    z_mk   = collect(z')
    img_mk = [RGBAf(red(img[r,c]), green(img[r,c]), blue(img[r,c]), 1.0f0)
              for c in 1:ncols, r in 1:nrows]

    # ── 5. GLMakie figure ────────────────────────────────────────────────────
    az = deg2rad(Float64(camera_angle[1]))
    el = deg2rad(Float64(camera_angle[2]))

    fig = Figure(backgroundcolor = background, size = (820, 780))
    ax  = Axis3(fig[1, 1];
        xlabel          = "Re(z)",
        ylabel          = "Im(z)",
        zlabel          = "",
        title           = "3D basin relief  (bands=$(Int(round(bands))), σ=$(smooth_z))" * (inverted ? "  [inverted]" : ""),
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
        backgroundcolor = (:black, 0.0),
        limits          = (Float64(xinterval[1]), Float64(xinterval[2]),
                           Float64(yinterval[1]), Float64(yinterval[2]),
                           0.0, Float64(relief_scale)))

    # ── 6. Surface with shading ─────────────────────────────────────────────
    surface!(ax, xs, ys, z_mk;
        color       = img_mk,
        shading     = true,
        interpolate = true)

    show_plot  && display(fig)
    isempty(savepath) || save(savepath, fig; px_per_unit = px_per_unit)
    return fig
end


# ─────────────────────────────────────────────────────────────────────────────
# §16  basin_sphere_portrait_3d  — S² basin portrait with radial relief
# ─────────────────────────────────────────────────────────────────────────────

function RiemannBasins.basin_sphere_portrait_3d(
        coefficientlistnum::Vector{T},
        coefficientlistden::Vector{T};
        ntheta        ::Int              = 300,
        nphi          ::Int              = 600,
        aproxprecision::Int              = 8,
        iterprecision ::Int              = 8,
        iter          ::Int              = 200,
        bands         ::Float64          = 6.0,
        relief_scale  ::Float64          = 0.12,
        ridge_fraction::Float64          = 0.1,
        smooth_z      ::Float64          = 2.0,
        inverted      ::Bool             = false,
        gseed         ::RGB{Float64}     = RGB{Float64}(0.2, 0.6, 0.9),
        colorblind    ::Bool             = false,
        background                       = :black,
        show_plot     ::Bool             = true,
        savepath      ::String           = "",
        px_per_unit   ::Int              = 2,
        camera_angle  ::Tuple{Real,Real} = (30.0, 20.0),
        cam_dist      ::Float64          = 2.2) where {T<:Number}
    # camera_angle = (azimuth_deg, elevation_deg)
    # azimuth  : angle in the XY-plane from the +X axis (degrees)
    # elevation: angle above (+) or below (−) the equatorial plane (degrees)
    #   elevation =  90 → camera above north pole (z=∞)
    #   elevation =   0 → equatorial view
    #   elevation = −90 → camera below south pole (z=0)

    # ── 1. Build rational-map internals ──────────────────────────────────────────
    cnum  = ComplexF64.(coefficientlistnum)
    cden  = ComplexF64.(coefficientlistden)
    hpair = RiemannBasins._pair_of_functions(cnum, cden)

    _, _, finitelist = RiemannBasins.fixed_points(cnum, cden, aproxprecision)
    fp_list = Vector{Tuple{Complex{Float64},Complex{Float64}}}(finitelist)
    push!(fp_list, (complex(1.0), complex(0.0)))   # ∞ always last
    nfixed = length(fp_list)

    # ── 2. Colour palette ───────────────────────────────────────────────────────
    cmap  = RiemannBasins._build_colormap(nfixed, gseed; colorblind = colorblind)

    # ── 3. Regular (θ,φ) parametric grid on S² ─────────────────────────────────
    # θ ∈ (ε, π−ε)  colatitude from north pole
    # φ ∈ [0, 2π]   longitude — nphi+1 columns; last == first to close the seam
    thetas = range(1e-4, π - 1e-4, length = ntheta)
    phis   = range(0.0,  2π,        length = nphi + 1)

    X = Float64[sin(θ) * cos(φ) for θ in thetas, φ in phis]
    Y = Float64[sin(θ) * sin(φ) for θ in thetas, φ in phis]
    Z = Float64[cos(θ)           for θ in thetas, φ in phis]

    # ── 4. Classify each grid point ─────────────────────────────────────────────
    basin_idx = Matrix{Int}(undef, ntheta, nphi + 1)
    iter_mat  = Matrix{Int}(undef, ntheta, nphi + 1)
    resid_mat = Matrix{Float64}(undef, ntheta, nphi + 1)
    Threads.@threads for j in 1:nphi
        for i in 1:ntheta
            x3 = Z[i, j]
            if x3 >= 1.0 - 1e-10
                pt_hom = (complex(1.0), complex(0.0))
            else
                x1, x2 = X[i, j], Y[i, j]
                pt_hom = (complex(x1, x2) / (1.0 - x3), complex(1.0))
            end
            pos, n, residual = RiemannBasins._position_iteration(
                fp_list, aproxprecision, iterprecision, pt_hom, hpair, iter)
            basin_idx[i, j] = max(0, pos)
            iter_mat[i, j]  = n
            resid_mat[i, j] = residual
        end
    end
    basin_idx[:, nphi + 1] = basin_idx[:, 1]
    iter_mat[:, nphi + 1]  = iter_mat[:, 1]
    resid_mat[:, nphi + 1] = resid_mat[:, 1]

    # ── 5. Colour matrix + smooth iteration field ────────────────────────────────
    tol_c   = 1.0 / 10^aproxprecision
    maxiter = Float64(max(1, maximum(iter_mat)))
    img_mk  = Matrix{RGBAf}(undef, ntheta, nphi + 1)
    t_mat   = zeros(Float64, ntheta, nphi + 1)
    for j in 1:nphi+1, i in 1:ntheta
        k   = basin_idx[i, j]
        nit = Float64(iter_mat[i, j])
        if k > 0 && resid_mat[i, j] > 0
            rv        = max(resid_mat[i, j], 1e-15)
            tv        = max(tol_c, 1e-15)
            log_ratio = log(log(tv) / log(rv)) / log(2.0)
            t_smooth  = clamp(nit + log_ratio, 0.0, maxiter)
        else
            t_smooth = nit
        end
        t_mat[i, j] = clamp(t_smooth / maxiter, 0.0, 1.0)

        c = cmap[k <= 0 ? 1 : min(k + 1, length(cmap))]
        if k > 0 && bands > 0.0
            c = RiemannBasins._modulate_lab(c, t_smooth, maxiter; bands = bands)
        end
        img_mk[i, j] = RGBAf(Float32(red(c)), Float32(green(c)), Float32(blue(c)), 1.0f0)
    end

    # ── 5b. Gaussian smoothing of t_mat on the (θ,φ) grid ───────────────────────
    # θ (rows): clamp  —  φ (cols): mod1 wrap (periódico)
    if smooth_z > 0.0
        r_ker = max(1, ceil(Int, 3smooth_z))
        ker   = [exp(-i^2 / (2smooth_z^2)) for i in -r_ker:r_ker]
        ker ./= sum(ker)
        nr, nc = ntheta, nphi + 1
        tmp = copy(t_mat)
        for i in 1:nr, j in 1:nc    # horizontal (φ) — periódico
            s = 0.0
            for di in -r_ker:r_ker
                jj = mod1(j + di, nphi)
                s += ker[di + r_ker + 1] * t_mat[i, jj]
            end
            tmp[i, j] = s
        end
        for i in 1:nr, j in 1:nc    # vertical (θ) — clamp
            s = 0.0
            for di in -r_ker:r_ker
                s += ker[di + r_ker + 1] * tmp[clamp(i + di, 1, nr), j]
            end
            t_mat[i, j] = s
        end
    end

    # ── 5c. Desplazamiento radial: r(θ,φ) = 1 + relief_scale · f(t_norm) ───────
    # inverted = false : atractores sobresalen (1−t alto cerca del atractor)
    # inverted = true  : frontera fractal sobresale (t alto en la frontera)
    X_r = copy(X)
    Y_r = copy(Y)
    Z_r = copy(Z)
    if relief_scale > 0.0
        for j in 1:nphi+1, i in 1:ntheta
            tn    = t_mat[i, j]
            f     = inverted ? tn : (1.0 - tn)
            phase = 2π * bands * tn
            disp  = f * (1.0 - ridge_fraction) + f * ridge_fraction * (0.5 + 0.5 * cos(phase))
            r     = 1.0 + relief_scale * disp
            X_r[i, j] = r * X[i, j]
            Y_r[i, j] = r * Y[i, j]
            Z_r[i, j] = r * Z[i, j]
        end
    end

    # ── 6. GLMakie figure — LScene (no axes; sphere speaks for itself) ───────────
    fig = Figure(backgroundcolor = background, size = (820, 820))
    ax  = LScene(fig[1, 1]; show_axis = false)

    surface!(ax, X_r, Y_r, Z_r;
        color       = img_mk,
        shading     = relief_scale > 0.0,
        interpolate = true)

    # ── 7. Set camera position from camera_angle ─────────────────────────────────
    az = deg2rad(Float64(camera_angle[1]))
    el = deg2rad(Float64(camera_angle[2]))
    cam_x = cam_dist * cos(el) * cos(az)
    cam_y = cam_dist * cos(el) * sin(az)
    cam_z = cam_dist * sin(el)
    Makie.update_cam!(ax.scene, Makie.Vec3f(cam_x, cam_y, cam_z),
                      Makie.Vec3f(0, 0, 0), Makie.Vec3f(0, 0, 1))

    show_plot  && display(fig)
    isempty(savepath) || save(savepath, fig; px_per_unit = px_per_unit)
    return fig
end

end  # module RiemannBasinsGLMakieExt
