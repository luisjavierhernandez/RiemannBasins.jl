# ─────────────────────────────────────────────────────────────────────────────
# relief_plane.jl
# Genera docs/figures/relief_plane.png
# Relieve 3D sobre el plano para Newton z³ − 1  (cuencas en [-1.5,1.5]²)
#
# Uso:  julia --project relief_plane.jl
# ─────────────────────────────────────────────────────────────────────────────

using RiemannBasins
using GLMakie
using Colors

# Newton z³ − 1 :  f(z) = (2z³ + 1) / (3z²)
# Coeficientes en orden ascendente: [a₀, a₁, ..., aₙ]
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 0z + 0z² + 2z³
den = complex([0.0, 0.0, 3.0])        # 0 + 0z + 3z²

fig = basin_relief_3d_makie(num, den;
    xinterval    = (-1.5, 1.5),
    yinterval    = (-1.5, 1.5),
    inverted     = true,       # frontera fractal = pico, atractor = valle
    bands        = 10.0,
    smooth_z     = 4.0,
    relief_scale = 0.01,
    camera_angle = (30, 70),   # (azimut°, elevación°)
    px_per_unit  = 2,
    show_plot    = false,
    savepath     = joinpath(@__DIR__, "relief_plane.png"))

@info "Guardado: docs/figures/relief_plane.png"
