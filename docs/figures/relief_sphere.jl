# ─────────────────────────────────────────────────────────────────────────────
# relief_sphere.jl
# Genera docs/figures/relief_sphere.png
# Relieve 3D sobre la esfera de Riemann para Newton z³ − 1
#
# Uso:  julia --project relief_sphere.jl
# ─────────────────────────────────────────────────────────────────────────────

using RiemannBasins
using GLMakie
using Colors

# Newton z³ − 1 :  f(z) = (2z³ + 1) / (3z²)
# Coeficientes en orden ascendente: [a₀, a₁, ..., aₙ]
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 0z + 0z² + 2z³
den = complex([0.0, 0.0, 3.0])        # 0 + 0z + 3z²

fig = basin_sphere_portrait_3d(num, den;
    ntheta       = 500,
    nphi         = 1000,
    inverted     = true,       # frontera fractal = protuberancia
    relief_scale = 0.12,
    smooth_z     = 2.0,
    bands        = 6.0,
    camera_angle = (30.0, -50.0),  # elevación negativa → polo sur visible
    cam_dist     = 2.2,            # distancia cámara: menor = esfera más grande
    px_per_unit  = 2,
    show_plot    = false,
    savepath     = joinpath(@__DIR__, "relief_sphere.png"))

@info "Guardado: docs/figures/relief_sphere.png"
