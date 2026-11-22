# gen_sphere_figures.jl
# Genera los retratos de cuencas sobre la esfera S² (GLMakie) para
# Newton y Halley aplicados a z³−1.
#
# Uso:   julia --project scripts/gen_sphere_figures.jl
# Salida en docs/articles/01_pcycle_basins/figures/
#   sphere_newton_z3.png   — Newton z³−1
#   sphere_halley_z3.png   — Halley z³−1
#
# Requisito: GLMakie instalado en el entorno del proyecto.
# En macOS se necesita entorno gráfico (no headless).

using GLMakie, RiemannBasins, Colors

outdir = joinpath(@__DIR__, "..", "docs", "articles", "01_pcycle_basins", "figures")
mkpath(outdir)

# Parámetros comunes de calidad
# size=(820,820), px_per_unit=4  →  3280×3280 px guardados  (~300 dpi en 11")
NTHETA    = 400
NPHI      = 800
ITER      = 400
BANDS     = 6.0
PX        = 4

println("=== Figura S1: Newton para z³−1 ===")
# N(z) = (2z³+1)/(3z²)   →   num=[1,0,0,2],  den=[0,0,3,0]
num_new = complex([1.0, 0.0, 0.0, 2.0])
den_new = complex([0.0, 0.0, 3.0, 0.0])

fig_new = basin_sphere_portrait_3d(num_new, den_new;
    ntheta      = NTHETA,
    nphi        = NPHI,
    iter        = ITER,
    bands       = BANDS,
    show_plot   = false,
    px_per_unit = PX,
    savepath    = joinpath(outdir, "sphere_newton_z3.png"))
println("  → sphere_newton_z3.png")

println("=== Figura S2: Halley para z³−1 ===")
# H(z) = (z⁴+2z)/(2z³+1)   →   num=[0,2,0,0,1],  den=[1,0,0,2]
num_hal = complex([0.0, 2.0, 0.0, 0.0, 1.0])
den_hal = complex([1.0, 0.0, 0.0, 2.0])

fig_hal = basin_sphere_portrait_3d(num_hal, den_hal;
    ntheta      = NTHETA,
    nphi        = NPHI,
    iter        = ITER,
    bands       = BANDS,
    show_plot   = false,
    px_per_unit = PX,
    savepath    = joinpath(outdir, "sphere_halley_z3.png"))
println("  → sphere_halley_z3.png")

println("\nFiguras esféricas guardadas en:")
println("  ", joinpath(outdir, "sphere_newton_z3.png"))
println("  ", joinpath(outdir, "sphere_halley_z3.png"))
