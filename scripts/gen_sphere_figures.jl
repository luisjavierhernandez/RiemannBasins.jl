# gen_sphere_figures.jl
# Genera los retratos de cuencas sobre la esfera S² (GLMakie) para
# Newton z³−1 y el mapa racional f(z)=(z²−1)/(z²+1) iterado dos veces.
#
# Uso:   julia --project scripts/gen_sphere_figures.jl
# Salida en docs/articles/01_pcycle_basins/figures/
#   sphere_newton_z3.png        — Newton z³−1  (tres cuencas simétricas)
#   sphere_rational_2cycle.png  — f²(z)=-2z²/(z⁴+1)  (dos cuencas, esfera entera)
#
# Requisito: GLMakie instalado en el entorno del proyecto.
# En macOS se necesita entorno gráfico (no headless).

using GLMakie, RiemannBasins, Colors

outdir = joinpath(@__DIR__, "..", "docs", "articles", "01_pcycle_basins", "figures")
mkpath(outdir)

# Parámetros comunes de calidad
# size=(820,820), px_per_unit=4  →  3280×3280 px guardados  (~300 dpi en 11")
NTHETA = 400
NPHI   = 800
ITER   = 400
BANDS  = 6.0
PX     = 4

# ─── Figura 1: Newton para z³−1 ─────────────────────────────────────────────
# N(z) = (2z³+1)/(3z²)   →   num=[1,0,0,2],  den=[0,0,3,0]
# Tres cuencas de igual tamaño (H_norm≈1).  Cámara a 25° sobre el ecuador,
# ligeramente inclinada para ver la simetría Z/3 en perspectiva.
println("=== Figura 1: Newton para z³−1 ===")
num_new = complex([1.0, 0.0, 0.0, 2.0])
den_new = complex([0.0, 0.0, 3.0, 0.0])

fig_new = basin_sphere_portrait_3d(num_new, den_new;
    ntheta       = NTHETA,
    nphi         = NPHI,
    iter         = ITER,
    bands        = BANDS,
    show_plot    = false,
    px_per_unit  = PX,
    camera_angle = (30.0, 25.0),   # az=30°, el=25°: simetría 3-fold visible
    savepath     = joinpath(outdir, "sphere_newton_z3.png"))
println("  → sphere_newton_z3.png")

# ─── Figura 2: Mapa racional f²(z) = −2z²/(z⁴+1) ───────────────────────────
# f(z) = (z²−1)/(z²+1) tiene el 2-ciclo {0,−1}.
# f²(z) = −2z²/(z⁴+1) tiene 0 y −1 como puntos fijos atractores.
# Las dos cuencas cubren toda la 2-esfera (conjunto de Julia de medida cero).
# Simetría: f²(1/z) = f²(z)  ⟹  el polo norte (z=∞) y el polo sur (z=0)
# tienen el MISMO color (ambos en la cuenca de 0).
# Cámara hacia el polo sur (z=0) para que esa cuenca se vea como el "centro".
println("=== Figura 2: Mapa racional f²(z) = −2z²/(z⁴+1) ===")
num_rat2 = complex([0.0, 0.0, -2.0])            # −2z²
den_rat2 = complex([1.0, 0.0, 0.0, 0.0, 1.0])  # z⁴+1

fig_rat = basin_sphere_portrait_3d(num_rat2, den_rat2;
    ntheta       = NTHETA,
    nphi         = NPHI,
    iter         = ITER,
    bands        = BANDS,
    show_plot    = false,
    px_per_unit  = PX,
    camera_angle = (45.0, -70.0),  # az=45°, el=−70°: polo sur (z=0) bien visible
    savepath     = joinpath(outdir, "sphere_rational_2cycle.png"))
println("  → sphere_rational_2cycle.png")

println("\nFiguras esféricas guardadas en:")
println("  ", joinpath(outdir, "sphere_newton_z3.png"))
println("  ", joinpath(outdir, "sphere_rational_2cycle.png"))
