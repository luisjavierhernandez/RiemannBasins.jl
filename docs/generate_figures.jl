# generate_figures.jl
# Genera las figuras que faltan para borrador_gaceta.tex
# Ejecutar desde el directorio raíz del paquete:
#   julia --project docs/generate_figures.jl
#
# Requisitos: RiemannBasins.jl compilado (Plots + GR backend).
# La figura relief_makie.png debe generarse manualmente con GLMakie.

using RiemannBasins
using Plots
gr()

figdir = joinpath(@__DIR__, "figures")
mkpath(figdir)

println("=== Generando figuras para borrador_gaceta.tex ===\n")

# ─────────────────────────────────────────────────────────────────────────────
# Fig. 1: Newton para z^4 - 1 en las dos cartas
#   N(z) = (3z^4 + 1) / (4z^3)
#   num: [3, 0, 0, 0, 1]  (term. const. primero: 3·z^0 + 0·z + 0·z^2 + 0·z^3 + 1·z^4)
#         Wait: N(z)= z - (z^4-1)/(4z^3) = (4z^4 - z^4 + 1)/(4z^3) = (3z^4+1)/(4z^3)
#   Numerador:  1 + 0z + 0z^2 + 0z^3 + 3z^4  → [1.0, 0.0, 0.0, 0.0, 3.0]
#   Denominador: 0 + 4z + 0z^2 + 0z^3 + ...  → [0.0, 0.0, 0.0, 4.0, 0.0]
#   (convenio: coeficiente del término z^k está en posición k+1)
# ─────────────────────────────────────────────────────────────────────────────
println("1. Newton z^4-1: dos cartas …")
num_z4 = complex([1.0, 0.0, 0.0, 0.0, 3.0])   # 1 + 3z^4
den_z4 = complex([0.0, 0.0, 0.0, 4.0, 0.0])   # 4z^3

savepath_z4 = joinpath(figdir, "newton_z4_two_charts.png")
two_chart_portrait(num_z4, den_z4;
    bands      = 8.0,
    show_plot  = false,
    savepath   = savepath_z4)
println("   Guardado en $savepath_z4")

# ─────────────────────────────────────────────────────────────────────────────
# Fig. 2: Método de Halley para z^3 - 1 en las dos cartas
#   H(z) = (z^4 + 2z) / (2z^3 + 1)
#   Numerador:  0 + 2z + 0z^2 + 0z^3 + 1z^4  → [0.0, 2.0, 0.0, 0.0, 1.0]
#   Denominador: 1 + 0z + 0z^2 + 2z^3 + ...  → [1.0, 0.0, 0.0, 2.0, 0.0]
# ─────────────────────────────────────────────────────────────────────────────
println("2. Halley z^3-1: dos cartas …")
num_h = complex([0.0, 2.0, 0.0, 0.0, 1.0])   # 2z + z^4
den_h = complex([1.0, 0.0, 0.0, 2.0, 0.0])   # 1 + 2z^3

savepath_h = joinpath(figdir, "halley_z3_two_charts.png")
two_chart_portrait(num_h, den_h;
    bands      = 8.0,
    show_plot  = false,
    savepath   = savepath_h)
println("   Guardado en $savepath_h")

# ─────────────────────────────────────────────────────────────────────────────
# Fig. 3: Nota para la figura GLMakie (debe generarse manualmente)
# ─────────────────────────────────────────────────────────────────────────────
println("""
3. Figura relief_makie.png (requiere ventana nativa OpenGL):
   Ejecutar en el REPL de Julia con GLMakie:

   using RiemannBasins
   # Newton z^3-1: N(z) = (2z^3+1)/(3z^2)
   num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z^3
   den = complex([0.0, 0.0, 3.0, 0.0])   # 3z^2
   basin_relief_3d_makie(num, den; bands=50.0, relief_scale=0.01, show_plot=true)
   # Hacer captura de pantalla y guardar como docs/figures/relief_makie.png
""")

println("=== Figuras completadas ===")
