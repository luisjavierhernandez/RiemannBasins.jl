# gen_parabolic_cycles.jl  (v3 — ventanas simétricas + orden cíclico del p-ciclo)
#
# Uso:   julia scripts/gen_parabolic_cycles.jl
#
# Figuras generadas en docs/articles/02_fatou_hopf/figures/:
#
#   cycle_parabolic_p2.png      f(z)=z²-3/4  : punto fijo parabólico (z*=-1/2, mult=-1)
#   cycle_parabolic_p4.png      f(z)=z²-5/4  : 2-ciclo parabólico, con arcos cíclicos
#   cycle_attract_basilica.png  f(z)=z²-1    : 2-ciclo superatractor, con arcos cíclicos
#   cycle_birkhoff_comparison.png: B_N parabólico vs atractor (hasta N=15000)

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors, Plots, Printf

outdir = joinpath(@__DIR__, "..", "docs", "articles", "02_fatou_hopf", "figures")
mkpath(outdir)

# ── Convención coeficientes: cnum[k+1] = coef de z^k·t^(d−k), cte primero ──

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Utilidades para superponer los puntos del ciclo sobre el portrait
#
# basin_portrait construye la cuadrícula con _rectangle que itera
#   im_val in d:-tol:c   →  fila 1 = Im_max (d),  fila nrows = Im_min (c)
# Plots.jl muestra la fila 1 arriba → el eje Im ya es correcto.
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Coordenada compleja z → posición píxel (px, py) del plot de basin_portrait.
function _z2px(z_re::Float64, z_im::Float64,
               a::Float64, b::Float64,
               c::Float64, d::Float64,
               expp::Int)
    tol   = 1.0 / 2^expp
    ncols = round(Int, (b - a) / tol) + 1
    nrows = round(Int, (d - c) / tol) + 1
    px = 1.0 + (z_re - a) / (b - a) * (ncols - 1)
    py = 1.0 + (d - z_im) / (d - c) * (nrows - 1)   # d (Im_max) → fila 1 (arriba)
    return px, py
end

# Semicírculo entre dos puntos reales (z1_re,0) y (z2_re,0) en coordenadas píxel.
#   above=true  →  Im ≥ 0,  camino de z1 → z2  (flecha al final = en z2)
#   above=false →  Im ≤ 0,  camino de z2 → z1  (flecha al final = en z1)
function _half_arc(z1_re::Float64, z2_re::Float64,
                   a::Float64, b::Float64,
                   c::Float64, d::Float64,
                   expp::Int;
                   above::Bool = true, n::Int = 70)
    center_re = (z1_re + z2_re) / 2.0
    radius    = abs(z2_re - z1_re) / 2.0
    t = above ? range(0.0, Float64(π), length=n) :
                range(Float64(π), 2Float64(π), length=n)
    px = Vector{Float64}(undef, n)
    py = Vector{Float64}(undef, n)
    for (i, ti) in enumerate(t)
        r  = center_re + radius * cos(ti)
        im = radius    * sin(ti)
        px[i], py[i] = _z2px(r, im, a, b, c, d, expp)
    end
    return px, py
end

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Mapa 1:  f(z) = z² − 3/4
#   Punto fijo parabólico: z* = −1/2,  mult. f′(z*) = 2·(−1/2) = −1
#   En la frontera cardioide / bulbo p=2 del conjunto de Mandelbrot.
#   Ventana:  [−1.8, 1.8] × [−1.0, 1.0]  (simétrica, z*=−1/2 casi centrado)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
println("=== Mapa 1:  f(z) = z²−3/4  ===")
cv1   = -3.0/4.0
cnum1 = complex([cv1, 0.0, 1.0])
cden1 = complex([1.0, 0.0, 0.0])

res1 = critical_orbit_analysis(cnum1, cden1; detect_parabolic=true, N_birkhoff=1200)
for r in res1
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

x1lo, x1hi, y1lo, y1hi = -1.8, 1.8, -1.0, 1.0
expp = 8

plt1 = basin_portrait(cnum1, cden1;
    xinterval = (x1lo, x1hi),
    yinterval = (y1lo, y1hi),
    iter = 600, bands = 10.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath  = "")

# Marcar el punto fijo z* = −1/2 con estrella
px_s, py_s = _z2px(-0.5, 0.0, x1lo, x1hi, y1lo, y1hi, expp)
scatter!(plt1, [px_s], [py_s];
    marker = :star5, markersize = 14,
    color  = :white, markerstrokecolor = :black, markerstrokewidth = 1.5,
    label  = "")

# Bucle de auto-mapa (270° de círculo) para indicar z* → z*
loop_r = 24.0                                     # radio en píxeles
lcx    = px_s                                     # centro del bucle: mismo Re
lcy    = py_s - loop_r                            # por encima del marcador (py menor = arriba)
t_loop = range(Float64(π)/2, Float64(π)/2 + 1.5Float64(π), length=55)   # 270°
loop_px = lcx .+ loop_r .* cos.(t_loop)
loop_py = lcy .+ loop_r .* sin.(t_loop)
plot!(plt1, loop_px, loop_py;
    color = :white, linewidth = 2.0, label = "", arrow = :arrow)

# Etiqueta al lado del bucle
annotate!(plt1, px_s + loop_r + 8, py_s - loop_r,
    text("z*=−1/2", :white, :left, 10))

savefig(plt1, joinpath(outdir, "cycle_parabolic_p2.png"))
println("  → cycle_parabolic_p2.png\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Mapa 2:  f(z) = z² − 5/4
#   2-ciclo parabólico: z₁=(−1+√2)/2 ≈ 0.207,  z₂=(−1−√2)/2 ≈ −1.207
#   mult. de f²  en z₁:  4(c+1) = 4·(−1/4) = −1
#   En la frontera bulbo p=2 / bulbo p=4 del conjunto de Mandelbrot.
#   Radio del semicírculo: √2/2 ≈ 0.707  <  0.8 = y2hi  → cabe en la ventana.
#   Ventana:  [−1.8, 1.8] × [−0.8, 0.8]
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
println("=== Mapa 2:  f(z) = z²−5/4  ===")
cv2   = -5.0/4.0
cnum2 = complex([cv2, 0.0, 1.0])
cden2 = complex([1.0, 0.0, 0.0])

res2 = critical_orbit_analysis(cnum2, cden2; detect_parabolic=true, N_birkhoff=1200)
for r in res2
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

x2lo, x2hi, y2lo, y2hi = -1.8, 1.8, -0.8, 0.8

# Puntos del 2-ciclo (soluciones de z²+z+(c+1)=0, c=−5/4):
z1_re = (-1.0 + sqrt(2.0)) / 2.0   # ≈  0.2071  (f: z1 → z2)
z2_re = (-1.0 - sqrt(2.0)) / 2.0   # ≈ −1.2071  (f: z2 → z1)

plt2 = basin_portrait(cnum2, cden2;
    xinterval = (x2lo, x2hi),
    yinterval = (y2lo, y2hi),
    iter = 600, bands = 10.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath  = "")

# Marcadores en z₁ y z₂
px1, py1 = _z2px(z1_re, 0.0, x2lo, x2hi, y2lo, y2hi, expp)
px2, py2 = _z2px(z2_re, 0.0, x2lo, x2hi, y2lo, y2hi, expp)
scatter!(plt2, [px1, px2], [py1, py2];
    marker = :circle, markersize = 11,
    color  = :white, markerstrokecolor = :black, markerstrokewidth = 1.5,
    label  = "")
annotate!(plt2, px1 + 12, py1 - 20, text("z₁", :white, 11, :bold))
annotate!(plt2, px2 - 12, py2 - 20, text("z₂", :white, 11, :bold, :right))

# Arcos semicirculares con flecha:  z₁→z₂ (arriba)  y  z₂→z₁ (abajo)
pxu, pyu = _half_arc(z1_re, z2_re, x2lo, x2hi, y2lo, y2hi, expp; above=true)
plot!(plt2, pxu, pyu; color = :white, linewidth = 2.0, label = "", arrow = :arrow)
pxd, pyd = _half_arc(z1_re, z2_re, x2lo, x2hi, y2lo, y2hi, expp; above=false)
plot!(plt2, pxd, pyd; color = :white, linewidth = 2.0, label = "", arrow = :arrow)

savefig(plt2, joinpath(outdir, "cycle_parabolic_p4.png"))
println("  → cycle_parabolic_p4.png\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Control:  f(z) = z² − 1  (Basílica)
#   2-ciclo superatractor: 0 → −1 → 0,  mult. = 0
#   Ventana:  [−1.8, 1.8] × [−0.9, 0.9]  (simétrica, z→−z)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
println("=== Control:  f(z) = z²−1  (Basílica)  ===")
cnum3 = complex([-1.0, 0.0, 1.0])
cden3 = complex([ 1.0, 0.0, 0.0])

res3 = critical_orbit_analysis(cnum3, cden3; detect_parabolic=true)
for r in res3
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

x3lo, x3hi, y3lo, y3hi = -1.8, 1.8, -0.9, 0.9

# 2-ciclo: 0 → −1 → 0  (el punto crítico 0 ES el primer punto del ciclo)
zb1_re, zb2_re = 0.0, -1.0

plt3 = basin_portrait(cnum3, cden3;
    xinterval = (x3lo, x3hi),
    yinterval = (y3lo, y3hi),
    iter = 400, bands = 8.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath  = "")

px_b1, py_b1 = _z2px(zb1_re, 0.0, x3lo, x3hi, y3lo, y3hi, expp)
px_b2, py_b2 = _z2px(zb2_re, 0.0, x3lo, x3hi, y3lo, y3hi, expp)
scatter!(plt3, [px_b1, px_b2], [py_b1, py_b2];
    marker = :circle, markersize = 11,
    color  = :white, markerstrokecolor = :black, markerstrokewidth = 1.5,
    label  = "")
annotate!(plt3, px_b1 + 12, py_b1 - 20, text("0",  :white, 11, :bold))
annotate!(plt3, px_b2 + 12, py_b2 - 20, text("−1", :white, 11, :bold))

# Arcos cíclicos: 0 → −1 (arriba) y −1 → 0 (abajo)
pxu3, pyu3 = _half_arc(zb1_re, zb2_re, x3lo, x3hi, y3lo, y3hi, expp; above=true)
plot!(plt3, pxu3, pyu3; color = :white, linewidth = 2.0, label = "", arrow = :arrow)
pxd3, pyd3 = _half_arc(zb1_re, zb2_re, x3lo, x3hi, y3lo, y3hi, expp; above=false)
plot!(plt3, pxd3, pyd3; color = :white, linewidth = 2.0, label = "", arrow = :arrow)

savefig(plt3, joinpath(outdir, "cycle_attract_basilica.png"))
println("  → cycle_attract_basilica.png\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gráfico comparativo de promedios de Birkhoff (hasta N=15000)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
println("=== Gráfico de promedios de Birkhoff ===")

c_att   = -1.0/2.0
cnum_a  = complex([c_att, 0.0, 1.0])
cden_a  = complex([1.0,   0.0, 0.0])
hpair_a = RiemannBasins._pair_of_functions(cnum_a, cden_a)
hpair1  = RiemannBasins._pair_of_functions(cnum1, cden1)
hpair2  = RiemannBasins._pair_of_functions(cnum2, cden2)

cpt = RiemannBasins._homogeneous_normalize((0.0+0.0im, 1.0+0.0im))

Ns = vcat(collect(100:100:2000), collect(3000:500:15000))

B_p2  = RiemannBasins._birkhoff_partial_averages(hpair1, cnum1, cden1, 2, cpt, Ns)
B_p4  = RiemannBasins._birkhoff_partial_averages(hpair2, cnum2, cden2, 2, cpt, Ns)
B_att = RiemannBasins._birkhoff_partial_averages(hpair_a, cnum_a, cden_a, 2, cpt, Ns)

A_p2      = B_p2[1] * 100.0 / log(100.0)
theory_p2 = [A_p2 * log(n) / n for n in Ns]

println("  B_N parabólico z²-3/4  en N=100:   $(B_p2[1])")
println("  B_N parabólico z²-3/4  en N=15000: $(B_p2[end])")
println("  B_N parabólico z²-5/4  en N=15000: $(B_p4[end])")
println("  B_N atractor   z²-1/2  en N=15000: $(B_att[end])")

ylim_lo = min(minimum(B_att) - 0.02, minimum(B_p2) - 0.01, -0.08)

p_bk = plot(Ns, B_p2;
    label      = "f(z) = z^2 - 3/4  (parabólico, p=1)",
    color      = :steelblue, linewidth = 2,
    xlabel     = "N",
    ylabel     = "B_N  =  (1/N) sum log f#(z_k)",
    title      = "Promedios de Birkhoff: parabólico vs atractor",
    legend     = :bottomright,
    size       = (820, 500),
    ylims      = (ylim_lo, 0.01),
    framestyle = :box)
plot!(p_bk, Ns, B_p4;
    label = "f(z) = z^2 - 5/4  (parabólico, p=2)",
    color = :darkorange, linewidth = 2)
plot!(p_bk, Ns, B_att;
    label = "f(z) = z^2 - 1/2  (atractor)",
    color = :seagreen, linewidth = 2)
plot!(p_bk, Ns, theory_p2;
    label = "A*log(N)/N  (teorico)",
    color = :steelblue, linestyle = :dash, linewidth = 1.5)
hline!(p_bk, [0.0]; color = :black, linewidth = 1, linestyle = :dot, label = "")

savefig(p_bk, joinpath(outdir, "cycle_birkhoff_comparison.png"))
println("  → cycle_birkhoff_comparison.png\n")
println("Hecho.")
