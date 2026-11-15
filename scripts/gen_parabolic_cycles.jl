# gen_parabolic_cycles.jl  (v2 — ventanas y Birkhoff corregidos)
#
# Uso: julia --project scripts/gen_parabolic_cycles.jl
#
# Figuras generadas en docs/articles/02_fatou_hopf/figures/:
#
#   cycle_parabolic_p2.png      f(z)=z²-3/4  : punto fijo parabólico (mult. -1 en f²)
#   cycle_parabolic_p4.png      f(z)=z²-5/4  : 2-ciclo parabólico (mult. -1)
#   cycle_attract_basilica.png  f(z)=z²-1    : 2-ciclo superatractor (control)
#   cycle_birkhoff_comparison.png: B_N para parabólico vs atractor (hasta N=15000)

using RiemannBasins, Colors, Plots, Printf

outdir = joinpath(@__DIR__, "..", "docs", "articles", "02_fatou_hopf", "figures")
mkpath(outdir)

# ── Convención: cnum[k+1] = coef de z^k·t^(d−k), cte primero ──────────────

# ============================================================
# Mapa 1:  f(z) = z² - 3/4
#   Punto fijo z = -1/2  con mult. f'(-1/2) = -1
#   → parabólico para f²,  frontera cardioide / bulbo p=2
#
#   Ventana centrada en z = -1/2:
#     Real: [-1.5,  0.5],  Im: [-1.0, 1.0]
# ============================================================
println("=== Mapa 1:  f(z) = z²-3/4  ===")
c1    = -3/4
cnum1 = complex([ c1, 0.0, 1.0])
cden1 = complex([1.0, 0.0, 0.0])

res1 = critical_orbit_analysis(cnum1, cden1; detect_parabolic=true, N_birkhoff=1200)
for r in res1
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

basin_portrait(cnum1, cden1;
    xinterval = (-1.5,  0.5),
    yinterval = (-1.0,  1.0),
    iter = 600, bands = 10.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath = joinpath(outdir, "cycle_parabolic_p2.png"))
println("  → cycle_parabolic_p2.png\n")

# ============================================================
# Mapa 2:  f(z) = z² - 5/4
#   2-ciclo {z₁≈ 0.207, z₂≈ -1.207},  mult. 4(c+1) = -1
#   → 2-ciclo parabólico,  frontera bulbo p=2 / bulbo p=4
#
#   Ventana centrada entre los dos puntos del ciclo:
#     midpoint = (-1.207 + 0.207)/2 ≈ -0.5
#     Real: [-2.0,  0.8],  Im: [-0.8, 0.8]
# ============================================================
println("=== Mapa 2:  f(z) = z²-5/4  ===")
c2    = -5/4
cnum2 = complex([ c2, 0.0, 1.0])
cden2 = complex([1.0, 0.0, 0.0])

res2 = critical_orbit_analysis(cnum2, cden2; detect_parabolic=true, N_birkhoff=1200)
for r in res2
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

basin_portrait(cnum2, cden2;
    xinterval = (-2.0,  0.8),
    yinterval = (-0.8,  0.8),
    iter = 600, bands = 10.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath = joinpath(outdir, "cycle_parabolic_p4.png"))
println("  → cycle_parabolic_p4.png\n")

# ============================================================
# Control:  f(z) = z² - 1  (Basílica)
#   2-ciclo {0, -1},  mult. = 0  (superatractor)
#
#   Ventana centrada entre 0 y -1:
#     Real: [-2.0,  0.8],  Im: [-0.9, 0.9]
# ============================================================
println("=== Control:  f(z) = z²-1  (Basílica)  ===")
cnum3 = complex([-1.0, 0.0, 1.0])
cden3 = complex([ 1.0, 0.0, 0.0])

res3 = critical_orbit_analysis(cnum3, cden3; detect_parabolic=true)
for r in res3
    @printf("  crit z≈%.4f+%.4fi  →  %-14s  c=%.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

basin_portrait(cnum3, cden3;
    xinterval = (-1.8,  1.8),
    yinterval = (-0.9,  0.9),
    iter = 400, bands = 8.0,
    colorstrategy = "positionplusiteration",
    show_plot = false,
    savepath = joinpath(outdir, "cycle_attract_basilica.png"))
println("  → cycle_attract_basilica.png\n")

# ============================================================
# Gráfico comparativo de promedios de Birkhoff (hasta N=15000)
#
#   Parabólico z²-3/4:  B_N → 0 como log(N)/N
#   Parabólico z²-5/4:  B_N → 0 como log(N)/N  (ciclo de orden 2)
#   Atractor   z²-1/2:  B_N → cte < 0  (convergencia geométrica)
#
#   Etiquetas en texto plano para evitar problemas de renderizado LaTeX.
# ============================================================
println("=== Gráfico de promedios de Birkhoff ===")

# Atractor puro: c = -1/2  (dentro de la cardioide principal)
#   Punto fijo z* ≈ -0.366,  mult. ≈ -0.732
c_att   = -1/2
cnum_a  = complex([c_att, 0.0, 1.0])
cden_a  = complex([1.0,   0.0, 0.0])
hpair_a = RiemannBasins._pair_of_functions(cnum_a, cden_a)
hpair1  = RiemannBasins._pair_of_functions(cnum1, cden1)
hpair2  = RiemannBasins._pair_of_functions(cnum2, cden2)

# Punto crítico z = 0 (en coordenadas homogéneas)
cpt = RiemannBasins._homogeneous_normalize((0.0+0.0im, 1.0+0.0im))

# Checkpoints hasta N = 15000:
#   densos al principio (capturan el transitorio) y más espaciados al final
Ns = vcat(collect(100:100:2000), collect(3000:500:15000))

B_p2  = RiemannBasins._birkhoff_partial_averages(hpair1, cnum1, cden1, 2, cpt, Ns)
B_p4  = RiemannBasins._birkhoff_partial_averages(hpair2, cnum2, cden2, 2, cpt, Ns)
B_att = RiemannBasins._birkhoff_partial_averages(hpair_a, cnum_a, cden_a, 2, cpt, Ns)

# Ajuste teórico A·log(N)/N calibrado en N=100
A_p2  = B_p2[1]  * 100.0 / log(100.0)
theory_p2 = [A_p2 * log(n) / n for n in Ns]

println("  B_N parabólico z²-3/4  en N=100:   $(B_p2[1])")
println("  B_N parabólico z²-3/4  en N=15000: $(B_p2[end])")
println("  B_N parabólico z²-5/4  en N=15000: $(B_p4[end])")
println("  B_N atractor   z²-1/2  en N=15000: $(B_att[end])")

# Rango vertical: desde el mínimo de B_att (que ya convergió) hasta 0.01
ylim_lo = min(minimum(B_att) - 0.02, minimum(B_p2) - 0.01)
ylim_lo = min(ylim_lo, -0.08)   # al menos -0.08 para ver la escala

p_bk = plot(Ns, B_p2;
    label     = "f(z) = z^2 - 3/4  (parabólico, p=1)",
    color     = :steelblue, linewidth = 2,
    xlabel    = "N",
    ylabel    = "B_N  =  (1/N) sum log f#(z_k)",
    title     = "Promedios de Birkhoff: parabólico vs atractor",
    legend    = :bottomright,
    size      = (820, 500),
    ylims     = (ylim_lo, 0.01),
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
