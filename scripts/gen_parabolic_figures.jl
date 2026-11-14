# gen_parabolic_figures.jl
#
# Genera retratos de cuencas para ejemplos parabólicos, con detección
# automática vía promedios de Birkhoff (Birkhoff average acceleration).
#
# Uso: julia --project scripts/gen_parabolic_figures.jl
#
# Figuras generadas en docs/articles/02_fatou_hopf/figures/:
#   parabolic_z_plus_z2.png       — f(z) = z + z²  (punto fijo parabólico en 0)
#   parabolic_z2_plus_quarter.png — f(z) = z² + 1/4 (cúspide de Mandelbrot, c = 1/4)
#   parabolic_birkhoff_plot.png   — gráfico de promedios de Birkhoff para ambos casos
#
# Control (atractores): Newton z³−1 y la basílica z²−1

using RiemannBasins, Colors, Plots, Printf

outdir = joinpath(@__DIR__, "..", "docs", "articles", "02_fatou_hopf", "figures")
mkpath(outdir)

# ── Convención de coeficientes ──────────────────────────────────────────────
# cnum[k+1] = coeficiente de z^k · t^(d−k)  (grado d, cte primero)
# Para f(z) = p(z)/q(z), cnum = coefs de p, cden = coefs de q (mismo grado d).

# ============================================================
# Ejemplo 1:  f(z) = z + z²   (parabólico en z = 0)
# ============================================================
println("=== Ejemplo 1:  f(z) = z + z²  ===")
# p(z) = z + z²  →  d=2, cnum = [0, 1, 1],  cden = [1, 0, 0]
cnum1 = complex([0.0, 1.0, 1.0])
cden1 = complex([1.0, 0.0, 0.0])

# Análisis de órbitas críticas con detección parabólica
res1 = critical_orbit_analysis(cnum1, cden1; detect_parabolic = true, N_birkhoff = 800)
println("  Puntos críticos:")
for r in res1
    @printf("    [%.4f+%.4fi : 1]  →  %s  (c = %.6f)\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

println("  Generando retrato de cuencas …")
basin_portrait(
    cnum1, cden1;
    xinterval     = (-2.0, 1.5),
    yinterval     = (-1.5, 1.5),
    iter          = 400,
    aproxprecision = 8,
    colorstrategy = "positionplusiteration",
    bands         = 8.0,
    show_plot     = false,
    savepath      = joinpath(outdir, "parabolic_z_plus_z2.png"))
println("  → parabolic_z_plus_z2.png")

# ============================================================
# Ejemplo 2:  f(z) = z² + 1/4  (cúspide del conjunto de Mandelbrot)
# ============================================================
println("\n=== Ejemplo 2:  f(z) = z² + 1/4  ===")
# p(z) = z² + 1/4  →  d=2, cnum = [0.25, 0, 1],  cden = [1, 0, 0]
cnum2 = complex([0.25, 0.0, 1.0])
cden2 = complex([1.0,  0.0, 0.0])

res2 = critical_orbit_analysis(cnum2, cden2; detect_parabolic = true, N_birkhoff = 800)
println("  Puntos críticos:")
for r in res2
    @printf("    [%.4f+%.4fi : 1]  →  %s  (c = %.6f)\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

println("  Generando retrato de cuencas …")
basin_portrait(
    cnum2, cden2;
    xinterval     = (-2.0, 1.0),
    yinterval     = (-1.5, 1.5),
    iter          = 500,
    aproxprecision = 8,
    colorstrategy = "positionplusiteration",
    bands         = 8.0,
    show_plot     = false,
    savepath      = joinpath(outdir, "parabolic_z2_plus_quarter.png"))
println("  → parabolic_z2_plus_quarter.png")

# ============================================================
# Comparación: Newton z³ − 1  (atractor puro, control)
# ============================================================
println("\n=== Control:  Newton(z³ − 1) = (2z³ + 1) / (3z²)  ===")
# f(z) = (2z³ + 1) / (3z²)  →  d=3
# p(z) = 2z³ + 1  →  cnum = [1, 0, 0, 2]
# q(z) = 3z²      →  cden = [0, 0, 3, 0]
cnum_N = complex([1.0, 0.0, 0.0, 2.0])
cden_N = complex([0.0, 0.0, 3.0, 0.0])

resN = critical_orbit_analysis(cnum_N, cden_N; detect_parabolic = true, N_birkhoff = 800)
println("  Puntos críticos:")
for r in resN
    @printf("    [%.4f+%.4fi : 1]  →  %s  (c = %.6f)\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        string(r.label), r.lyapunov_constant)
end

# ============================================================
# Gráfico de promedios de Birkhoff (ilustración del artículo)
# ============================================================
println("\n=== Gráfico comparativo de promedios de Birkhoff ===")
hpair1 = RiemannBasins._pair_of_functions(cnum1, cden1)
hpair2 = RiemannBasins._pair_of_functions(cnum2, cden2)

# Órbita del punto crítico −1/2 para Ejemplo 1
cpt1 = RiemannBasins._homogeneous_normalize((-0.5 + 0.0im, 1.0 + 0.0im))
# Órbita del punto crítico 0 para Ejemplo 2
cpt2 = RiemannBasins._homogeneous_normalize((0.0 + 0.0im, 1.0 + 0.0im))

Ns  = collect(50:50:2000)
B1  = RiemannBasins._birkhoff_partial_averages(hpair1, cnum1, cden1, 2, cpt1, Ns)
B2  = RiemannBasins._birkhoff_partial_averages(hpair2, cnum2, cden2, 2, cpt2, Ns)

# Curva teórica: A·log(N)/N para ajuste visual
A1_fit = B1[1] * 50 / log(50)      # calibrado en N=50
A2_fit = B2[1] * 50 / log(50)
theory1 = [A1_fit * log(n) / n for n in Ns]
theory2 = [A2_fit * log(n) / n for n in Ns]

p_bk = plot(Ns, B1;
    label     = raw"$f(z)=z+z^2$",
    color     = :steelblue,
    linewidth = 2,
    xlabel    = raw"$N$",
    ylabel    = raw"$B_N = \frac{1}{N}\sum_{k=0}^{N-1}\log f^\#(z_k)$",
    title     = "Promedios de Birkhoff de log(f⁺) a lo largo de la órbita crítica",
    legend    = :bottomright,
    size      = (800, 480))

plot!(p_bk, Ns, B2;
    label     = raw"$f(z)=z^2+\tfrac{1}{4}$",
    color     = :firebrick,
    linewidth = 2)

plot!(p_bk, Ns, theory1;
    label     = raw"$A\log(N)/N$ (ajuste teórico)",
    color     = :steelblue,
    linestyle = :dash,
    linewidth = 1.5)

plot!(p_bk, Ns, theory2;
    label     = "",
    color     = :firebrick,
    linestyle = :dash,
    linewidth = 1.5)

hline!(p_bk, [0.0]; color = :black, linewidth = 1, label = "")

savefig(p_bk, joinpath(outdir, "parabolic_birkhoff_plot.png"))
println("  → parabolic_birkhoff_plot.png")

println("\nTodo completado.")
