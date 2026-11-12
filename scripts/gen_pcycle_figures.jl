# gen_pcycle_figures.jl
# Genera las figuras para draft_pcycle_basins.tex
# Uso: julia --project scripts/gen_pcycle_figures.jl
#
# Figuras generadas en docs/figures/:
#   pcycle_basilica_2cycle.png   — basílica z²−1: cuencas del 2-ciclo
#   pcycle_newton_z3.png         — Newton z³−1: dos cartas (puntos fijos)
#   pcycle_z2_fixed.png          — f(z)=z²: cuenca del punto fijo superatractor en 0
#   pcycle_airplane_3cycle.png   — avión z²+c: cuencas del 3-ciclo superatractor

using RiemannBasins, Colors

outdir = joinpath(@__DIR__, "..", "docs", "figures")
mkpath(outdir)

println("=== Figura 1: basílica z²−1  (cuencas del 2-ciclo) ===")
# f(z) = z² − 1   →   num = [−1, 0, 1],  den = [1]
num_bas = complex([-1.0,  0.0, 1.0])
den_bas = complex([ 1.0])

basin_portrait_pcycle(
    num_bas, den_bas, 2;
    xinterval     = (-1.8, 1.8),
    yinterval     = (-1.8, 1.8),
    iter          = 300,
    colorstrategy = "positionplusiteration",
    show_plot     = false,
    savepath      = joinpath(outdir, "pcycle_basilica_2cycle.png"))
println("  → pcycle_basilica_2cycle.png")

println("=== Figura 1b: basílica z²−1  (subcuencas LCHab) ===")
basin_portrait_pcycle(
    num_bas, den_bas, 2;
    xinterval     = (-1.8, 1.8),
    yinterval     = (-1.8, 1.8),
    iter          = 300,
    colorstrategy       = "pcycle_lch",
    hue_band            = 120.0,
    show_cyclic_palette = true,
    show_plot           = false,
    savepath            = joinpath(outdir, "pcycle_basilica_2cycle_lch.png"))
println("  → pcycle_basilica_2cycle_lch.png")

println("=== Figura 2: Newton para z³−1  (dos cartas) ===")
# N(z) = (2z³+1)/(3z²)   →   num = [1,0,0,2],  den = [0,0,3,0]
num_new = complex([1.0, 0.0, 0.0, 2.0])
den_new = complex([0.0, 0.0, 3.0, 0.0])

two_chart_portrait(
    num_new, den_new;
    xinterval = (-1.8, 1.8),
    yinterval = (-1.8, 1.8),
    iter      = 300,
    show_plot = false,
    savepath  = joinpath(outdir, "pcycle_newton_z3.png"))
println("  → pcycle_newton_z3.png")

println("=== Figura 3: f(z) = z²  (cuenca del punto fijo superatractor en 0) ===")
# f(z) = z²   →   num = [0,0,1],  den = [1]
num_z2 = complex([0.0, 0.0, 1.0])
den_z2 = complex([1.0])

basin_portrait(
    num_z2, den_z2;
    xinterval = (-1.5, 1.5),
    yinterval = (-1.5, 1.5),
    iter      = 300,
    show_plot = false,
    savepath  = joinpath(outdir, "pcycle_z2_fixed.png"))
println("  → pcycle_z2_fixed.png")

println("=== Figura 4: avión z²+c  (cuencas del 3-ciclo) ===")
# c es la raíz real de c³+2c²+c+1=0  ≈ −1.7549
# Se calcula con Newton en pocas iteraciones:
c_air = let x = -1.75
    for _ in 1:60; x -= (x^3 + 2x^2 + x + 1) / (3x^2 + 4x + 1) end; x
end
println("  c_airplane = ", c_air)
# f(z) = z² + c   →   num = [c, 0, 1],  den = [1]
num_air = complex([c_air, 0.0, 1.0])
den_air = complex([1.0])

basin_portrait_pcycle(
    num_air, den_air, 3;
    xinterval     = (-2.0, 2.0),
    yinterval     = (-0.75, 0.75),
    iter          = 400,
    colorstrategy = "positionplusiteration",
    show_plot     = false,
    savepath      = joinpath(outdir, "pcycle_airplane_3cycle.png"))
println("  → pcycle_airplane_3cycle.png")


println("=== Figura 5: avión z²+c  (cuencas del 3-ciclo, coloración LCHab) ===")
# Same map as Fig 4, but rendered with the perceptual pcycle_lch strategy
basin_portrait_pcycle(
    num_air, den_air, 3;
    xinterval     = (-2.0, 2.0),
    yinterval     = (-0.75, 0.75),
    iter          = 400,
    colorstrategy       = "pcycle_lch",
    hue_band            = 120.0,
    show_cyclic_palette = true,
    show_plot           = false,
    savepath            = joinpath(outdir, "pcycle_airplane_3cycle_lch.png"))
println("  → pcycle_airplane_3cycle_lch.png")



println("\n=== Figura 7: centro de período 5  (5-ciclo superatractor, LCHab) ===")
# c5: centro del componente hiperbólico de período 5 en el conjunto de Mandelbrot
# (ángulo 1/5 en la cardioide principal)
c5 = let x = complex(-0.5024, 0.5629)
    for _ in 1:100
        z = 0+0im; dz = 0+0im
        for _ in 1:5
            z, dz = z^2 + x, 2z*dz + 1
        end
        x -= z/dz
    end; x
end
println("  c5 = ", c5)   # ≈ -0.50434...+0.56277...i
num_p5 = ComplexF64.([c5, 0.0, 1.0])
den_p5 = ComplexF64.([1.0])

basin_portrait_pcycle(
    num_p5, den_p5, 5;
    xinterval     = (-2.0,  1.1),
    yinterval     = (-1.2,  1.9),
    iter          = 500,
    colorstrategy       = "pcycle_lch",
    hue_band            = 150.0,
    show_cyclic_palette = true,
    show_plot           = false,
    savepath            = joinpath(outdir, "pcycle_period5_lch.png"))
println("  → pcycle_period5_lch.png")

println("\nListo. Figuras en docs/figures/")

println("=== Figura 6: conejo de Douady z²+c  (3-ciclo superatractor, LCHab) ===")
# c es la raíz compleja de c³+2c²+c+1=0  ≈ −0.1226+0.7449i
# (misma ecuación que el avión, pero raíz compleja)
c_rabbit = let x = complex(-0.12, 0.75)
    for _ in 1:80
        x -= (x^3 + 2x^2 + x + 1) / (3x^2 + 4x + 1)
    end; x
end
println("  c_rabbit = ", c_rabbit)
num_rab = complex([c_rabbit, 0.0, 1.0])
den_rab = complex([1.0])

basin_portrait_pcycle(
    num_rab, den_rab, 3;
    xinterval     = (-1.5, 1.5),
    yinterval     = (-1.5, 1.5),
    iter          = 400,
    colorstrategy       = "pcycle_lch",
    hue_band            = 120.0,
    show_cyclic_palette = true,
    show_plot           = false,
    savepath            = joinpath(outdir, "pcycle_rabbit_3cycle_lch.png"))
println("  → pcycle_rabbit_3cycle_lch.png")
