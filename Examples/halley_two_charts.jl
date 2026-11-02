# halley_two_charts.jl
# Side-by-side basin portraits in both stereographic charts
# for Halley's method applied to  z³ − 1
#
# (Replaces halley_sphere.jl — sphere visualization has been removed
#  from RiemannBasins v1.0.0.  The two-chart portrait shows the full
#  Riemann sphere P¹(ℂ) through both stereographic coordinate charts.)

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors

# Halley's method for z³ − 1:  f(z) = (2z⁴ + z) / (2z³ + 1)
num = complex([0.0, 2.0, 0.0, 0.0, 1.0])   # 2z + z⁴
den = complex([1.0, 0.0, 0.0, 2.0, 0.0])   # 1 + 2z³

plt = two_chart_portrait(num, den;
    expprecision = 8,
    colorblind   = false,
    savepath     = joinpath(@__DIR__, "halley_two_charts.png"))

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points: $(length(finitelist))")
println("Saved: halley_two_charts.png")
