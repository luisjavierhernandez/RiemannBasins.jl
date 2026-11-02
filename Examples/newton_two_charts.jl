# newton_two_charts.jl
# Side-by-side basin portraits in both stereographic charts
# for Newton's method applied to  z³ − 1

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors

# Newton's method for z³ − 1:  f(z) = (2z³ + 1) / (3z²)
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z³
den = complex([0.0, 0.0, 3.0, 0.0])   # 3z²

# Standard palette
plt = two_chart_portrait(num, den;
    expprecision = 8,
    colorblind   = false,
    savepath     = joinpath(@__DIR__, "newton_two_charts.png"))

# Colorblind-friendly version
plt_cb = two_chart_portrait(num, den;
    expprecision = 8,
    colorblind   = true,
    savepath     = joinpath(@__DIR__, "newton_two_charts_cb.png"))

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points ($(length(finitelist))):")
for (k, fp) in enumerate(finitelist)
    println("  fp[$k] = $fp")
end
println("Saved: newton_two_charts.png  and  newton_two_charts_cb.png")
