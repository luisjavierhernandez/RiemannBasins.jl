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
    expprecision = 9,   # Chart 2: 512×512 (vs 256 con exp=8 → borroso)
    colorblind   = false,
    savepath     = joinpath(@__DIR__, "halley_two_charts.png"))

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points: $(length(finitelist))")
println("Saved: halley_two_charts.png")

# Copy to article figures directory
artdir = joinpath(@__DIR__, "..", "docs", "articles",
                  "00_gaceta_cuencas_atractoras", "figures")
cp(joinpath(@__DIR__, "halley_two_charts.png"),
   joinpath(artdir, "halley_z3_two_charts.png"); force=true)
println("Copied to: $artdir/halley_z3_two_charts.png")
