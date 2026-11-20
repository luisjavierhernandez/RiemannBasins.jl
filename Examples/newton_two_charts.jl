# newton_two_charts.jl
# Side-by-side basin portraits in both stereographic charts
# for Newton's method applied to  z³ − 1

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors

# Newton's method for z³ − 1:  f(z) = (2z³ + 1) / (3z²)
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z³
den = complex([0.0, 0.0, 3.0, 0.0])   # 3z²

# expprecision=9: Chart 1 → 1536×1536 grid, Chart 2 → 512×512
#   (with expprecision=8 Chart 2 was only 256×256, visually blurry in the right panel)
outdir = @__DIR__
artdir = joinpath(@__DIR__, "..", "docs", "articles",
                  "00_gaceta_cuencas_atractoras", "figures")

# Standard palette
plt = two_chart_portrait(num, den;
    expprecision = 9,
    bands        = 8.0,
    colorblind   = false,
    show_plot    = false,
    savepath     = joinpath(outdir, "newton_two_charts.png"))
println("Saved: newton_two_charts.png")

# Copy to article figures directory
cp(joinpath(outdir, "newton_two_charts.png"),
   joinpath(artdir, "newton_two_charts.png"); force=true)
println("Copied to: $artdir/newton_two_charts.png")

# Colorblind-friendly version
plt_cb = two_chart_portrait(num, den;
    expprecision = 9,
    bands        = 8.0,
    colorblind   = true,
    show_plot    = false,
    savepath     = joinpath(outdir, "newton_two_charts_cb.png"))
println("Saved: newton_two_charts_cb.png")

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points ($(length(finitelist))):")
for (k, fp) in enumerate(finitelist)
    println("  fp[$k] = $fp")
end
