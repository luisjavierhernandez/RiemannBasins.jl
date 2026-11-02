# relief_3d.jl
# 3-D fractal relief portraits of basins of attraction on the Riemann sphere
#
# The surface looks nearly flat (2-D with subtle embossing):
#   • colour  → attractor identity (different hue per basin)
#   • height  → smooth iteration count with cyclic banding (fractal ridges)
#
# relief_scale controls the total z range (default 0.03, very low).
# camera_angle=(30,82) is almost top-down so the relief reads as texture.
#
# Run from the repository root:
#   julia --project Examples/relief_3d.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors

# ── Example 1: Newton z³−1, low relief (default settings) ───────────────────
num3 = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z³
den3 = complex([0.0, 0.0, 3.0, 0.0])   # 3z²

println("▶  Example 1 — Newton z³−1  (low relief, near top-down)")
basin_relief_3d(num3, den3;
    expprecision  = 8,
    iter          = 200,
    bands         = 8.0,
    relief_scale  = 0.03,     # very low — looks almost 2-D
    ridge_fraction= 0.15,     # ridges use 15% of each colour zone
    camera_angle  = (30, 82), # near top-down (82° elevation)
    savepath      = joinpath(@__DIR__, "relief_newton_z3.png"))
println("   Saved: relief_newton_z3.png")

# ── Example 2: slightly more relief to appreciate the 3-D effect ────────────
println("▶  Example 2 — Newton z³−1  (moderate relief, oblique view)")
basin_relief_3d(num3, den3;
    expprecision  = 8,
    iter          = 200,
    bands         = 8.0,
    relief_scale  = 0.15,
    ridge_fraction= 0.20,
    camera_angle  = (40, 65),
    savepath      = joinpath(@__DIR__, "relief_newton_z3_oblique.png"))
println("   Saved: relief_newton_z3_oblique.png")

# ── Example 3: colorblind-safe palette ──────────────────────────────────────
println("▶  Example 3 — Newton z³−1  (Okabe-Ito / tab10, low relief)")
basin_relief_3d(num3, den3;
    expprecision  = 8,
    iter          = 200,
    bands         = 8.0,
    relief_scale  = 0.03,
    camera_angle  = (30, 82),
    colorblind    = true,
    savepath      = joinpath(@__DIR__, "relief_newton_z3_cb.png"))
println("   Saved: relief_newton_z3_cb.png")

# ── Example 4: Newton z⁴−1, four basins ─────────────────────────────────────
num4 = complex([1.0, 0.0, 0.0, 0.0, 3.0])
den4 = complex([0.0, 0.0, 0.0, 4.0, 0.0])

println("▶  Example 4 — Newton z⁴−1  (four basins, low relief)")
basin_relief_3d(num4, den4;
    expprecision  = 8,
    iter          = 200,
    bands         = 8.0,
    relief_scale  = 0.03,
    camera_angle  = (30, 82),
    savepath      = joinpath(@__DIR__, "relief_newton_z4.png"))
println("   Saved: relief_newton_z4.png")

println("\nDone.")
