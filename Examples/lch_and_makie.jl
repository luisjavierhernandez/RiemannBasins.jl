# Examples/lch_and_makie.jl
# Demonstrates the two aesthetic improvements added in v1.1:
#   1. LCHab colour arcs in basin_portrait  (perceptually uniform hues)
#   2. basin_relief_3d_makie                (GLMakie, real OpenGL lighting)
#
# Run from the package root:
#   julia --project Examples/lch_and_makie.jl

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using RiemannBasins, Colors, Plots

# Newton's method for z^3 - 1
num3 = complex([1.0, 0.0, 0.0, 2.0])   # 2z^3 + 1
den3 = complex([0.0, 0.0, 3.0, 0.0])   # 3z^2

# ── 1. LCHab palette — basin_portrait ────────────────────────────────────────
# Hues distributed as a uniform arc in LCHab space:
# all basins have the same perceived brightness and saturation.

basin_portrait(num3, den3;
    iter = 200, bands = 50.0,
    savepath = "portrait_lch_default.png")

# Different seed colours rotate the whole palette while keeping uniformity
basin_portrait(num3, den3;
    gseed    = RGB{Float64}(0.9, 0.3, 0.1),   # warm red seed
    iter     = 200, bands = 50.0,
    savepath = "portrait_lch_warm.png")

basin_portrait(num3, den3;
    gseed    = RGB{Float64}(0.1, 0.7, 0.4),   # green seed
    iter     = 200, bands = 50.0,
    savepath = "portrait_lch_green.png")

# ── 2. GLMakie 3D relief ─────────────────────────────────────────────────────
# Requires GLMakie to be installed in the environment.
# The extension loads automatically when you `using GLMakie`.

using GLMakie   # activates RiemannBasinsGLMakieExt

# Standard view — defaults: bands=50.0, relief_scale=0.01
basin_relief_3d_makie(num3, den3; iter=200)

# Oblique view — shows the stepped mesas and ridge texture
basin_relief_3d_makie(num3, den3;
    iter         = 200,
    camera_angle = (40, 55),
    savepath     = "relief_makie_oblique.png")

# Newton z^4 - 1  (four basins, four-fold symmetry)
num4 = complex([1.0, 0.0, 0.0, 0.0, 3.0])   # 3z^4 + 1
den4 = complex([0.0, 0.0, 0.0, 4.0, 0.0])   # 4z^3

basin_relief_3d_makie(num4, den4;
    iter         = 200,
    camera_angle = (30, 75),
    savepath     = "relief_makie_z4.png")
