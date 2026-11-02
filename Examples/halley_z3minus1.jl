# Activate the RiemannBasins project before loading the package
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# ============================================================
#  Halley's method for  z³ − 1
#
#  g(z) = z · (z³ + 2) / (2z³ + 1)    (degree-4 rational)
#
#  Numerator:   z⁴ + 2z  →  [0, 2, 0, 0, 1]
#  Denominator: 2z³ + 1  →  [1, 0, 0, 2, 0]
# ============================================================
using RiemannBasins

num = complex([0.0, 2.0, 0.0, 0.0, 1.0])
den = complex([1.0, 0.0, 0.0, 2.0, 0.0])

plt = basin_portrait(
    num, den;
    expprecision  = 8,
    colorstrategy = "positionplusiteration")

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points ($(length(finitelist))):")
for (k, fp) in enumerate(finitelist); println("  fp[$k] = $fp"); end
