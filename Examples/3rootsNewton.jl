# Activate the RiemannBasins project before loading the package
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# ============================================================
#  Newton's method for  z³ − 1
#
#  f(z) = (2z³ + 1) / (3z²)
#
#  Fixed points (roots of z³−1): 1,  e^{2πi/3},  e^{4πi/3}
# ============================================================
using RiemannBasins
using Colors

# Coefficient lists  (constant term first, as in Polynomials.jl)
num = complex([1.0, 0.0, 0.0, 2.0])   #  1 + 0·z + 0·z² + 2·z³
den = complex([0.0, 0.0, 3.0, 0.0])   #  0 + 0·z + 3·z² + 0·z³

plt = basin_portrait(
    num, den;
    expprecision  = 8,          # 256 points per unit side
    colorstrategy = "positionplusiteration")

superfixed, infixed, finitelist = fixed_points(num, den)
println("Fixed points found ($(length(finitelist))):")
for (k, fp) in enumerate(finitelist)
    println("  fp[$k] = $fp")
end
