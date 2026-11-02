# newton_animation.jl
# GIF animation: Newton's method for z³ − c  as c varies in [0.5, 2.0]

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using RiemannBasins, Colors

# For each value of c, compute Newton's method for z³ − c:
#   f(z) = (2z³ + c) / (3z²)
make_coeffs(c) = (complex([c, 0.0, 0.0, 2.0]),
                  complex([0.0, 0.0, 3.0, 0.0]))

anim = animate_basins(make_coeffs, range(0.5, 2.0, length=30);
    expprecision = 7,      # 128×128 grid — faster for animation
    fps          = 6,
    colorblind   = false,
    filename     = joinpath(@__DIR__, "newton_c_animation.gif"))

println("GIF saved to Examples/newton_c_animation.gif")
