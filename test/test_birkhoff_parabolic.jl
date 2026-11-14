# test/test_birkhoff_parabolic.jl
#
# Tests and illustrations for the Birkhoff average acceleration tools
# introduced in src/RiemannBasins.jl for parabolic basin detection.
#
# Canonical parabolic example: f(z) = z/(1 - z²)   (degree 2 rational map)
#   Written homogeneously:  F(z,t) = zt,  G(z,t) = t² - z²
#   Critical points: roots of Wronskian W = F_z·G - F·G_z
#   The fixed point z=0 is parabolic with multiplier f'(0) = 1.
#
# Classic textbook example: f(z) = z + z² (polynomial, parabolic at 0)
#   Homogeneous: F = z(z+t) = z²+zt,  G = t²
#   i.e. numerator coefficients [0, 0, 1, 1] (const→high), denominator [0, 0, 1] padded
#   Wait — rational form: f(z)=z+z²  means p(z)=z+z², q(z)=1
#   Homogeneous P(z,t)=zt+z²=z(t+z), Q(z,t)=t²
#   cnum = [0, 1, 1] (coefficients of z⁰t², z¹t¹, z²t⁰ for deg-2 numerator? )
#   Let's use the convention in RiemannBasins: cnum[k+1] = coeff of z^k * t^(d-k)
#   For p(z) = z + z²: d=2, coefficients = [0, 1, 1]  (z⁰t², z¹t¹, z²t⁰)
#   For q(z) = 1:       pad to degree 2: [1, 0, 0]
#
# Run with: julia --project test/test_birkhoff_parabolic.jl

include("../src/RiemannBasins.jl")
using .RiemannBasins
using LinearAlgebra

println("=" ^ 60)
println("Birkhoff average acceleration — parabolic basin tests")
println("=" ^ 60)

# ----------------------------------------------------------------
# Helper: build hpair from coefficient lists (same as internal code)
# ----------------------------------------------------------------
function make_hpair(cnum, cden)
    return RiemannBasins._pair_of_functions(
        ComplexF64.(cnum), ComplexF64.(cden))
end

# ----------------------------------------------------------------
# Example 1:  f(z) = z + z²   (parabolic fixed point at z = 0)
# ----------------------------------------------------------------
#   p(z) = z + z²,  q(z) = 1  (degree 2 rational map)
#   cnum = [0, 1, 1],  cden = [1, 0, 0]  (constant-first convention)
println("\n--- Example 1:  f(z) = z + z²   ---")
cnum1 = ComplexF64.([0, 1, 1])
cden1 = ComplexF64.([1, 0, 0])
d1    = 2
hpair1 = RiemannBasins._pair_of_functions(cnum1, cden1)

# Starting point inside the parabolic basin of 0
start1 = (complex(0.1 + 0.05im), complex(1.0 + 0.0im))
start1 = RiemannBasins._homogeneous_normalize(start1[1], start1[2])

cps    = [100, 200, 400, 800]
B1     = RiemannBasins._birkhoff_partial_averages(hpair1, cnum1, cden1, d1,
                                                   start1, cps)
println("  Birkhoff averages at N = $cps:")
for (n, b) in zip(cps, B1)
    @printf("    B[%4d] = %+.6f\n", n, b)
end

B1_accel = RiemannBasins._aitken_accelerate(B1)
println("  Aitken-accelerated (last value): $(isempty(B1_accel) ? "N/A" : B1_accel[end])")

# Rate check: |B[N/4]| / |B[N]| should be > 1 for parabolic
ratio1 = abs(B1[2]) / abs(B1[end])
println("  Ratio |B[N/4]|/|B[N]| = $(@sprintf("%.3f", ratio1))  (> 1.3 → parabolic)")

# ----------------------------------------------------------------
# Example 2:  f(z) = z²   (superattracting fixed point at 0)
# ----------------------------------------------------------------
println("\n--- Example 2:  f(z) = z²   ---")
cnum2 = ComplexF64.([0, 0, 1])   # z²
cden2 = ComplexF64.([1, 0, 0])   # 1
d2    = 2
hpair2 = RiemannBasins._pair_of_functions(cnum2, cden2)

start2 = (complex(0.3 + 0.1im), complex(1.0 + 0.0im))
start2 = RiemannBasins._homogeneous_normalize(start2[1], start2[2])

B2 = RiemannBasins._birkhoff_partial_averages(hpair2, cnum2, cden2, d2, start2, cps)
println("  Birkhoff averages at N = $cps:")
for (n, b) in zip(cps, B2)
    @printf("    B[%4d] = %+.6f\n", n, b)
end
B2_accel = RiemannBasins._aitken_accelerate(B2)
println("  Aitken-accelerated (last value): $(isempty(B2_accel) ? "N/A" : B2_accel[end])")
ratio2 = abs(B2[2]) / abs(B2[end])
println("  Ratio |B[N/4]|/|B[N]| = $(@sprintf("%.3f", ratio2))  (≈ 1 → attracting/superattracting)")

# ----------------------------------------------------------------
# Example 3: critical_orbit_analysis on f(z) = z + z²
# ----------------------------------------------------------------
println("\n--- Example 3:  critical_orbit_analysis on f(z) = z + z²  ---")
res1 = critical_orbit_analysis(cnum1, cden1; detect_parabolic = true, N_birkhoff = 800)
for r in res1
    @printf("  critical_pt ≈ [%.4f+%.4fi : %.4f+%.4fi]  label = %-16s  c = %.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        real(r.critical_point[2]), imag(r.critical_point[2]),
        string(r.label), r.lyapunov_constant)
end

# ----------------------------------------------------------------
# Example 4: critical_orbit_analysis on f(z) = z²  (superattracting)
# ----------------------------------------------------------------
println("\n--- Example 4:  critical_orbit_analysis on f(z) = z²  ---")
res2 = critical_orbit_analysis(cnum2, cden2; detect_parabolic = true)
for r in res2
    @printf("  critical_pt ≈ [%.4f+%.4fi : %.4f+%.4fi]  label = %-16s  c = %.6f\n",
        real(r.critical_point[1]), imag(r.critical_point[1]),
        real(r.critical_point[2]), imag(r.critical_point[2]),
        string(r.label), r.lyapunov_constant)
end

println("\nDone.")
