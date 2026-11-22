# measure_sphere_areas.jl
# Computa areas esfericas (en estereorradianes) de las cuencas de atraccion
# para cada ejemplo del articulo draft_pcycle_basins.tex
# Uso: julia --project scripts/measure_sphere_areas.jl
#      (ejecutar desde 2026_RiemannBasins/)

using RiemannBasins
using Printf

# Total sphere = 4π sr
TOTAL_SR = 4π

function balance_metrics(areas::Dict{Int,Float64})
    vals = collect(values(areas))
    filter!(v -> v > 0, vals)
    fracs = vals ./ TOTAL_SR
    H = -sum(f * log(f) for f in fracs if f > 0)
    H_norm = length(fracs) > 1 ? H / log(length(fracs)) : 0.0
    minmax = minimum(fracs) / maximum(fracs)
    return H_norm, minmax
end

function print_result(name, areas)
    println("\n── $name ──")
    total = sum(values(areas))
    for (k, v) in sort(collect(areas))
        pct = 100v / TOTAL_SR
        @printf("  Cuenca %2d : %6.3f sr  (%5.1f%%)\n", k, v, pct)
    end
    @printf("  Total clasificado : %.3f sr  (%.1f%% de 4π)\n", total, 100total/TOTAL_SR)
    H, r = balance_metrics(areas)
    @printf("  Entropía norm.    : %.3f\n", H)
    @printf("  Ratio min/max     : %.3f\n", r)
    return H, r
end

results = Dict{String, Tuple{Float64,Float64}}()

println("=== Calculando cuencas esféricas ===")

# ── 1. Basílica  f(z)=z²−1,  2-ciclo {0,−1}
println("1/8 Basílica (z²−1, p=2)...")
num_bas = complex([-1.0, 0.0, 1.0]);  den_bas = complex([1.0])
m = basin_measure_pcycle(num_bas, den_bas, 2; iter=400, aproxprecision=8)
H, r = print_result("Basílica  z²−1  p=2", m)
results["Basílica"] = (H, r)

# ── 2. Avión  f(z)=z²+c_air, 3-ciclo
println("2/8 Avión (z²+c_air, p=3)...")
c_air = let x = -1.75
    for _ in 1:60; x -= (x^3 + 2x^2 + x + 1) / (3x^2 + 4x + 1) end; x
end
num_air = complex([c_air, 0.0, 1.0]);  den_air = complex([1.0])
m = basin_measure_pcycle(num_air, den_air, 3; iter=400, aproxprecision=8)
H, r = print_result("Avión  z²+c_air  p=3", m)
results["Avión"] = (H, r)

# ── 3. Conejo  f(z)=z²+c_rab, 3-ciclo
println("3/8 Conejo (z²+c_rab, p=3)...")
c_rabbit = let x = complex(-0.12, 0.75)
    for _ in 1:80; x -= (x^3 + 2x^2 + x + 1) / (3x^2 + 4x + 1) end; x
end
num_rab = complex([c_rabbit, 0.0, 1.0]);  den_rab = complex([1.0])
m = basin_measure_pcycle(num_rab, den_rab, 3; iter=400, aproxprecision=8)
H, r = print_result("Conejo  z²+c_rab  p=3", m)
results["Conejo"] = (H, r)

# ── 4. Período 5  f(z)=z²+c5
println("4/8 Período-5 (z²+c5, p=5)...")
c5 = let x = complex(-0.5024, 0.5629)
    for _ in 1:100
        z = 0+0im; dz = 0+0im
        for _ in 1:5; z, dz = z^2 + x, 2z*dz + 1 end
        x -= z/dz
    end; x
end
num_p5 = ComplexF64.([c5, 0.0, 1.0]);  den_p5 = ComplexF64.([1.0])
m = basin_measure_pcycle(num_p5, den_p5, 5; iter=500, aproxprecision=8)
H, r = print_result("Período-5  z²+c5  p=5", m)
results["Período-5"] = (H, r)

# ── 5. Newton  N(z)=(2z³+1)/(3z²),  puntos fijos
println("5/8 Newton z³−1 (p=1)...")
num_new = complex([1.0, 0.0, 0.0, 2.0]);  den_new = complex([0.0, 0.0, 3.0, 0.0])
m = basin_measure_pcycle(num_new, den_new, 1; iter=400, aproxprecision=8)
H, r = print_result("Newton z³−1  p=1", m)
results["Newton"] = (H, r)

# ── 6. Halley  H(z)=(z⁴+2z)/(2z³+1),  puntos fijos
println("6/8 Halley z³−1 (p=1)...")
num_hal = complex([0.0, 2.0, 0.0, 0.0, 1.0]);  den_hal = complex([1.0, 0.0, 0.0, 2.0])
m = basin_measure_pcycle(num_hal, den_hal, 1; iter=400, aproxprecision=8)
H, r = print_result("Halley z³−1  p=1", m)
results["Halley"] = (H, r)

# ── 7. Racional (z²−1)/(z²+1),  2-ciclo {0,−1}
println("7/8 Racional (z²−1)/(z²+1) (p=2)...")
num_rat = complex([-1.0, 0.0, 1.0]);  den_rat = complex([1.0, 0.0, 1.0])
m = basin_measure_pcycle(num_rat, den_rat, 2; iter=500, aproxprecision=8)
H, r = print_result("Racional (z²−1)/(z²+1)  p=2", m)
results["Racional"] = (H, r)

# ── 8. z²  (punto fijo superatractor en 0)
println("8/8 f(z)=z² (p=1)...")
num_z2 = complex([0.0, 0.0, 1.0]);  den_z2 = complex([1.0])
m = basin_measure_pcycle(num_z2, den_z2, 1; iter=300, aproxprecision=8)
H, r = print_result("z²  p=1", m)
results["z²"] = (H, r)

# ── Resumen ordenado por entropía
println("\n\n========== RANKING por entropía normalizada (mayor → más equilibrado) ==========")
@printf("%-22s  %s  %s\n", "Ejemplo", "Entropía", "min/max")
for (name, (H, r)) in sort(collect(results), by=x->-x[2][1])
    @printf("%-22s   %.3f      %.3f\n", name, H, r)
end
println("\nNota: Entropía=1 ⟺ todas las cuencas igual de grandes")
println("       min/max=1 ⟺ cuenca más pequeña = cuenca más grande")
println("       Total esfera = 4π ≈ $(round(TOTAL_SR, digits=4)) sr")
