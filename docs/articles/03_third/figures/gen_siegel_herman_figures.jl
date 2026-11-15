# gen_siegel_herman_figures.jl
#
# Genera las figuras del artículo "Siegel disks, Herman rings and the Hopf fibration"
# (docs/articles/03_third/draft_siegel_herman.tex).
#
# Uso:
#   cd docs/articles/03_third/figures
#   julia --project=../../../../ gen_siegel_herman_figures.jl
#
# Las figuras se guardan como PNG en el directorio actual.
# Requiere: RiemannBasins, Plots (incluidos en Project.toml).

using RiemannBasins
using Plots
using LinearAlgebra

# ── Parámetros globales de salida ──────────────────────────────────────────
const DPI   = 300
const NPIX  = 600          # resolución cuadrada de cada imagen
const OUTDIR = @__DIR__    # directorio donde está este script

# ── Helpers de coloreado ──────────────────────────────────────────────────

"""Convierte un ángulo en [0,2π) a un color HSV."""
function angle_to_color(θ)
    h = mod(θ / (2π), 1.0)
    return HSV(360h, 0.85, 0.95)
end

"""Paleta para clasificación de Fatou: color por etiqueta."""
function label_color(label)
    label == :superattracting && return colorant"#1a1a2e"   # azul muy oscuro
    label == :attracting       && return colorant"#16213e"   # azul oscuro
    label == :parabolic        && return colorant"#f7a844"   # naranja
    label == :siegel           && return colorant"#2ecc71"   # verde
    label == :herman           && return colorant"#e74c3c"   # rojo
    return colorant"#ecf0f1"                                  # blanco: indeterminate/Julia
end

# ══════════════════════════════════════════════════════════════════════════
# Figura 1: Retrato de cuencas de z²+c con disco de Siegel de razón dorada
# ══════════════════════════════════════════════════════════════════════════
#
# El disco de Siegel de razón dorada es el componente de Fatou de
#     f_c(z) = z² + c₀,   c₀ = λ(2−λ)/4,   λ = exp(2πiρ),  ρ = (√5−1)/2
# centrado en el punto fijo indeferente p = λ/2.
# (ρ satisface la condición de Brjuno porque es de tipo acotado.)

function figure_siegel_disk()
    println("Figura 1: disco de Siegel (razón dorada)…")

    ρ  = (sqrt(5) - 1) / 2        # razón dorada ≈ 0.6180
    λ  = exp(2π * im * ρ)         # multiplicador del punto fijo
    c₀ = λ * (2 - λ) / 4          # parámetro de la familia cuadrática
    p₀ = λ / 2                    # punto fijo indeferente

    println("  c₀ ≈ ", round(c₀, digits=5))
    println("  p₀ ≈ ", round(p₀, digits=5))
    println("  |λ| = ", abs(λ))

    # Coeficientes homogéneos: f(z) = z² + c₀  → num = [c₀, 0, 1], den = [1]
    cnum = [c₀, 0.0 + 0.0im, 1.0 + 0.0im]
    cden = [1.0 + 0.0im]

    # Verificar clasificación del punto crítico z=0
    res = critical_orbit_analysis(cnum, cden; max_iter=800, N_birkhoff=1200)
    for r in res
        z_crit = r.critical_point[2] ≈ 0 ? Inf : r.critical_point[1] / r.critical_point[2]
        println("  Crítico ≈ ", round(z_crit, digits=4),
                "  →  ", r.label,
                "  ρ ≈ ", round(r.rotation_number, digits=6))
    end

    # Cuadrícula de píxeles en el plano z
    xmin, xmax = -1.5, 1.5
    ymin, ymax = -1.5, 1.5
    xs = range(xmin, xmax; length=NPIX)
    ys = range(ymin, ymax; length=NPIX)

    img = Matrix{Float64}(undef, NPIX, NPIX)   # tiempo de escape (NaN = no escapa)

    max_it = 300
    R2     = 4.0^2   # radio de escape al cuadrado

    for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        z = complex(x, y)
        escaped = false
        for n in 1:max_it
            z = z^2 + c₀
            if abs2(z) > R2
                img[i, j] = n / max_it
                escaped = true
                break
            end
        end
        escaped || (img[i, j] = NaN)
    end

    # Colorear: puntos escapados por tiempo suavizado, Fatou interior en verde
    p = heatmap(xs, ys, img';
        color=:inferno, aspect_ratio=1,
        colorbar=false, axis=nothing, border=:none,
        title="Disco de Siegel — ρ = (√5−1)/2",
        titlefontsize=10, dpi=DPI, size=(NPIX, NPIX))

    # Marcar el punto fijo y el centro del disco
    scatter!(p, [real(p₀)], [imag(p₀)];
        marker=:circle, markersize=5, color=:lime, label="p₀")
    scatter!(p, [0.0], [0.0];
        marker=:star5, markersize=5, color=:yellow, label="crit")

    outfile = joinpath(OUTDIR, "fig1_siegel_disk.png")
    savefig(p, outfile)
    println("  Guardada: ", outfile)
    return p
end

# ══════════════════════════════════════════════════════════════════════════
# Figura 2: Coloreado por número de rotación estimado en la familia z²+c
# ══════════════════════════════════════════════════════════════════════════
#
# Recorre los parámetros c en la cardioide principal del conjunto de Mandelbrot
# y colorea por la razón de rotación estimada del punto crítico.
# La cardioide principal: c = e^{iθ}/2 − e^{2iθ}/4,  θ ∈ [0, 2π).

function figure_rotation_field()
    println("Figura 2: campo de número de rotación…")

    N_param = 300          # número de parámetros en la cardioide
    thetas  = range(0, 2π; length=N_param + 1)[1:end-1]

    rhos    = Float64[]
    cs_re   = Float64[]
    cs_im   = Float64[]

    for θ in thetas
        λ = exp(im * θ)
        c = λ / 2 - λ^2 / 4    # punto en la cardioide para multiplicador λ
        push!(cs_re, real(c))
        push!(cs_im, imag(c))

        cnum = [c, 0.0 + 0.0im, 1.0 + 0.0im]
        cden = [1.0 + 0.0im]

        if abs(abs(λ) - 1.0) < 1e-5   # solo borde: |λ|=1
            res = critical_orbit_analysis(cnum, cden;
                      max_iter=600, N_birkhoff=800, detect_parabolic=true)
            ρ_est = isempty(res) ? NaN : res[1].rotation_number
        else
            ρ_est = mod(angle(λ) / (2π), 1.0)   # rotación exacta en cardioide
        end

        push!(rhos, ρ_est)
    end

    # Scatter coloreado por ρ (hue), tamaño constante
    hue_vals = mod.(rhos, 1.0)

    p = scatter(cs_re, cs_im;
        marker_z=hue_vals, color=:hsv,
        markersize=4, markerstrokewidth=0,
        aspect_ratio=1, legend=false,
        title="Número de rotación a lo largo de la cardioide",
        xlabel="Re(c)", ylabel="Im(c)",
        titlefontsize=9, dpi=DPI, size=(NPIX, NPIX))

    outfile = joinpath(OUTDIR, "fig2_rotation_field.png")
    savefig(p, outfile)
    println("  Guardada: ", outfile)
    return p
end

# ══════════════════════════════════════════════════════════════════════════
# Figura 3: Anillo de Herman — familia de Blaschke B_a(z)
# ══════════════════════════════════════════════════════════════════════════
#
# B_a(z) = exp(2πiρ) · z · (z−a)/(1−a·z),   a ∈ ℝ, 0 < a < 1
#
# La circunferencia unitaria |z|=1 es invariante (es un producto de Blaschke).
# Para a suficientemente grande (a→1⁻) la restricción f|_{S¹} es una difeo-
# morfismo del círculo con número de rotación ρ.  En la esfera de Riemann,
# el anillo de Herman es el componente de Fatou que contiene S¹.
#
# Como función racional de grado 2:
#   numerador  = exp(2πiρ) · [z² − az]  =  exp(2πiρ)·(-a·z + z²)
#   denominador= [1 − az]
# En coeficientes homogéneos (constante-primero):
#   cnum = [0, -a·exp(2πiρ), exp(2πiρ)]
#   cden = [1, -a]

function figure_herman_ring()
    println("Figura 3: anillo de Herman (Blaschke)…")

    ρ   = (sqrt(5) - 1) / 2        # razón dorada
    λ   = exp(2π * im * ρ)
    a   = 0.7                      # parámetro de Blaschke ∈ (0,1)

    # B_a(z) = λ · z(z−a)/(1−az)
    # num(z) = λ(z²−az) = λ·(-a)·z + λ·z²  → coefs [0, -aλ, λ]
    # den(z) = 1 − az               → coefs [1, -a]
    cnum = [complex(0.0), -a * λ, λ]
    cden = [complex(1.0), complex(-a)]

    println("  B_a con a=", a, "  ρ=", round(ρ, digits=6))

    # Verificar clasificación
    res = critical_orbit_analysis(cnum, cden; max_iter=600, N_birkhoff=1000)
    for r in res
        t = r.critical_point[2]
        z_crit = abs(t) < 1e-10 ? Inf : r.critical_point[1] / t
        println("  Crítico ≈ ", round(z_crit, digits=4),
                "  →  ", r.label,
                "  ρ ≈ ", round(r.rotation_number, digits=6))
    end

    # Cuadrícula de píxeles — vista en el plano afín con |z| ≤ 2.5
    xmin, xmax = -2.5, 2.5
    ymin, ymax = -2.5, 2.5
    xs = range(xmin, xmax; length=NPIX)
    ys = range(ymin, ymax; length=NPIX)

    # Colorear por destino del punto: escapa a ∞ (oscuro) / se queda (claro)
    img    = fill(0.0, NPIX, NPIX)
    max_it = 200
    hpair  = RiemannBasins._pair_of_functions(ComplexF64.(cnum), ComplexF64.(cden))

    for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        z   = complex(x, y)
        hpt = abs(z) >= 1 ? (z, complex(1.0)) : (complex(1.0), 1 / z)
        hpt = (Complex{Float64}(hpt[1]), Complex{Float64}(hpt[2]))
        norm_max = max(abs(hpt[1]), abs(hpt[2]))
        norm_max < 1e-15 && continue
        hpt = (hpt[1] / norm_max, hpt[2] / norm_max)

        escaped  = false
        converge = false
        for n in 1:max_it
            hpt = RiemannBasins._rational_function(hpair, hpt)
            rr  = abs(hpt[2]) < 1e-10 ? Inf : abs(hpt[1] / hpt[2])
            if rr > 1e3        # escapa a ∞
                img[i, j] = n / max_it
                escaped = true
                break
            end
            if rr < 1e-6       # colapsa a 0
                img[i, j] = -n / max_it
                converge = true
                break
            end
        end
        (!escaped && !converge) && (img[i, j] = NaN)  # anillo de Herman ≈ NaN
    end

    # Mapa de colores: escapa→azul oscuro, converge→rojo, anillo→blanco
    cmap = cgrad([:navy, :white, :firebrick])
    p = heatmap(xs, ys, img';
        color=cmap, aspect_ratio=1,
        colorbar=false, axis=nothing, border=:none,
        title="Anillo de Herman — B_{0.7}(z),  ρ = (√5−1)/2",
        titlefontsize=9, dpi=DPI, size=(NPIX, NPIX))

    # Dibujar la circunferencia unitaria
    θs = range(0, 2π; length=300)
    plot!(p, cos.(θs), sin.(θs); color=:lime, linewidth=1.5, label="|z|=1")

    outfile = joinpath(OUTDIR, "fig3_herman_ring.png")
    savefig(p, outfile)
    println("  Guardada: ", outfile)
    return p
end

# ══════════════════════════════════════════════════════════════════════════
# Figura 4: Órbita cuasi-periódica en el disco de Siegel — toro sólido en S³
# ══════════════════════════════════════════════════════════════════════════
#
# Visualiza la órbita del punto crítico z=0 bajo f_c₀ (razón dorada)
# coloreada por la fase cumulativa de rotación.

function figure_quasiperiodic_orbit()
    println("Figura 4: órbita cuasi-periódica en el disco de Siegel…")

    ρ  = (sqrt(5) - 1) / 2
    λ  = exp(2π * im * ρ)
    c₀ = λ * (2 - λ) / 4

    cnum = [c₀, 0.0 + 0.0im, 1.0 + 0.0im]
    cden = [1.0 + 0.0im]
    hpair = RiemannBasins._pair_of_functions(ComplexF64.(cnum), ComplexF64.(cden))

    # Órbita de z=0 (punto crítico)
    N_orb = 2000
    orbit = Vector{Complex{Float64}}(undef, N_orb)
    hpt   = (complex(0.0), complex(1.0))   # z=0 en homogéneo
    for k in 1:N_orb
        hpt      = RiemannBasins._rational_function(hpair, hpt)
        orbit[k] = hpt[2] ≈ 0 ? Inf + 0im : hpt[1] / hpt[2]
    end

    # Fase acumulada (ángulo respecto al centro del disco)
    p₀      = λ / 2
    angles  = mod.(angle.(orbit .- p₀), 2π)
    colors  = angles ./ (2π)

    p = scatter(real.(orbit), imag.(orbit);
        marker_z=colors, color=:hsv, markersize=1.5,
        markerstrokewidth=0, legend=false, aspect_ratio=1,
        title="Órbita cuasi-periódica en D²×S¹ (disco de Siegel, N=$(N_orb))",
        xlabel="Re(z)", ylabel="Im(z)",
        titlefontsize=9, dpi=DPI, size=(NPIX, NPIX))

    # Marcar el punto fijo
    scatter!(p, [real(p₀)], [imag(p₀)]; marker=:star5,
        markersize=8, color=:white, label="p₀")

    outfile = joinpath(OUTDIR, "fig4_quasiperiodic_orbit.png")
    savefig(p, outfile)
    println("  Guardada: ", outfile)
    return p
end

# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
println("=" ^ 60)
println("Generando figuras para draft_siegel_herman.tex")
println("=" ^ 60)

figure_siegel_disk()
figure_rotation_field()
figure_herman_ring()
figure_quasiperiodic_orbit()

println()
println("Listo. Figuras guardadas en: ", OUTDIR)
