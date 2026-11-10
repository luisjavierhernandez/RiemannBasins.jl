using Test
using RiemannBasins
using Colors

@testset "RiemannBasins unit tests" begin

    # ── sphere_bijection ─────────────────────────────────────────────────────
    @testset "sphere_bijection" begin
        # [0:1] → south pole (0, 0, -1)
        @test sphere_bijection((0.0+0im, 1.0+0im)) ≈ [0.0, 0.0, -1.0]

        # [1:0] → north pole (0, 0, 1)
        @test sphere_bijection((1.0+0im, 0.0+0im)) ≈ [0.0, 0.0, 1.0]

        # super-zero [0:0] → origin sentinel
        @test sphere_bijection((0.0+0im, 0.0+0im)) == [0.0, 0.0, 0.0]

        # [1:1] → equatorial point (1, 0, 0)
        @test sphere_bijection((1.0+0im, 1.0+0im)) ≈ [1.0, 0.0, 0.0] atol=1e-12

        # [im:1] → equatorial point (0, 1, 0)
        @test sphere_bijection((1.0im, 1.0+0im)) ≈ [0.0, 1.0, 0.0] atol=1e-12

        # image always lies on S²: |φ(p)|² = 1 for any non-super-zero p
        for p in [(1.0+2.0im, 1.0+0im), (3.0+0im, 1.0+im), (1.0+0im, 2.0+3.0im)]
            v = sphere_bijection(p)
            @test sum(v.^2) ≈ 1.0 atol=1e-12
        end
    end

    # ── chordal_metric ───────────────────────────────────────────────────────
    @testset "chordal_metric" begin
        north = (1.0+0im, 0.0+0im)
        south = (0.0+0im, 1.0+0im)
        p     = (1.0+0im, 1.0+0im)
        q     = (2.0+0im, 1.0+0im)

        # distance from a point to itself is 0
        @test chordal_metric(p, p) ≈ 0.0 atol=1e-14

        # symmetry
        @test chordal_metric(p, q) ≈ chordal_metric(q, p) atol=1e-14

        # antipodal points (north–south poles) are distance 2 apart
        @test chordal_metric(north, south) ≈ 2.0 atol=1e-12

        # non-negativity
        @test chordal_metric(p, q) ≥ 0.0

        # triangle inequality: d(p, r) ≤ d(p, q) + d(q, r)
        r = (1.0+1.0im, 1.0+0im)
        @test chordal_metric(p, r) ≤ chordal_metric(p, q) + chordal_metric(q, r) + 1e-12
    end

    # ── fixed_points ─────────────────────────────────────────────────────────
    @testset "fixed_points – Newton z³−1" begin
        # f(z) = (2z³+1)/(3z²)
        num = complex([1.0, 0.0, 0.0, 2.0])
        den = complex([0.0, 0.0, 3.0, 0.0])
        superfixed, infixed, finitelist = fixed_points(num, den)

        # three finite fixed points
        @test length(finitelist) == 3

        # ∞ is a fixed point of this map
        @test infixed == true

        # (superfixed depends on the code's homogeneous evaluation at [0:0])
        # the finite fixed points are the cube roots of unity: 1, e^{2πi/3}, e^{4πi/3}
        roots_of_unity = [
            complex(1.0),
            exp(2π*im/3),
            exp(4π*im/3)
        ]
        fp_values = [fp[1] / fp[2] for fp in finitelist]   # affine coords
        for ω in roots_of_unity
            @test any(abs(z - ω) < 1e-6 for z in fp_values)
        end
    end

    @testset "fixed_points – Newton z²−1" begin
        # f(z) = (z²+1)/(2z)  →  fixed points: z=1 and z=−1
        num = complex([1.0, 0.0, 1.0])   # 1 + z²
        den = complex([0.0, 2.0, 0.0])   # 2z
        superfixed, infixed, finitelist = fixed_points(num, den)

        # two finite fixed points
        @test length(finitelist) == 2

        # ∞ is a fixed point (degree num = degree den = 2)
        @test infixed == true

        # the fixed points are ±1
        fp_values = [fp[1] / fp[2] for fp in finitelist]
        @test any(abs(z - 1.0) < 1e-6 for z in fp_values)
        @test any(abs(z + 1.0) < 1e-6 for z in fp_values)
    end

    # ── okabe_ito_colors ─────────────────────────────────────────────────────
    @testset "okabe_ito_colors" begin
        oi = okabe_ito_colors()

        # returns exactly 8 colours
        @test length(oi) == 8

        # all entries are RGB{Float64}
        @test all(c isa RGB{Float64} for c in oi)

        # all channel values in [0, 1]
        for c in oi
            @test 0.0 ≤ red(c)   ≤ 1.0
            @test 0.0 ≤ green(c) ≤ 1.0
            @test 0.0 ≤ blue(c)  ≤ 1.0
        end

        # all colours are distinct
        @test length(unique(oi)) == 8
    end

    # ── basin_portrait ───────────────────────────────────────────────────────
    @testset "basin_portrait – smoke tests" begin
        num = complex([1.0, 0.0, 0.0, 2.0])
        den = complex([0.0, 0.0, 3.0, 0.0])
        kw  = (expprecision=4, iter=20, show_plot=false)

        @test basin_portrait(num, den; kw..., colorstrategy="positionplusiteration") !== nothing
        @test basin_portrait(num, den; kw..., colorstrategy="positionfixedpoints")   !== nothing
        @test basin_portrait(num, den; kw..., colorstrategy="iteration")             !== nothing
        @test basin_portrait(num, den; kw..., colorblind=true)                       !== nothing
    end

    # ── two_chart_portrait ───────────────────────────────────────────────────
    @testset "two_chart_portrait – smoke test" begin
        num = complex([1.0, 0.0, 0.0, 2.0])
        den = complex([0.0, 0.0, 3.0, 0.0])
        plt = two_chart_portrait(num, den;
            expprecision = 4,
            iter         = 20,
            show_plot    = false)
        @test plt !== nothing
    end

    # ── basin_relief_3d ──────────────────────────────────────────────────────
    @testset "basin_relief_3d – smoke tests" begin
        num = complex([1.0, 0.0, 0.0, 2.0])
        den = complex([0.0, 0.0, 3.0, 0.0])
        kw  = (expprecision=4, iter=20, show_plot=false)

        @test basin_relief_3d(num, den; kw...)                    !== nothing
        @test basin_relief_3d(num, den; kw..., colorblind=true)   !== nothing
        @test basin_relief_3d(num, den; kw..., bands=12.0)        !== nothing
        @test basin_relief_3d(num, den; kw..., inverted=true)     !== nothing

        # inverted surface must have its maximum at the Julia-set boundary:
        # non-converging points (pos ≤ 0) are assigned z = relief_scale,
        # while basin-centre points have z close to 0.
        # We cannot easily inspect the z matrix here, but a second call
        # with inverted=true should simply not error and return a figure.
        @test basin_relief_3d(num, den; kw..., inverted=true, bands=6.0) !== nothing
    end

    # ── basin_portrait_pcycle ────────────────────────────────────────────────
    @testset "basin_portrait_pcycle – smoke tests" begin
        # Newton z³−1: f(z) = (2z³+1)/(3z²) — three fixed points, no 2-cycles
        num = complex([1.0, 0.0, 0.0, 2.0])
        den = complex([0.0, 0.0, 3.0, 0.0])
        kw  = (expprecision=4, iter=20, show_plot=false)

        # p=1 should reproduce a portrait with the same fixed points as basin_portrait
        @test basin_portrait_pcycle(num, den, 1; kw...) !== nothing

        # p=2 finds fixed points of f², groups them, should not error
        @test basin_portrait_pcycle(num, den, 2; kw...) !== nothing

        # colorstrategy and colorblind flags
        @test basin_portrait_pcycle(num, den, 1; kw...,
            colorstrategy="positionfixedpoints") !== nothing
        @test basin_portrait_pcycle(num, den, 2; kw...,
            colorblind=true) !== nothing
    end

    # ── §14 Critical orbit analysis (Fatou's theorem) ───────────────────────
    @testset "§14 critical_points – Möbius degree 1" begin
        # f(z) = (z+1)/z  → W(z) = p′q − pq′ = 1·z − (z+1)·1 = −1 (constant)
        # 2d−2 = 0 finite critical points expected.
        num_m = complex([1.0, 1.0])
        den_m = complex([0.0, 1.0])
        cp = critical_points(num_m, den_m)
        @test isempty(cp)
    end

    @testset "§14 critical_points – f(z)=z²" begin
        # W(z) = 2z → 1 finite critical point at z=0
        num_z2 = complex([0.0, 0.0, 1.0])
        den_z2 = complex([1.0])
        cp = critical_points(num_z2, den_z2)
        @test length(cp) == 1
        @test abs(cp[1][1]) < 1e-5    # z-component ≈ 0
    end

    @testset "§14 critical_points – basílica z²−1" begin
        # W(z) = 2z → 1 critical point at z=0
        num_b = complex([-1.0, 0.0, 1.0])
        den_b = complex([1.0])
        cp = critical_points(num_b, den_b)
        @test length(cp) == 1
        @test abs(cp[1][1]) < 1e-5    # z-component ≈ 0
    end

    @testset "§14 critical_points – Newton z³−1" begin
        # f(z) = (2z³+1)/(3z²); W = 6z(z³−1) → roots z=0, 1, e^{2πi/3}, e^{4πi/3}
        num_n = complex([1.0, 0.0, 0.0, 2.0])
        den_n = complex([0.0, 0.0, 3.0, 0.0])
        cp = critical_points(num_n, den_n)
        @test length(cp) == 4
        # Recover affine values (second coord is t, normalised to 1 in first chart)
        aff = [p[1] / p[2] for p in cp]
        @test any(abs(z)         < 1e-5 for z in aff)      # z = 0
        @test any(abs(z - 1)     < 1e-5 for z in aff)      # z = 1
        ω = exp(2π * im / 3)
        @test any(abs(z - ω)     < 1e-5 for z in aff)      # z = e^{2πi/3}
        @test any(abs(z - conj(ω)) < 1e-5 for z in aff)    # z = e^{4πi/3}
    end

    @testset "§14 critical_orbit_analysis – f(z)=z² superattracting ∞" begin
        # z=0 is a superattracting fixed point: f(0)=0 and f′(0)=0
        num_z2 = complex([0.0, 0.0, 1.0])
        den_z2 = complex([1.0])
        res = critical_orbit_analysis(num_z2, den_z2)
        @test length(res) == 1
        r = res[1]
        @test r.label           == :superattracting
        @test r.period          == 1
        @test r.lyapunov_constant < 1e-8
        @test length(r.cycle_points) == 1
        @test abs(r.cycle_points[1][1]) < 1e-5   # cycle point ≈ z=0
    end

    @testset "§14 critical_orbit_analysis – basílica z²−1 superattracting 2-cycle" begin
        # The unique critical point z=0 lies IN the 2-cycle {0,−1}.
        # Because f′(0)=0, the Lyapunov constant c=0 → :superattracting.
        num_b = complex([-1.0, 0.0, 1.0])
        den_b = complex([1.0])
        res = critical_orbit_analysis(num_b, den_b)
        @test length(res) == 1
        r = res[1]
        @test r.label            == :superattracting
        @test r.period           == 2
        @test r.lyapunov_constant < 1e-8
        @test length(r.cycle_points) == 2
    end

    @testset "§14 critical_orbit_analysis – Newton z³−1" begin
        # Three critical points (cube roots of unity) land on superattracting
        # fixed points of Newton's method (the roots of z³−1).
        # The fourth critical point z=0 maps to ∞ in one step; ∞ is a repelling
        # fixed point (multiplier 3/2), so that orbit is :indeterminate.
        num_n = complex([1.0, 0.0, 0.0, 2.0])
        den_n = complex([0.0, 0.0, 3.0, 0.0])
        res = critical_orbit_analysis(num_n, den_n)
        @test length(res) == 4
        @test count(r -> r.label == :superattracting && r.period == 1, res) == 3
        @test count(r -> r.label == :indeterminate,                     res) == 1
    end

    @testset "§14 max_attracting_period" begin
        # Möbius (degree 1): no critical points → no bound → 0
        @test max_attracting_period(complex([1.0, 1.0]), complex([0.0, 1.0])) == 0

        # f(z)=z²: superattracting fixed point at 0 → P = 1
        @test max_attracting_period(complex([0.0, 0.0, 1.0]), complex([1.0])) == 1

        # Basílica z²−1: superattracting 2-cycle {0,−1} → P = 2
        @test max_attracting_period(complex([-1.0, 0.0, 1.0]), complex([1.0])) == 2

        # Newton z³−1: attracting fixed points (period 1) → P = 1
        @test max_attracting_period(
            complex([1.0, 0.0, 0.0, 2.0]),
            complex([0.0, 0.0, 3.0, 0.0])) == 1
    end

    @testset "§14 print_critical_orbit_summary – smoke test" begin
        # Verifies the function runs, prints, and returns the same results
        # as critical_orbit_analysis without crashing.
        num_b = complex([-1.0, 0.0, 1.0])
        den_b = complex([1.0])
        res = print_critical_orbit_summary(num_b, den_b)
        @test length(res) == 1
        @test res[1].period == 2
        @test res[1].label  == :superattracting
    end

end
