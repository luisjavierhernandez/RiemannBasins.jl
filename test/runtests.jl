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

end
