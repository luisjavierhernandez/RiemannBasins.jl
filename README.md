# RiemannBasins.jl

**Visualisation of attraction basins of rational maps on the Riemann sphere**

`RiemannBasins.jl` is a Julia package for computing and visualising the basins
of attraction of rational maps on the projective line **P¹(ℂ)** — the Riemann
sphere. Unlike tools that work only in the affine complex plane ℂ, this package
uses **homogeneous coordinates** throughout, so the point at infinity and its
basin are always visible.

Key features:

- Full projective arithmetic on P¹(ℂ); an augmented super-zero `[0:0]`
  absorbs numerical indeterminations without exceptions.
- **Period-p cycle basins** via symbolic polynomial composition of f^p.
- **Perceptually uniform colorisation** in the LCHab colour space: every basin
  has identical perceived lightness, and convergence speed is encoded by
  lightness modulation (banding).
- **Colorblind-safe palette** (Okabe–Ito 2008) available with `colorblind=true`.
- 2-D basin portraits in one or both stereographic charts.
- 3-D relief surfaces (Plots.jl/GR and GLMakie).
- Animated GIFs over one-parameter families.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/luisjavierhernandez/RiemannBasins.jl")
```

---

## Quick start

```julia
using RiemannBasins

# Newton's method for z³ − 1:  f(z) = (2z³ + 1) / (3z²)
# Coefficients in ascending order [a₀, a₁, …, aₙ]
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z³
den = complex([0.0, 0.0, 3.0, 0.0])   # 3z²

basin_portrait(num, den; bands = 8.0)              # standard affine chart
two_chart_portrait(num, den; bands = 8.0)          # both stereographic charts
basin_portrait_pcycle(num, den, 2; bands = 8.0)    # basins of period-2 cycles
basin_relief_3d_makie(num, den; relief_scale = 0.01)  # 3-D relief (GLMakie)
```

Add `colorblind = true` to any function for the Okabe–Ito accessible palette.

---

## Main functions

| Function | Description |
|---|---|
| `basin_portrait` | Basin portrait in the standard affine chart |
| `two_chart_portrait` | Simultaneous portrait in both stereographic charts |
| `basin_portrait_pcycle` | Basins of an attracting cycle of period p |
| `basin_relief_3d_makie` | 3-D relief surface via GLMakie |
| `basin_relief_3d` | 3-D relief surface via Plots.jl/GR |
| `animate_basins` | Animated GIF over a one-parameter family |
| `fixed_points` | Fixed points of f in P¹(ℂ) |
| `chordal_metric` | Chordal (spherical) distance between two points |
| `okabe_ito_colors` | Okabe–Ito palette for up to eight categories |

---

## Examples

The `Examples/` folder contains ready-to-run scripts:

| Script | Content |
|---|---|
| `3rootsNewton.jl` | Newton's method for z³ − 1 |
| `newton_two_charts.jl` | Two-chart portrait |
| `halley_z3minus1.jl` | Halley's method for z³ − 1 |
| `halley_two_charts.jl` | Halley two-chart portrait |
| `relief_3d.jl` | 3-D relief surface |
| `newton_animation.jl` | Animated basin over a parameter family |
| `lch_and_makie.jl` | LCHab colorisation and GLMakie demo |

---

## Documentation

A detailed article (in Spanish) describing the mathematics and the package is
available as `docs/borrador_gaceta.tex` (submitted to *La Gaceta de la Real
Sociedad Matemática Española*).

A software paper describing the package is available as `paper.md`
(prepared for submission to the *Journal of Open Source Software*).

---

## Citation

If you use `RiemannBasins.jl` in your research, please cite:

```bibtex
@software{RiemannBasins2026,
  author = {Hern{\'a}ndez Paricio, Luis Javier},
  title  = {{RiemannBasins.jl}: Visualisation of Attraction Basins of
             Rational Maps on the {Riemann} Sphere},
  year   = {2026},
  url    = {https://github.com/luisjavierhernandez/RiemannBasins.jl}
}
```

---

## License

MIT © 2026 Luis Javier Hernández Paricio — see [LICENSE](LICENSE).
