---
title: 'RiemannBasins.jl: Visualisation of attraction basins of rational maps on the Riemann sphere'
tags:
  - Julia
  - complex dynamics
  - Riemann sphere
  - rational maps
  - basins of attraction
  - periodic cycles
  - colorblind-safe visualisation
authors:
  - name: Luis Javier Hernández Paricio
    orcid: 0000-0003-4528-7781
    affiliation: 1
affiliations:
  - name: Departamento de Matemáticas y Computación, Universidad de La Rioja, Logroño, Spain
    index: 1
date: 2 November 2026
bibliography: paper.bib
---

# Summary

The iteration of rational maps is one of the most visually rich subjects in
mathematics.  Starting from any point $z_0$ in the complex plane, repeated
application of a rational map $f$ produces an orbit that either converges to
an attracting cycle, diverges to infinity, or wanders indefinitely near a
fractal boundary known as the Julia set [@Milnor2006; @Beardon1991].  Colouring each starting point by the
attractor its orbit reaches produces the *basin portrait* of $f$.

A fundamental limitation of standard tools for this visualisation is that they
work exclusively in the affine complex plane $\mathbb{C}$, treating the point
at infinity as a special case or ignoring it entirely.  This is mathematically
incomplete: a rational map $f$ of degree $d$ extends to a holomorphic self-map
of the Riemann sphere $\mathbb{P}^1(\mathbb{C})$, where infinity is a perfectly
ordinary point.  By the Lefschetz fixed-point theorem, $f$ has exactly $d+1$
fixed points in $\mathbb{P}^1(\mathbb{C})$ — some of which necessarily lie at
infinity — and their basins are invisible when working only in $\mathbb{C}$.

`RiemannBasins.jl` is a Julia [@Julia2017] package that visualises attraction basins of
rational maps directly on $\mathbb{P}^1(\mathbb{C})$, using homogeneous
coordinates throughout.  It makes the full global dynamics visible — including
the basins at infinity — and produces figures that are perceptually honest and
accessible to people with colour-vision deficiency.

# Statement of Need

## Mathematical motivation

Working in $\mathbb{P}^1(\mathbb{C})$ rather than $\mathbb{C}$ has concrete
consequences.  For Newton's method applied to $z^3-1$, for example, the
degree-3 map has four fixed points in $\mathbb{P}^1(\mathbb{C})$: the three
cube roots of unity (superattracting) and $\infty$ (repelling).  The basin of
the repeller at $\infty$ is the boundary of all three attracting basins and is
a key feature of the global dynamics.  Any software that works only in
$\mathbb{C}$ either misses this point entirely or treats it as a numerical
infinity that causes exceptions.

`RiemannBasins.jl` avoids this by representing every point as a homogeneous
pair $[z_0:z_1]\in\mathbb{P}^1(\mathbb{C})$ and extending arithmetic to an
*augmented sphere*
$\mathbb{P}^{1+}(\mathbb{C})=\mathbb{P}^1(\mathbb{C})\cup\{[0:0]\}$,
where $[0:0]$ is a computational *super-zero* that absorbs numerical
indeterminations (such as $0/0$ or $\infty/\infty$) without raising exceptions.
Iteration loops proceed through all points of the sphere, including $\infty$,
without special-casing.

## Gap in existing software

General-purpose dynamical-systems packages such as `DynamicalSystems.jl`
[@Datseris2023] do not specialise in rational maps on
$\mathbb{P}^1(\mathbb{C})$.  Ad-hoc scripts for Newton fractals in Python or
MATLAB typically work in $\mathbb{C}$, use arbitrary or aesthetically motivated
colour assignments, and offer no accessibility features.  `RiemannBasins.jl`
fills this gap by combining rigorous projective arithmetic, perceptual
colourimetry, and inclusive visualisation in a single documented and tested
package.

# Mathematics

## Homogeneous coordinates and iteration

A rational map $f(z)=p(z)/q(z)$ with $d=\max(\deg p,\deg q)$ is lifted to
$\mathbb{P}^1(\mathbb{C})$ by homogenising $p$ and $q$ to degree $d$:

$$F([z_0:z_1]) = \bigl[p_0(z_0,z_1):q_0(z_0,z_1)\bigr].$$

Coefficients are stored in ascending order $[a_0,a_1,\dots,a_d]$, consistent
with the `Polynomials.jl` convention, so that $[a_0,\dots,a_d]$ represents
$a_0+a_1z+\cdots+a_dz^d$.  At each pixel the chordal (spherical) metric is used
to decide convergence, making the tolerance independent of the chart.

## Period-$p$ cycle basins

For attracting cycles of period $p\geq 2$, the function $f^p$ is required.
`RiemannBasins.jl` computes it by symbolic polynomial composition: if
$f(z)=p(z)/q(z)$ and $g(z)=r(z)/s(z)$, the numerator of $f\circ g$ is

$$\operatorname{num}(f\circ g)=\sum_{k=0}^{d_f}a_k\,r^k\,s^{d_f-k},$$

where $a_0,\dots,a_{d_f}$ are the numerator coefficients of $f$.  Fixed points
of $f^p$ are then grouped into orbits under $f$ via the successor map
$\zeta_i\mapsto f(\zeta_i)$, so that all elements of the same $p$-cycle receive
a single colour in the basin portrait.  This *Option-A colouring* gives an
immediately readable picture even for high-period cycles.

## Perceptually uniform colorisation

Colours are constructed in the LCHab perceptual colour space [@Fairchild2013],
distributing hue angles uniformly so that every basin has identical perceived
lightness and chroma.  Convergence speed is encoded by modulating the lightness
channel (banding), revealing the fractal structure of the Julia set without
distorting the colour assignment.  An alternative palette implementing the
Okabe–Ito recommendations [@OkabeIto2008] is available via the
`colorblind=true` keyword; it guarantees legibility for observers with
deuteranopia or protanopia.

# Functionality

The main public functions are listed in Table 1.

| Function | Description |
|---|---|
| `basin_portrait` | Basin portrait in the standard affine chart |
| `two_chart_portrait` | Simultaneous portrait in both stereographic charts |
| `basin_portrait_pcycle` | Basins of an attracting cycle of period $p$ |
| `basin_relief_3d_makie` | 3-D relief surface via GLMakie (colour $\perp$ elevation) |
| `basin_relief_3d` | 3-D relief surface via Plots.jl/GR |
| `animate_basins` | Animated GIF over a one-parameter family |
| `fixed_points` | Exact fixed points of $f$ in $\mathbb{P}^1(\mathbb{C})$ |
| `chordal_metric` | Chordal (spherical) distance between two points |
| `okabe_ito_colors` | Okabe–Ito palette for up to eight categories |

Table: Main exported functions of `RiemannBasins.jl`.

A minimal working example:

```julia
using RiemannBasins

# Newton's method for z³ − 1:  f(z) = (2z³+1)/(3z²)
# Coefficients in ascending order [a₀, a₁, …, aₙ]
num = complex([1.0, 0.0, 0.0, 2.0])   # 1 + 2z³
den = complex([0.0, 0.0, 3.0, 0.0])   # 3z²

basin_portrait(num, den; bands = 8.0)             # affine chart
two_chart_portrait(num, den; bands = 8.0)         # both charts
basin_portrait_pcycle(num, den, 2; bands = 8.0)   # period-2 basins
basin_relief_3d_makie(num, den; relief_scale = 0.01) # 3-D relief
```

The `bands` parameter controls the number of visible banding rings around each
basin; higher values reveal finer structure near the Julia set.  The
`relief_scale` parameter in the 3-D functions controls the elevation amplitude:
low values foreground colour information, high values foreground the fractal
topography.  All functions accept `colorblind=true` to switch to the
Okabe–Ito palette.

# Development history

The package was first prototyped in 2019 as an exploratory tool accompanying
research on the topology of iteration spaces of rational maps
[@HernandezParicio2019].
The current version [@RiemannBasins2026] is a complete rewrite that adds homogeneous-coordinate
arithmetic, the LCHab and Okabe–Ito colour engines, 3-D relief surfaces via
GLMakie, and symbolic computation of period-$p$ iterates.

# Acknowledgements

The author thanks the Universidad de La Rioja for institutional support and the
Julia community for the libraries on which this package builds: `Polynomials.jl`,
`Colors.jl`, `Plots.jl`, and `GLMakie.jl`.

# References
