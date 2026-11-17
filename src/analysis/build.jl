# Building and editing an `Orbit` from Keplerian elements, and two-body propagation: the
# port of `try_keplerian`, `try_keplerian_apsis_radii`, `try_keplerian_mean_anomaly`, the
# `with_*` and `set_*` calls and `at_epoch` of src/astro/orbit.rs of ANISE (b3aaadc).
# Angles are in rad here (ANISE: degrees). Everything is generic in the number type; what
# decides does so on `primal`, in the conversions of brouwer.jl.

# The orbit of six Keplerian elements with the epoch, center, frame, GM and radii given:
# ANISE `try_keplerian`. The six elements get the one number type.
function _orbit_from_kep(sma, ecc, inc, raan, aop, ta, et, center, frame, gm, radii)
    T = promote_type(
        typeof(sma), typeof(ecc), typeof(inc), typeof(raan), typeof(aop),
        typeof(ta), Float64
    )
    state = _kep_to_cart(T(sma), T(ecc), T(inc), T(raan), T(aop), T(ta), gm)
    return Orbit(state, et, center, frame, gm, radii)
end

# ANISE `try_keplerian_apsis_radii`: the semi-major axis and eccentricity of the two radii.
function _apsides_to_sma_ecc(periapsis, apoapsis)
    primal(apoapsis) > eps(Float64) || throw(
        SpiceError(
            :degeneratecase,
            "undefined: the radius of apoapsis is $apoapsis km, not above $(eps(Float64))"
        )
    )
    primal(periapsis) > eps(Float64) || throw(
        SpiceError(
            :degeneratecase,
            "undefined: the radius of periapsis is $periapsis km, not above $(eps(Float64))"
        )
    )
    sma = (apoapsis + periapsis) / 2
    return sma, apoapsis / sma - 1
end

"""
    keplerian_orbit([ctx,] et, center, frame; shape..., inc, raan, aop, anomaly[, gm]) -> Orbit

Return the [`Orbit`](@ref) of Keplerian elements about `center` in `frame` at `et`
(seconds past J2000 TDB); lengths in km, angles in rad. The counterpart of
[`conics_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/conics_c.html); `Orbit::try_keplerian`,
`try_keplerian_apsis_radii` and `try_keplerian_mean_anomaly` of ANISE.

- `shape`: `sma` and `ecc`, or the radii `periapsis` and `apoapsis`;
- `anomaly`: `ta` (true) or `ma` (mean);
- `gm` (km³/s²): defaults to `BODYnnn_GM` of the kernel pool, and then the radii are
  `BODYnnn_RADII` when present; with `gm` given the radii are NaN, and only names of `center` and
  `frame` are looked up in the kernels.

As in ANISE, a negative `ecc` and a semi-major axis of the wrong sign for `ecc` change
sign. The elements and the epoch may be dual numbers.

Throws `SpiceError`:
- `:badarg` unless exactly one shape pair and exactly one of `ta` and `ma` are given;
- `:degeneratecase` for a parabola (`|ecc - 1| < 1e-11`), a true anomaly beyond the
  asymptotes of a hyperbola, an infinite radius, a zero semi-parameter, an apsis radius
  not above `eps(Float64)`, or a mean anomaly that does not converge;
- `:kernelvarnotfound` if `gm` is not given and the pool has no GM for `center`.

# Examples
```jldoctest
julia> o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, inc=0.5, raan=1.0, aop=2.0,
                           ta=0.3, gm=398600.4415);

julia> round.(collect(values(keplerian(o))), digits=9)
6-element Vector{Float64}:
 8000.0
    0.1
    0.5
    1.0
    2.0
    0.3
```
"""
function keplerian_orbit(
        ks::KernelSet, et::Real, center::NameOrID, frame::NameOrID;
        sma = nothing, ecc = nothing, periapsis = nothing, apoapsis = nothing,
        inc, raan, aop, ta = nothing, ma = nothing, gm = nothing
    )
    shape = (sma !== nothing) + (ecc !== nothing) + (periapsis !== nothing) +
        (apoapsis !== nothing)
    (
        (sma !== nothing && ecc !== nothing && shape == 2) ||
            (periapsis !== nothing && apoapsis !== nothing && shape == 2)
    ) ||
        throw(
        SpiceError(
            :badarg, "keplerian_orbit: give exactly one of the pairs " *
                "sma and ecc, periapsis and apoapsis"
        )
    )
    (ta === nothing) == (ma === nothing) && throw(
        SpiceError(
            :badarg,
            "keplerian_orbit: give exactly one of ta and ma"
        )
    )
    c, f = body_code(ks, center), frame_id(ks, frame)
    if gm === nothing
        gm = bodvcd(ks, c, "GM")[1]
        radii = bodfnd(ks, c, "RADII") ? SVector{3, Float64}(bodvcd(ks, c, "RADII")) : _NO_RADII
    else
        radii = _NO_RADII
    end
    a, e = sma !== nothing ? (sma, ecc) : _apsides_to_sma_ecc(periapsis, apoapsis)
    # `try_keplerian_mean_anomaly`: the true anomaly of the mean one, with the given ecc
    nu = ta !== nothing ? ta : _ma_to_ta(promote(ma, e)...)
    return _orbit_from_kep(a, e, inc, raan, aop, nu, et, c, f, gm, radii)
end
keplerian_orbit(ctx::SpiceContext, et::Real, center::NameOrID, frame::NameOrID; kw...) =
    keplerian_orbit(kernelset(ctx), et, center, frame; kw...)
keplerian_orbit(et::Real, center::NameOrID, frame::NameOrID; kw...) =
    keplerian_orbit(default_context(), et, center, frame; kw...)

const _KEP_NAMES = (:sma, :ecc, :inc, :raan, :aop, :ta)

"""
    with_element(o, name, value) -> Orbit

Return a copy of `o` with the Keplerian element `name` (`:sma`, `:ecc`, `:inc`, `:raan`,
`:aop` or `:ta`) set to `value` (km, dimensionless or rad) and the other five, the epoch,
center, frame, GM and radii unchanged; `with_sma_km`,
`with_ecc`, `with_inc_deg` and the like of the `Orbit` of ANISE. This call has no CSPICE
counterpart.

The orbit is rebuilt as by [`keplerian_orbit`](@ref), with its sign rules; the angles are
not wrapped.

Throws `SpiceError(:badarg)` for another `name`, and `SpiceError(:degeneratecase)` where
an element of `o` or the rebuild is undefined: for a circular orbit, which has no true
anomaly, every `name` but `:ta` throws.

# Examples
```jldoctest
julia> o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, inc=0.5, raan=1.0, aop=2.0,
                           ta=0.3, gm=398600.4415);

julia> round(with_element(o, :sma, 9000.0).sma, digits=6)
9000.0
```
"""
function with_element(o::Orbit, name::Symbol, value::Real)
    name in _KEP_NAMES || throw(
        SpiceError(
            :badarg,
            "with_element: no Keplerian element :$name; the names are " * join(_KEP_NAMES, ", ")
        )
    )
    # ANISE evaluates the other five elements in this order before the conversion
    vals = map(n -> n === name ? value : element(o, n), _KEP_NAMES)
    return _orbit_from_kep(vals..., o.et, o.center, o.frame, o.gm, o.radii)
end

"""
    with_apsides(o, periapsis, apoapsis) -> Orbit

Return a copy of `o` with the apsis radii `periapsis` and `apoapsis` (km), and the
inclination, `raan`, `aop` and true anomaly unchanged;
`Orbit::with_apoapsis_periapsis_km` of ANISE, whose arguments come in the other order. This call has no CSPICE counterpart.

Throws `SpiceError(:degeneratecase)` for a radius not above `eps(Float64)`, for a circular
`o`, and where [`keplerian_orbit`](@ref) does.

# Examples
```jldoctest
julia> o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, inc=0.5, raan=1.0, aop=2.0,
                           ta=0.3, gm=398600.4415);

julia> p = with_apsides(o, 7000.0, 9000.0); round.((p.periapsis, p.apoapsis), digits=6)
(7000.0, 9000.0)
```
"""
function with_apsides(o::Orbit, periapsis::Real, apoapsis::Real)
    inc, raan, aop, ta = o.inc, o.raan, o.aop, o.ta
    sma, ecc = _apsides_to_sma_ecc(periapsis, apoapsis)
    return _orbit_from_kep(sma, ecc, inc, raan, aop, ta, o.et, o.center, o.frame, o.gm, o.radii)
end

"""
    two_body(o, et2) -> Orbit

Return `o` at the epoch `et2` (seconds past J2000 TDB, earlier or later, a dual number
too) by Kepler motion, without perturbations. The counterpart of
[`prop2b_c`](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/prop2b_c.html);
`Orbit::at_epoch` of ANISE, except for some hyperbolas ([deviations from ANISE](@ref "Deviations from ANISE")).

Throws `SpiceError(:degeneratecase)` for a circular orbit (`ecc < 1e-11`), which has no
mean anomaly, where another equinoctial element of `o` is undefined, and where
[`keplerian_orbit`](@ref) does.

# Examples
```jldoctest
julia> o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, inc=0.5, raan=1.0, aop=2.0,
                           ta=0.3, gm=398600.4415);

julia> t = two_body(o, 3000.0); t.et, round(t.sma, digits=6)
(3000.0, 8000.0)
```
"""
function two_body(o::Orbit{T}, et2::Real) where {T}
    et2 isa Float64 && o.et isa Float64 && et2 == o.et && return o
    a, h, k = _el_equinoctial_a(o), _el_equinoctial_h(o), _el_equinoctial_k(o)
    p, q, lambda = _el_equinoctial_p(o), _el_equinoctial_q(o), _el_equinoctial_lambda(o)
    n = _el_mean_motion(o)
    lambda2 = lambda + n * (et2 - o.et)
    sma, ecc, inc, raan, aop, ma = _equinoctial_to_kep(a, h, k, p, q, lambda2)
    # ANISE takes the mean anomaly as `lambda - atan2(h, k)`, which differs from the mean
    # anomaly of `o` by a multiple of 2 pi when raan + aop leaves (-pi, pi]. That does not
    # matter for an ellipse, whose motion has the period 2 pi, but a hyperbola has none, so
    # here its mean anomaly is the one of `o` advanced by n dt.
    primal(ecc) > 1 && (ma = _el_ma(o) + n * (et2 - o.et))
    ta = _ma_to_ta(promote(ma, ecc)...)
    return _orbit_from_kep(sma, ecc, inc, raan, aop, ta, et2, o.center, o.frame, o.gm, o.radii)
end
