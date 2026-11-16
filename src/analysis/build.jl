# Building and editing an `Orbit` from Keplerian elements, and two-body propagation: the
# port of `try_keplerian`, `try_keplerian_apsis_radii`, `try_keplerian_mean_anomaly`, the
# `with_*` and `set_*` calls and `at_epoch` of src/astro/orbit.rs of ANISE (b3aaadc).
# Angles are in rad here (ANISE: degrees). Everything is generic in the number type; what
# decides does so on `primal`, in the conversions of brouwer.jl.

# The orbit of six Keplerian elements with the epoch, center, frame, GM and radii given:
# ANISE `try_keplerian`. The six elements get the one number type.
function _orbit_from_kep(sma, ecc, inc, raan, aop, ta, et, center, frame, gm, radii)
    T = promote_type(typeof(sma), typeof(ecc), typeof(inc), typeof(raan), typeof(aop),
        typeof(ta), Float64)
    state = _kep_to_cart(T(sma), T(ecc), T(inc), T(raan), T(aop), T(ta), gm)
    return Orbit(state, et, center, frame, gm, radii)
end

# ANISE `try_keplerian_apsis_radii`: the semi-major axis and eccentricity of the two radii.
function _apsides_to_sma_ecc(periapsis, apoapsis)
    primal(apoapsis) > eps(Float64) || throw(SpiceError(:degeneratecase,
        "undefined: the radius of apoapsis is $apoapsis km, not above $(eps(Float64))"))
    primal(periapsis) > eps(Float64) || throw(SpiceError(:degeneratecase,
        "undefined: the radius of periapsis is $periapsis km, not above $(eps(Float64))"))
    sma = (apoapsis + periapsis) / 2
    return sma, apoapsis / sma - 1
end

"""
    keplerian_orbit([ctx,] et, center, frame; sma, ecc, inc, raan, aop, ta[, gm]) -> Orbit
    keplerian_orbit([ctx,] et, center, frame; periapsis, apoapsis, inc, raan, aop, ma[, gm]) -> Orbit

The [`Orbit`](@ref) of Keplerian elements about `center` in `frame` at the epoch `et`
(seconds past J2000 TDB). Lengths are in km and angles in rad. No CSPICE call has this
interface; `conics` returns the state of conic elements about a body of given GM, and
this call matches `Orbit::try_keplerian`, `try_keplerian_apsis_radii` and
`try_keplerian_mean_anomaly` of ANISE.

The shape is given by `sma` and `ecc` or by `periapsis` and `apoapsis` (radii), and the
position on the orbit by the true anomaly `ta` or the mean anomaly `ma`; `inc`, `raan`
and `aop` are required. A negative `ecc` changes its sign, and a semi-major axis of the
wrong sign for the eccentricity (a positive one for `ecc > 1`, a negative one for
`ecc < 1`) changes its sign, as in ANISE.

`gm` (km³/s²) defaults to the `BODYnnn_GM` of the kernel pool, and then so do the radii
of `center`; with `gm` given, the radii are NaN (so the altitudes, `:latitude` and
`:height` of the orbit throw `SpiceError(:badarg)`) and the call reads no kernels when
`center` and `frame` are NAIF IDs. With a `KernelSet` or a
[`SpiceContext`](@ref) first the kernels are those of the snapshot or of the context;
without, those of the default context. The elements and the epoch may be dual numbers.

Throws `SpiceError(:badarg)` unless the shape is given by exactly one of the two pairs
(`sma` and `ecc`, or `periapsis` and `apoapsis`) and the position by exactly one of `ta`
and `ma`. It throws `SpiceError(:degeneratecase)`
where ANISE fails: a parabola (`|ecc - 1| < 1e-11`), a true anomaly beyond the asymptotes
of a hyperbola, an infinite radius, a zero semi-parameter, a radius of apsis not above
`eps(Float64)`, or a mean anomaly that does not converge; and
`SpiceError(:kernelvarnotfound)` if `gm` is not given and the pool has no GM for `center`.

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
function keplerian_orbit(ks::KernelSet, et::Real, center::NameOrID, frame::NameOrID;
        sma=nothing, ecc=nothing, periapsis=nothing, apoapsis=nothing,
        inc, raan, aop, ta=nothing, ma=nothing, gm=nothing)
    shape = (sma !== nothing) + (ecc !== nothing) + (periapsis !== nothing) +
        (apoapsis !== nothing)
    ((sma !== nothing && ecc !== nothing && shape == 2) ||
        (periapsis !== nothing && apoapsis !== nothing && shape == 2)) ||
        throw(SpiceError(:badarg, "keplerian_orbit: give exactly one of the pairs " *
            "sma and ecc, periapsis and apoapsis"))
    (ta === nothing) == (ma === nothing) && throw(SpiceError(:badarg,
        "keplerian_orbit: give exactly one of ta and ma"))
    c, f = body_code(ks, center), frame_id(ks, frame)
    if gm === nothing
        gm = bodvcd(ks, c, "GM")[1]
        radii = bodfnd(ks, c, "RADII") ? SVector{3,Float64}(bodvcd(ks, c, "RADII")) : _NO_RADII
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

A copy of the [`Orbit`](@ref) `o` with the Keplerian element `name` (one of `:sma`,
`:ecc`, `:inc`, `:raan`, `:aop`, `:ta`) set to `value` (km or rad) and the other five
unchanged. The epoch, center, frame, GM and radii stay. `value` is in km (`:sma`), rad
or dimensionless (`:ecc`). No CSPICE call has this interface; `oscltx` and `conics`
convert between the state and the elements, and this call matches `Orbit::with_sma_km`,
`with_ecc`, `with_inc_deg`, `with_raan_deg`, `with_aop_deg` and `with_ta_deg` of ANISE.

The orbit is rebuilt as by [`keplerian_orbit`](@ref), so its sign rules for `sma` and
`ecc` apply, and the angles are not wrapped. A circular orbit has no true anomaly, so
every name but `:ta` throws for it, as in ANISE.

Throws `SpiceError(:badarg)` for another `name`, and `SpiceError(:degeneratecase)` where
an element of `o` or the conversion is undefined (see [`element`](@ref) and
[`keplerian_orbit`](@ref)).

# Examples
```jldoctest
julia> o = keplerian_orbit(0.0, 399, 1; sma=8000.0, ecc=0.1, inc=0.5, raan=1.0, aop=2.0,
                           ta=0.3, gm=398600.4415);

julia> round(with_element(o, :sma, 9000.0).sma, digits=6)
9000.0
```
"""
function with_element(o::Orbit, name::Symbol, value::Real)
    name in _KEP_NAMES || throw(SpiceError(:badarg,
        "with_element: no Keplerian element :$name; the names are " * join(_KEP_NAMES, ", ")))
    # ANISE evaluates the other five elements in this order before the conversion
    vals = map(n -> n === name ? value : element(o, n), _KEP_NAMES)
    return _orbit_from_kep(vals..., o.et, o.center, o.frame, o.gm, o.radii)
end

"""
    with_apsides(o, periapsis, apoapsis) -> Orbit

A copy of the [`Orbit`](@ref) `o` with the radii of periapsis and apoapsis (km) given,
and its inclination, `raan`, `aop` and true anomaly unchanged. No CSPICE call has this
interface; it matches `Orbit::with_apoapsis_periapsis_km` of ANISE (whose arguments are
in the order apoapsis, periapsis).

Throws `SpiceError(:degeneratecase)` for a radius not above `eps(Float64)`, for a circular
orbit `o` (no true anomaly), and where [`keplerian_orbit`](@ref) does.

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

The [`Orbit`](@ref) `o` advanced to the epoch `et2` by two-body motion: the equinoctial
mean longitude moves by the mean motion times `et2 - o.et`, and the other equinoctial
elements stay. No force model is applied; `et2` may be earlier or later than `o.et`,
and may be a dual number, and then so is the state. No CSPICE call has this interface;
`prop2b` propagates a state by two-body motion, and this call matches `Orbit::at_epoch`
of ANISE. For a hyperbola it differs from ANISE; see
[Coverage and validation](@ref "Beyond CSPICE").

Throws `SpiceError(:degeneratecase)` where an equinoctial element is undefined (a circular
orbit has no mean anomaly, see [`element`](@ref)) or the conversion fails (see
[`keplerian_orbit`](@ref)).

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
