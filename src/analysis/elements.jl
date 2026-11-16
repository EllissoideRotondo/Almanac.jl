# The scalar quantities of an `Orbit`. `element` is one `if` chain over the private
# `_el_<name>(o)`; `ELEMENT_NAMES` is the one ordered list of the names, in the order of
# the table of the design (later groups of elements are inserted in that order).

const ELEMENT_NAMES = (:x, :y, :z, :vx, :vy, :vz, :rmag, :vmag,
    :energy, :hmag, :hx, :hy, :hz, :light_time)

"""
    elements() -> NTuple{N,Symbol}

The names that [`element`](@ref) accepts, in the order of the manual. This call has no
CSPICE counterpart.
"""
elements() = ELEMENT_NAMES

"""
    element(o::Orbit, name::Symbol) -> T

One scalar quantity of the orbit `o`, in km, s and rad and their combinations; the
property `o.name` returns the same. This call has no CSPICE counterpart; the names are
those of `StateElement` of ANISE (`element(o, :rmag)` is `rmag_km`, `:hmag` is `hmag`,
`:energy` is `energy_km2_s2`, `:light_time` is `light_time`). The light time is the
distance divided by [`clight`](@ref), in s.

The vectors `position`, `velocity` and `hvec` are properties only: `element` throws for
them, so that a name known only at run time never gives a union of scalar and vector.

Throws `SpiceError(:badarg)` if `name` is not in [`elements`](@ref).

# Examples
```jldoctest
julia> o = Orbit([3.0, 4.0, 0.0, 0.0, 0.0, 2.0], 0.0, 399, 1, 4.0e5);

julia> element(o, :rmag), element(o, :hmag)
(5.0, 10.0)
```
"""
function element(o::Orbit, name::Symbol)
    name === :x && return _el_x(o)
    name === :y && return _el_y(o)
    name === :z && return _el_z(o)
    name === :vx && return _el_vx(o)
    name === :vy && return _el_vy(o)
    name === :vz && return _el_vz(o)
    name === :rmag && return _el_rmag(o)
    name === :vmag && return _el_vmag(o)
    name === :energy && return _el_energy(o)
    name === :hmag && return _el_hmag(o)
    name === :hx && return _el_hx(o)
    name === :hy && return _el_hy(o)
    name === :hz && return _el_hz(o)
    name === :light_time && return _el_light_time(o)
    throw(SpiceError(:badarg, "element: unknown name :$name; the names are " *
        join(elements(), ", ")))
end

"Throw the `:degeneratecase` error of `name`, with `what` for the reason."
_degenerate(name::Symbol, what::AbstractString) =
    throw(SpiceError(:degeneratecase, "element :$name is undefined: $what"))

_hvec(o::Orbit) = cross(o.position, o.velocity)

_el_x(o::Orbit) = getfield(o, :state)[1]
_el_y(o::Orbit) = getfield(o, :state)[2]
_el_z(o::Orbit) = getfield(o, :state)[3]
_el_vx(o::Orbit) = getfield(o, :state)[4]
_el_vy(o::Orbit) = getfield(o, :state)[5]
_el_vz(o::Orbit) = getfield(o, :state)[6]
_el_rmag(o::Orbit) = norm(o.position)
_el_vmag(o::Orbit) = norm(o.velocity)
_el_energy(o::Orbit) = _el_vmag(o)^2 / 2 - getfield(o, :gm) / _el_rmag(o)
_el_hmag(o::Orbit) = norm(_hvec(o))
_el_hx(o::Orbit) = _hvec(o)[1]
_el_hy(o::Orbit) = _hvec(o)[2]
_el_hz(o::Orbit) = _hvec(o)[3]
_el_light_time(o::Orbit) = _el_rmag(o) / clight()
