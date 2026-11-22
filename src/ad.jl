# Derivative rules for automatic differentiation packages other than ForwardDiff.
#
# Every quantity that depends on the epoch reaches it through a few internal functions,
# the "derivative points": `spk_jet`, `frame_rotation`, `frame_jet`, `abcorr_jet`,
# `abcorr_position` and `abcorr_state`. The derivative of each with respect to the
# epoch is known analytically: it is the next order of the same jet, as in `spkder` and
# `pxder`. The package extensions give these derivatives to the AD packages, so that no
# package differentiates the evaluators themselves:
#
# - packages that transform the code (Enzyme, Mooncake, Zygote and the other users of
#   ChainRulesCore) see a `Float64` epoch and call a rule at a derivative point. The rule
#   takes the value and the derivative from `_ad_rule`.
# - packages with a number type of their own (ReverseDiff, Tracker) pass that number as
#   the epoch. A method of each derivative point for that type evaluates the jet at the
#   plain epoch and assembles the result from it with `_ad_lift`, as the ForwardDiff
#   extension does for a `Dual`.
#
# The Orbit layer adds one derivative point that is not a function of the epoch: the
# Kepler solver `_ma_to_ta(M, e)` (`src/analysis/brouwer.jl`), with a method of `_ad_rule`
# there that returns the true anomaly and its partials with respect to M and e. The rules
# of ChainRulesCore, Enzyme and Mooncake take it; the tracked numbers do not need it.
#
# Two root solvers run on plain numbers and are constants for the AD packages: the
# bisection of `nearpt` (`_m_nearpt_plain`, for the multiplier of the near point) and
# that of `prop2b` (`_prop2b_root`, for the universal variable). The derivative of the
# root is that of the implicit function, from Newton steps at the converged root in the
# number type of the arguments (`_m_nearpt`, `_prop2b`). A `Float64` takes one step,
# through `_ad_pin(a, p)` (`src/spk_types/numeric.jl`), whose value is `p` and whose
# derivative, in the rules of ChainRulesCore, Enzyme and Mooncake, is that of `a`; so
# the value is the plain one bit for bit, and a reverse-mode package that differentiates
# the code after the root gets the derivative of the root. `_ad_pin` also gives the
# angles of an orbit (:ta, :aop, :raan) the value of an arc cosine and the derivatives
# of an arc tangent (`_AD_CONSTANT_ORBIT` holds those arc cosines).
#
# The derivative is that of the function that `spkder` and `pxder` differentiate. For an
# SPK type with a separate velocity series (type 3), the derivative of the position is
# the stored velocity. For a frame chain with a CK, dynamic or switch frame, the
# derivative of a rotation is the exact derivative, and not the rate that `sxform`
# returns.

"""
Functions on the way from a public call to a derivative point whose results do not
depend on the epoch: lookups by name or code (`bodvcd` and `bodfnd` for the GM and the
radii of an `Orbit`), the snapshot of a context, and the parsing of corrections and method
strings. The extensions mark them as constant, so that no AD
package looks into them (they hold strings, dictionaries, atomics and `try` blocks). The
root solvers `_m_nearpt_plain` and `_prop2b_root` are constants too: their callers take
the derivative of the root by a Newton step (see the header of this file).
"""
const _AD_CONSTANT = (
    body_code, frame_id, kernelset, timedefaults, dsktolerances,
    default_context, Aberration, bodvcd, bodfnd,
    _spk_abcorr, _geo_method, _geo_abcorr, _geo_body, _geo_namfrm, _geo_frinfo,
    _geo_fixfid, _geo_radii, _geo_shape, _lsk, _unitim_code,
    _m_nearpt_plain, _prop2b_root,
)

"The type to dispatch on for the callable `f`: `Type{f}` for a type, `typeof(f)` otherwise."
_ad_typeof(f) = f isa Type ? Type{f} : typeof(f)

"""
The derivative points, each with the types of its arguments and the position of the
epoch among them (`Real` there). The extensions for ChainRulesCore, Enzyme and Mooncake
generate their rules from this table; Mooncake takes the epoch as `Float64`.
"""
const _AD_POINTS = (
    (spk_jet, (KernelSet, Int32, Real, Int32, Int32, Val), 3),
    (frame_rotation, (KernelSet, Int32, Int32, Real), 4),
    (frame_jet, (KernelSet, Int32, Int32, Real, Val), 4),
    (abcorr_jet, (KernelSet, Int32, Real, Int32, Aberration, Int32, Val), 3),
    (abcorr_position, (KernelSet, Int32, Real, Int32, Aberration, Int32), 3),
    (abcorr_state, (KernelSet, Int32, Real, Int32, Aberration, Int32), 3),
)

# --- value and first derivative at an epoch ---------------------------------------------

"""
    _ad_rule(f, args...) -> (y, dy)

Return the value `y = f(args...)` of a derivative point and its derivative `dy` with
respect to the epoch, with the same structure as `y` (a tuple of vectors or matrices, a
matrix, or a tuple of a vector and the light time). The epoch can be any real number
that the point accepts, so a rule also works when the epoch is a `ForwardDiff.Dual`
(forward mode over a reverse-mode rule). No allocation.
"""
function _ad_rule(
        ::typeof(spk_jet), ks::KernelSet, targ::Int32, et::Real, ref::Int32,
        obs::Int32, ::Val{N}
    ) where {N}
    jet = spk_jet(ks, targ, et, ref, obs, Val(N + 1))
    return spk_jet(ks, targ, et, ref, obs, Val(N)), Base.tail(jet)
end

function _ad_rule(::typeof(frame_rotation), ks::KernelSet, from::Int32, to::Int32, et::Real)
    jet = _exact_frame_jet(ks, from, to, et, Val(1))
    return frame_rotation(ks, from, to, et), jet[2]
end

function _ad_rule(
        ::typeof(frame_jet), ks::KernelSet, from::Int32, to::Int32, et::Real,
        ::Val{N}
    ) where {N}
    jet = _exact_frame_jet(ks, from, to, et, Val(N + 1))
    return frame_jet(ks, from, to, et, Val(N)), Base.tail(jet)
end

function _ad_rule(
        ::typeof(abcorr_jet), ks::KernelSet, targ::Int32, et::Real, ref::Int32,
        corr::Aberration, obs::Int32, ::Val{N}
    ) where {N}
    jet = abcorr_jet(ks, targ, et, ref, corr, obs, Val(N + 1))
    return abcorr_jet(ks, targ, et, ref, corr, obs, Val(N)), Base.tail(jet)
end

function _ad_rule(
        ::typeof(abcorr_position), ks::KernelSet, targ::Int32, et::Real,
        ref::Int32, corr::Aberration, obs::Int32
    )
    r, lt = abcorr_position(ks, targ, et, ref, corr, obs)
    dr = abcorr_jet(ks, targ, et, ref, corr, obs, Val(1))[2]
    return (r, lt), (dr, _ad_lt_rate(r, dr))
end

function _ad_rule(
        ::typeof(abcorr_state), ks::KernelSet, targ::Int32, et::Real,
        ref::Int32, corr::Aberration, obs::Int32
    )
    s, lt = abcorr_state(ks, targ, et, ref, corr, obs)
    jet = abcorr_jet(ks, targ, et, ref, corr, obs, Val(2))
    r = s[SVector(1, 2, 3)]
    return (s, lt), (vcat(jet[2], jet[3]), _ad_lt_rate(r, jet[2]))
end

"The derivative `dy` of a rule times the tangent `s`, over the structure of `dy`."
_ad_scale(dy::Tuple, s) = map(d -> _ad_scale(d, s), dy)
_ad_scale(dy, s) = dy * s

"Derivative of the light time `|r| / c` for the position `r` with derivative `dr`."
@inline function _ad_lt_rate(r::Vec3, dr::Vec3)
    n = _norm3(r)
    return primal(n) == 0 ? zero(n) : _m_dot(r, dr) / (n * CLIGHT)
end

# --- results for a number type that carries derivatives ---------------------------------

"""
    _ad_lift(k, jet, x) -> value

Return order `k` of the jet at the epoch `x`, a number that carries derivatives, where
`jet[m+1]` is order `m` at `primal(x)`. The result is the Taylor polynomial of order `K =
_ad_order(typeof(x))` around the plain epoch, evaluated at `x`, so its value is order
`k` and its derivative parts of every order up to `K` are exact. Requires
`length(jet) > k + K`.
"""
@inline function _ad_lift(k::Int, jet::Tuple, x::T) where {T <: Real}
    K = _ad_order(T)
    d = x - primal(x)       # the derivative parts of `x`, with plain value 0
    return jet[k + 1] + _ad_terms(k, jet, d, Val(K))
end

@inline _ad_terms(k::Int, jet::Tuple, d, ::Val{K}) where {K} =
    sum(ntuple(j -> jet[k + j + 1] * (d^j / factorial(j)), Val(K)))

"Orders 0 to N at the epoch `x` from the jet `jet` of order `N + _ad_order(typeof(x))`."
@inline _ad_lift_jet(jet::Tuple, x::Real, ::Val{N}) where {N} =
    ntuple(i -> _ad_lift(i - 1, jet, x), Val(N + 1))

"""
Light time `|r| / c` of a position with derivative parts; zero for the zero vector. The
callers give it the plain value of the light time of the plain call (`_with_primal`),
which CSPICE takes from the norm of the position before the frame change.
"""
@inline function _ad_lt(r::Vec3{T}) where {T}
    iszero(primal(r[1])) && iszero(primal(r[2])) && iszero(primal(r[3])) && return zero(T)
    return sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]) / CLIGHT
end

# The methods below are what the methods of the derivative points for the number type
# of an extension call (see `@_ad_epoch_methods`).

function _ad_spk_jet(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, obs::Int32,
        v::Val{N}
    ) where {T <: Real, N}
    jet = spk_jet(ks, targ, primal(et), ref, obs, Val(N + _ad_order(T)))
    return _ad_lift_jet(jet, et, v)
end

function _ad_frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::T) where {T <: Real}
    jet = rotation_jet(ks, from, to, primal(et), Val(_ad_order(T)))
    return _ad_lift(0, jet, et)
end

# The value of order 1 of `frame_jet` is the rate of `sxform` at the plain epoch (built
# from the angular velocity of a CK, the blocks of `zzdynfrm_` and `tisbod_`), as for the
# rules of `_ad_rule`, and not the derivative of the rotation; the derivative parts of
# both orders are the next orders of the exact jet, as for a corrected state. The other
# orders of `frame_jet` at a plain epoch are those of `rotation_jet`.
function _ad_frame_jet(
        ks::KernelSet, from::Int32, to::Int32, et::T,
        v::Val{N}
    ) where {T <: Real, N}
    t = primal(et)
    jet = rotation_jet(ks, from, to, t, Val(N + _ad_order(T)))
    N == 1 || return _ad_lift_jet(jet, et, v)
    r0, dr0 = frame_jet(ks, from, to, t, v)
    return (
        _ad_lift(0, Base.setindex(jet, r0, 1), et),
        _ad_lift(1, Base.setindex(jet, dr0, 2), et),
    )
end

# `_exact_frame_jet` (`pxder`) keeps the exact jet: every order is a derivative of the
# rotation.
function _ad_exact_frame_jet(
        ks::KernelSet, from::Int32, to::Int32, et::T,
        v::Val{N}
    ) where {T <: Real, N}
    jet = rotation_jet(ks, from, to, primal(et), Val(N + _ad_order(T)))
    return _ad_lift_jet(jet, et, v)
end

function _ad_abcorr_jet(
        ks::KernelSet, targ::Int32, et::T, ref::Int32, corr::Aberration,
        obs::Int32, v::Val{N}
    ) where {T <: Real, N}
    jet = abcorr_jet(ks, targ, primal(et), ref, corr, obs, Val(N + _ad_order(T)))
    return _ad_lift_jet(jet, et, v)
end

# The value of the corrected position (state) is that of the plain call, as for a
# `Dual`; see the extension for ForwardDiff.
function _ad_abcorr_position(
        ks::KernelSet, targ::Int32, et::T, ref::Int32,
        corr::Aberration, obs::Int32
    ) where {T <: Real}
    t = primal(et)
    jet = abcorr_jet(ks, targ, t, ref, corr, obs, Val(_ad_order(T)))
    r0, lt0 = abcorr_position(ks, targ, t, ref, corr, obs)
    r = _ad_lift(0, Base.setindex(jet, r0, 1), et)
    return r, _with_primal(_ad_lt(r), lt0)
end

function _ad_abcorr_state(
        ks::KernelSet, targ::Int32, et::T, ref::Int32,
        corr::Aberration, obs::Int32
    ) where {T <: Real}
    t = primal(et)
    jet = abcorr_jet(ks, targ, t, ref, corr, obs, Val(_ad_order(T) + 1))
    s0, lt0 = abcorr_state(ks, targ, t, ref, corr, obs)
    r = _ad_lift(0, Base.setindex(jet, s0[SVector(1, 2, 3)], 1), et)
    v = _ad_lift(1, Base.setindex(jet, s0[SVector(4, 5, 6)], 2), et)
    return vcat(r, v), _with_primal(_ad_lt(r), lt0)
end

"""
    @_ad_epoch_methods T

Define the methods of the derivative points for an epoch of the number type `T` of an
extension (`ForwardDiff.Dual`, the `TrackedReal` of ReverseDiff and Tracker). Each one
takes the jet at the plain epoch and assembles the result with `_ad_lift`. The method of
`frame_jet` for `Val{0}` resolves the ambiguity with `frame_jet(..., et::Real, ::Val{0})`.
`_exact_frame_jet` (the jet of `pxder`) gets a method too, because its generic method
calls `frame_jet`, whose order 1 is the rate of `sxform`.
"""
macro _ad_epoch_methods(T)
    A = Almanac
    return esc(
        quote
            $A.spk_jet(
                ks::$KernelSet, targ::Int32, et::$T, ref::Int32, obs::Int32,
                v::Val
            ) = $A._ad_spk_jet(ks, targ, et, ref, obs, v)
            $A.frame_rotation(ks::$KernelSet, from::Int32, to::Int32, et::$T) =
                $A._ad_frame_rotation(ks, from, to, et)
            $A.frame_jet(ks::$KernelSet, from::Int32, to::Int32, et::$T, v::Val) =
                $A._ad_frame_jet(ks, from, to, et, v)
            $A.frame_jet(ks::$KernelSet, from::Int32, to::Int32, et::$T, v::Val{0}) =
                $A._ad_frame_jet(ks, from, to, et, v)
            $A._exact_frame_jet(ks::$KernelSet, from::Int32, to::Int32, et::$T, v::Val) =
                $A._ad_exact_frame_jet(ks, from, to, et, v)
            $A.abcorr_jet(
                ks::$KernelSet, targ::Int32, et::$T, ref::Int32,
                corr::$Aberration, obs::Int32, v::Val
            ) = $A._ad_abcorr_jet(ks, targ, et, ref, corr, obs, v)
            $A.abcorr_position(
                ks::$KernelSet, targ::Int32, et::$T, ref::Int32,
                corr::$Aberration, obs::Int32
            ) = $A._ad_abcorr_position(ks, targ, et, ref, corr, obs)
            $A.abcorr_state(
                ks::$KernelSet, targ::Int32, et::$T, ref::Int32,
                corr::$Aberration, obs::Int32
            ) = $A._ad_abcorr_state(ks, targ, et, ref, corr, obs)
        end
    )
end
