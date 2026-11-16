# Opt-in comparison of the Orbit calls of Almanac.jl with the ANISE Python package.
#
#   julia --startup-file=no --project=test/anise test/anise/runtests.jl
#
# The group is not in `GROUPS` of test/runtests.jl: it needs Python, which CondaPkg installs
# on the first run (network), and the `anise` wheel pinned in CondaPkg.toml.
#
# The reference of the port is the source of ANISE main at commit b3aaadc. The pin is the newest
# release on PyPI, 1.0.0b0. Where 1.0.0b0 differs from main by definition the quantity is in
# `KNOWN_DIFFERENT_IN_PIN.pin`: it is compared and reported, but not tested. The set is emptied
# when the pin moves to the first release that contains main's definitions. The deliberate
# deviations of Almanac.jl (docs/src/explanation/coverage.md, "Beyond CSPICE") stay in
# `KNOWN_DIFFERENT_IN_PIN.deliberate`.
#
# Both libraries get the same Cartesian state, epoch and GM; the radii are given to the ANISE
# `Frame` as an `Ellipsoid`, so no ANISE data file is needed. A difference is
# |a - b| / max(|a|, |b|, 1), where for an angle (rad) |a - b| is the circular difference; where
# one library has no value (a thrown error, or NaN) and the other has, the state is counted as a
# definedness mismatch, and where neither has one the state counts as "both undefined". The
# range convention of each angle ((-pi, pi] or [0, 2 pi]) is checked separately, because the
# circular difference cannot see it. The tolerance of the conditioned quantities is in units of
# their error model: angles of near-circular orbits in eps/ecc, rates in eps/(2 ms), the
# rounding error of the difference quotient of ANISE.

using Test, Random, Printf, StaticArrays, LinearAlgebra
using Almanac
using PythonCall

const AA = pyimport("anise.astro")
const AT = pyimport("anise.time")

const FULL = get(ENV, "ALMANAC_TEST_FULL", "0") == "1"
const N = FULL ? 20000 : 2000                     # states per class
# the variables and defaults of test/helpers.jl (not included: it loads the CSPICE oracle)
const SEED = parse(Int, get(ENV, "ALMANAC_TEST_SEED", "20261011"))
const T0 = time()

const ANISE_VERSION = pyconvert(String, pyimport("importlib.metadata").version("anise"))

# ---------------------------------------------------------------------------------------------
# What is known to differ
# ---------------------------------------------------------------------------------------------

# Evidence for the entries of `pin`: the sdist of 1.0.0b0, src/astro/orbit.rs, compared with
# the source of main at b3aaadc (diff of orbit.rs: ta check in try_keplerian, fpa, vinf_*, TNW).
#   https://files.pythonhosted.org/packages/3f/68/054a14788661f4db8803f4c4cd3fdee5ffbf3304efb3cadc041d829e519f/anise-1.0.0b0.tar.gz
#   sha256 b01ab498ba5f0949c55cdd43acf7ecb56a231df3ed49db870a536f704507ecf4
# (ma_deg and hyperbolic_anomaly_deg are identical in both, so they are tested.)
const KNOWN_DIFFERENT_IN_PIN = (
    # (a) 1.0.0b0 differs from main by definition (checked in the sdist of 1.0.0b0, orbit.rs).
    pin = Dict(
        "el:fpa" => "fpa_deg: 1.0.0b0 has cos_fpa = 1 + e cos(nu)/denom, main (1 + e cos(nu))/denom",
        "vinf_periapsis" => "1.0.0b0 divides mu by vmag^2, main by C3",
        "vinf_turn_angle" => "1.0.0b0 uses vmag^2 for C3, main C3",
        "*:hypin" => "1.0.0b0 rejects an inbound hyperbolic true anomaly in try_keplerian (between_0_360, so only the outbound leg passes); main accepts it. Applies to keplerian_orbit, with_element, with_apsides and two_body results that are inbound hyperbolas",
        "lf:tnw" => "no TNW frame in 1.0.0b0 (not in the Python API: no dcm_from_tnw_to_inertial, no tnw_difference)",
        "bplane" => "no b-plane calls in the Python API of 1.0.0b0 (compared with GMAT in test_analysis_orbit.jl)",
    ),
    # (b) deliberate deviations of Almanac.jl, which stay after the pin moves.
    deliberate = Dict(
        "two_body:hyp_wrap" => "R11 (also the keys ending in :hyp_wrap: the rates of a local frame, which ANISE differences over at_epoch): ANISE recovers M as lambda - atan2(h, k), off by 2 pi j for a hyperbola when raan + aop > pi; Almanac advances M_h directly",
    ),
)

# ---------------------------------------------------------------------------------------------
# Bodies and sampling
# ---------------------------------------------------------------------------------------------

struct Body
    id::Int
    gm::Float64
    radii::NTuple{3,Float64}
    frame::Py          # ANISE frame, same GM and radii
end
function Body(id, gm, radii)
    shape = AA.Ellipsoid(radii[1], radii[3], radii[2])
    return Body(id, gm, radii, AA.Frame(id, 1, gm, shape))
end
const EARTH = Body(399, 398600.4415, (6378.1363, 6378.1363, 6356.7517))
const MARS = Body(499, 42828.37, (3396.19, 3396.19, 3376.2))

struct Sample
    cls::Symbol
    body::Body
    kep::NTuple{6,Float64}      # sma, ecc, inc, raan, aop, ta (km, rad)
    et::Float64
    o::Orbit{Float64}
    a::Py                       # ANISE orbit
end

const CLASSES = (:elliptic, :circular, :equatorial, :retrograde, :hyp_out, :hyp_in)

anise_epoch(et) = AT.Epoch.from_tdb_seconds(et)
anise_orbit(s::AbstractVector, et, body::Body) =
    AA.Orbit.from_cartesian(s[1], s[2], s[3], s[4], s[5], s[6], anise_epoch(et), body.frame)
almanac_orbit(s, et, body::Body) = Orbit(SVector{6,Float64}(s), et, body.id, 1, body.gm, collect(body.radii))

function make_sample(rng, cls, body)
    u(lo, hi) = lo + (hi - lo) * rand(rng)
    R = (body.radii[1] + body.radii[2]) / 2
    twopi = 2pi
    inc = u(0.05, pi - 0.05)
    ecc = u(0.001, 0.9)
    rp = u(1.05, 6.0) * R
    if cls === :circular
        ecc = 10.0^u(-9, -6)
    elseif cls === :equatorial
        inc = 10.0^u(-7, -3)
        rand(rng) < 0.5 && (inc = pi - inc)
        ecc = u(0.001, 0.5)
    elseif cls === :retrograde
        inc = u(pi / 2 + 0.05, pi - 0.05)
    elseif cls === :hyp_out || cls === :hyp_in
        ecc = u(1.05, 4.0)
        rp = u(1.05, 10.0) * R
    end
    sma = rp / (1 - ecc)
    raan, aop = u(0, twopi), u(0, twopi)
    if ecc > 1
        tinf = acos(-1 / ecc)
        ta = u(0.02, 0.95 * tinf)
        cls === :hyp_in && (ta = -ta)
    else
        ta = u(0, twopi)
    end
    et = Float64(rand(rng, -10^9:10^9))
    ko = keplerian_orbit(et, body.id, 1; sma, ecc, inc, raan, aop, ta, gm=body.gm)
    st = collect(ko.state)
    return Sample(cls, body, (sma, ecc, inc, raan, aop, ta), et,
        almanac_orbit(st, et, body), anise_orbit(st, et, body))
end

# ---------------------------------------------------------------------------------------------
# Accumulators
# ---------------------------------------------------------------------------------------------

mutable struct Acc
    n::Int                  # both defined
    both_undef::Int
    mismatch::Int           # exactly one defined
    maxd::Float64
    worst::String
    mexample::String
    maxn::Float64           # the tested value: the difference in units of its error model
    worstn::String
    normalised::Bool
    range::Int              # angle outside its range convention
    rexample::String
end
Acc() = Acc(0, 0, 0, 0.0, "", "", 0.0, "", false, 0, "")
const ACCS = Dict{String,Acc}()
acc(key) = get!(ACCS, key, Acc())

describe(s::Sample) = @sprintf("%s %d et=%.0f state=%s", s.cls, s.body.id, s.et, repr(collect(s.o.state)))

function record!(key, d, s, extra=""; nv=d)
    a = acc(key)
    a.n += 1
    if d > a.maxd || a.worst == ""
        a.maxd = d
        a.worst = describe(s) * extra
    end
    if nv > a.maxn || a.worstn == ""
        a.maxn = nv
        a.worstn = describe(s) * extra
    end
    nv === d || (a.normalised = true)
    return a
end
# x, y: the value of Almanac and of ANISE, `nothing` where undefined. `dist(x, y)` is the difference.
function pair!(key, x, y, s, dist; extra="", nfun=nothing)
    a = acc(key)
    if x === nothing && y === nothing
        a.both_undef += 1
    elseif x === nothing || y === nothing
        a.mismatch += 1
        a.mexample == "" && (a.mexample = (x === nothing ? "Almanac undefined, ANISE " * string(y) :
            "ANISE undefined, Almanac " * string(x)) * " | " * describe(s) * extra)
    else
        d = dist(x, y)
        record!(key, d, s, extra; nv=nfun === nothing ? d : nfun(x, y))
    end
    return nothing
end

function scalar_dist(x, y; angle=false)
    (isinf(x) || isinf(y)) && return x == y ? 0.0 : Inf
    d = angle ? abs(mod(x - y + pi, 2pi) - pi) : abs(x - y)
    return d / max(abs(x), abs(y), 1.0)
end
outcome!(key, x, y, s; angle=false, extra="", nfun=nothing) =
    pair!(key, x, y, s, (p, q) -> scalar_dist(p, q; angle), extra=extra, nfun=nfun)

# Error model of the angles: they come from the arc cosine of a normalised dot product, whose
# rounding error eps is amplified by 1/sin(angle), and the vectors it uses (the eccentricity
# vector, the node vector) have a relative error eps/ecc and eps/sin(inc). The difference is
# measured in units of eps (1/ecc + 1/sin inc) / max(|sin angle|, 1e-8), the angle being the true
# anomaly (ta, ma, ea), the argument of periapsis (aop), or the larger of the two amplifiers
# (the rebuilt states). ecc is that of the sample, at most 1.
function cond_scale(s::Sample, which::Symbol)
    f(x) = 1 / max(abs(sin(x)), 1e-8)
    base = 1 / min(s.o.ecc, 1.0) + 1 / max(sin(s.o.inc), 1e-8)
    fac = which === :ta ? f(s.o.ta) : which === :aop ? f(s.o.aop) : max(f(s.o.ta), f(s.o.aop))
    return eps() * base * fac
end
const COND_KEYS = Dict("el:ta" => :ta, "el:ma" => :ta, "el:ea" => :ta, "el:aop" => :aop,
    "with_element:ta" => :both, "with_element:aop" => :both, "with_apsides" => :both)
# Error model of a rate: ANISE differences two at_epoch states 1 ms apart. at_epoch solves the
# Kepler equation with a Newton iteration that stops when two iterates differ by MA_EPSILON =
# 1e-12 (src/astro/utils.rs) and returns the older one, so each state carries an error up to
# 1e-12 rad in the anomaly and the rate an error up to 1e-12/(2 ms) = 5e-10 (absolute). The
# difference of a rate is measured in units of that.
const RATE_UNIT = 2e-3 / 1e-12
is_rate(key) = occursin(r"^lf:\w+:Rdot", key) || occursin(r"^rel:\w+:vel", key)

# Suffix of the key of the quantities that ANISE differences over at_epoch (rates): an
# outbound hyperbola with raan + aop <= pi is as good as an ellipse, :hyp_wrap is R11, and an
# inbound leg is rejected by 1.0.0b0.
function hyp_suffix(s::Sample)
    s.o.ecc > 1 || return ""
    s.o.ta > pi && return ":hypin"
    return s.o.raan + s.o.aop > pi ? ":hyp_wrap" : ":hyp"
end

# range convention of the angle elements (both ends inclusive, with a slack of 1e-12 rad)
const RANGES = Dict(:inc => (0.0, pi), :raan => (0.0, 2pi), :aop => (0.0, 2pi), :ta => (0.0, 2pi),
    :ea => (-pi, pi), :ha => (0.0, 2pi), :aol => (0.0, 2pi), :tlong => (0.0, 2pi),
    :right_ascension => (0.0, 2pi), :declination => (-pi / 2, pi / 2),
    :latitude => (-pi / 2, pi / 2), :longitude => (-pi, pi),
    :brouwer_inc => (0.0, pi), :brouwer_raan => (0.0, 2pi), :brouwer_aop => (0.0, 2pi),
    :brouwer_ma => (0.0, 2pi), :ma => (0.0, 2pi))   # :ma of an ellipse only
function check_range!(name, x, y, s)
    haskey(RANGES, name) || return
    name === :ma && s.o.ecc > 1 && return
    lo, hi = RANGES[name]
    for (who, v) in (("Almanac", x), ("ANISE", y))
        (v === nothing || lo - 1e-12 <= v <= hi + 1e-12) && continue
        a = acc("el:$name")
        a.range += 1
        a.rexample == "" && (a.rexample = "$who $v outside [$lo, $hi] | " * describe(s))
    end
end

# Python errors seen, by message (ANISE raises for an undefined quantity)
const PYERRORS = Dict{String,Int}()
function pyerr!(e)
    e isa PyException || return false
    # ANISE raises the type `Exception` (a PhysicsError); any other type (an AttributeError of a
    # misspelled method, a TypeError) is a failure of this file, not an undefined quantity
    pyis(e.t, pybuiltins.Exception) || return false
    k = first(split(sprint(showerror, e), '\n'))
    k = replace(k, r"-?\d+(\.\d+)?(e[-+]?\d+)?" => "<x>")
    k = length(k) > 110 ? k[1:110] : k
    PYERRORS[k] = get(PYERRORS, k, 0) + 1
    return true
end

# Almanac side: a value, or nothing where it throws :degeneratecase (or the quantity is NaN)
function almanac_val(f)
    try
        v = f()
        return (v isa Real && isfinite(v) || (v isa Real && isinf(v))) ? Float64(v) : nothing
    catch e
        e isa Almanac.SpiceError && return nothing
        rethrow()
    end
end
# ANISE side: a value, or nothing where it raises (or the value is NaN)
function anise_val(f)
    try
        v = f()
        v = v isa Py ? pyconvert(Float64, v) : Float64(v)
        return isnan(v) ? nothing : v
    catch e
        pyerr!(e) && return nothing
        rethrow()
    end
end

# ---------------------------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------------------------

# name, ANISE accessor, kind of accessor (:attr, :call, :dur), unit (:one, :deg, :degs), measure
const ELEMENT_TABLE = [
    (:x, "x_km", :attr, :one, :lin), (:y, "y_km", :attr, :one, :lin), (:z, "z_km", :attr, :one, :lin),
    (:vx, "vx_km_s", :attr, :one, :lin), (:vy, "vy_km_s", :attr, :one, :lin), (:vz, "vz_km_s", :attr, :one, :lin),
    (:rmag, "rmag_km", :call, :one, :lin), (:vmag, "vmag_km_s", :call, :one, :lin),
    (:sma, "sma_km", :call, :one, :lin), (:ecc, "ecc", :call, :one, :lin),
    (:inc, "inc_deg", :call, :deg, :ang), (:raan, "raan_deg", :call, :deg, :ang),
    (:aop, "aop_deg", :call, :deg, :ang), (:ta, "ta_deg", :call, :deg, :ang),
    (:ma, "ma_deg", :call, :deg, :ang), (:ea, "ea_deg", :call, :deg, :ang),
    (:ha, "hyperbolic_anomaly_deg", :call, :deg, :ang),
    (:aol, "aol_deg", :call, :deg, :ang), (:tlong, "tlong_deg", :call, :deg, :ang),
    (:period, "period", :dur, :one, :lin), (:mean_motion, "mean_motion_deg_s", :call, :degs, :lin),
    (:ta_dot, "ta_dot_deg_s", :call, :degs, :lin),
    (:energy, "energy_km2_s2", :call, :one, :lin), (:c3, "c3_km2_s2", :call, :one, :lin),
    (:fpa, "fpa_deg", :call, :deg, :ang), (:semi_parameter, "semi_parameter_km", :call, :one, :lin),
    (:semi_minor_axis, "semi_minor_axis_km", :call, :one, :lin),
    (:hmag, "hmag", :call, :one, :lin), (:hx, "hx", :call, :one, :lin), (:hy, "hy", :call, :one, :lin),
    (:hz, "hz", :call, :one, :lin),
    (:periapsis, "periapsis_km", :call, :one, :lin), (:apoapsis, "apoapsis_km", :call, :one, :lin),
    (:periapsis_altitude, "periapsis_altitude_km", :call, :one, :lin),
    (:apoapsis_altitude, "apoapsis_altitude_km", :call, :one, :lin),
    (:sma_altitude, "sma_altitude_km", :call, :one, :lin), (:altitude, "altitude_km", :call, :one, :lin),
    (:right_ascension, "right_ascension_deg", :call, :deg, :ang),
    (:declination, "declination_deg", :call, :deg, :ang),
    (:velocity_declination, "velocity_declination_deg", :call, :deg, :ang),
    (:latitude, "latitude_deg", :call, :deg, :ang), (:longitude, "longitude_deg", :call, :deg, :ang),
    (:height, "height_km", :call, :one, :lin), (:light_time, "light_time", :dur, :one, :lin),
    (:equinoctial_a, "equinoctial_a_km", :call, :one, :lin),
    (:equinoctial_h, "equinoctial_h", :call, :one, :lin), (:equinoctial_k, "equinoctial_k", :call, :one, :lin),
    (:equinoctial_p, "equinoctial_p", :call, :one, :lin), (:equinoctial_q, "equinoctial_q", :call, :one, :lin),
    (:equinoctial_lambda, "equinoctial_lambda_mean_deg", :call, :deg, :lin),
    (:brouwer_sma, "sma_brouwer_short_km", :call, :one, :lin),
    (:brouwer_ecc, "ecc_brouwer_short", :call, :one, :lin),
    (:brouwer_inc, "inc_brouwer_short_deg", :call, :deg, :ang),
    (:brouwer_raan, "raan_brouwer_short_deg", :call, :deg, :ang),
    (:brouwer_aop, "aop_brouwer_short_deg", :call, :deg, :ang),
    (:brouwer_ma, "ma_brouwer_short_deg", :call, :deg, :ang),
]
# the b-plane names have no counterpart in the Python API of 1.0.0b0
const NOT_IN_PIN = [:bplane_bt, :bplane_br, :bplane_mag, :bplane_angle, :bplane_ltof]
@assert Set(first.(ELEMENT_TABLE)) ∪ Set(NOT_IN_PIN) == Set(elements())

function anise_element(a::Py, (name, pyname, kind, unit, _))
    v = if kind === :attr
        pygetattr(a, pyname)
    elseif kind === :call
        pygetattr(a, pyname)()
    else
        pygetattr(a, pyname)().to_seconds()
    end
    x = pyconvert(Float64, v)
    return unit === :deg ? deg2rad(x) : unit === :degs ? deg2rad(x) : x
end

function compare_elements!(s::Sample)
    for row in ELEMENT_TABLE
        name, _, _, _, measure = row
        x = almanac_val(() -> element(s.o, name))
        y = anise_val(() -> anise_element(s.a, row))
        angle = measure === :ang
        # the hyperbolic mean anomaly is unbounded, not an angle modulo 2 pi
        name === :ma && s.o.ecc > 1 && (angle = false)
        check_range!(name, x, y, s)
        outcome!("el:$name", x, y, s; angle,
            nfun=haskey(COND_KEYS, "el:$name") ? (p, q) -> scalar_dist(p, q; angle) / cond_scale(s, COND_KEYS["el:$name"]) : nothing)
    end
end

# ---------------------------------------------------------------------------------------------
# Local frames and relative states
# ---------------------------------------------------------------------------------------------

const ANISE_LF = (ric="ric", rcn="rcn", vnc="vnc")
pymat(m) = pyconvert(Matrix{Float64}, m)

function maxrel(A, B, floor)
    return maximum(abs.(A .- B)) / max(maximum(abs.(A)), maximum(abs.(B)), floor)
end

function compare_local_frames!(s::Sample)
    for kind in (:ric, :rcn, :vnc)
        pyname = getfield(ANISE_LF, kind)
        # Almanac
        R = try
            Matrix(local_frame(s.o, kind))
        catch e
            e isa Almanac.SpiceError ? nothing : rethrow()
        end
        S = try
            local_frame_state(s.o, kind)
        catch e
            e isa Almanac.SpiceError ? nothing : rethrow()
        end
        Rd = S === nothing ? nothing : Matrix(S[4:6, 1:3])
        # ANISE
        d3 = try
            pymat(pygetattr(s.a, "dcm3x3_from_$(pyname)_to_inertial")().rot_mat)
        catch e
            pyerr!(e) ? nothing : rethrow()
        end
        d6 = nothing
        try
            d = pygetattr(s.a, "dcm_from_$(pyname)_to_inertial")()
            dt = d.rot_mat_dt
            d6 = pyis(dt, pybuiltins.None) ? nothing : pymat(dt)
        catch e
            pyerr!(e) || rethrow()
        end
        d3 !== nothing && any(isnan, d3) && (d3 = nothing)
        d6 !== nothing && any(isnan, d6) && (d6 = nothing)
        pair!("lf:$kind:R", R, d3, s, (p, q) -> maximum(abs.(p .- q)))
        pair!("lf:$kind:Rdot" * hyp_suffix(s), Rd, d6, s, (p, q) -> maxrel(p, q, 1e-12),
            nfun=(p, q) -> maximum(abs.(p .- q)) * RATE_UNIT)
    end
end

function compare_relative!(rng, s::Sample)
    dr = [(rand(rng) - 0.5) * 2 * 10.0^rand(rng, 0:2) for _ in 1:3]
    dv = [(rand(rng) - 0.5) * 2 * 10.0^rand(rng, -4:-2) for _ in 1:3]
    st = collect(s.o.state) .+ [dr; dv]
    o2 = almanac_orbit(st, s.et, s.body)
    a2 = anise_orbit(st, s.et, s.body)
    for (kind, fname) in ((:ric, "ric_difference_unchecked"), (:rcn, "ric_difference_unchecked"),
                          (:vnc, "vnc_difference_unchecked"))
        x = try
            relative_state(o2, s.o, kind)
        catch e
            e isa Almanac.SpiceError ? nothing : rethrow()
        end
        y = try
            pyconvert(Vector{Float64}, pygetattr(s.a, fname)(a2).cartesian_pos_vel())
        catch e
            pyerr!(e) ? nothing : rethrow()
        end
        y !== nothing && any(isnan, y) && (y = nothing)
        for (part, r) in ((:pos, 1:3), (:vel, 4:6))
            fl = part === :pos ? 1.0 : 1e-3
            vel = part === :vel
            pair!("rel:$kind:$part" * (vel ? hyp_suffix(s) : ""), x, y, s,
                (p, q) -> norm(p[r] - q[r]) / max(norm(p[r]), norm(q[r]), fl);
                nfun=vel ? (p, q) -> norm(p[r] - q[r]) * RATE_UNIT / norm(dr) : nothing)
        end
    end
end

# ---------------------------------------------------------------------------------------------
# Construction, editing and propagation
# ---------------------------------------------------------------------------------------------

# ANISE state of an orbit, or nothing where the call raises
function anise_state(f)
    try
        v = pyconvert(Vector{Float64}, f().cartesian_pos_vel())
        return any(isnan, v) ? nothing : v
    catch e
        pyerr!(e) && return nothing
        rethrow()
    end
end
# compares two states (position and velocity vectors separately, relative to their norms)
function compare_states!(key, xo, y, s; extra="")
    pair!(key, xo === nothing ? nothing : collect(xo.state), y, s, (x, y) ->
        max(norm(x[1:3] - y[1:3]) / max(norm(x[1:3]), norm(y[1:3]), 1.0),
            norm(x[4:6] - y[4:6]) / max(norm(x[4:6]), norm(y[4:6]), 1e-3)); extra=extra,
        nfun=haskey(COND_KEYS, key) ? (p, q) -> max(norm(p[1:3] - q[1:3]) / max(norm(p[1:3]), norm(q[1:3]), 1.0),
            norm(p[4:6] - q[4:6]) / max(norm(p[4:6]), norm(q[4:6]), 1e-3)) / cond_scale(s, COND_KEYS[key]) : nothing)
end

# the key gets the suffix ":hypin" where the Almanac result is an inbound hyperbola
function key_for(base, xo)
    xo === nothing && return base
    return (xo.ecc > 1 && xo.ta > pi) ? base * ":hypin" : base
end

function almanac_try(f)
    try
        return f()
    catch e
        e isa Almanac.SpiceError && return nothing
        rethrow()
    end
end

function compare_construction!(rng, s::Sample)
    sma, ecc, inc, raan, aop, ta = s.kep
    et, body = s.et, s.body
    ep = anise_epoch(et)
    deg(x) = rad2deg(x)
    # from Keplerian elements with the true anomaly
    xo = almanac_try(() -> keplerian_orbit(et, body.id, 1; sma, ecc, inc, raan, aop, ta, gm=body.gm))
    y = anise_state(() -> AA.Orbit.from_keplerian(sma, ecc, deg(inc), deg(raan), deg(aop), deg(ta), ep, body.frame))
    compare_states!(key_for("keplerian_orbit", xo), xo, y, s)
    # with the mean anomaly (the mean anomaly of the sample)
    ma = almanac_val(() -> s.o.ma)
    if ma !== nothing
        xo = almanac_try(() -> keplerian_orbit(et, body.id, 1; sma, ecc, inc, raan, aop, ma, gm=body.gm))
        y = anise_state(() -> AA.Orbit.from_keplerian_mean_anomaly(sma, ecc, deg(inc), deg(raan), deg(aop), deg(ma), ep, body.frame))
        compare_states!(key_for("keplerian_orbit_ma", xo), xo, y, s)
    end
    # with the apsides (ellipses)
    if ecc < 1
        rp, ra = sma * (1 - ecc), sma * (1 + ecc)
        xo = almanac_try(() -> keplerian_orbit(et, body.id, 1; periapsis=rp, apoapsis=ra, inc, raan, aop, ta, gm=body.gm))
        y = anise_state(() -> AA.Orbit.from_keplerian_apsis_radii(ra, rp, deg(inc), deg(raan), deg(aop), deg(ta), ep, body.frame))
        compare_states!(key_for("keplerian_orbit_apsides", xo), xo, y, s)
    end
    # with_element: the six names
    u(lo, hi) = lo + (hi - lo) * rand(rng)
    values = (sma=sma * u(0.8, 1.25), ecc=clamp(ecc * u(0.8, 1.2), 0.0005, ecc < 1 ? 0.95 : 6.0),
              inc=u(0.05, pi - 0.05), raan=u(0, 2pi), aop=u(0, 2pi), ta=ecc > 1 ? u(0.02, 0.5) : u(0, 2pi))
    pynames = (sma="with_sma_km", ecc="with_ecc", inc="with_inc_deg", raan="with_raan_deg",
               aop="with_aop_deg", ta="with_ta_deg")
    for name in (:sma, :ecc, :inc, :raan, :aop, :ta)
        v = getfield(values, name)
        # keep the shape of the orbit: a new sma of a hyperbola stays negative
        name === :sma && ecc > 1 && (v = -abs(v))
        xo = almanac_try(() -> with_element(s.o, name, v))
        pv = name in (:inc, :raan, :aop, :ta) ? deg(v) : v
        y = anise_state(() -> pygetattr(s.a, getfield(pynames, name))(pv))
        compare_states!(key_for("with_element:$name", xo), xo, y, s)
    end
    # with_apsides: ellipses only
    if ecc < 1
        rp = u(1.05, 6.0) * (body.radii[1] + body.radii[2]) / 2
        ra = rp * u(1.0001, 2.5)
        xo = almanac_try(() -> with_apsides(s.o, rp, ra))
        y = anise_state(() -> s.a.with_apoapsis_periapsis_km(ra, rp))
        compare_states!(key_for("with_apsides", xo), xo, y, s)
    end
    # two_body: whole milliseconds, so that both libraries get the same interval
    dt = ecc < 1 ? round(u(-1.5, 1.5) * 2pi * sqrt(abs(sma)^3 / body.gm), digits=3) : round(u(-5000, 5000), digits=3)
    xo = almanac_try(() -> two_body(s.o, et + dt))
    y = anise_state(() -> s.a.at_epoch(anise_epoch(et + dt)))
    base = "two_body"
    if ecc > 1
        # R11: ANISE's mean anomaly is off by 2 pi j where raan + aop leaves (-pi, pi]
        wrap = s.o.raan + s.o.aop > pi
        base = wrap ? "two_body:hyp_wrap" : "two_body:hyp"
    end
    compare_states!(key_for(base, xo), xo, y, s; extra=" dt=$dt")
end

# ---------------------------------------------------------------------------------------------
# v-infinity
# ---------------------------------------------------------------------------------------------

function compare_vinf!(rng, s::Sample)
    s.o.ecc > 1 || return
    u(lo, hi) = lo + (hi - lo) * rand(rng)
    turn = u(10.0, 150.0)
    x = almanac_val(() -> vinf_periapsis(s.o, deg2rad(turn)))
    y = anise_val(() -> s.a.vinf_periapsis_km(turn))
    outcome!("vinf_periapsis", x, y, s; extra=" turn=$turn deg")
    rp = u(1.05, 10.0) * (s.body.radii[1] + s.body.radii[2]) / 2
    x = almanac_val(() -> vinf_turn_angle(s.o, rp))
    y = anise_val(() -> deg2rad(pyconvert(Float64, s.a.vinf_turn_angle_deg(rp))))
    outcome!("vinf_turn_angle", x, y, s; angle=true, extra=" rp=$rp")
end

# ---------------------------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------------------------

function run_all()
    rng = Xoshiro(SEED)
    samples = Sample[]
    for cls in CLASSES, i in 1:N
        body = iseven(i) ? EARTH : MARS
        s = try
            make_sample(rng, cls, body)
        catch e
            e isa Almanac.SpiceError ? nothing : rethrow()
        end
        s === nothing && continue
        push!(samples, s)
        compare_elements!(s)
        compare_local_frames!(s)
        compare_relative!(rng, s)
        compare_construction!(rng, s)
        compare_vinf!(rng, s)
    end
    return samples
end

println("ANISE ", ANISE_VERSION, ", ", N, " states per class, seed ", SEED, FULL ? " (exhaustive)" : "")
const SAMPLES = run_all()
for (k, n) in sort(collect(PYERRORS); by=last, rev=true)
    println("ANISE raised ", n, "x: ", k)
end
println("run time of the comparison: ", round(time() - T0, digits=1), " s")

# ---------------------------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------------------------

function known_label(key)
    P, D = KNOWN_DIFFERENT_IN_PIN.pin, KNOWN_DIFFERENT_IN_PIN.deliberate
    endswith(key, ":hypin") && haskey(P, "*:hypin") && return ("pin", P["*:hypin"])
    endswith(key, ":hyp_wrap:hypin") && return ("pin", P["*:hypin"])
    haskey(P, key) && return ("pin", P[key])
    endswith(key, ":hyp_wrap") && return ("deliberate", D["two_body:hyp_wrap"])
    haskey(D, key) && return ("deliberate", D[key])
    return nothing
end

println()
println("Notice: the quantities below marked 'not tested' are reported only. [pin] differ between")
println("1.0.0b0 and ANISE main by definition (empty the set when the pin moves); [deliberate] are")
println("deviations of Almanac.jl that stay. Not compared at all: TNW frame, b-plane (absent in 1.0.0b0).")
println()
@printf("%-34s %7s %9s %9s %11s\n", "quantity", "n", "undef", "mismatch", "largest")
for key in sort(collect(keys(ACCS)))
    a = ACCS[key]
    lab = known_label(key)
    @printf("%-34s %7d %9d %9d %11.3e %s\n", key, a.n, a.both_undef, a.mismatch, a.maxd,
        lab === nothing ? "" : "[$(lab[1]), not tested]")
    a.normalised && @printf("    in units of %s: %.3f at: %s\n",
        is_rate(key) ? "1e-12/(2 ms)" : "the angle model", a.maxn, a.worstn)
    a.range > 0 && println("    RANGE: ", a.range, "x, e.g. ", a.rexample)
    a.worst != "" && println("    largest at: ", a.worst)
    a.mexample != "" && println("    mismatch:   ", a.mexample)
end

# Tolerances, measured (seed 20261011; the default tier, 2000 states per class, and the exhaustive
# tier, 20000 per class), each set to a round figure above the largest measured value, with the
# reason of its size.
const TOL = Dict{String,Float64}(
    # Angles of near-circular orbits, in units of eps/ecc (see cond_units): C_ANGLE
    # Brouwer iteration: the worst state is a near-equatorial one (inc ~ 1e-7), where the
    # angles are ill-conditioned; the largest difference is 3.5e-10 (inc)
    "el:brouwer_sma" => 1e-8, "el:brouwer_ecc" => 1e-8, "el:brouwer_inc" => 1e-8,
    "el:brouwer_raan" => 1e-8, "el:brouwer_aop" => 1e-8, "el:brouwer_ma" => 1e-8,
    # Duration of ANISE has a resolution of 1 ns, measured 1.0e-9 s
    "el:light_time" => 1e-8,
    # two_body: the interval is whole milliseconds in both libraries, but the mean motion
    # of ANISE goes through degrees; the worst state is near-equatorial; measured 3.6e-9
    "two_body" => 1e-8, "two_body:hyp" => 1e-9,
)
# C_ANGLE: the angles in units of the model of cond_scale; largest measured 36 (default tier) and
# 66 (exhaustive tier, el:ma, whose Kepler-equation step has an amplifier beyond the model);
# the limit is 100. C_RATE: the rates in units of 1e-12/(2 ms), the error that the stopping
# rule of the Kepler iteration of ANISE (|E_n - E_n+1| < MA_EPSILON = 1e-12, the older iterate
# returned) gives to each of the two states of the difference quotient; largest measured 1.45
# (both tiers, ellipses) and 0.02 (outbound hyperbolas); the limit is 5.
const C_ANGLE = 100.0
const C_RATE = 5.0
function tol(key)
    haskey(TOL, key) && return TOL[key]
    haskey(COND_KEYS, key) && return C_ANGLE
    is_rate(key) && return C_RATE
    # Other quantities: the largest difference is 7.7e-12 (tlong), the angle sums and
    # equinoctial elements, which differ by the last bits of sin/cos in different libms;
    # the default is about one decade above (largest 8.4e-12 in the exhaustive tier, 2e-13 in the
    # default tier).
    return 1e-10
end

@testset "ANISE $(ANISE_VERSION)" begin
    for key in sort(collect(keys(ACCS)))
        known_label(key) === nothing || continue
        a = ACCS[key]
        @testset "$key" begin
            @test a.mismatch == 0
            @test a.range == 0
            @test a.maxn <= tol(key)
        end
    end
end

# ---------------------------------------------------------------------------------------------
# Pinned values: ANISE results for 20 states, for test/test_analysis_orbit.jl once the pin is a
# release that contains the definitions of main
# ---------------------------------------------------------------------------------------------

function pinned_values()
    out = String[]
    pick = Sample[]
    for cls in (:elliptic, :circular, :equatorial, :retrograde, :hyp_out)
        n = 0
        for s in SAMPLES
            s.cls === cls || continue
            push!(pick, s)
            (n += 1) == 4 && break
        end
    end
    push!(out, "# ANISE $(ANISE_VERSION), seed $SEED; elements in km, rad, s (angles converted from degrees)")
    push!(out, "# Use a tolerance of 1e-7 for :ta, :aop, :ma, :ea of the near-circular states (ecc >= 1e-9):")
    push!(out, "# they are good to ~4e-8 (1/ecc amplification), the other values to 1e-10.")
    push!(out, "const ANISE_PINNED = [")
    skip = Set([:fpa, :ha])   # fpa differs in the pin; ha is unbounded or errors for ellipses
    for s in pick
        vals = String[]
        for row in ELEMENT_TABLE
            row[1] in skip && continue
            y = anise_val(() -> anise_element(s.a, row))
            y === nothing && continue
            push!(vals, ":$(row[1]) => $(repr(y))")
        end
        for kind in (:ric, :vnc)
            d3 = try
                pymat(pygetattr(s.a, "dcm3x3_from_$(getfield(ANISE_LF, kind))_to_inertial")().rot_mat)
            catch
                nothing
            end
            d3 === nothing || push!(vals, ":$(kind)_R => $(repr(vec(d3)))")
        end
        push!(out, "    (body = $(s.body.id), et = $(repr(s.et)), state = $(repr(collect(s.o.state))),")
        push!(out, "     values = Dict($(join(vals, ", ")))),")
    end
    push!(out, "]")
    return join(out, "\n")
end

println()
println(pinned_values())
println()
println("total time: ", round(time() - T0, digits=1), " s")
