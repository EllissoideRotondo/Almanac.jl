# A test SPK file of every SPK data type other than 2: `write_type_kernel` writes one
# segment per case with the segment writers of CSPICE and returns the cases as
# `TypeCase`s. Each case has its own target body, relative to the Earth (399).
#
# Used by the group `spktypes` and by `benchmark/ephemerides.jl`. The module that
# includes this file must have the names of `SPICE`, `LinearAlgebra` and `Random`.

const LIBCSPICE = SPICE.libcspice
const CInt = Cint                  # `SpiceInt`, and `integer`, `logical`, `ftnlen` of f2c

"Bound of the relative difference between the states of Almanac and CSPICE for a case."
const TYPE_RTOL = 5.0e-13

const EARTH = 399
const GM_EARTH = 398600.4418

"""
One generated segment. `epochs` are the special epochs (knots, record bounds) inside
`[t0, t1]`; `h` is the step of the finite difference of the velocity; `rtol` is the
bound of the state comparison with CSPICE.
"""
struct TypeCase
    label::String
    dtype::Int
    body::Int
    t0::Float64
    t1::Float64
    epochs::Vector{Float64}
    h::Float64
    rtol::Float64
end

"Raise the pending CSPICE error, if there is one."
cspice_check() = SPICE.handleerror()

# --- reference trajectories -------------------------------------------------------------

"State (6-vector) at `dt` seconds after `state0` on the two-body orbit about the Earth."
orbit_state(state0, dt) = collect(SPICE.prop2b(GM_EARTH, state0, dt))

"Two-body acceleration at the position `r`."
orbit_acc(r) = -GM_EARTH .* r ./ norm(r)^3

"Ellipse: perigee radius 7000 km, eccentricity 0.2, inclined; period about 2.3 hours."
const ELLIPSE0 = [7000.0, 0.0, 0.0, 0.0, 7.0 * cosd(40.0) * 1.1783, 7.0 * sind(40.0) * 1.1783]
"Hyperbola: perigee radius 8000 km, speed above the escape speed."
const HYPERBOLA0 = [8000.0, 0.0, 0.0, 0.0, 11.0, 4.0]

"6 x n matrix of the states of the orbit from `state0` at `epochs` (time zero is `epochs[1]`)."
function state_matrix(state0, epochs)
    m = Matrix{Float64}(undef, 6, length(epochs))
    for (j, t) in enumerate(epochs)
        m[:, j] = orbit_state(state0, t - epochs[1])
    end
    return m
end

"12 x n matrix of type 18 subtype 0 packets: position, velocity, velocity, acceleration."
function packets12(states)
    p = Matrix{Float64}(undef, 12, size(states, 2))
    for j in axes(states, 2)
        p[1:6, j] = states[:, j]
        p[7:9, j] = states[4:6, j]
        p[10:12, j] = orbit_acc(states[1:3, j])
    end
    return p
end

"`n` ascending epochs from `t0` with steps of `step` seconds that vary by ±30 percent."
function uneven_epochs(rng, n, t0, step)
    out = Vector{Float64}(undef, n)
    out[1] = t0
    for i in 2:n
        out[i] = out[i - 1] + step * (0.7 + 0.6 * rand(rng))
    end
    return out
end

# --- segment writers ------------------------------------------------------------------

function write_type5!(h, body, frame, states, epochs)
    ccall(
        (:spkw05_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, Ref{Cdouble}, Ref{Cdouble},
        ), h, body, EARTH, frame, epochs[1],
        epochs[end], "type 5", GM_EARTH, length(epochs), states, epochs
    )
    return cspice_check()
end

"Types 8 and 12: states at equal steps."
function write_equal!(h, dtype, body, frame, degree, states, epoch1, step)
    n = size(states, 2)
    last = epoch1 + (n - 1) * step
    types = (
        CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring, CInt, CInt, Ref{Cdouble},
        Cdouble, Cdouble,
    )
    if dtype == 8
        ccall(
            (:spkw08_c, LIBCSPICE), Cvoid, (
                CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble,
            ), h, body, EARTH, frame,
            epoch1, last, "type 8", degree, n, states, epoch1, step
        )
    else
        ccall(
            (:spkw12_c, LIBCSPICE), Cvoid, (
                CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble,
            ), h, body, EARTH, frame,
            epoch1, last, "type 12", degree, n, states, epoch1, step
        )
    end
    cspice_check()
    return last
end

"Types 9 and 13: states at arbitrary epochs."
function write_unequal!(h, dtype, body, frame, degree, states, epochs)
    n = length(epochs)
    if dtype == 9
        ccall(
            (:spkw09_c, LIBCSPICE), Cvoid, (
                CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble},
            ), h, body, EARTH, frame,
            epochs[1], epochs[end], "type 9", degree, n, states, epochs
        )
    else
        ccall(
            (:spkw13_c, LIBCSPICE), Cvoid, (
                CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble},
            ), h, body, EARTH, frame,
            epochs[1], epochs[end], "type 13", degree, n, states, epochs
        )
    end
    return cspice_check()
end

function write_type18!(h, subtyp, body, frame, degree, packets, epochs)
    ccall(
        (:spkw18_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
            Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble},
        ), h, subtyp, body, EARTH, frame,
        epochs[1], epochs[end], "type 18", degree, length(epochs), packets, epochs
    )
    return cspice_check()
end

"""
Type 19 through the Fortran-interface routine. `minis` is a vector of
`(subtype, degree, packets, epochs)`; `bounds` are the interval bounds.
"""
function write_type19!(h, body, frame, minis, bounds, later::Bool)
    @assert length(bounds) == length(minis) + 1
    nintvl = Ref{CInt}(length(minis))
    npkts = CInt[length(m[4]) for m in minis]
    subtps = CInt[m[1] for m in minis]
    degres = CInt[m[2] for m in minis]
    packts = reduce(vcat, [vec(m[3]) for m in minis])
    epochs = reduce(vcat, [m[4] for m in minis])
    seg = "type 19"
    ccall(
        (:spkw19_, LIBCSPICE), CInt, (
            Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
            Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{CInt}, Ref{CInt}, Ref{CInt},
            Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{CInt}, CInt, CInt,
        ),
        h, body, EARTH, frame, bounds[1], bounds[end], seg, nintvl, npkts, subtps, degres,
        packts, epochs, bounds, CInt(later), ncodeunits(frame), ncodeunits(seg)
    )
    return cspice_check()
end

"Types 1 and 21 through the Fortran-interface routines; `dlines` has one line per column."
function write_mda!(h, dtype, body, frame, first, dlines, epochs)
    seg = "type $dtype"
    n = Ref{CInt}(length(epochs))
    if dtype == 1
        @assert size(dlines, 1) == 71
        ccall(
            (:spkw01_, LIBCSPICE), CInt, (
                Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
                Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{Cdouble}, Ref{Cdouble},
                CInt, CInt,
            ), h, body, EARTH, frame, first, epochs[end], seg, n, dlines, epochs,
            ncodeunits(frame), ncodeunits(seg)
        )
    else
        ccall(
            (:spkw21_, LIBCSPICE), CInt, (
                Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
                Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{CInt}, Ref{Cdouble},
                Ref{Cdouble}, CInt, CInt,
            ), h, body, EARTH, frame, first, epochs[end], seg, n,
            size(dlines, 1), dlines, epochs, ncodeunits(frame), ncodeunits(seg)
        )
    end
    return cspice_check()
end

"""
Difference lines of order `m` (15 for type 1) with final epochs `epochs`: reference
epoch inside the interval of the line, step sizes of the size of the interval, a
reference state of an orbit and a random difference table that decreases with the
order, as the table of an integrator does.
"""
function difference_lines(rng, m, t0, epochs)
    dl = zeros(4m + 11, length(epochs))
    prev = t0
    for (i, tf) in enumerate(epochs)
        span = tf - prev
        dl[1, i] = prev + span * (0.3 + 0.4 * rand(rng))                 # TL
        dl[2:(m + 1), i] = span .* (0.5 .+ rand(rng, m))                     # G
        s = orbit_state(ELLIPSE0, dl[1, i] - t0)
        dl[(m + 2):2:(m + 6), i] = s[1:3]
        dl[(m + 3):2:(m + 7), i] = s[4:6]
        kq = rand(rng, (m - 5):(m - 1), 3)      # KQMAX1 = max(KQ) + 1 is not above m
        for c in 1:3, j in 1:m
            dl[m + 7 + (c - 1) * m + j, i] = 1.0e-3 * randn(rng) / 2.0^j             # DT(j, c)
        end
        dl[4m + 8, i] = maximum(kq) + 1                                    # KQMAX1
        dl[(4m + 9):(4m + 11), i] = kq
        prev = tf
    end
    return dl
end

"Chebyshev coefficients that decrease with the degree; `scale` is the size of term 0."
cheb_coefs(rng, ncoef, scale) = [scale * randn(rng) / 3.0^(k - 1) for k in 1:ncoef]

function write_type14!(rng, h, body, frame, degree, bounds)
    n = length(bounds) - 1
    ncoef = degree + 1
    recs = Matrix{Float64}(undef, 2 + 6 * ncoef, n)
    for i in 1:n
        recs[1, i] = (bounds[i] + bounds[i + 1]) / 2
        recs[2, i] = (bounds[i + 1] - bounds[i]) / 2
        recs[3:end, i] = vcat((cheb_coefs(rng, ncoef, c <= 3 ? 7000.0 : 7.0) for c in 1:6)...)
    end
    ccall(
        (:spk14b_c, LIBCSPICE), Cvoid, (
            CInt, Cstring, CInt, CInt, Cstring, Cdouble, Cdouble,
            CInt,
        ), h, "type 14", body, EARTH, frame, bounds[1], bounds[end], degree
    )
    cspice_check()
    # Two calls, to cover the addition of records to a segment that is open.
    half = n ÷ 2
    ccall(
        (:spk14a_c, LIBCSPICE), Cvoid, (CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, half,
        recs[:, 1:half], bounds[1:half]
    )
    cspice_check()
    ccall(
        (:spk14a_c, LIBCSPICE), Cvoid, (CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, n - half,
        recs[:, (half + 1):n], bounds[(half + 1):n]
    )
    cspice_check()
    ccall((:spk14e_c, LIBCSPICE), Cvoid, (CInt,), h)
    return cspice_check()
end

"Type 20: `n` records of `intlen` days from the Julian date `initjd + initfr`; returns the coverage."
function write_type20!(rng, h, body, frame, polydg, n, intlen, dscale, tscale, initjd, initfr)
    recs = Matrix{Float64}(undef, 3 * (polydg + 2), n)
    for i in 1:n, c in 1:3
        rows = (c - 1) * (polydg + 2) .+ (1:(polydg + 2))
        recs[rows, i] = vcat(cheb_coefs(rng, polydg + 1, 0.01), 0.5 * randn(rng))
    end
    first = (initjd - 2451545.0 + initfr) * 86400.0
    last = first + n * intlen * 86400.0
    ccall(
        (:spkw20_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble,
        ), h, body,
        EARTH, frame, first, last, "type 20", intlen, n, polydg, recs, dscale, tscale, initjd,
        initfr
    )
    cspice_check()
    return first, last
end

"Type 3: `n` records of `intlen` seconds from `t0`; the velocity series is not the derivative of the position series."
function write_type3!(rng, h, body, frame, degree, n, t0, intlen)
    @assert n >= 1 && intlen > 0
    ncoef = degree + 1
    recs = Matrix{Float64}(undef, 6 * ncoef, n)
    for i in 1:n
        recs[:, i] = vcat((cheb_coefs(rng, ncoef, c <= 3 ? 7000.0 : 7.0) for c in 1:6)...)
    end
    last = t0 + n * intlen
    ccall(
        (:spkw03_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, CInt, Ref{Cdouble}, Cdouble,
        ), h, body, EARTH, frame, t0, last, "type 3",
        intlen, n, degree, recs, t0
    )
    cspice_check()
    return last
end

function write_type15!(h, body, frame, t0, t1, epoch, tp, pa, p, ecc, j2flg, pv, j2, radius)
    ccall(
        (:spkw15_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
            Cdouble, Cdouble, Cdouble,
        ), h, body, EARTH, frame, t0, t1, "type 15", epoch, tp, pa,
        p, ecc, Float64(j2flg), pv, GM_EARTH, j2, radius
    )
    return cspice_check()
end

function write_type17!(h, body, frame, t0, t1, epoch, eqel, rapol, decpol)
    ccall(
        (:spkw17_c, LIBCSPICE), Cvoid, (
            CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, Ref{Cdouble}, Cdouble, Cdouble,
        ), h, body, EARTH, frame, t0, t1, "type 17",
        epoch, eqel, rapol, decpol
    )
    return cspice_check()
end

# --- the test kernel ------------------------------------------------------------------

"""
    write_type_kernel(path, rng) -> Vector{TypeCase}

Write the SPK file `path` with one segment per case and return the cases.
Postcondition: the file is closed and every case has `t0 < t1`.
"""
function write_type_kernel(path::AbstractString, rng::AbstractRNG)
    cases = TypeCase[]
    h = SPICE.spkopn(path, "Almanac test kernel", 0)
    body = -9000
    add!(label, dtype, t0, t1, epochs, hfd; rtol = TYPE_RTOL) = push!(
        cases,
        TypeCase(label, dtype, body, t0, t1, [t for t in epochs if t0 <= t <= t1], hfd, rtol)
    )
    J = "J2000"

    # Type 5: more than 100 states, so that the epoch directory is in use.
    ep = uneven_epochs(rng, 230, 1000.0, 60.0)
    body -= 1; write_type5!(h, body, J, state_matrix(ELLIPSE0, ep), ep)
    add!("type 5, ellipse", 5, ep[1], ep[end], ep, 1.0e-2)
    ep = uneven_epochs(rng, 40, -5000.0, 3000.0)
    body -= 1; write_type5!(h, body, J, state_matrix(HYPERBOLA0, ep), ep)
    add!("type 5, hyperbola", 5, ep[1], ep[end], ep, 1.0)
    ep = [0.0, 40000.0]
    body -= 1; write_type5!(h, body, J, state_matrix(ELLIPSE0, ep), ep)
    add!("type 5, two states, many revolutions", 5, ep[1], ep[end], ep, 1.0e-2)

    # Types 8 and 12: even and odd window sizes.
    for (dtype, degree) in ((8, 7), (8, 4), (8, 1), (12, 7), (12, 5), (12, 1))
        n, t0, step = 150, 2000.0, 60.0
        ep = [t0 + (i - 1) * step for i in 1:n]
        body -= 1
        t1 = write_equal!(h, dtype, body, J, degree, state_matrix(ELLIPSE0, ep), t0, step)
        # Knots, and the middles between knots, where an odd window moves.
        add!("type $dtype, degree $degree", dtype, t0, t1, vcat(ep, ep[1:(end - 1)] .+ step / 2), 1.0e-3)
    end

    # Types 9 and 13, with the largest degree and one segment in another frame.
    for (dtype, degree, frame) in (
            (9, 9, J), (9, 6, J), (9, 27, J), (9, 3, "ECLIPJ2000"),
            (13, 7, J), (13, 5, J), (13, 27, J), (13, 3, "ECLIPJ2000"),
        )
        ep = uneven_epochs(rng, 230, 500.0, 60.0)
        body -= 1; write_unequal!(h, dtype, body, frame, degree, state_matrix(ELLIPSE0, ep), ep)
        mids = (ep[1:(end - 1)] .+ ep[2:end]) ./ 2
        add!("type $dtype, degree $degree, $frame", dtype, ep[1], ep[end], vcat(ep, mids), 1.0e-3)
    end
    # The smallest segment that the writer accepts: as many states as the window.
    ep = uneven_epochs(rng, 4, 0.0, 60.0)
    body -= 1; write_unequal!(h, 9, body, J, 3, state_matrix(ELLIPSE0, ep), ep)
    add!("type 9, degree 3, 4 states", 9, ep[1], ep[end], ep, 1.0e-3)

    # Type 18: both subtypes, a small and the largest window.
    for (subtyp, degree) in ((0, 7), (0, 15), (1, 7), (1, 15))
        ep = uneven_epochs(rng, 230, 500.0, 60.0)
        st = state_matrix(ELLIPSE0, ep)
        body -= 1; write_type18!(h, subtyp, body, J, degree, subtyp == 0 ? packets12(st) : st, ep)
        add!("type 18, subtype $subtyp, degree $degree", 18, ep[1], ep[end], ep, 1.0e-3)
    end

    # Type 19: mini-segments of the three subtypes, with packets outside the interval
    # bounds (padding) and one mini-segment with a packet directory; both values of
    # the selection flag; and a segment with more than 100 intervals.
    for later in (true, false)
        minis = Any[]
        bounds = Float64[0.0]
        for (subtyp, degree, npkt) in ((0, 7, 30), (1, 7, 130), (2, 11, 25), (1, 3, 5), (2, 3, 2))
            ep = uneven_epochs(rng, npkt, bounds[end] - (npkt > 6 ? 150.0 : 0.0), 60.0)
            push!(bounds, npkt > 6 ? ep[end - 3] + 10.0 : ep[end])
            st = state_matrix(ELLIPSE0, ep)
            push!(minis, (subtyp, degree, subtyp == 0 ? packets12(st) : st, ep))
        end
        body -= 1; write_type19!(h, body, J, minis, bounds, later)
        knots = vcat(bounds, (m[4] for m in minis)...)
        add!(
            "type 19, mixed subtypes, select " * (later ? "later" : "earlier"), 19, bounds[1],
            bounds[end], knots, 1.0e-3
        )
    end
    minis = Any[]
    bounds = Float64[0.0]
    for i in 1:130
        ep = [bounds[end] + 50.0 * (j - 1) for j in 1:6]
        push!(bounds, ep[end])
        push!(minis, (isodd(i) ? 1 : 2, 3, state_matrix(ELLIPSE0, ep), ep))
    end
    body -= 1; write_type19!(h, body, J, minis, bounds, true)
    add!("type 19, 130 intervals", 19, bounds[1], bounds[end], bounds, 1.0e-3)

    # Types 1 and 21: more than 100 difference lines.
    for (dtype, m) in ((1, 15), (21, 15), (21, 25))
        t0 = 100.0
        ep = uneven_epochs(rng, 130, t0 + 80.0, 80.0)
        body -= 1; write_mda!(h, dtype, body, J, t0, difference_lines(rng, m, t0, ep), ep)
        add!("type $dtype, order $m", dtype, t0, ep[end], ep, 1.0e-3)
    end

    # Type 14: records of different lengths, more than 100 of them.
    bounds = uneven_epochs(rng, 131, -3000.0, 900.0)
    body -= 1; write_type14!(rng, h, body, J, 9, bounds)
    add!("type 14, degree 9", 14, bounds[1], bounds[end], bounds, 1.0e-2)
    bounds = uneven_epochs(rng, 4, 0.0, 900.0)
    body -= 1; write_type14!(rng, h, body, "ECLIPJ2000", 2, bounds)
    add!("type 14, degree 2, ECLIPJ2000", 14, bounds[1], bounds[end], bounds, 1.0e-2)

    # Type 20: scaled units (au and days) and a start that is a fraction of a day.
    for (polydg, n, intlen, dscale, tscale) in (
            (8, 40, 0.5, 149597870.7, 86400.0),
            (1, 3, 2.0, 1.0, 1.0), (12, 7, 1.25, 1000.0, 3600.0),
        )
        body -= 1
        t0, t1 = write_type20!(rng, h, body, J, polydg, n, intlen, dscale, tscale, 2451545.0, 0.25)
        add!("type 20, degree $polydg", 20, t0, t1, [t0 + k * intlen * 86400.0 for k in 0:n], 1.0)
    end

    # Type 15: the four values of the J2 flag, many revolutions on both sides of the
    # periapsis epoch, and a hyperbola (no J2 drift).
    tp = [sind(30.0) * sind(50.0), -cosd(30.0) * sind(50.0), cosd(50.0)]     # orbit normal
    pa = normalize(cross([0.0, 0.0, 1.0], tp))                               # in the plane
    pv = [0.0, 0.0, 1.0]
    for (j2flg, ecc, p) in (
            (0, 0.2, 9000.0), (1, 0.2, 9000.0), (2, 0.05, 8000.0),
            (3, 0.2, 9000.0), (0, 0.6, 12000.0), (0, 1.5, 20000.0), (0, 0.0, 7500.0),
        )
        epoch = 1234.5
        t0, t1 = ecc > 1 ? (epoch - 4.0e4, epoch + 6.0e4) : (epoch - 8.0e4, epoch + 1.2e5)
        body -= 1
        write_type15!(h, body, J, t0, t1, epoch, tp, pa, p, ecc, j2flg, pv, 1.0826e-3, 6378.14)
        add!("type 15, flag $j2flg, ecc $ecc", 15, t0, t1, [epoch], 1.0e-2)
    end

    # Type 17: an ellipse about a pole that is not the frame pole, and a near-circle.
    for (hh, kk, a) in ((0.12, -0.09, 12000.0), (1.0e-4, 2.0e-4, 42164.0), (0.6, 0.55, 30000.0))
        epoch = -777.0
        n = sqrt(GM_EARTH / a^3)
        eqel = [a, hh, kk, 0.7, 0.2, -0.15, 3.0e-7, n * (1 + 1.0e-3), -5.0e-7]
        body -= 1
        write_type17!(h, body, J, epoch - 1.0e5, epoch + 2.0e5, epoch, eqel, deg2rad(20.0), deg2rad(70.0))
        add!("type 17, h $hh, k $kk", 17, epoch - 1.0e5, epoch + 2.0e5, [epoch], 1.0e-2)
    end

    # Type 3: position and velocity series apart. It is the last case, so that the cases
    # before it keep their random data.
    t0 = -2000.0
    body -= 1
    t1 = write_type3!(rng, h, body, J, 9, 12, t0, 900.0)
    add!("type 3, degree 9", 3, t0, t1, [t0 + k * 900.0 for k in 0:12], 1.0e-2)

    SPICE.spkcls(h)
    @assert all(c -> c.t0 < c.t1, cases)
    return cases
end
