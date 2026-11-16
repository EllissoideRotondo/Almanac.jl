# The reference cases of the tests of ANISE (tests/astro/orbit.rs at b3aaadc), which
# ANISE checks against GMAT: a state, the GM, and the elements it must give.
# `GMAT_CASES` is a vector of `(name, state, gm, expected)`; `expected` maps an element
# name to its value in rad (or km, km^2/s^2 ...), converted with `deg2rad` from the
# degrees that ANISE prints, digits kept. Cases named `*_kep` use the Cartesian state
# that ANISE prints for its `Orbit::keplerian(...)` of the same test, which stands in
# for the element-based constructor (`keplerian_orbit` is a later task).
# Angles are listed in `GMAT_ANGLES`: they are compared modulo 2π, because ANISE reports
# a node of 360 degrees where CSPICE and Almanac give 0, and `:ea` in (-π, π].

const GMAT_GM = 398_600.441_5
const GMAT_ANGLES = (:inc, :raan, :aop, :ta, :tlong, :ea, :ma, :fpa)

# Cartesian state of the classical elements (km, rad), perifocal frame rotated by
# Rz(raan) Rx(inc) Rz(aop); used for the hyperbolic and the flight path angle cases,
# whose states ANISE builds with `Orbit::keplerian`.
function gmat_kep_state(sma, ecc, inc, raan, aop, ta, gm=GMAT_GM)
    p = sma * (1 - ecc^2)
    r = p / (1 + ecc * cos(ta))
    sv = sqrt(gm / p)
    pq = SVector(r * cos(ta), r * sin(ta), 0.0)
    vq = SVector(-sv * sin(ta), sv * (ecc + cos(ta)), 0.0)
    rz(a) = @SMatrix [cos(a) -sin(a) 0.0; sin(a) cos(a) 0.0; 0.0 0.0 1.0]
    rx(a) = @SMatrix [1.0 0.0 0.0; 0.0 cos(a) -sin(a); 0.0 sin(a) cos(a)]
    m = rz(raan) * rx(inc) * rz(aop)
    return vcat(m * pq, m * vq)
end

const GMAT_CASES = Any[]
_d(x) = deg2rad(x)

# val_state_def_circ_inc
push!(GMAT_CASES, (name="circ_inc", gm=GMAT_GM,
    state=SVector(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0),
    expected=Dict(:energy => -25.842_247_282_849_137, :hx => 35_065.806_679_607_005,
        :hy => 35_065.806_679_607_005, :hz => 24_796.292_541_9, :sma => 7_712.186_117_895_043,
        :ecc => 0.000_999_582_831_432_052_5, :inc => _d(63.434_003_407_751_14),
        :raan => _d(135.0), :aop => _d(90.0), :ta => _d(0.0), :tlong => _d(225.0),
        :ea => _d(0.0), :ma => _d(0.0), :apoapsis => 7_719.895_086_731_299,
        :periapsis => 7_704.477_149_058_786, :semi_parameter => 7_712.178_412_142_147,
        :period => 6_740.269_063_641)))
push!(GMAT_CASES, (name="circ_inc_kep", gm=GMAT_GM,
    state=SVector(8_057.976_452_202_976, -0.196_740_370_290_888_9, 1_475.383_214_274_138,
        -0.166_470_488_584_076_31, 6.913_868_638_275_646_5, 0.910_157_981_443_279_1),
    expected=Dict(:sma => 8_191.929_999_999_999, :ecc => 1.000_000_000_388_51e-06,
        :inc => _d(12.849_999_999_999_987), :raan => _d(306.614),
        :aop => _d(314.189_999_994_618_1), :ta => _d(99.887_700_005_381_9),
        :energy => -24.328_848_116_377_95, :period => 7_378.877_993_955,
        :hx => -10_200.784_799_426_574, :hy => -7_579.639_346_783_497,
        :hz => 55_711.757_929_384_25, :tlong => _d(0.691_700_000_000_082_6),
        :ea => _d(99.887_643_560_656_85), :ma => _d(99.887_587_115_926_96),
        :apoapsis => 8_191.938_191_930_002, :periapsis => 8_191.921_808_069_997)))

# val_state_def_elliptical
push!(GMAT_CASES, (name="elliptical", gm=GMAT_GM,
    state=SVector(5_946.673_548_288_958, 1_656.154_606_023_661, 2_259.012_129_598_249,
        -3.098_683_050_943_824, 4.579_534_132_135_011, 6.246_541_551_539_432),
    expected=Dict(:energy => -25.842_247_282_849_144, :period => 6_740.269_063_641,
        :hx => 0.015_409_898_034_704_383, :hy => -44_146.106_010_690_01,
        :hz => 32_364.892_694_481_765, :sma => 7_712.186_117_895_041,
        :ecc => 0.158_999_999_999_999_95, :inc => _d(53.753_69),
        :raan => _d(1.998_632_864_211_17e-05), :aop => _d(359.787_880_000_004),
        :ta => _d(25.434_003_407_751_188), :tlong => _d(25.221_903_394_083_824),
        :ea => _d(21.763_052_882_584_79), :ma => _d(18.385_336_330_516_39),
        :apoapsis => 8_938.423_710_640_353, :periapsis => 6_485.948_525_149_73,
        :semi_parameter => 7_517.214_340_648_537)))
push!(GMAT_CASES, (name="elliptical_kep", gm=GMAT_GM,
    state=SVector(8_087.161_618_048_522_5, -0.197_452_943_772_520_73, 1_480.726_901_246_883,
        -0.000_168_592_186_843_952_16, 6.886_845_792_370_852, 0.936_931_260_302_891_8),
    expected=Dict(:sma => 8_191.930_000_000_003, :ecc => 0.024_500_000_000_000_348,
        :inc => _d(12.850_000_000_000_016), :raan => _d(306.614),
        :aop => _d(314.190_000_000_000_4), :ta => _d(99.887_699_999_999_58),
        :energy => -24.328_848_116_377_94, :period => 7_378.877_993_955,
        :hx => -10_197.722_829_337_885, :hy => -7_577.364_166_057_776,
        :hz => 55_695.034_928_191_49, :tlong => _d(0.691_699_999_999_855_2),
        :ea => _d(98.501_748_370_880_22), :ma => _d(97.113_427_049_323_43),
        :apoapsis => 8_392.632_285_000_007, :periapsis => 7_991.227_715_000_001,
        :semi_parameter => 8_187.012_794_017_503)))

# val_state_def_circ_eq
push!(GMAT_CASES, (name="circ_eq", gm=GMAT_GM,
    state=SVector(-38_892.724_449_149_02, 16_830.384_772_891_86, 0.722_659_929_135_562_2,
        -1.218_008_333_846_6, -2.814_651_172_605_98, 1.140_294_223_185_661e-5),
    expected=Dict(:energy => -4.702_902_670_552_006, :period => 86_820.776_152_981,
        :hx => 2.225_951_522_241_969_5, :hy => -0.436_714_326_090_944_6,
        :hz => 129_969.001_391_865_75, :sma => 42_378.129_999_999_98,
        :ecc => 9.999_999_809_555_511e-9, :inc => _d(0.001_000_000_401_564_538_6),
        :raan => _d(78.9), :aop => _d(65.399_999_847_186_78),
        :ta => _d(12.300_000_152_813_197), :tlong => _d(156.599_999_999_999_97),
        :ea => _d(12.300_000_030_755_777), :ma => _d(12.299_999_908_698_359),
        :apoapsis => 42_378.130_423_781_27, :periapsis => 42_378.129_576_218_69,
        :semi_parameter => 42_378.129_999_999_976)))
push!(GMAT_CASES, (name="circ_eq_kep", gm=GMAT_GM,
    state=SVector(18_190.717_357_886_37, -118.107_162_539_218_69, 0.000_253_845_647_633_053_35,
        0.030_396_440_130_264_88, 4.680_909_107_924_576, 4.907_089_816_726_583e-8),
    expected=Dict(:sma => 18_191.098_000_000_013, :ecc => 9.999_999_997_416_087e-7,
        :inc => _d(1.207_418_269_725_733_3e-6), :raan => _d(306.543),
        :aop => _d(314.320_000_025_403_66), :ta => _d(98.764_999_974_596_28),
        :energy => -10.955_920_349_063_035, :period => 24_417.396_242_566,
        :hx => -0.001_194_024_028_558_358_7, :hy => -0.000_884_918_835_027_750_6,
        :hz => 85_152.684_597_507_06, :tlong => _d(359.627_999_999_999_93),
        :ea => _d(98.764_943_347_932_57), :ma => _d(98.764_886_721_264_56),
        :apoapsis => 18_191.116_191_098_008, :periapsis => 18_191.079_808_902_017)))

# val_state_def_equatorial
push!(GMAT_CASES, (name="equatorial", gm=GMAT_GM,
    state=SVector(-7273.338970882, 253.990592670, 0.022164861,
        -0.258285289, -7.396322922, -0.000645451),
    expected=Dict(:sma => 7278.136188379306, :ecc => 4.99846643158263e-05,
        :inc => _d(0.005000000478594339), :raan => _d(360.0),
        :aop => _d(177.9999736473912), :ta => _d(2.650_826_247_094_554e-5))))

# hyperbolic_inbound_true_anomaly_is_accepted: sma = -5000 km, ecc = 2 (nu_inf = 120 deg);
# the outbound state of Orbit::keplerian(sma, ecc, 10, 20, 30, 30) and its mirror image
# with the velocity reversed, whose true anomaly ANISE reports as 330 deg.
let (sma, ecc) = (-5000.0, 2.0), s = gmat_kep_state(-5000.0, 2.0, _d(10), _d(20), _d(30), _d(30))
    push!(GMAT_CASES, (name="hyperbolic_outbound", gm=GMAT_GM, state=s,
        expected=Dict(:sma => sma, :ecc => ecc, :inc => _d(10.0), :raan => _d(20.0),
            :aop => _d(30.0), :ta => _d(30.0))))
    push!(GMAT_CASES, (name="hyperbolic_inbound", gm=GMAT_GM,
        state=SVector(s[1], s[2], s[3], -s[4], -s[5], -s[6]),
        expected=Dict(:sma => sma, :ecc => ecc, :ta => _d(330.0))))
end

# val_fpa_matches_geometric_definition: fpa = asin(r.v / (|r||v|)) for sma = 20000 km,
# inc = raan = aop = 0.
for ecc in (0.1, 0.5, 0.9), ta in (0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0)
    s = gmat_kep_state(20_000.0, ecc, 0.0, 0.0, 0.0, _d(ta))
    r, v = s[1:3], s[4:6]
    push!(GMAT_CASES, (name="fpa_ecc$(ecc)_ta$(Int(ta))", gm=GMAT_GM, state=s,
        expected=Dict(:fpa => asin(dot(r, v) / (norm(r) * norm(v))))))
end
