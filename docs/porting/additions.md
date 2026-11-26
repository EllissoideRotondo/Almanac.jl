# Beyond CSPICE

These features have no counterpart in CSPICE:

- **Derivatives of any order** of states and rotations: [`spkder`](@ref),
  [`pxder`](@ref). See [Time derivatives](@ref).
- **Automatic differentiation** through ForwardDiff.jl ([dual-number
  epochs](@ref "Differentiate with ForwardDiff")) and through Enzyme.jl, Mooncake.jl, Zygote.jl,
  ReverseDiff.jl and Tracker.jl ([other AD packages](@ref "Differentiate with other AD packages")).
- **Queries resolved once**, for many epochs: [`state_query`](@ref),
  [`rotation_query`](@ref). See [Contexts and threads](@ref).
- **An aberration correction parsed once**: [`Aberration`](@ref). After
  `ab = Aberration("LT+S")`, write `spkezr("MOON", et, "J2000", ab, "EARTH")`.
- **Tables** for inner loops and GPUs: [`ephemeris_table`](@ref),
  [`rotation_table`](@ref). See [Ephemeris and rotation tables](@ref).
- **Bodies and frames defined by Julia functions**: [`add_ephemeris!`](@ref),
  [`add_frame!`](@ref). See [Custom bodies and frames](@ref).
- **Kernels built in memory**, with no file: [`load_samples!`](@ref),
  [`load_function!`](@ref), and [`spkops`](@ref), [`ckops`](@ref), [`pckops`](@ref) with
  [`furnsh`](@ref) of the handle. See [Kernels in memory](@ref).
- **Scratch files kept in memory**: [`dasops`](@ref), [`ekops`](@ref). CSPICE writes a
  temporary file instead.
- **Text kernels in JSON**: [`furnsh`](@ref) reads a file named `*.json` as a text kernel.

## Text kernels in JSON

A file named `*.json` is read as a text kernel written as one JSON object, in the
format of [jsonspice](https://github.com/jacobwilliams/jsonspice), which
[PEPPER](https://github.com/jacobwilliams/pepper) also reads. [`furnsh`](@ref),
[`ldpool`](@ref), [`unload`](@ref), [`kinfo`](@ref) and meta-kernels accept it;
CSPICE does not. The extension decides the format: a JSON file with another extension
is read as a classic text kernel, and loads nothing, because it has no `\begindata`.

```json
// earth.json
{
    "BODY399_RADII": [6378.1366, 6378.1366, 6356.7519],
    "+NAIF_BODY_NAME": "SPOCK",
    "+NAIF_BODY_CODE": -1000,
    "EPOCH": "@2026-NOV-11/12:00",
}
```

- Each member is one assignment. A key `"+NAME"` appends to the variable, like `+=`.
- Numbers, `true` (1) and `false` (0) are numeric values. Strings are string values,
  and a string that starts with `@` is a time.
- An array holds numbers and times, or strings.
- `//` line comments and trailing commas are allowed; any other JSON5 syntax is an
  error.

JSON numbers are read correctly rounded, while the numbers of a text kernel are read
as CSPICE reads them, which can be a few units in the last place off, so one value
written in the two formats can give pools that differ in the last bits.
