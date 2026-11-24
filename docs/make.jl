# Build the manual: julia --project=docs docs/make.jl
using Documenter, Almanac

DocMeta.setdocmeta!(Almanac, :DocTestSetup, :(using Almanac); recursive = true)

makedocs(;
    sitename = "Almanac.jl",
    modules = [Almanac],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        collapselevel = 1,
        size_threshold = 600_000, size_threshold_warn = 300_000,
        search_size_threshold_warn = 800_000,  # 571 docstrings; about 150 KiB gzipped
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "tutorials/contexts.md",
            "tutorials/derivatives.md",
            "tutorials/tables.md",
            "tutorials/providers.md",
            "tutorials/writing.md",
            "tutorials/orbits.md",
            "tutorials/events.md",
        ],
        "Porting from CSPICE" => [
            "manual/conventions.md",
            "manual/coverage.md",
            "manual/additions.md",
            "manual/cspice_defects.md",
        ],
        "Validation" => [
            "manual/numerics.md",
            "manual/validation.md",
        ],
        "API reference" => [
            "reference/kernels.md",
            "reference/time.md",
            "reference/states.md",
            "reference/frames.md",
            "reference/geometry.md",
            "reference/events.md",
            "reference/orbits.md",
            "Mathematics" => [
                "reference/math/vectors.md",
                "reference/math/rotations.md",
                "reference/math/coordinates.md",
                "reference/math/ellipsoids.md",
                "reference/math/conics.md",
            ],
            "Files" => [
                "reference/files/daf.md",
                "reference/files/writers.md",
                "reference/files/das.md",
                "reference/files/ek.md",
            ],
        ],
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/Orbitology/Almanac.jl.git", devbranch = "main")
