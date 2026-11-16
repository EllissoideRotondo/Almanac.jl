# Group 51: JSON kernels (`src/jsonkernel.jl`).
#
# A JSON kernel is a text kernel in the format of jsonspice, which PEPPER also reads.
# The tests write each kernel twice, as a text kernel and as JSON, and compare the
# pools that the two give; the rest checks the errors of the parser and that a JSON
# kernel loads, unloads and lists files as a text kernel does. No kernel file of NAIF
# is needed.

@testset "JSON kernels" begin
    dir = mktempdir()
    put(name, text) = (p = joinpath(dir, name); write(p, text); p)
    value(ctx, n) = try
        gdpool(ctx, n)
    catch
        gcpool(ctx, n)
    end
    pool(ctx) = Dict{String,Any}(n => value(ctx, n) for n in gnpool(ctx, "*"))

    @testset "the same pool as the text kernel" begin
        tk = put("a.tpc", """
            KPL/PCK
            \\begindata
            BODY399_RADII   = ( 6378.1366 6378.1366 6356.7519 )
            NAIF_BODY_NAME += 'SPOCK'
            NAIF_BODY_CODE += -1000
            EPOCH           = @2026-NOV-11/12:00
            FLAGS           = ( 1 0 )
            WORDS           = ( 'ab' 'it''s' )
            ONE             = 7
            \\begintext
            """)
        js = put("a.json", """
            \ufeff// a JSON kernel
            {
                "BODY399_RADII": [6378.1366, 6378.1366, 6356.7519],
                "+NAIF_BODY_NAME": "SPOCK",
                "+NAIF_BODY_CODE": -1000,
                "EPOCH": "@2026-NOV-11/12:00",   // a time
                "FLAGS": [true, false],
                "WORDS": ["ab", "it's"],
                "ONE": 7,
            }
            """)
        a, b = SpiceContext(), SpiceContext()
        furnsh(a, tk)
        furnsh(b, js)
        @test pool(a) == pool(b)
        @test bodn2c(b, "SPOCK") == -1000
        @test kinfo(b, js)[1] == "TEXT"
        unload(b, js)
        @test isempty(gnpool(b, "*"))
    end

    @testset "append, order and ldpool" begin
        p = put("b.json", """{"X": 1, "+X": [2, 3], "Y": "a", "Y": "b", "+Z": 4}""")
        ctx = SpiceContext()
        ldpool(ctx, p)
        @test gdpool(ctx, "X") == [1.0, 2.0, 3.0]
        @test gcpool(ctx, "Y") == ["b"]
        @test gdpool(ctx, "Z") == [4.0]
    end

    @testset "escapes and numbers" begin
        bs = '\\'
        p = put("c.json", "{\"S\": \"a$(bs)\"b$(bs)$(bs)c$(bs)/d$(bs)u00e9$(bs)ud83d$(bs)ude00\", " *
                          "\"N\": [-0.5e-3, 1E2, 0, 12.25]}")
        ctx = SpiceContext()
        furnsh(ctx, p)
        @test gcpool(ctx, "S") == ["a\"b\\c/dé😀"]
        @test gdpool(ctx, "N") == [-0.5e-3, 100.0, 0.0, 12.25]
    end

    @testset "a JSON meta-kernel" begin
        put("m1.json", """{"A": 1}""")
        put("m2.tpc", "\\begindata\nB = 2\n\\begintext\n")
        m = put("meta.json", """
            {"PATH_VALUES": ["$(escape_string(dir))"], "PATH_SYMBOLS": ["K"],
             "KERNELS_TO_LOAD": ["\$K/m1.json", "\$K/m2.tpc"]}
            """)
        ctx = SpiceContext()
        furnsh(ctx, m)
        @test gdpool(ctx, "A") == [1.0] && gdpool(ctx, "B") == [2.0]
        @test ktotal(ctx) == 3
        @test kinfo(ctx, m)[1] == "META"
        @test !expool(ctx, "KERNELS_TO_LOAD")
    end

    @testset "errors leave the context unchanged" begin
        ctx = SpiceContext()
        furnsh(ctx, put("ok.json", """{"KEEP": 1}"""))
        bad = [
            "[1]",                                  # not an object
            """{"A": 1""",                          # not closed
            """{"A": 1} x""",                       # text after the object
            """{"A" 1}""",                          # no colon
            """{"A": 1 "B": 2}""",                  # no comma
            """{"A": null}""",                      # null
            """{"A": {"B": 1}}""",                  # nested object
            """{"A": [[1]]}""",                     # nested array
            """{"A": []}""",                        # empty array
            """{"A": [1, "x"]}""",                  # mixed
            """{"A": ""}""",                        # empty string
            """{"A": "$("x"^81)"}""",               # longer than the pool keeps
            """{"": 1}""",                          # empty name
            """{"A B": 1}""",                       # blank in a name
            """{"$("N"^33)": 1}""",                 # name too long
            """{"A": 01}""",                        # leading zero
            """{"A": 1.}""",                        # no fraction digit
            """{"A": "\\x"}""",                     # bad escape
            """{"A": "@"}""",                       # no time
            """/* c */ {"A": 1}""",                 # block comment
            """{"A": 1, "B": 2,, }""",              # two commas
        ]
        for (k, text) in enumerate(bad)
            p = put("bad$k.json", text)
            err = try
                furnsh(ctx, p)
                nothing
            catch e
                e
            end
            @test err isa SpiceError && err.kind === :badkernel
            @test gnpool(ctx, "*") == ["KEEP"]
        end
        @test ktotal(ctx) == 1
    end

    @testset "the name decides, not the content" begin
        ctx = SpiceContext()
        @test_throws SpiceError furnsh(ctx, put("text.json", "\\begindata\nA = 1\n\\begintext\n"))
        furnsh(ctx, put("json.tpc", """{"A": 1}"""))     # no \\begindata: nothing loads
        @test isempty(gnpool(ctx, "*"))
        @test Almanac.is_json_kernel("X.JSON") && !Almanac.is_json_kernel("x.json.tpc")
    end
end
