# Group: static checks of the package code, with no kernel. Aqua: method ambiguities,
# unbound type parameters, undefined exports, stale dependencies, compat bounds, type
# piracy, persistent tasks. ExplicitImports: every imported name is explicit, owned by
# the module it comes from, and public. The extensions are checked with the package, as
# far as the AD backends are loaded; the public-name checks cover the names taken from
# Base and the dependencies, since an extension uses the internals of Almanac and of its
# backend by design.

import Aqua
import ExplicitImports as EI

group("static checks: Aqua and ExplicitImports") do
    @testset "Aqua" begin
        Aqua.test_all(Almanac)
    end
    @testset "ExplicitImports" begin
        @test EI.check_no_implicit_imports(Almanac) === nothing
        @test EI.check_no_stale_explicit_imports(Almanac) === nothing
        @test EI.check_all_explicit_imports_via_owners(Almanac) === nothing
        deps = (Base, Core, Almanac.StaticArrays, Almanac.LinearAlgebra, Almanac.Mmap)
        @test EI.check_all_explicit_imports_are_public(Almanac; from = deps) === nothing
        @test EI.check_all_qualified_accesses_via_owners(Almanac) === nothing
        @test EI.check_no_self_qualified_accesses(Almanac) === nothing
        # Names that Base does not declare public, each kept for a reason:
        # `broadcastable`, the documented interface of broadcasting; `RefValue`, the
        # concrete type of `Ref(x)` for a field without boxing; `setindex`, the copy of
        # a tuple with one element replaced; `Filesystem.open`, the only open with
        # O_EXCL (a file is created only if it does not exist). Julia 1.10 has no
        # `public` declaration, so names such as `Base.tail` count as internal there.
        @test VERSION < v"1.11" || EI.check_all_qualified_accesses_are_public(
            Almanac; from = deps,
            ignore = (:broadcastable, :RefValue, :setindex, :open)
        ) === nothing
    end
end
