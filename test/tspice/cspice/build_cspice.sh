#!/usr/bin/env bash
# Build the recording program of tspice_c, NAIF's test suite of CSPICE N0067, for
# record.sh (TSPICE_SOURCE=cspice).
#
#     test/tspice/cspice/build_cspice.sh OUTDIR
#
# Downloads NAIF's archive (below) into OUTDIR/tspice_c.tar unless it is there, checks
# its SHA-256, extracts it into OUTDIR/tspice_c, writes the wrappers and the driver with
# instrument_c.py into OUTDIR/gen, compiles CSPICE, its test utilities and the test
# families with NAIF's flags (gcc), redirects the calls of the test program to the wrappers
# (objcopy), and links OUTDIR/tspice_c_rec. The
# last line of the output is the path of the program. Nothing here is committed to
# Almanac.jl. Environment: JOBS (parallel compilations, default 1).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
url=https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/PC_Linux_64bit/tspice_c.tar
sha=373f3213cabab2717d21fc0110da703a713184d6cf73454372d0ae5a18c4adde
jobs=${JOBS:-1}
tar="$out/tspice_c.tar"
src="$out/tspice_c"

if ! echo "$sha  $tar" | sha256sum -c --quiet 2>/dev/null; then
    echo "downloading $url" >&2
    curl -fsSL --retry 3 -o "$tar.part" "$url"
    mv "$tar.part" "$tar"
    echo "$sha  $tar" | sha256sum -c --quiet || { echo "the archive does not have the SHA-256 $sha" >&2; exit 1; }
fi
if [ ! -f "$src/.extracted" ] || [ "$src/.extracted" -ot "$tar" ]; then
    rm -rf "$src"
    tar -xf "$tar" -C "$out"
    touch "$src/.extracted"
fi

python3 "$here/instrument_c.py" "$src" "$here/../routines.txt" "$out/gen" 2> "$out/gen.log" ||
    { cat "$out/gen.log" >&2; exit 1; }
grep -v '^mark ' "$out/gen.log" >&2 || true
cp "$here/tspice_log.c" "$here/tspice_log.h" "$out/gen/"

# NAIF's flags (mkprodct.csh of each directory); -w: NAIF's sources are not ours to fix
cflags="-c -ansi -m64 -O2 -DNON_UNIX_STDIO"
compile() {   # DIR NAME: the .c files of DIR into OUTDIR/obj/NAME, then OUTDIR/obj/NAME.a
    local dir=$1 name=$2 o="$out/obj/$2"
    mkdir -p "$o"
    find "$dir" -maxdepth 1 -name '*.c' -print0 |
        xargs -0 -P "$jobs" -I{} bash -c 'f={}; b=$(basename "$f" .c); [ "$0/$b.o" -nt "$f" ] ||
            gcc '"$cflags"' -w -I"$1" -I"$2" -o "$0/$b.o" "$f"' "$o" "$dir" "$src/include"
    rm -f "$out/obj/$name.a"
    ar crs "$out/obj/$name.a" "$o"/*.o
}
echo "compiling CSPICE and the test families with $jobs jobs" >&2
compile "$src/src/cspice" cspice
compile "$src/src/csupport" csupport
compile "$src/src/tutils_c" tutils_c
compile "$src/src/tspice_c" tspice_c
# the generated files, warnings on
mkdir -p "$out/obj/gen"
for f in driver tlog_wrap tspice_log; do
    inc="-I$out/gen -I$src/include"
    [ "$f" = driver ] && inc="$inc -I$src/src/tutils_c"
    gcc $cflags -Wall $inc -o "$out/obj/gen/$f.o" "$out/gen/$f.c"
done
# the calls of the test program (the families, NAIF's test and support utilities) go to
# the wrappers; CSPICE's calls between its own routines do not
for lib in tspice_c tutils_c csupport; do
    objcopy --redefine-syms="$out/gen/wrap.syms" "$out/obj/$lib.a" "$out/obj/$lib.wrapped.a"
done
gcc -m64 -o "$out/tspice_c_rec" "$out"/obj/gen/{driver,tlog_wrap,tspice_log}.o \
    -Wl,--start-group "$out"/obj/{tspice_c,tutils_c,csupport}.wrapped.a "$out/obj/cspice.a" -Wl,--end-group -lm
echo "$out/tspice_c_rec"
