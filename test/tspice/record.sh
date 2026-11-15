#!/usr/bin/env bash
# Record the calls that NAIF's test suite (tspice) makes to SPICELIB, for replay.jl.
#
#     test/tspice/record.sh OUTDIR [FAMILY ...]
#
# Clones PEPPER (https://github.com/jacobwilliams/pepper, which carries tspice and a
# SPICELIB interface on which tspice passes) at the commit below into OUTDIR/pepper,
# instruments that interface with instrument.py for the routines of routines.txt,
# builds it with fpm (gfortran 16, as PEPPER requires), and runs each family in its own
# directory with the log OUTDIR/logs/F_XXX.log. Nothing here is committed to Almanac.jl:
# the tspice sources stay in the clone, and the logs and the copies of the kernels that
# the families write stay in OUTDIR.
#
# Environment: PEPPER_COMMIT (default below), JOBS (parallel families, default nproc),
# TIMEOUT (seconds per family, default 1800).
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
shift
commit=${PEPPER_COMMIT:-454d6ba}
jobs=${JOBS:-$(nproc)}
timeout=${TIMEOUT:-1800}

if [ ! -d "$out/pepper/.git" ]; then
    git clone -q https://github.com/jacobwilliams/pepper.git "$out/pepper"
fi
git -C "$out/pepper" checkout -q "$commit"
git -C "$out/pepper" checkout -q -- src/legacy        # undo an earlier instrumentation
rm -f "$out/pepper/src/legacy/tspice_log.f90"
python3 "$here/instrument.py" "$out/pepper" "$here/routines.txt"
(cd "$out/pepper/tspice" && fpm build --profile release)
exe="$out/pepper/tspice/$(cd "$out/pepper/tspice" && fpm run --profile release --runner echo 2>/dev/null | tail -1)"

mkdir -p "$out/logs" "$out/work"
if [ $# -gt 0 ]; then
    families=$(printf '%s\n' "$@" | tr 'a-z' 'A-Z' | sed 's/^\(F_\)\{0,1\}/F_/')
else
    families=$("$exe" --list)
fi

run_family() {
    local fam=$1 wd="$out/work/$1"
    rm -rf "$wd" "$out/logs/$fam.log" "$out/logs/$fam.log.files"
    mkdir -p "$wd"
    (cd "$wd" && TSPICE_LOG="$out/logs/$fam.log" timeout "$timeout" "$exe" "$fam" > output.txt 2>&1) \
        && echo "$fam passed" || echo "$fam FAILED (see $wd/output.txt)"
}
export -f run_family
export out exe timeout
printf '%s\n' $families | xargs -P "$jobs" -I{} bash -c 'run_family {}'
