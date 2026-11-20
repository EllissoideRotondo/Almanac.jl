#!/usr/bin/env bash
# Test of the C logger (tspice_log.c): build test_log.c with stub hooks and f2c types,
# run it, and compare the log it writes with the expected text.
#
#     bash test/tspice/cspice/test_log.sh
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
t=$(mktemp -d)
trap 'rm -rf "$t"' EXIT

# the four f2c types of CSPICE on 64-bit Linux (the archive's f2c.h)
cat > "$t/f2c.h" <<'EOF'
typedef int integer;
typedef int logical;
typedef int ftnlen;
typedef double doublereal;
EOF

gcc -ansi -m64 -O2 -Wall -I"$t" -I"$here" -o "$t/test_log" "$here/tspice_log.c" "$here/test_log.c" -lm
log="$t/F_TEST.log"
TSPICE_LOG="$log" TSPICE_MAXPURE=1 TSPICE_MAXBIN=1000000000 "$t/test_log"

# the binary file of the long array: FNV-1a over the bytes of the doubles 1 to 2001
h=$(python3 -c '
import struct
h = 0xcbf29ce484222325
for b in b"".join(struct.pack("<d", k) for k in range(1, 2002)):
    h = ((h ^ b) * 0x100000001b3) & 0xffffffffffffffff
print("%016X" % h)')

expected="@ VNORM 1 F
<v1 d 3
 1.00000000000000000E+000
 2.00000000000000000E+000
 3.00000000000000000E+000
>result d 1
 1.00000000000000000E+000
! F
@ BODN2C 1 F
<name c 1
|EARTH
>code i 1
399
>found l 1
T
! F
@ CLPOOL 1 F
! F
@ BODN2C 1 F
<name c 1
|EARTH
>code i 1
399
>found l 1
T
! F
@ TKVRSN 1 F
<item c 1
%410a42
B $log.files/bin/$h.d.gz
<x d -1
! F"
# the end of a call that did not fail: "! F " with a trailing blank, as in Fortran
expected=$(echo "$expected" | sed 's/^! F$/! F /')
if [ "$(cat "$log")" = "$expected" ]; then
    echo "ok   log"
else
    echo "FAIL log"
    diff <(echo "$expected") "$log" || true
    exit 1
fi
[ -f "$log.files/bin/$h.d.gz" ] && echo "ok   binary file" || { echo "FAIL binary file"; ls -R "$log.files"; exit 1; }
cmp <(gzip -dc "$log.files/bin/$h.d.gz") <(python3 -c '
import struct, sys
sys.stdout.buffer.write(b"".join(struct.pack("<d", k) for k in range(1, 2002)))') \
    && echo "ok   binary values" || { echo "FAIL binary values"; exit 1; }
