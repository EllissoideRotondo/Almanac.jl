"""Instrument PEPPER's SPICELIB interface (src/legacy) so that NAIF's test suite
(tspice) records every call it makes: the routine, its inputs, its outputs and the
error it signals. The record is replayed against Almanac.jl by replay.jl.

    python instrument.py PEPPER_DIR ROUTINES_FILE

PEPPER_DIR is a copy of the PEPPER repository; its src/legacy files are rewritten in
place and src/legacy/tspice_log.f90 is added. ROUTINES_FILE lists the routines to
record, one name per line (`#` starts a comment); replay.jl writes it from the routines
it knows. A routine is skipped, with a note on standard error, when its wrapper or the
intents of its arguments cannot be read, or when an array argument has no length that
the log can use.

Nothing NAIF wrote is copied by this script; it changes PEPPER's own wrappers.
"""

import os
import re
import sys

DECL = re.compile(
    r"^\s*(?P<type>real\(wp\)|integer\(ip\)|integer|logical|character\s*\(\s*[^)]*\)|double precision)"
    r"\s*(?:,\s*(?P<attrs>[^:]*))?::\s*(?P<rest>.*)$", re.I)
HEAD = re.compile(
    r"^\s*(?:recursive\s+|pure\s+|elemental\s+|module\s+)*"
    r"(?:(?P<ftype>real\(wp\)|integer\(ip\)|integer|logical|character\s*\(\s*[^)]*\))\s+)?"
    r"(?P<kind>subroutine|function)\s+(?P<name>[a-z0-9_]+)\s*\((?P<args>[^)]*)\)?", re.I)

# Routines whose file arguments are copied when the call is recorded (the files of a
# family change and disappear while it runs): the snapshot of the working directory is
# taken before the call.
FILE_ROUTINES = {"furnsh", "ldpool", "spkcov", "spkobj", "ckcov", "ckobj", "pckcov",
                 "dskobj", "dsksrf", "getfat", "spklef", "pcklof", "cklpf", "eklef"}


# Lengths that the arguments do not give: SPICELIB's documented sizes.
FIXED_LENGTHS = {
    "oscltx": {"elts": "20"},
    "limbpt": {"npts": "max(0, ncuts)", "points": "3*sum(npts(1:max(0, ncuts)))",
               "epochs": "sum(npts(1:max(0, ncuts)))", "tangts": "3*sum(npts(1:max(0, ncuts)))"},
    "termpt": {"npts": "max(0, ncuts)", "points": "3*sum(npts(1:max(0, ncuts)))",
               "epochs": "sum(npts(1:max(0, ncuts)))", "trmvcs": "3*sum(npts(1:max(0, ncuts)))"},
    # N counts characters (DASADC, DASUDC) or records (SPKW10, CKW05), not array elements
    "dasadc": {"data": "merge((n + epos - bpos) / max(1, epos - bpos + 1), 0, epos >= bpos .and. n > 0)"},
    "dasudc": {"data": "merge((last - first + 1 + epos - bpos) / max(1, epos - bpos + 1), 0, "
                        "epos >= bpos .and. last >= first)"},
    "spkw10": {"consts": "8", "elems": "max(0, 10*n)"},
    "ckw05": {"packts": "max(0, n*merge(8, merge(4, merge(14, 7, subtyp == 2), subtyp == 1), subtyp == 0))"},
}


def joined_lines(text):
    """Fortran lines with `&` continuations joined, and comments removed."""
    out, buf = [], ""
    for line in text.splitlines():
        code = line.split("!", 1)[0].rstrip()
        if buf:
            code = code.lstrip()
            if code.startswith("&"):
                code = code[1:]
        if code.endswith("&"):
            buf += code[:-1] + " "
            continue
        out.append(buf + code)
        buf = ""
    return out


def split_names(rest):
    """`a ( 3 ), b, c(n, *)` -> [('a', '3'), ('b', None), ('c', 'n, *')]."""
    out, depth, cur = [], 0, ""
    for ch in rest:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    res = []
    for item in out:
        item = item.strip()
        if not item or "=" in item:
            continue
        m = re.match(r"([a-z0-9_]+)\s*(?:\((.*)\))?$", item, re.I)
        if m:
            res.append((m.group(1).lower(), m.group(2)))
    return res


def library_intents(src):
    """{routine: {arg: intent}} from the declarations of PEPPER's modules."""
    table = {}
    for fn in sorted(os.listdir(src)):
        if not fn.endswith(".f90"):
            continue
        name, intents = None, {}
        for line in joined_lines(open(os.path.join(src, fn), encoding="latin-1").read()):
            h = HEAD.match(line)
            if h and not re.match(r"^\s*end\b", line, re.I):
                name, intents = h.group("name").lower(), {}
                continue
            if name and re.match(r"^\s*end\s+(subroutine|function)\b", line, re.I):
                if intents:
                    table.setdefault(name, []).append(intents)
                name = None
                continue
            if name:
                d = DECL.match(line)
                if d and d.group("attrs") and "intent" in d.group("attrs").lower():
                    it = re.search(r"intent\s*\(\s*(in\s*out|inout|in|out)\s*\)", d.group("attrs"), re.I)
                    for arg, _ in split_names(d.group("rest")):
                        intents[arg] = it.group(1).lower().replace(" ", "")
    return table


def kind_of(ftype):
    t = ftype.lower().replace(" ", "")
    if t.startswith("real") or t.startswith("double"):
        return "d"
    if t.startswith("integer"):
        return "i"
    if t.startswith("logical"):
        return "l"
    return "c"


def length_expr(dims, args, decls):
    """The number of elements of an array argument, as a Fortran expression, or None."""
    if dims is None:
        return "1"
    parts = [p.strip().lower() for p in dims.split(",")]
    for k, p in enumerate(parts):                    # explicit bounds `lo:hi`
        if ":" in p:
            lo, hi = (x.strip() for x in p.split(":", 1))
            if hi == "*" or not re.match(r"^-?\d+$", lo):    # cells and windows: not recorded
                return None
            parts[k] = f"({hi})-({lo})+1"
    if parts[-1] == "*":
        # SPICELIB's convention: the count of an assumed-size array is the argument N
        # (or NVALS, NV, NSIZE), which is an input for an input array and an output
        # for an output array.
        for cand in ("n", "dim", "nvals", "nv", "nsize", "size", "nvars", "npts"):
            if cand in args and decls.get(cand, ("", None)) == ("i", None):
                lead = "*".join(f"({p})" for p in parts[:-1]) or "1"
                return f"max(0, ({lead})*({cand}))"
        return None
    return "*".join(f"({p})" for p in parts)


def instrument_file(path, wanted, intents, done):
    text = open(path, encoding="latin-1").read()
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        h = HEAD.match(line.split("!", 1)[0])
        if not h or re.match(r"^\s*end\b", line, re.I) or h.group("name").lower() not in wanted:
            out.append(line)
            i += 1
            continue
        name = h.group("name").lower()
        # The whole routine, to its end line.
        j = i
        while not re.match(rf"^\s*end\s+(subroutine|function)\s+{name}\b", lines[j], re.I):
            j += 1
        body = lines[i:j + 1]
        res = instrument_routine(name, h, body, intents.get(name, []))
        if isinstance(res, str):
            print(f"skip {name}: {res}", file=sys.stderr)
            out.extend(body)
        else:
            out.extend(res)
            done.add(name)
        i = j + 1
    open(path, "w", encoding="latin-1").write("\n".join(out))


def instrument_routine(name, h, body, cands):
    code = joined_lines("\n".join(body))
    head = HEAD.match(code[0])
    args = [a.strip().lower() for a in (head.group("args") or "").split(",") if a.strip()]
    decls = {}
    for line in code[1:]:
        d = DECL.match(line)
        if d:
            for arg, dims in split_names(d.group("rest")):
                decls[arg] = (kind_of(d.group("type")), dims)
    if any(a not in decls for a in args):
        return "an argument has no declaration in the wrapper"
    intent = {} if not args else None
    for c in cands:
        if all(a in c for a in args):
            intent = c
            break
    if intent is None:
        return "no library procedure declares the intents of all its arguments"
    if any(re.match(r"^\s*(return|go\s*to|entry)\b", l, re.I) for l in code[1:]):
        return "the wrapper has a RETURN, GO TO or ENTRY"
    lens = {}
    fixed = FIXED_LENGTHS.get(name, {})
    for a in args:
        if a in fixed:
            lens[a] = fixed[a]
            continue
        k, dims = decls[a]
        n = length_expr(dims, args, decls)
        if n is None:
            return f"the length of argument {a} ({dims}) is not known"
        lens[a] = n
    # Insert after the last declaration (the first executable line follows).
    last_decl = 0
    for k, line in enumerate(body):
        s = line.split("!", 1)[0].strip().lower()
        if k == 0 or not s:
            continue
        if ("::" in s or s.startswith("use ") or s.startswith("implicit") or
                s.startswith("parameter") or s.startswith("&") or body[k - 1].rstrip().endswith("&")):
            last_decl = k
    pre = ["      call tlog_begin('" + name.upper() + "')"]
    if name in FILE_ROUTINES:
        pre.append("      call tlog_snapshot()")
    post = []
    for a in args:
        k, _ = decls[a]
        it = intent[a]
        n = lens[a]
        rec = f"tlog_{k}('{a}', {a}, {n})" if n != "1" else f"tlog_{k}1('{a}', {a})"
        if it in ("in", "inout"):
            pre.append(f"      call {rec.replace(chr(39) + a + chr(39), chr(39) + '<' + a + chr(39), 1)}")
        if it in ("out", "inout"):
            post.append(f"      call {rec.replace(chr(39) + a + chr(39), chr(39) + '>' + a + chr(39), 1)}")
    if head.group("kind").lower() == "function":
        fk = kind_of(head.group("ftype") or "real(wp)")
        post.append(f"      call tlog_{fk}1('>result', {name})")
    post.append("      call tlog_end()")
    # `use tspice_log` goes after the first `use` line.
    new = list(body)
    for k, line in enumerate(new):
        if re.match(r"^\s*(use|implicit)\b", line, re.I):
            new.insert(k, "      use tspice_log")
            last_decl += 1
            break
    new = new[:last_decl + 1] + pre + new[last_decl + 1:-1] + post + [new[-1]]
    return new


LOGMOD = r'''!> Call log for the tspice record (written by Almanac.jl test/tspice/instrument.py).
!> The log goes to the file named by the environment variable TSPICE_LOG; without it
!> nothing is written. Format, one item per line:
!>   @ NAME DEPTH FAILED   a call (DEPTH > 1: made by a routine that is being logged)
!>   <ARG K N / >ARG K N   an input / output argument of kind K (d, i, l, c), N values
!>                         (N = -1: more than 2000, not logged),
!>                         one per following line (strings: their text, trailing blanks
!>                         removed, after a `|`)
!>   S DIR                 the files of the working directory, copied to DIR
!>   ! FAILED SHORT        the end of the call: the error status and short message
module tspice_log
   use pepper_kinds, only: ip, wp
   use pepper_instance_module, only: pepper_context
   use pepper_legacy_context, only: legacy_context
   implicit none
   private
   public :: tlog_begin, tlog_end, tlog_snapshot, tlog_d, tlog_d1, tlog_i, tlog_i1, &
      tlog_l, tlog_l1, tlog_c, tlog_c1
   integer, save :: unit = 0
   logical, save :: opened = .false.
   logical, save :: checked = .false.
   integer, save :: depth = 0
   integer, save :: nsnap = 0
   integer, parameter :: maxlog = 2000   ! longer arrays are logged as N = -1, not replayed
   character(1024), save :: logpath = ' '
contains
   logical function on()
      integer :: n, stat
      if (.not. checked) then
         checked = .true.
         call get_environment_variable('TSPICE_LOG', logpath, n, stat)
         if (stat == 0 .and. n > 0) then
            open(newunit=unit, file=trim(logpath), status='replace', action='write')
            opened = .true.
         end if
      end if
      on = opened
   end function on

   subroutine tlog_begin(name)
      character(*), intent(in) :: name
      type(pepper_context), pointer :: ctx
      if (.not. on()) return
      depth = depth + 1
      ctx => legacy_context()
      write(unit, '(a,1x,a,1x,i0,1x,l1)') '@', name, depth, ctx%failed()
   end subroutine tlog_begin

   subroutine tlog_end()
      type(pepper_context), pointer :: ctx
      character(64) :: msg
      if (.not. on()) return
      ctx => legacy_context()
      msg = ' '
      if (ctx%failed()) call ctx%getmsg('SHORT', msg)
      write(unit, '(a,1x,l1,1x,a)') '!', ctx%failed(), trim(msg)
      flush(unit)
      depth = depth - 1
   end subroutine tlog_end

   subroutine tlog_snapshot()
      character(1100) :: dir
      if (.not. on()) return
      nsnap = nsnap + 1
      write(dir, '(a,a,i0)') trim(logpath), '.files/', nsnap
      call execute_command_line('mkdir -p "' // trim(dir) // '" && find . -maxdepth 1 -type f ' // &
         '-exec cp -l {} "' // trim(dir) // '/" \; 2>/dev/null')
      write(unit, '(a,1x,a)') 'S', trim(dir)
   end subroutine tlog_snapshot

   subroutine tlog_d(label, x, n)
      character(*), intent(in) :: label
      real(wp), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      integer :: k
      if (.not. on()) return
      if (n > maxlog) then
         write(unit, '(a,1x,a,1x,i0)') label, 'd', -1
         return
      end if
      write(unit, '(a,1x,a,1x,i0)') label, 'd', n
      do k = 1, n
         write(unit, '(es25.17e3)') x(k)
      end do
   end subroutine tlog_d
   subroutine tlog_d1(label, x)
      character(*), intent(in) :: label
      real(wp), intent(in) :: x
      real(wp) :: v(1)
      v(1) = x
      call tlog_d(label, v, 1_ip)
   end subroutine tlog_d1

   subroutine tlog_i(label, x, n)
      character(*), intent(in) :: label
      integer(ip), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      integer :: k
      if (.not. on()) return
      if (n > maxlog) then
         write(unit, '(a,1x,a,1x,i0)') label, 'i', -1
         return
      end if
      write(unit, '(a,1x,a,1x,i0)') label, 'i', n
      do k = 1, n
         write(unit, '(i0)') x(k)
      end do
   end subroutine tlog_i
   subroutine tlog_i1(label, x)
      character(*), intent(in) :: label
      integer(ip), intent(in) :: x
      integer(ip) :: v(1)
      v(1) = x
      call tlog_i(label, v, 1_ip)
   end subroutine tlog_i1

   subroutine tlog_l(label, x, n)
      character(*), intent(in) :: label
      logical, intent(in) :: x(*)
      integer(ip), intent(in) :: n
      integer :: k
      if (.not. on()) return
      if (n > maxlog) then
         write(unit, '(a,1x,a,1x,i0)') label, 'l', -1
         return
      end if
      write(unit, '(a,1x,a,1x,i0)') label, 'l', n
      do k = 1, n
         write(unit, '(l1)') x(k)
      end do
   end subroutine tlog_l
   subroutine tlog_l1(label, x)
      character(*), intent(in) :: label
      logical, intent(in) :: x
      logical :: v(1)
      v(1) = x
      call tlog_l(label, v, 1_ip)
   end subroutine tlog_l1

   subroutine tlog_c(label, x, n)
      character(*), intent(in) :: label
      character(*), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      integer :: k
      if (.not. on()) return
      if (n > maxlog) then
         write(unit, '(a,1x,a,1x,i0)') label, 'c', -1
         return
      end if
      write(unit, '(a,1x,a,1x,i0)') label, 'c', n
      do k = 1, n
         write(unit, '(a,a)') '|', trim(x(k))
      end do
   end subroutine tlog_c
   subroutine tlog_c1(label, x)
      character(*), intent(in) :: label
      character(*), intent(in) :: x
      if (.not. on()) return
      write(unit, '(a,1x,a,1x,i0)') label, 'c', 1
      write(unit, '(a,a)') '|', trim(x)
   end subroutine tlog_c1
end module tspice_log
'''


def main():
    pepper, routines = sys.argv[1], sys.argv[2]
    wanted = set()
    for line in open(routines):
        line = line.split("#", 1)[0].strip().lower()
        if line:
            wanted.add(line)
    src = os.path.join(pepper, "src")
    intents = library_intents(src)
    legacy = os.path.join(src, "legacy")
    done = set()
    for fn in sorted(os.listdir(legacy)):
        if fn.endswith(".f90"):
            instrument_file(os.path.join(legacy, fn), wanted, intents, done)
    open(os.path.join(legacy, "tspice_log.f90"), "w").write(LOGMOD)
    missing = sorted(wanted - done)
    print(f"instrumented {len(done)} routines; not instrumented: {len(missing)}", file=sys.stderr)
    if missing:
        print("  " + " ".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main()
