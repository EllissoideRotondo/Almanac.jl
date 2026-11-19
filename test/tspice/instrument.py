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
# Routines without state: their outputs depend on their inputs alone (no kernel, pool,
# file, time default or saved variable), so a repeat of a call adds nothing to the replay.
PURE = set("""
axisar azlrec b1900 b1950 cgv2el chbder chbigr chbint chbval clight conics convrt cyllat
cylrec cylsph dafps dafus dazldr dcyldr dgeodr diags2 dlatdr dnearp dp2hx dpr drdazl drdcyl
drdgeo drdlat drdsph dsphdr ducrss dvcrss dvdot dvhat dvnorm dvsep edlimb ednmpt edpnt
el2cgv eqncpv eul2m eul2xf georec gr2jul halfpi hrmesp hrmint hx2dp ident inedpl inelpl
inrypl invert invort invstm isrot j1900 j1950 j2000 j2100 jul2gr jyear latcyl latrec latsph
lgresp lgrind lgrint m2eul m2q mtxm mtxmg mtxv mtxvg mxm mxmg mxmt mxmtg mxv mxvg nearpt
npedln npelpt nplnpt nvc2pl nvp2pl oscelt oscltx pjelpl pl2nvc pl2nvp pl2psv pltar pltexp
pltnp pltnrm pltvol polyds prop2b psv2pl q2m qderiv qdq2av qxq radrec rav2xf raxisa recazl
reccyl recgeo reclat recrad recsph rotate rotmat rotvec rpd rquad saelgv sharpr sphcyl
sphlat sphrec spd stelab stlabx stmp03 surfnm surfpt surfpv twopi twovec twovxf tyear ucrss
unorm unormg vadd vaddg vcrss vdist vdistg vdot vdotg vequ vequg vhat vhatg vlcom vlcom3
vlcomg vminug vminus vnorm vnormg vpack vperp vprjp vprjpi vproj vprojg vrel vrelg vrotv
vscl vsclg vsep vsepg vsub vsubg vtmv vtmvg vupack vzero vzerog xf2eul xf2rav xpose xposeg
""".split())

# Queries: they read the loaded kernels, the pool or the time settings and change none
# of them (and keep no position, as the DAF and DLA searches do).
QUERY = set("""
azlcpo bodc2n bodc2s bodfnd bodn2c bods2c bodvcd bodvrd ccifrm cidfrm ckfrot ckfxfm ckgp
ckgpav cnmfrm deltet dtpool et2lst et2utc etcal frinfo frmnam gcpool gdpool getfov getfvn
gipool gnpool illumf illumg ilumin latsrf limbpt ltime namfrm occult pgrrec phaseq pxform
pxfrm2 recpgr scdecd sce2c sce2s sce2t scencd scfmt scs2e sct2e sctiks sctype sincpt
spkacs spkapo spkaps spkcpo spkcpt spkcvo spkcvt spkez spkezp spkezr spkgeo spkgps spkltc
spkpos spkssb srfc2s srfcss srfnrm srfs2c srfscc str2et subpnt subslr sxform tangpt termpt
tipbod tisbod tkfram trgsep unitim utc2et xfmsta timout tparse tpictr dpgrdr drdpgr
kinfo kdata ktotal
""".split())

FILE_ROUTINES = {"furnsh", "ldpool", "spkcov", "spkobj", "ckcov", "ckobj", "pckcov",
                 "dskobj", "dsksrf", "getfat", "spklef", "pcklof", "cklpf", "eklef"}
# Opens of an existing file: a snapshot of that file (its first argument) at a call from
# the test program, not from inside SPICELIB, gives the replay SPICELIB's file where
# Almanac's copy misses writes. Only the one file: a family such as F_ZZDSKBSR opens each
# of thousands of files in its directory.
OPEN_ROUTINES = {"dafopr", "dafopw", "dasopr", "dasopw", "ekopr", "ekopw", "spkopa", "dskopa"}
# Closes: the file that a close from the test program (not from inside SPICELIB) leaves,
# if it was open for write, is copied after the close, so that the replay can compare
# the file that Almanac wrote with SPICELIB's.
CLOSE_ROUTINES = {"dafcls", "dascls", "spkcls", "ckcls", "pckcls", "ekcls", "dskcls"}
# Routines that write to an open file and are not recorded (no replay, or the length of
# an array is not known): a call is logged with a `W` line and its handle only, so that
# the replay knows that its copy of the file misses the write. Without the mark the file
# that Almanac closes would differ from SPICELIB's for a reason that is not Almanac's.
WRITE_MARKS = set("""
ckw04a ckw04b ckw04e ckw06 spkw19 spk14a pckw20 pck03a pck03b pck03e
dafarr dafra dafwda dafwdr dafwfr dafwcr dafws dafrs
daswfr dasacr dasrcr dasufs dasurc dasurd dasuri daswbr daswrc daswrd daswri dascud dassdr dasacu
sgbwfs sgbwvs sgwes sgwfpk sgwvpk ekacli ekacld ekaclc ekifld ekffld
spks01 spks02 spks03 spks05 spks08 spks09 spks10 spks12 spks13 spks14 spks15 spks17
spks18 spks19 spks20 spks21
""".split())
# Arguments logged by a routine of tspice_log instead of as an array: the summary of
# DAFBNA packs its integers into doubles, which the decimal log cannot carry; it is
# logged unpacked, as `<dc` and `<ic`.
# None: the argument is not logged (a workspace, or the record pointers of an EK fast load,
# which the replay keeps as the state that Almanac's EKIFLD returns).
CUSTOM_ARGS = {
    "dafbna": {"sum": "call tlog_dafsum(handle, sum)"},
    "dafws": {"sum": "call tlog_dafsum_cur(sum)"},
    "dafrs": {"sum": "call tlog_dafsum_cur(sum)"},
    "dskmi2": {"work": None},
    "ekaclc": {"rcptrs": None, "wkindx": None},
    "ekacld": {"rcptrs": None, "wkindx": None},
    "ekacli": {"rcptrs": None, "wkindx": None},
    "ekffld": {"rcptrs": None},
}
# Calls after the call, before its outputs: the notes from which later calls take their
# lengths (the rows and column declarations of an EK fast load, the degree of SPK14B).
CUSTOM_POST = {
    "ekifld": "call tlog_ek_note(handle, segno, ncols, nrows, cnames, decls)",
    "spk14b": "call tlog_spk14_note(handle, chbdeg)",
}


# Lengths that the arguments do not give: SPICELIB's documented sizes.
FIXED_LENGTHS = {
    "oscltx": {"elts": "20"},
    "limbpt": {"npts": "max(0, ncuts)", "points": "3*sum(npts(1:max(0, ncuts)))",
               "epochs": "sum(npts(1:max(0, ncuts)))", "tangts": "3*sum(npts(1:max(0, ncuts)))"},
    "termpt": {"npts": "max(0, ncuts)", "points": "3*sum(npts(1:max(0, ncuts)))",
               "epochs": "sum(npts(1:max(0, ncuts)))", "trmvcs": "3*sum(npts(1:max(0, ncuts)))"},
    # N counts characters (DASADC, DASUDC) or records (SPKW10, CKW05), not array elements
    "dasadc": {"data": "tlog_nsub(n, bpos, epos, len(data))"},
    "dasudc": {"data": "tlog_nsub(last - first + 1, bpos, epos, len(data))"},
    "spkw10": {"consts": "8", "elems": "max(0, 10*n)"},
    "pckw02": {"cdata": "max(0, 3*(polydg+1)*n)"},
    "pckw20": {"cdata": "max(0, 3*(polydg+2)*n)"},
    "spkw18": {"packts": "max(0, n*merge(12, 6, subtyp == 0))"},
    # arrays of assumed size whose length another argument gives
    "dafps": {"dc": "max(0, nd)", "ic": "max(0, ni)", "sum": "max(0, nd + (ni + 1)/2)"},
    "dafus": {"sum": "max(0, nd + (ni + 1)/2)", "dc": "max(0, nd)", "ic": "max(0, ni)"},
    "dafgda": {"data": "max(0, eaddr - baddr + 1)"},
    "dafgsr": {"data": "max(0, end - begin + 1)"},
    "dasrdd": {"data": "max(0, last - first + 1)"},
    "dasrdi": {"data": "max(0, last - first + 1)"},
    "dasudd": {"data": "max(0, last - first + 1)"},
    "dasudi": {"data": "max(0, last - first + 1)"},
    "dasrdc": {"data": "tlog_nsub(last - first + 1, bpos, epos, len(data))"},
    "dlabfs": {"dladsc": "8"},
    "dlabbs": {"dladsc": "8"},
    "dlafns": {"dladsc": "8", "nxtdsc": "8"},
    "dlafps": {"dladsc": "8", "prvdsc": "8"},
    "dskgd": {"dladsc": "8", "dskdsc": "24"},
    "dskn02": {"dladsc": "8"},
    "dskv02": {"dladsc": "8"},
    "dskp02": {"dladsc": "8"},
    "dskz02": {"dladsc": "8"},
    "dski02": {"dladsc": "8"},
    "dskd02": {"dladsc": "8"},
    "dskb02": {"dladsc": "8"},
    "dskx02": {"dladsc": "8"},
    "dskxsi": {"srflst": "max(0, nsurf)", "dladsc": "8", "dskdsc": "24", "dc": "max(0, maxd)",
               "ic": "max(0, maxi)"},
    "dskxv": {"srflst": "max(0, nsurf)"},
    "spksfs": {"descr": "5"},
    "spkpds": {"descr": "5"},
    "spkuds": {"descr": "5"},
    "ckgr02": {"descr": "5"},
    "ckgr03": {"descr": "5"},
    "cknr02": {"descr": "5"},
    "cknr03": {"descr": "5"},
    "ckw01": {"sclkdp": "max(0, nrec)", "quats": "max(0, 4*nrec)", "avvs": "max(0, 3*nrec)"},
    "ckw02": {"start": "max(0, nrec)", "stop": "max(0, nrec)", "quats": "max(0, 4*nrec)",
              "avvs": "max(0, 3*nrec)", "rates": "max(0, nrec)"},
    "ckw03": {"sclkdp": "max(0, nrec)", "quats": "max(0, 4*nrec)", "avvs": "max(0, 3*nrec)",
              "starts": "max(0, nints)"},
    # the integer spatial index: SICGRD - 1 + MAXCGR + the voxel pointers, the voxel-plate
    # list, the vertex pointers (NV) and the vertex-plate list (dsk02.inc)
    "dskw02": {"corpar": "10", "spaixd": "10", "spaixi":
               "merge(7 + 100000 + spaixi(5) + spaixi(6) + max(0, nv) + spaixi(7), huge(0_ip), "
               "min(spaixi(5), spaixi(6), spaixi(7)) >= 0 .and. max(spaixi(5), spaixi(6), spaixi(7), nv) < 50000000)"},
    "ekbseg": {"cnames": "max(0, ncols)", "decls": "max(0, ncols)"},
    "chbint": {"cp": "max(0, degp + 1)"},
    "chbval": {"cp": "max(0, degp + 1)"},
    "chbigr": {"cp": "max(0, degp + 1)"},
    "chbder": {"cp": "max(0, degp + 1)", "dpdxs": "max(0, nderiv + 1)"},
    "polyds": {"coeffs": "max(0, deg + 1)", "p": "max(0, nderiv + 1)"},
    "xposeg": {"matrix": "max(0, nrow*ncol)", "xposem": "max(0, nrow*ncol)"},
    "swpool": {"names": "max(0, nnames)"},
    "scpart": {"pstart": "max(0, nparts)", "pstop": "max(0, nparts)"},
    "vdistg": {"v1": "max(0, ndim)", "v2": "max(0, ndim)"},
    "vdotg": {"v1": "max(0, ndim)", "v2": "max(0, ndim)"},
    "vrelg": {"v1": "max(0, ndim)", "v2": "max(0, ndim)"},
    "vsepg": {"v1": "max(0, ndim)", "v2": "max(0, ndim)"},
    "vzerog": {"v": "max(0, ndim)"},
    # mini-segments: the packet size of each subtype (spk19.inc, ck06.inc)
    "spkw19": {"npkts": "max(0, nintvl)", "subtps": "max(0, nintvl)", "degres": "max(0, nintvl)",
               "ivlbds": "merge(nintvl + 1, 0, nintvl > 0)", "epochs": "tlog_psum(npkts, subtps, nintvl, 1, 1, 1, 1)",
               "packts": "tlog_psum(npkts, subtps, nintvl, 12, 6, 6, 0)"},
    "ckw06": {"npkts": "max(0, nmini)", "subtps": "max(0, nmini)", "degres": "max(0, nmini)",
              "rates": "max(0, nmini)", "ivlbds": "merge(nmini + 1, 0, nmini > 0)",
              "sclkdp": "tlog_psum(npkts, subtps, nmini, 1, 1, 1, 1)",
              "packts": "tlog_psum(npkts, subtps, nmini, 8, 4, 14, 7)"},
    # the EK fast load: the rows and the declarations that EKIFLD noted
    "ekifld": {"cnames": "max(0, ncols)", "decls": "max(0, ncols)", "rcptrs": "max(0, nrows)"},
    "ekaclc": {"cvals": "tlog_ek_nvals(handle, segno, column, entszs)",
               "entszs": "tlog_ek_nent(handle, segno, column)", "nlflgs": "tlog_ek_nnul(handle, segno, column)"},
    "ekacld": {"dvals": "tlog_ek_nvals(handle, segno, column, entszs)",
               "entszs": "tlog_ek_nent(handle, segno, column)", "nlflgs": "tlog_ek_nnul(handle, segno, column)"},
    "ekacli": {"ivals": "tlog_ek_nvals(handle, segno, column, entszs)",
               "entszs": "tlog_ek_nent(handle, segno, column)", "nlflgs": "tlog_ek_nnul(handle, segno, column)"},
    # records of 2 + 6 (CHBDEG + 1) doubles, the degree of SPK14B
    "spk14a": {"coeffs": "tlog_spk14_len(handle, ncsets)", "epochs": "max(0, ncsets)"},
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


# Assumed-size arrays counted by an argument other than the first count of the routine:
# the plates by NP (DSKRB2 and DSKMI2 also take NV), and the coordinate parameters,
# which have the fixed size NSYPAR.
NAMED_COUNTS = {"plates": "np", "vrtces": "nv"}
FIXED_SIZES = {"corpar": "10"}


def length_expr(dims, args, decls, arg=None):
    """The number of elements of an array argument, as a Fortran expression, or None."""
    if dims is None:
        return "1"
    if arg in FIXED_SIZES and dims.strip() == "*":
        return FIXED_SIZES[arg]
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
        named = [NAMED_COUNTS[arg]] if arg in NAMED_COUNTS else []
        for cand in named + ["n", "dim", "nvals", "nv", "nsize", "size", "nvars", "npts"]:
            if cand in args and decls.get(cand, ("", None)) == ("i", None):
                lead = "*".join(f"({p})" for p in parts[:-1]) or "1"
                return f"max(0, ({lead})*({cand}))"
        return None
    return "*".join(f"({p})" for p in parts)


# The GF searches whose windows are logged (cell_args). The other GF entry points take
# callbacks (GFEVNT, GFFOVE, GFOCCE, GFUDB, GFUDS), which the log cannot carry.
CELL_ROUTINES = set("gfdist gfilum gfoclt gfpa gfposc gfrfov gfrr gfsep gfsntc gfsubc gftfov".split())


def cell_args(name, args, decls):
    """{arg: record} for the cell arguments of a routine of CELL_ROUTINES: a window
    `(lbcell:*)` is logged as its size, its cardinality and its elements (CELL(-1),
    CELL(0) and the CARD elements that follow, at most 10**6), `record` being the call
    with `{label}` for the label; a workspace of windows `(lbcell:mw, nw)` is not logged
    (record None)."""
    if name not in CELL_ROUTINES:
        return {}
    cells = {}
    for a in args:
        k, dims = decls[a]
        if dims is None or not dims.strip().lower().startswith("lbcell"):
            continue
        if dims.replace(" ", "").lower() != "lbcell:*":
            cells[a] = None
        elif k == "d":
            cells[a] = (f"tlog_d('{{label}}', {a}(-1), "
                        f"2_ip + nint(max(0.0_wp, min({a}(0), {a}(-1), 1.0e6_wp)), ip))")
        elif k == "i":
            cells[a] = (f"tlog_i('{{label}}', {a}(-1), "
                        f"2_ip + max(0_ip, min({a}(0), {a}(-1), 1000000_ip)))")
    return cells


def instrument_file(path, wanted, intents, done):
    text = open(path, encoding="latin-1").read()
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        h = HEAD.match(line.split("!", 1)[0])
        if not h or re.match(r"^\s*end\b", line, re.I) or h.group("name").lower() not in wanted | WRITE_MARKS:
            out.append(line)
            i += 1
            continue
        name = h.group("name").lower()
        # The whole routine, to its end line.
        j = i
        while not re.match(rf"^\s*end\s+(subroutine|function)\s+{name}\b", lines[j], re.I):
            j += 1
        body = lines[i:j + 1]
        res = (instrument_routine(name, h, body, intents.get(name, [])) if name in wanted
               else "not in the routines file")
        if isinstance(res, str) and name in WRITE_MARKS:
            print(f"mark {name}: {res}", file=sys.stderr)
            res = mark_routine(name, h, body)
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
    custom = CUSTOM_ARGS.get(name, {})
    cells = cell_args(name, args, decls)
    for a in args:
        if a in custom or a in cells:
            lens[a] = None
            continue
        if a in fixed:
            lens[a] = fixed[a]
            continue
        k, dims = decls[a]
        n = length_expr(dims, args, decls, a)
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
    elif name in OPEN_ROUTINES:
        pre.append(f"      call tlog_snapshot_file({args[0]})")
    elif name in CLOSE_ROUTINES:
        pre.append(f"      call tlog_close_name({args[0]})")
    post = ["      call tlog_mark()"]
    if name in CLOSE_ROUTINES:
        post.append("      call tlog_close_snapshot()")
    if name in CUSTOM_POST:
        post.append("      " + CUSTOM_POST[name])
    for a in args:
        k, _ = decls[a]
        it = intent[a]
        n = lens[a]
        if a in custom:
            custom[a] is None or pre.append("      " + custom[a])
            continue
        if a in cells:
            if cells[a] is not None and it in ("in", "inout"):
                pre.append("      call " + cells[a].format(label="<" + a))
            if cells[a] is not None and it in ("out", "inout"):
                post.append("      call " + cells[a].format(label=">" + a))
            continue
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


# Integer arguments of a marked writer that the replay maps to its files.
MARK_HANDLES = {"handle", "inhan", "han", "handl", "newh"}


def mark_routine(name, h, body):
    """The wrapper of a writer that is not recorded, with a mark: the call, a `W` line and
    the handles that it takes, with no other argument."""
    code = joined_lines("\n".join(body))
    head = HEAD.match(code[0])
    args = [a.strip().lower() for a in (head.group("args") or "").split(",") if a.strip()]
    decls = {}
    for line in code[1:]:
        d = DECL.match(line)
        if d:
            for arg, dims in split_names(d.group("rest")):
                decls[arg] = (kind_of(d.group("type")), dims)
    if any(re.match(r"^\s*(return|go\s*to|entry)\b", l, re.I) for l in code[1:]):
        return "the wrapper has a RETURN, GO TO or ENTRY"
    last_decl = 0
    for k, line in enumerate(body):
        s = line.split("!", 1)[0].strip().lower()
        if k == 0 or not s:
            continue
        if ("::" in s or s.startswith("use ") or s.startswith("implicit") or
                s.startswith("parameter") or s.startswith("&") or body[k - 1].rstrip().endswith("&")):
            last_decl = k
    pre = ["      call tlog_begin('" + name.upper() + "')", "      call tlog_writer()"]
    for a in args:
        if a in MARK_HANDLES and decls.get(a) == ("i", None):
            pre.append(f"      call tlog_i1('<{a}', {a})")
    post = ["      call tlog_mark()", "      call tlog_end()"]
    new = list(body)
    for k, line in enumerate(new):
        if re.match(r"^\s*(use|implicit)\b", line, re.I):
            new.insert(k, "      use tspice_log")
            last_decl += 1
            break
    return new[:last_decl + 1] + pre + new[last_decl + 1:-1] + post + [new[-1]]


LOGMOD = r'''!> Call log for the tspice record (written by Almanac.jl test/tspice/instrument.py).
!> The log goes to the file named by the environment variable TSPICE_LOG (a file or a
!> FIFO that record.sh compresses); without it nothing is written. Format, one item per
!> line:
!>   @ NAME DEPTH FAILED   a call (DEPTH > 1: made by a routine that is being logged)
!>   <ARG K N / >ARG K N   an input / output argument of kind K (d, i, l, c), N values
!>                         (N = -1: more than 2000, not logged),
!>                         one per following line (strings: their text, trailing blanks
!>                         removed, after a `|`)
!>   S DIR                 the files of the working directory, linked into DIR (a load),
!>                         or the file of an open (before it) or of a close (after it),
!>                         copied into DIR; for a close DIR has no file past the cap
!>   B PATH                the values of the next argument, of N = -1, of a call
!>                         from the test program, in the binary file PATH (gzip; doubles
!>                         or 4-byte integers, native order; strings of length L: `cL`);
!>                         an input or output, within TSPICE_MAXBIN bytes
!>   W                     a call of a writer that is not recorded: only its handles
!>                         follow (WRITE_MARKS of instrument.py)
!>   ! FAILED SHORT        the end of the call: the error status and short message
!> A call is kept in memory until it ends. A call of a routine without state (the list
!> PURE of instrument.py) is not written if its name, error status before the call and
!> inputs are those of an earlier call of the family; a call of a query (QUERY), which
!> reads the kernels or the settings but changes nothing, if they are those of an earlier
!> call with no call of another routine between them. Once the family has made
!> TSPICE_MAXPURE (default 20000) distinct calls of a routine without state, its other
!> calls are not written either (queries have no such limit). TSPICE_ALL=1 writes all.
module tspice_log
   use, intrinsic :: iso_fortran_env, only: int64
   use pepper_kinds, only: ip, wp
   use pepper_instance_module, only: pepper_context, pepper_zzddhnfo => zzddhnfo
   use pepper_legacy_context, only: legacy_context
   implicit none
   private
   public :: tlog_nsub, tlog_begin, tlog_end, tlog_mark, tlog_snapshot, tlog_snapshot_file, tlog_d, tlog_d1, &
      tlog_i, tlog_i1, tlog_l, tlog_l1, tlog_c, tlog_c1, tlog_close_name, tlog_close_snapshot, tlog_writer, &
      tlog_dafsum, tlog_dafsum_cur, tlog_psum, tlog_ek_note, tlog_ek_nvals, tlog_ek_nent, tlog_ek_nnul, &
      tlog_spk14_note, tlog_spk14_len
   integer, save :: unit = 0
   logical, save :: opened = .false.
   logical, save :: checked = .false.
   logical, save :: keepall = .false.
   integer, save :: depth = 0
   integer, save :: nsnap = 0
   integer, save :: nwritten = 0
   integer, parameter :: maxlog = 2000   ! longer arrays are logged as N = -1, not replayed
   character(1024), save :: logpath = ' '
   character(1024), save :: closename = ' '   ! the file of a close, if open for write
   logical, save :: bincheck = .false.
   integer(ip), save :: ekh(4000) = 0, eks(4000) = 0, eksz(4000) = 0, ekrows(4000) = 0
   logical, save :: eknul(4000) = .false.
   character(32), save :: ekcol(4000) = ' '
   integer, save :: nek = 0
   integer(ip), save :: s14h(64) = 0, s14d(64) = 0
   integer(int64), save :: binbytes = 0, maxbin = 268435456   ! binary files of long inputs
   character(32), save :: callname = ' '
   character(len=:), allocatable, save :: buf
   integer, save :: blen = 0, inend = 0
   integer(int64), allocatable, save :: seen(:)
   integer, save :: nseen = 0
   integer, save :: maxpure = 20000
   integer(int64), save :: state = 0           ! calls of routines that can change the state
   character(8), save :: capname(400) = ' '
   integer, save :: capcount(400) = 0, ncapped = 0
contains
   logical function on()
      integer :: n, stat
      character(8) :: all
      if (.not. checked) then
         checked = .true.
         call get_environment_variable('TSPICE_LOG', logpath, n, stat)
         if (stat == 0 .and. n > 0) then
            open(newunit=unit, file=trim(logpath), status='unknown', action='write', &
               access='stream', form='unformatted')
            opened = .true.
            allocate(character(len=1048576) :: buf)
            allocate(seen(0:65535)); seen = 0
            call get_environment_variable('TSPICE_ALL', all, n, stat)
            keepall = stat == 0 .and. all(1:1) == '1'
            call get_environment_variable('TSPICE_MAXPURE', all, n, stat)
            if (stat == 0 .and. n > 0) read(all, *) maxpure
         end if
      end if
      on = opened
   end function on

   subroutine emit(s)
      character(*), intent(in) :: s
      character(len=:), allocatable :: tmp
      integer :: need
      need = blen + len(s) + 1
      if (need > len(buf)) then
         allocate(character(len=max(2*len(buf), need)) :: tmp)
         tmp(1:blen) = buf(1:blen)
         call move_alloc(tmp, buf)
      end if
      buf(blen+1:blen+len(s)) = s
      buf(blen+len(s)+1:blen+len(s)+1) = achar(10)
      blen = need
   end subroutine emit

   subroutine tlog_begin(name)
      character(*), intent(in) :: name
      type(pepper_context), pointer :: ctx
      character(80) :: ln
      if (.not. on()) return
      depth = depth + 1
      if (depth == 1) then
         blen = 0
         inend = 0
         callname = name
      end if
      ctx => legacy_context()
      write(ln, '(a,1x,a,1x,i0,1x,l1)') '@', name, depth, ctx%failed()
      call emit(trim(ln))
   end subroutine tlog_begin

   !> The inputs of the call are logged: they end here.
   subroutine tlog_mark()
      if (.not. on()) return
      if (depth == 1) inend = blen
   end subroutine tlog_mark

   subroutine tlog_end()
      type(pepper_context), pointer :: ctx
      character(64) :: msg
      if (.not. on()) return
      ctx => legacy_context()
      msg = ' '
      if (ctx%failed()) call ctx%getmsg('SHORT', msg)
      call emit('! ' // merge('T', 'F', ctx%failed()) // ' ' // trim(msg))
      depth = depth - 1
      if (depth > 0) return
      if (.not. keepall .and. (is_pure(callname) .or. is_query(callname))) then
         if (inend == 0) inend = blen
         ! a query repeats only in the same state: since the last call of another routine
         if (.not. insert(ieor(fnv(buf(1:inend)), merge(state, 0_int64, is_query(callname))))) &
            return                                       ! a repeat
         ! a family that has made enough distinct calls of a routine without state
         if (is_pure(callname) .and. .not. below_cap(callname)) return
      else
         state = state + 1
      end if
      write(unit) buf(1:blen)
      nwritten = nwritten + 1
      if (mod(nwritten, 256) == 0) flush(unit)
   end subroutine tlog_end

   !> Count a new call of the routine without state NAME; false once the family has made
   !> maxpure distinct ones (TSPICE_MAXPURE, default 20000; 0: no limit).
   logical function below_cap(name)
      character(*), intent(in) :: name
      integer :: k
      do k = 1, ncapped
         if (capname(k) == name) exit
      end do
      if (k > ncapped) then
         ncapped = min(ncapped + 1, size(capname))
         k = ncapped
         capname(k) = name
         capcount(k) = 0
      end if
      capcount(k) = capcount(k) + 1
      below_cap = maxpure <= 0 .or. capcount(k) <= maxpure
   end function below_cap

   integer(int64) function fnv(s)
      character(*), intent(in) :: s
      integer :: k
      integer(int64), parameter :: prime = 1099511628211_int64
      fnv = -3750763034362895579_int64
      do k = 1, len(s)
         fnv = ieor(fnv, int(iachar(s(k:k)), int64)) * prime
      end do
      if (fnv == 0) fnv = 1
   end function fnv

   !> Add the key to the set of the inputs seen; false if it was there.
   recursive logical function insert(key) result(new)
      integer(int64), intent(in) :: key
      integer(int64), allocatable :: old(:)
      integer :: k, m
      m = size(seen)
      if (2*(nseen + 1) > m) then
         call move_alloc(seen, old)
         allocate(seen(0:2*m-1)); seen = 0
         nseen = 0
         do k = 0, m - 1
            if (old(k) /= 0) new = insert(old(k))
         end do
         m = size(seen)
      end if
      k = int(modulo(key, int(m, int64)))
      do while (seen(k) /= 0)
         if (seen(k) == key) then
            new = .false.
            return
         end if
         k = mod(k + 1, m)
      end do
      seen(k) = key
      nseen = nseen + 1
      new = .true.
   end function insert

   !> Link the files of the working directory into a new snapshot directory. If no file
   !> changed since the last snapshot (no newer file, the same count), the new one is a
   !> symbolic link to it, relative so that the recording can move: F_DSKKPR loads 5000
   !> files, one at a time.
   subroutine tlog_snapshot()
      character(1100) :: dir, prev
      if (.not. on()) return
      nsnap = nsnap + 1
      write(dir, '(a,a,i0)') trim(logpath), '.files/', nsnap
      write(prev, '(a,a,i0)') trim(logpath), '.files/', nsnap - 1
      call execute_command_line('d="' // trim(dir) // '"; p="' // trim(prev) // '"; ' // &
         'n=$(ls -A . | wc -l); if [ -f "$p/.stamp" ] && [ "$(cat "$p/.count")" = "$n" ] ' // &
         '&& [ -z "$(find . -maxdepth 1 -type f -newer "$p/.stamp" -print -quit)" ]; then ' // &
         'ln -s "$(basename "$(readlink -f "$p")")" "$d"; else mkdir -p "$d" && ln -f ./* "$d/" 2>/dev/null; ' // &
         'echo "$n" > "$d/.count"; touch "$d/.stamp"; fi')
      call emit('S ' // trim(dir))
   end subroutine tlog_snapshot

   !> Copy the file NAME into a new snapshot directory, at a call from the test program.
   !> A copy, not a link: the test may later open the file for write and change it. A
   !> file equal to its last copy (kept in .files/latest) is linked to that copy, which
   !> nothing changes. Past 256 MB of copies in a family the file is linked, not copied:
   !> F_ZZDSKBSR writes and opens 5000 DSK files of 4 MB, and never changes one after.
   subroutine tlog_snapshot_file(name)
      character(*), intent(in) :: name
      character(1100) :: dir
      if (.not. on() .or. depth /= 1) return
      nsnap = nsnap + 1
      write(dir, '(a,a,i0)') trim(logpath), '.files/', nsnap
      call execute_command_line('d="' // trim(dir) // '"; f="' // trim(name) // '"; ' // &
         'l="' // trim(logpath) // '.files/latest"; mkdir -p "$d" "$l"; b=$(basename "$f"); ' // &
         'if [ -f "$f" ]; then if cmp -s "$f" "$l/$b"; then ln "$l/$b" "$d/"; ' // &
         'else t=$(cat "$l/.bytes" 2>/dev/null || echo 0); if [ "$t" -gt 268435456 ]; then ' // &
         'ln -f "$f" "$d/"; else cp "$f" "$d/" && ln -f "$d/$b" "$l/$b" && ' // &
         'echo $((t + $(stat -c %s "$f"))) > "$l/.bytes"; fi; fi; fi 2>/dev/null; true')
      call emit('S ' // trim(dir))
   end subroutine tlog_snapshot_file

   !> Note the file of HANDLE before a close from the test program, if the file is open
   !> for write (or new): tlog_close_snapshot copies it after the close.
   subroutine tlog_close_name(handle)
      integer(ip), intent(in) :: handle
      type(pepper_context), pointer :: ctx
      character(1024) :: fname
      integer(ip) :: arc, bff, amh
      logical :: found
      closename = ' '
      if (.not. on() .or. depth /= 1) return
      ctx => legacy_context()
      if (ctx%failed()) return
      found = .false.
      call pepper_zzddhnfo(ctx, handle, fname, arc, bff, amh, found)
      if (found .and. (amh == 2 .or. amh == 4)) closename = fname
   end subroutine tlog_close_name

   !> Copy the file of the close into a new snapshot directory, as tlog_snapshot_file
   !> does, but past the cap of 256 MB of copies nothing is copied or linked (the
   !> directory stays empty): a family may write, close and delete thousands of files,
   !> and a link would keep each of them on the disk.
   subroutine tlog_close_snapshot()
      character(1100) :: dir
      if (.not. on() .or. depth /= 1) return
      if (closename == ' ') return
      nsnap = nsnap + 1
      write(dir, '(a,a,i0)') trim(logpath), '.files/', nsnap
      call execute_command_line('d="' // trim(dir) // '"; f="' // trim(closename) // '"; ' // &
         'l="' // trim(logpath) // '.files/latest"; mkdir -p "$d" "$l"; b=$(basename "$f"); ' // &
         'if [ -f "$f" ]; then if cmp -s "$f" "$l/$b"; then ln "$l/$b" "$d/"; ' // &
         'else t=$(cat "$l/.bytes" 2>/dev/null || echo 0); if [ "$t" -le 268435456 ]; then ' // &
         'cp "$f" "$d/" && ln -f "$d/$b" "$l/$b" && ' // &
         'echo $((t + $(stat -c %s "$f"))) > "$l/.bytes"; fi; fi; fi 2>/dev/null; true')
      call emit('S ' // trim(dir))
      closename = ' '
   end subroutine tlog_close_snapshot

   !> The mark of a writer that is not recorded (WRITE_MARKS of instrument.py).
   subroutine tlog_writer()
      if (.not. on()) return
      call emit('W')
   end subroutine tlog_writer

   !> The summary SUM of DAFBNA as its double and integer components, `<dc` and `<ic`:
   !> the integers are packed into doubles, whose bits the decimal log does not keep.
   !> No components for a handle that is not of a DAF (DAFBNA then signals the error).
   subroutine tlog_dafsum(handle, sum)
      integer(ip), intent(in) :: handle
      real(wp), intent(in) :: sum(*)
      type(pepper_context), pointer :: ctx
      character(1024) :: fname
      integer(ip) :: arc, bff, amh, nd, ni
      logical :: found
      real(wp) :: dc(125)
      integer(ip) :: ic(250)
      if (.not. on()) return
      ctx => legacy_context()
      nd = 0
      ni = 0
      found = .false.
      if (.not. ctx%failed()) call pepper_zzddhnfo(ctx, handle, fname, arc, bff, amh, found)
      if (found .and. arc == 1) call ctx%dafhsf(handle, nd, ni)
      if (ctx%failed() .or. nd < 0 .or. nd > 124 .or. ni < 2 .or. ni > 250) then
         nd = 0
         ni = 0
      end if
      if (nd + ni > 0) call ctx%dafus(sum, nd, ni, dc, ic)
      call tlog_d('<dc', dc, nd)
      call tlog_i('<ic', ic, ni)
   end subroutine tlog_dafsum

   subroutine head(label, kind, n)
      character(*), intent(in) :: label
      character(1), intent(in) :: kind
      integer, intent(in) :: n
      character(80) :: ln
      write(ln, '(a,1x,a,1x,i0)') label, kind, n
      call emit(trim(ln))
   end subroutine head

   !> True if the array LABEL of N values of SIZE bytes, too long for the log, goes
   !> to a binary file: an argument of a call from the test program, of at most 50 million
   !> values, while the family has written less than TSPICE_MAXBIN bytes (default 256 MB)
   !> of such files.
   logical function bin_ok(label, n, size)
      character(*), intent(in) :: label
      integer(ip), intent(in) :: n
      integer, intent(in) :: size
      integer :: l, stat
      character(24) :: v
      if (.not. bincheck) then
         bincheck = .true.
         call get_environment_variable('TSPICE_MAXBIN', v, l, stat)
         if (stat == 0 .and. l > 0) read(v, *) maxbin
         call execute_command_line('mkdir -p "' // trim(logpath) // '.files/bin"')
      end if
      bin_ok = depth == 1 .and. n <= 50000000 .and. &
         binbytes + int(size, int64) * n <= maxbin
   end function bin_ok

   !> The path of the binary file of values whose hash is H, of kind K; true if it exists.
   logical function bin_path(h, k, path)
      integer(int64), intent(in) :: h
      character(*), intent(in) :: k
      character(*), intent(out) :: path
      write(path, '(a,a,z16.16,a,a)') trim(logpath), '.files/bin/', h, '.', k
      inquire(file=trim(path) // '.gz', exist=bin_path)
   end function bin_path

   !> FNV-1a of the NB low bytes of W, added to the hash H. Byte by byte: a hash of whole
   !> words would keep a change of the sign bit in the sign bit alone, and two sign
   !> changes (a quaternion and its negative) would cancel.
   subroutine fnv_word(h, w, nb)
      integer(int64), intent(inout) :: h
      integer(int64), intent(in) :: w
      integer, intent(in) :: nb
      integer :: b
      do b = 0, nb - 1
         h = ieor(h, iand(ishft(w, -8*b), 255_int64)) * 1099511628211_int64
      end do
   end subroutine fnv_word

   !> Write the values to a binary file (native byte order, gzip), named by their hash so
   !> that a repeat is stored once, and log its path in a `B` line before the argument.
   subroutine bin_d(x, n)
      real(wp), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      character(1200) :: path
      integer(int64) :: h
      integer :: k, u
      h = -3750763034362895579_int64
      do k = 1, n
         call fnv_word(h, transfer(x(k), 0_int64), 8)
      end do
      if (.not. bin_path(h, 'd', path)) then
         open(newunit=u, file=trim(path), access='stream', form='unformatted', status='replace')
         write(u) x(1:n)
         close(u)
         call execute_command_line('gzip -1 -f "' // trim(path) // '"')
         binbytes = binbytes + 8_int64 * n
      end if
      call emit('B ' // trim(path) // '.gz')
   end subroutine bin_d
   subroutine bin_i(x, n)
      integer(ip), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      character(1200) :: path
      integer(int64) :: h
      integer :: k, u
      h = -3750763034362895579_int64
      do k = 1, n
         call fnv_word(h, int(x(k), int64), 4)
      end do
      if (.not. bin_path(h, 'i4', path)) then
         open(newunit=u, file=trim(path), access='stream', form='unformatted', status='replace')
         write(u) x(1:n)
         close(u)
         call execute_command_line('gzip -1 -f "' // trim(path) // '"')
         binbytes = binbytes + 4_int64 * n
      end if
      call emit('B ' // trim(path) // '.gz')
   end subroutine bin_i

   !> The strings to a binary file: N strings of their declared length L, the kind `cL`.
   subroutine bin_c(x, n)
      character(*), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      character(1200) :: path
      character(16) :: kind
      integer(int64) :: h
      integer :: k, j, u
      h = -3750763034362895579_int64
      do k = 1, n
         do j = 1, len(x(1))
            call fnv_word(h, int(iachar(x(k)(j:j)), int64), 1)
         end do
      end do
      write(kind, '(a,i0)') 'c', len(x(1))
      if (.not. bin_path(h, trim(kind), path)) then
         open(newunit=u, file=trim(path), access='stream', form='unformatted', status='replace')
         write(u) x(1:n)
         close(u)
         call execute_command_line('gzip -1 -f "' // trim(path) // '"')
         binbytes = binbytes + int(len(x(1)), int64) * n
      end if
      call emit('B ' // trim(path) // '.gz')
   end subroutine bin_c

   !> Sum over the K = 1:N mini-segments of NPKTS(K) times the size of the subtype
   !> SUBTPS(K) (subtypes 0 to 3: S0 to S3); huge, which logs the array as not held, for
   !> a subtype outside or a count past a million.
   integer(ip) function tlog_psum(npkts, subtps, n, s0, s1, s2, s3)
      integer(ip), intent(in) :: npkts(*), subtps(*), n
      integer, intent(in) :: s0, s1, s2, s3
      integer :: k, sz
      integer(int64) :: t
      tlog_psum = huge(0_ip)
      if (n > 1000000) return
      t = 0
      do k = 1, n
         select case (subtps(k))
         case (0); sz = s0
         case (1); sz = s1
         case (2); sz = s2
         case (3); sz = s3
         case default; sz = 0
         end select
         if (sz == 0) return
         t = t + int(max(0_ip, npkts(k)), int64) * sz
      end do
      if (t < huge(0_ip)) tlog_psum = int(t, ip)
   end function tlog_psum

   !> The columns of the EK segments that EKIFLD began: handle, segment, column name, entry
   !> size (-1: variable), whether nulls are allowed, and the rows of the segment.
   subroutine tlog_ek_note(handle, segno, ncols, nrows, cnames, decls)
      integer(ip), intent(in) :: handle, segno, ncols, nrows
      character(*), intent(in) :: cnames(*), decls(*)
      type(pepper_context), pointer :: ctx
      character(512) :: u
      integer :: k, j, i, l, ios, sz
      if (.not. on()) return
      ctx => legacy_context()
      if (ctx%failed()) return
      do k = 1, ncols
         u = ' '
         l = 0
         do j = 1, len_trim(decls(k))
            if (decls(k)(j:j) == ' ' .or. l >= len(u)) cycle
            l = l + 1
            u(l:l) = decls(k)(j:j)
            if (u(l:l) >= 'a' .and. u(l:l) <= 'z') u(l:l) = achar(iachar(u(l:l)) - 32)
         end do
         sz = 1
         i = index(u, 'SIZE=')
         if (i > 0) then
            j = index(u(i+5:), ',')
            if (j == 0) j = len_trim(u(i+5:)) + 1
            if (u(i+5:i+4+j-1) == 'VARIABLE') then
               sz = -1
            else
               read(u(i+5:i+4+j-1), *, iostat=ios) sz
               if (ios /= 0) sz = 1
            end if
         end if
         nek = mod(nek, size(ekh)) + 1
         ekh(nek) = handle
         eks(nek) = segno
         ekcol(nek) = upper32(cnames(k))
         eksz(nek) = sz
         eknul(nek) = index(u, 'NULLS_OK=TRUE') > 0
         ekrows(nek) = nrows
      end do
   end subroutine tlog_ek_note

   function upper32(x) result(u)
      character(*), intent(in) :: x
      character(32) :: u
      integer :: j
      u = adjustl(x)
      do j = 1, 32
         if (u(j:j) >= 'a' .and. u(j:j) <= 'z') u(j:j) = achar(iachar(u(j:j)) - 32)
      end do
   end function upper32

   !> The note of a column of an EK fast load, or 0.
   integer function ek_find(handle, segno, column)
      integer(ip), intent(in) :: handle, segno
      character(*), intent(in) :: column
      character(32) :: c
      integer :: j, k
      c = upper32(column)
      ek_find = 0
      do j = 0, size(ekh) - 1             ! the latest note first
         k = modulo(nek - 1 - j, size(ekh)) + 1
         if (ekh(k) == handle .and. eks(k) == segno .and. ekcol(k) == c .and. ekrows(k) > 0) then
            ek_find = k
            return
         end if
      end do
   end function ek_find

   !> The lengths of the arrays of EKACLC, EKACLD, EKACLI: the values (the entry sizes of
   !> the rows of a variable column, else rows times size), the entry sizes (variable
   !> columns), the null flags (columns with nulls); huge for a column without a note.
   integer(ip) function tlog_ek_nvals(handle, segno, column, entszs)
      integer(ip), intent(in) :: handle, segno, entszs(*)
      character(*), intent(in) :: column
      integer :: k, r
      integer(int64) :: t
      tlog_ek_nvals = huge(0_ip)
      k = ek_find(handle, segno, column)
      if (k == 0) return
      if (eksz(k) >= 0) then
         t = int(eksz(k), int64) * ekrows(k)
      else
         t = 0
         do r = 1, ekrows(k)
            t = t + max(0_ip, entszs(r))
         end do
      end if
      if (t < huge(0_ip)) tlog_ek_nvals = int(t, ip)
   end function tlog_ek_nvals
   integer(ip) function tlog_ek_nent(handle, segno, column)
      integer(ip), intent(in) :: handle, segno
      character(*), intent(in) :: column
      integer :: k
      tlog_ek_nent = huge(0_ip)
      k = ek_find(handle, segno, column)
      if (k == 0) return
      tlog_ek_nent = merge(ekrows(k), 0_ip, eksz(k) < 0)
   end function tlog_ek_nent
   integer(ip) function tlog_ek_nnul(handle, segno, column)
      integer(ip), intent(in) :: handle, segno
      character(*), intent(in) :: column
      integer :: k
      tlog_ek_nnul = huge(0_ip)
      k = ek_find(handle, segno, column)
      if (k == 0) return
      tlog_ek_nnul = merge(ekrows(k), 0_ip, eknul(k))
   end function tlog_ek_nnul

   !> The degree of the type 14 segment that SPK14B began on HANDLE, for SPK14A.
   subroutine tlog_spk14_note(handle, chbdeg)
      integer(ip), intent(in) :: handle, chbdeg
      type(pepper_context), pointer :: ctx
      integer :: k
      if (.not. on()) return
      ctx => legacy_context()
      if (ctx%failed()) return
      do k = 1, size(s14h)
         if (s14h(k) == handle .or. s14h(k) == 0) exit
      end do
      if (k > size(s14h)) k = 1
      s14h(k) = handle
      s14d(k) = chbdeg
   end subroutine tlog_spk14_note
   integer(ip) function tlog_spk14_len(handle, ncsets)
      integer(ip), intent(in) :: handle, ncsets
      integer :: k
      tlog_spk14_len = huge(0_ip)
      do k = 1, size(s14h)
         if (s14h(k) == handle .and. handle /= 0) then
            if (s14d(k) >= 0 .and. ncsets < 1000000) tlog_spk14_len = max(0_ip, ncsets) * (2 + 6*(s14d(k) + 1))
            return
         end if
      end do
   end function tlog_spk14_len

   !> The summary of DAFWS or DAFRS, for the DAF of the current search (DAFGH).
   subroutine tlog_dafsum_cur(sum)
      real(wp), intent(in) :: sum(*)
      type(pepper_context), pointer :: ctx
      integer(ip) :: h
      if (.not. on()) return
      ctx => legacy_context()
      h = 0
      if (.not. ctx%failed()) then
         call ctx%dafgh(h)
         ! no search: the call itself signals; this lookup leaves no error behind
         if (ctx%failed()) then
            call ctx%reset()
            h = 0
         end if
      end if
      call tlog_dafsum(h, sum)
   end subroutine tlog_dafsum_cur

   subroutine tlog_d(label, x, n)
      character(*), intent(in) :: label
      real(wp), intent(in) :: x(*)
      integer(ip), intent(in) :: n
      integer :: k
      character(32) :: ln
      if (.not. on()) return
      if (n > maxlog) then
         if (bin_ok(label, n, 8)) call bin_d(x, n)
         call head(label, 'd', -1)
         return
      end if
      call head(label, 'd', int(n))
      do k = 1, n
         write(ln, '(es25.17e3)') x(k)
         call emit(trim(ln))
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
      character(24) :: ln
      if (.not. on()) return
      if (n > maxlog) then
         if (bin_ok(label, n, 4)) call bin_i(x, n)
         call head(label, 'i', -1)
         return
      end if
      call head(label, 'i', int(n))
      do k = 1, n
         write(ln, '(i0)') x(k)
         call emit(trim(ln))
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
         call head(label, 'l', -1)
         return
      end if
      call head(label, 'l', int(n))
      do k = 1, n
         call emit(merge('T', 'F', x(k)))
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
         if (bin_ok(label, n, len(x(1)))) call bin_c(x, n)
         call head(label, 'c', -1)
         return
      end if
      call head(label, 'c', int(n))
      do k = 1, n
         call emit(tlog_str(x(k)))
      end do
   end subroutine tlog_c
   subroutine tlog_c1(label, x)
      character(*), intent(in) :: label
      character(*), intent(in) :: x
      if (.not. on()) return
      call head(label, 'c', 1)
      call emit(tlog_str(x))
   end subroutine tlog_c1
   !> A string as one line of the log: `|` and the string without its trailing blanks, or,
   !> if it holds a character outside printable ASCII (a line feed would end the line),
   !> `%` and the hexadecimal codes of its characters.
   function tlog_str(x) result(line)
      character(*), intent(in) :: x
      character(len=:), allocatable :: line
      character(len=*), parameter :: hexd = '0123456789abcdef'
      integer :: k, c, n
      n = len_trim(x)
      if (all([(iachar(x(k:k)) >= 32 .and. iachar(x(k:k)) <= 126, k = 1, n)])) then
         line = '|' // x(1:n)
         return
      end if
      allocate(character(len=1 + 2*n) :: line)
      line(1:1) = '%'
      do k = 1, n
         c = iachar(x(k:k))
         line(2*k:2*k) = hexd(c/16 + 1:c/16 + 1)
         line(2*k + 1:2*k + 1) = hexd(mod(c, 16) + 1:mod(c, 16) + 1)
      end do
   end function tlog_str
   !> Number of strings that NCHARS characters fill at the positions BPOS:EPOS of strings
   !> of length SLEN; huge(0), which logs the array as not held, for positions outside.
   integer(ip) function tlog_nsub(nchars, bpos, epos, slen)
      integer(ip), intent(in) :: nchars, bpos, epos
      integer, intent(in) :: slen
      if (bpos < 1 .or. epos > slen) then
         tlog_nsub = huge(0_ip)
      else if (epos < bpos .or. nchars < 1) then
         tlog_nsub = 0
      else
         tlog_nsub = (nchars + epos - bpos) / (epos - bpos + 1)
      end if
   end function tlog_nsub
@PURE@
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
    cases = ", &\n            ".join(", ".join(f"'{n.upper()}'" for n in sorted(PURE)[k:k + 6])
                                     for k in range(0, len(PURE), 6))
    pure = ("   logical function is_pure(name)\n      character(*), intent(in) :: name\n"
            "      select case (trim(name))\n      case (" + cases + ")\n"
            "         is_pure = .true.\n      case default\n         is_pure = .false.\n"
            "      end select\n   end function is_pure")
    qcases = ", &\n            ".join(", ".join(f"'{n.upper()}'" for n in sorted(QUERY)[k:k + 6])
                                      for k in range(0, len(QUERY), 6))
    pure += ("\n   logical function is_query(name)\n      character(*), intent(in) :: name\n"
             "      select case (trim(name))\n      case (" + qcases + ")\n"
             "         is_query = .true.\n      case default\n         is_query = .false.\n"
             "      end select\n   end function is_query")
    open(os.path.join(legacy, "tspice_log.f90"), "w").write(LOGMOD.replace("@PURE@", pure))
    missing = sorted(wanted - done)
    print(f"instrumented {len(done)} routines; not instrumented: {len(missing)}", file=sys.stderr)
    if missing:
        print("  " + " ".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main()
