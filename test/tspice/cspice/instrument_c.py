"""Instrument CSPICE's f2c routines so that NAIF's test suite of CSPICE (tspice_c)
records the calls that its f2c-translated families make, into the log that tspice_log.c
writes (its format is described in tspice_log.h).

    python3 instrument_c.py TSPICE_C_DIR ROUTINES_FILE OUT_DIR

TSPICE_C_DIR is the extracted tspice_c archive of NAIF. OUT_DIR gets:

- tlog_wrap.c: for each recorded routine R, `__wrap_R_`, which logs the call and calls
  `R_`; the hooks that tspice_log.c calls; the lists PURE and QUERY of tables.py;
- wrap.syms: the lines `R_ __wrap_R_`, for `objcopy --redefine-syms` on the objects of
  the test program (the families and NAIF's test and support utilities), so that their
  calls, and only theirs, go through the wrappers: CSPICE's calls between its own
  routines are not logged;
- driver.c: the main program of tspice_c, which runs one f2c family by name
  (`--list` lists them).

The C prototypes come from TSPICE_C_DIR/include/SpiceZfc.h; the kinds, directions and
dimensions of the arguments from signatures.tsv. A routine that cannot be recorded is
reported on standard error as `skip NAME: reason`. Nothing NAIF wrote is copied: the
driver is made from tspice_c.pgm at build time.
"""

import os
import re
import sys
from collections import namedtuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import signatures  # noqa: E402
import tables as T  # noqa: E402

Proto = namedtuple("Proto", "ret params character")
Generated = namedtuple("Generated", "code wrapped symbols skips marks")

PROTO = re.compile(r"^extern\s+(.*?)\s*\b(\w+)_\((.*)\);\s*$")
RENAME = re.compile(r"^\s*#\s*define\s+(\w+)_\s+(\w+_)\s*$")

# Lengths that the arguments do not give (SPICELIB's documented sizes), as C
# expressions: `{x}` is the C parameter of the Fortran argument x (a pointer: `*{n}` is the value of a scalar, `{x}[k-1]` is X(K)),
# `{x_len}` its string length.
M = "tlog_max0"
FIXED_LENGTHS_C = {
    "oscltx": {"elts": "20"},
    "ekssum": {a: "tlog_okn(*{ncols})" for a in ("cnames", "dtypes", "sizes", "strlns", "indexd", "nullok")},
    "limbpt": {"npts": f"{M}(*{{ncuts}})", "points": f"3*tlog_sum_i({{npts}}, {M}(*{{ncuts}}))",
               "epochs": f"tlog_sum_i({{npts}}, {M}(*{{ncuts}}))",
               "tangts": f"3*tlog_sum_i({{npts}}, {M}(*{{ncuts}}))"},
    "termpt": {"npts": f"{M}(*{{ncuts}})", "points": f"3*tlog_sum_i({{npts}}, {M}(*{{ncuts}}))",
               "epochs": f"tlog_sum_i({{npts}}, {M}(*{{ncuts}}))",
               "trmvcs": f"3*tlog_sum_i({{npts}}, {M}(*{{ncuts}}))"},
    "dasadc": {"data": "tlog_nsub(*{n}, *{bpos}, *{epos}, {data_len})"},
    "dasudc": {"data": "tlog_nsub(*{last} - *{first} + 1, *{bpos}, *{epos}, {data_len})"},
    "spkw10": {"consts": "8", "elems": f"{M}(10*(long)*{{n}})"},
    "pckw02": {"cdata": f"{M}(3*(*{{polydg}} + 1)*(long)*{{n}})"},
    "pckw20": {"cdata": f"{M}(3*(*{{polydg}} + 2)*(long)*{{n}})"},
    "spkw18": {"packts": f"{M}(*{{n}}*(long)(*{{subtyp}} == 0 ? 12 : 6))"},
    "dafps": {"dc": f"{M}(*{{nd}})", "ic": f"{M}(*{{ni}})", "sum": f"{M}(*{{nd}} + (*{{ni}} + 1)/2)"},
    "dafus": {"sum": f"{M}(*{{nd}} + (*{{ni}} + 1)/2)", "dc": f"{M}(*{{nd}})", "ic": f"{M}(*{{ni}})"},
    "dafgda": {"data": f"{M}(*{{eaddr}} - *{{baddr}} + 1)"},
    "dafgsr": {"data": f"{M}(*{{end}} - *{{begin}} + 1)"},
    "dasrdd": {"data": f"{M}(*{{last}} - *{{first}} + 1)"},
    "dasrdi": {"data": f"{M}(*{{last}} - *{{first}} + 1)"},
    "dasudd": {"data": f"{M}(*{{last}} - *{{first}} + 1)"},
    "dasudi": {"data": f"{M}(*{{last}} - *{{first}} + 1)"},
    "dasrdc": {"data": "tlog_nsub(*{last} - *{first} + 1, *{bpos}, *{epos}, {data_len})"},
    "dlabfs": {"dladsc": "8"},
    "dlabbs": {"dladsc": "8"},
    "dlafns": {"dladsc": "8", "nxtdsc": "8"},
    "dlafps": {"dladsc": "8", "prvdsc": "8"},
    "dskgd": {"dladsc": "8", "dskdsc": "24"},
    **{r: {"dladsc": "8"} for r in ("dskn02", "dskv02", "dskp02", "dskz02", "dski02", "dskd02", "dskb02",
                                    "dskx02")},
    "dskxsi": {"srflst": f"{M}(*{{nsurf}})", "dladsc": "8", "dskdsc": "24", "dc": f"{M}(*{{maxd}})",
               "ic": f"{M}(*{{maxi}})"},
    "dskxv": {"srflst": f"{M}(*{{nsurf}})", "vtxarr": f"3*{M}(*{{nrays}})", "dirarr": f"3*{M}(*{{nrays}})",
              "xptarr": f"3*{M}(*{{nrays}})", "fndarr": f"{M}(*{{nrays}})"},
    "spksfs": {"descr": "5"},
    "spkpds": {"descr": "5"},
    "spkuds": {"descr": "5"},
    "ckgr02": {"descr": "5", "record": "10"},
    "ckgr03": {"descr": "5", "record": "8"},
    "cknr02": {"descr": "5"},
    "cknr03": {"descr": "5"},
    "cknr04": {"descr": "5"},
    "cknm06": {"descr": "5"},
    "ckmp06": {"descr": "5"},
    # a type 6 record: the epoch, the subtype, the clock rate and the packet of the subtype
    "ckgr06": {"descr": "5", "record": "3 + (tlog_nint({record}[1]) == 0 ? 8 : tlog_nint({record}[1]) == 1 ? 4"
                                       " : tlog_nint({record}[1]) == 2 ? 14 : 7)"},
    "ckw04a": {"pktsiz": f"{M}(*{{npkts}})", "sclkdp": f"{M}(*{{npkts}})",
               "pktdat": f"{M}(tlog_sum_i({{pktsiz}}, {M}(*{{npkts}})))"},
    "ckw01": {"sclkdp": f"{M}(*{{nrec}})", "quats": f"{M}(4*(long)*{{nrec}})", "avvs": f"{M}(3*(long)*{{nrec}})"},
    "ckw02": {"start": f"{M}(*{{nrec}})", "stop": f"{M}(*{{nrec}})", "quats": f"{M}(4*(long)*{{nrec}})",
              "avvs": f"{M}(3*(long)*{{nrec}})", "rates": f"{M}(*{{nrec}})"},
    "ckw03": {"sclkdp": f"{M}(*{{nrec}})", "quats": f"{M}(4*(long)*{{nrec}})", "avvs": f"{M}(3*(long)*{{nrec}})",
              "starts": f"{M}(*{{nints}})"},
    # the integer spatial index: SICGRD - 1 + MAXCGR + the voxel pointers, the voxel-plate
    # list, the vertex pointers (NV) and the vertex-plate list (dsk02.inc)
    "dskw02": {"corpar": "10", "spaixd": "10", "spaixi":
               "(({spaixi}[4] >= 0 && {spaixi}[5] >= 0 && {spaixi}[6] >= 0 && {spaixi}[4] < 50000000"
               " && {spaixi}[5] < 50000000 && {spaixi}[6] < 50000000 && *{nv} < 50000000)"
               f" ? 7L + 100000 + {{spaixi}}[4] + {{spaixi}}[5] + {M}(*{{nv}}) + {{spaixi}}[6] : TLOG_HUGE)"},
    "ekbseg": {"cnames": f"{M}(*{{ncols}})", "decls": f"{M}(*{{ncols}})"},
    "chbint": {"cp": f"{M}(*{{degp}} + 1)"},
    "chbval": {"cp": f"{M}(*{{degp}} + 1)"},
    "chbigr": {"cp": f"{M}(*{{degp}} + 1)"},
    "chbder": {"cp": f"{M}(*{{degp}} + 1)", "dpdxs": f"{M}(*{{nderiv}} + 1)"},
    "polyds": {"coeffs": f"{M}(*{{deg}} + 1)", "p": f"{M}(*{{nderiv}} + 1)"},
    "xposeg": {"matrix": f"{M}(*{{nrow}}*(long)*{{ncol}})", "xposem": f"{M}(*{{nrow}}*(long)*{{ncol}})"},
    "swpool": {"names": f"{M}(*{{nnames}})"},
    "scpart": {"pstart": f"{M}(*{{nparts}})", "pstop": f"{M}(*{{nparts}})"},
    **{r: {"v1": f"{M}(*{{ndim}})", "v2": f"{M}(*{{ndim}})"} for r in ("vdistg", "vdotg", "vrelg", "vsepg")},
    "vzerog": {"v": f"{M}(*{{ndim}})"},
    # mini-segments: the packet size of each subtype (spk19.inc, ck06.inc)
    "spkw19": {"npkts": f"{M}(*{{nintvl}})", "subtps": f"{M}(*{{nintvl}})", "degres": f"{M}(*{{nintvl}})",
               "ivlbds": "(*{nintvl} > 0 ? *{nintvl} + 1 : 0)",
               "epochs": "tlog_psum({npkts}, {subtps}, *{nintvl}, 1, 1, 1, 1)",
               "packts": "tlog_psum({npkts}, {subtps}, *{nintvl}, 12, 6, 6, 0)"},
    "ckw06": {"npkts": f"{M}(*{{nmini}})", "subtps": f"{M}(*{{nmini}})", "degres": f"{M}(*{{nmini}})",
              "rates": f"{M}(*{{nmini}})", "ivlbds": "(*{nmini} > 0 ? *{nmini} + 1 : 0)",
              "sclkdp": "tlog_psum({npkts}, {subtps}, *{nmini}, 1, 1, 1, 1)",
              "packts": "tlog_psum({npkts}, {subtps}, *{nmini}, 8, 4, 14, 7)"},
    # the EK fast load: the rows and the declarations that EKIFLD noted
    "ekifld": {"cnames": f"{M}(*{{ncols}})", "decls": f"{M}(*{{ncols}})", "rcptrs": f"{M}(*{{nrows}})"},
    **{r: {vals: "tlog_ek_nvals(*{handle}, *{segno}, {column}, {column_len}, {entszs})",
           "entszs": "tlog_ek_nent(*{handle}, *{segno}, {column}, {column_len})",
           "nlflgs": "tlog_ek_nnul(*{handle}, *{segno}, {column}, {column_len})"}
       for r, vals in (("ekaclc", "cvals"), ("ekacld", "dvals"), ("ekacli", "ivals"))},
    # records of 2 + 6 (CHBDEG + 1) doubles, the degree of SPK14B
    "spk14a": {"coeffs": "tlog_spk14_len(*{handle}, *{ncsets})", "epochs": f"{M}(*{{ncsets}})"},
    "ckw05": {"packts": f"{M}(*{{n}}*(long)(*{{subtyp}} == 0 ? 8 : *{{subtyp}} == 1 ? 4"
                        " : *{subtyp} == 2 ? 14 : 7))"},
}
# Arguments logged by a routine of tspice_log.c instead of as an array: the summary of
# DAFBNA packs its integers into doubles, which the decimal log cannot carry; it is
# logged unpacked, as `<dc` and `<ic`.
# None: the argument is not logged (a workspace, or the record pointers of an EK fast load,
# which the replay keeps as the state that Almanac's EKIFLD returns).
CUSTOM_ARGS_C = {
    "dafbna": {"sum": "tlog_dafsum(*{handle}, {sum}, '<')"},
    "dafws": {"sum": "tlog_dafsum_cur({sum}, '<')"},
    "dafrs": {"sum": "tlog_dafsum_cur({sum}, '<')"},
    "dskmi2": {"work": None},
    "chbder": {"partdp": None},
    "ekaclc": {"rcptrs": None, "wkindx": None},
    "ekacld": {"rcptrs": None, "wkindx": None},
    "ekacli": {"rcptrs": None, "wkindx": None},
    "ekffld": {"rcptrs": None},
}
# Output arguments logged by a routine of tspice_log.c, after the call: the summary that
# DAFGS returns, unpacked as `>dc` and `>ic` like the inputs above (its length is that of
# the DAF of the current search, which no argument gives).
CUSTOM_OUT_C = {"dafgs": {"sum": "tlog_dafsum_cur({sum}, '>')"}}
# Calls after the call, before its outputs: the notes from which later calls take their
# lengths (the rows and column declarations of an EK fast load, the degree of SPK14B).
CUSTOM_POST_C = {
    "ekifld": "tlog_ek_note(*{handle}, *{segno}, *{ncols}, *{nrows}, {cnames}, {cnames_len}, {decls}, {decls_len})",
    "spk14b": "tlog_spk14_note(*{handle}, *{chbdeg})",
}
# The input window of WNVALD holds N endpoints that are not yet a window (its
# cardinality is not set): they are logged as N values after the size and the old
# cardinality.
CELL_INPUTS_C = {"wnvald": {"window": 'tlog_d("{label}", {window} + 4, 2 + tlog_clamp(*{n}))'}}

# The CSPICE routines that tspice_log.c calls through its hooks.
HOOKS = {
    "failed": ("logical tlog_hook_failed(void)", "return {f}();"),
    "getmsg": ("void tlog_hook_getmsg_short(char *out, ftnlen len)", '{f}((char *)"SHORT", out, (ftnlen)5, len);'),
    "zzddhnfo": ("void tlog_hook_zzddhnfo(integer h, char *fname, ftnlen len, integer *arc, integer *bff, "
                 "integer *amh, logical *found)", "{f}(&h, fname, arc, bff, amh, found, len);"),
    "dafhsf": ("void tlog_hook_dafhsf(integer h, integer *nd, integer *ni)", "{f}(&h, nd, ni);"),
    "dafus": ("void tlog_hook_dafus(const doublereal *sum, integer nd, integer ni, doublereal *dc, integer *ic)",
              "{f}((doublereal *)sum, &nd, &ni, dc, ic);"),
    "dafgh": ("void tlog_hook_dafgh(integer *h)", "{f}(h);"),
    "reset": ("void tlog_hook_reset(void)", "{f}();"),
}
KIND_FN = {"d": "tlog_d", "i": "tlog_i", "l": "tlog_l"}
RESULT_TYPE = {"d": "doublereal", "i": "integer", "l": "logical"}


def parse_prototypes(text):
    """{routine: Proto(return type, [(C type, C name)], is a CHARACTER function)}."""
    out = {}
    for line in text.splitlines():
        m = PROTO.match(line.strip())
        if not m:
            continue
        params = []
        body = m.group(3).strip()
        if body not in ("", "void"):
            for p in body.split(","):
                pm = re.match(r"^(.*?)(\w+)$", p.strip())
                params.append((pm.group(1).strip(), pm.group(2)))
        out[m.group(2)] = Proto(m.group(1).strip(), params, "Character" in m.group(1))
    return out


def parse_renames(text):
    """{routine: symbol} of SpiceZrnm.h, which renames a few f2c symbols (vnorm_ is
    zz_vnorm_): the wrappers take the symbol, the log the routine's name."""
    return {m.group(1): m.group(2) for m in map(RENAME.match, text.splitlines()) if m}


def decl(params):
    return ", ".join(t + n if t.endswith("*") else f"{t} {n}" for t, n in params) or "void"


def call_args(params):
    return ", ".join(n for _, n in params)


def c_names(proto, sig):
    """{fortran arg: C parameter} and {fortran arg + "_len": its ftnlen}, or a reason."""
    args = sig["args"]
    chars = [a for a, k, _, _ in args if k == "c"]
    if len(proto.params) != len(args) + len(chars):
        return None, "the prototype does not match the signature"
    names = {a[0]: proto.params[k][1] for k, a in enumerate(args)}
    for k, a in enumerate(chars):
        t, n = proto.params[len(args) + k]
        if t != "ftnlen":
            return None, "the prototype does not match the signature"
        names[a + "_len"] = n
    return names, None


def length_c(dims, args, kinds, arg, names):
    """The number of elements of an array argument, as a C expression, or None."""
    if not dims:
        return "1"
    d = dims.lower()
    if arg in T.FIXED_SIZES and d == "*":
        return T.FIXED_SIZES[arg]
    parts = d.split(",")
    for k, p in enumerate(parts):
        if ":" in p:
            lo, hi = p.split(":", 1)
            if hi == "*" or not re.match(r"^-?\d+$", lo):
                return None
            parts[k] = f"({hi})-({lo})+1"

    def tr(e):
        bad = []

        def sub(m):
            n = m.group(0)
            if n in args and kinds.get(n) == ("i", ""):
                return f"(*{names[n]})"
            bad.append(n)
            return n
        e = re.sub(r"[a-z_][a-z0-9_]*", sub, e)
        return None if bad else e

    if parts[-1] == "*":
        named = [T.NAMED_COUNTS[arg]] if arg in T.NAMED_COUNTS else []
        for cand in named + ["n", "dim", "nvals", "nv", "nsize", "size", "nvars", "npts"]:
            if cand in args and kinds.get(cand) == ("i", ""):
                lead = [tr(p) for p in parts[:-1]]
                if None in lead:
                    return None
                lead = "*".join(f"({p})" for p in lead) or "1"
                return f"tlog_max0(({lead})*(long)(*{names[cand]}))"
        return None
    if all(re.fullmatch(r"\d+", p) for p in parts):
        n = 1
        for p in parts:
            n *= int(p)
        return str(n)
    parts = [tr(p) for p in parts]
    if None in parts:
        return None
    return "*".join(f"({p})" for p in parts)


def cell_records(name, sig):
    """{arg: record with {label} and the C names, or None} of a routine of CELL_ROUTINES."""
    if name not in T.CELL_ROUTINES:
        return {}
    cells = {}
    for a, k, dims, _ in sig["args"]:
        if k == "p" or not dims.lower().startswith("lbcell"):
            continue
        if dims.replace(" ", "").lower() != "lbcell:*":
            cells[a] = None
        elif k == "d":
            cells[a] = f'tlog_d("{{label}}", {{{a}}} + 4, 2 + tlog_card_d({{{a}}}))'
        elif k == "i":
            cells[a] = f'tlog_i("{{label}}", {{{a}}} + 4, 2 + tlog_card_i({{{a}}}))'
    return cells


def record(kind, label, c, n, names, a):
    if kind == "c":
        return f'tlog_c("{label}", {c}, {names[a + "_len"]}, {n})'
    return f'{KIND_FN[kind]}("{label}", {c}, {n})'


def wrapper(sym, proto, pre, post):
    ret = "int" if proto.ret.startswith("int") else proto.ret
    args = call_args(proto.params)
    body = [f"extern {ret} {sym}({decl(proto.params)});",
            f"{ret} __wrap_{sym}({decl(proto.params)})", "{", f"    {ret} r;", "    tlog_enter();",
            "    if (tlog_depth() > 1 || !tlog_on()) {", f"        r = {sym}({args});",
            "        tlog_leave();", "        return r;", "    }"]
    body += [f"    {x};" for x in pre]
    body.append(f"    r = {sym}({args});")
    body += [f"    {x};" for x in post]
    body += ["    tlog_leave();", "    return r;", "}", ""]
    return "\n".join(body)


def instrument_routine(name, protos, table, sym):
    """The wrapper of a recorded routine, or a reason."""
    if name not in protos:
        return "no prototype in SpiceZfc.h"
    if name not in table:
        return "no signature in signatures.tsv"
    proto, sig = protos[name], table[name]
    if proto.character or sig["ftype"] == "c":
        return "a CHARACTER function"
    if not proto.ret.startswith("int") and proto.ret not in RESULT_TYPE.values():
        return f"a function of type {proto.ret}"
    if any("?" in (k, it) for _, k, _, it in sig["args"]):
        return "an argument without a kind or a direction"
    procs = {a for a, k, _, _ in sig["args"] if k == "p"}
    if procs and name not in T.PROC_ROUTINES:
        return "an argument is a procedure"
    names, why = c_names(proto, sig)
    if why:
        return why
    args = [a for a, _, _, _ in sig["args"]]
    kinds = {a: (k, dims) for a, k, dims, _ in sig["args"]}
    intent = {a: it for a, _, _, it in sig["args"]}
    fixed = FIXED_LENGTHS_C.get(name, {})
    custom = CUSTOM_ARGS_C.get(name, {})
    custom_out = CUSTOM_OUT_C.get(name, {})
    cells = cell_records(name, sig)
    lens = {}
    for a in args:
        if a in custom or a in custom_out or a in cells or a in procs:
            lens[a] = None
        elif a in fixed:
            lens[a] = fixed[a].format(**names)
        else:
            n = length_c(kinds[a][1], args, kinds, a, names)
            if n is None:
                return f"the length of argument {a} ({kinds[a][1]}) is not known"
            lens[a] = n
    first = names[args[0]] if args else None
    pre = [f'tlog_begin("{name.upper()}")']
    if name in T.FILE_ROUTINES:
        pre.append("tlog_snapshot()")
    elif name in T.OPEN_ROUTINES:
        pre.append(f"tlog_snapshot_file({first}, {names[args[0] + '_len']})")
    elif name in T.CLOSE_ROUTINES:
        pre.append(f"tlog_close_name(*{first})")
    for a in args:
        if a in T.MARK_HANDLES and kinds[a][0] == "i" and lens.get(a) == "1" and intent[a] in ("in", "inout"):
            pre.append(f"tlog_handle_file(*{names[a]})")
    post = ["tlog_mark()"]
    if name in T.CLOSE_ROUTINES:
        post.append("tlog_close_snapshot()")
    if name in CUSTOM_POST_C:
        post.append(CUSTOM_POST_C[name].format(**names))
    for a in args:
        if a in procs:
            continue
        k, it, n, c = kinds[a][0], intent[a], lens[a], names[a]
        if a in custom:
            if custom[a] is not None:
                pre.append(custom[a].format(**names))
            continue
        if a in custom_out:
            post.append(custom_out[a].format(**names))
            continue
        if a in cells:
            if cells[a] is not None and it in ("in", "inout"):
                rec = CELL_INPUTS_C.get(name, {}).get(a, cells[a])
                pre.append(rec.format(label="<" + a, **names))
            if cells[a] is not None and it in ("out", "inout"):
                post.append(cells[a].format(label=">" + a, **names))
            continue
        if it in ("in", "inout"):
            pre.append(record(k, "<" + a, c, n, names, a))
        if it in ("out", "inout"):
            post.append(record(k, ">" + a, c, n, names, a))
    if sig["ftype"]:
        post.append(f'{KIND_FN[sig["ftype"]]}(">result", &r, 1)')
    post.append("tlog_end()")
    return wrapper(sym, proto, pre, post)


def mark_routine(name, protos, table, sym):
    """The wrapper of a writer that is not recorded: a `W` line and its handles only."""
    if name not in protos or name not in table:
        return "no prototype or no signature"
    proto, sig = protos[name], table[name]
    names, why = c_names(proto, sig)
    if why:
        return why
    pre = [f'tlog_begin("{name.upper()}")', "tlog_writer()"]
    for a, k, dims, _ in sig["args"]:
        if a in T.MARK_HANDLES and k == "i" and dims == "":
            pre.append(f'tlog_i("<{a}", {names[a]}, 1)')
    return wrapper(sym, proto, pre, ["tlog_mark()", "tlog_end()"])


def name_list(fn, names):
    items = ",\n        ".join(", ".join(f'"{n.upper()}"' for n in names[k:k + 6]) for k in range(0, len(names), 6))
    return (f"int {fn}(const char *name)\n{{\n    static const char *const names[] = {{\n        {items}\n    }};\n"
            f"    return bsearch(&name, names, sizeof names / sizeof names[0], sizeof names[0], tlog_cmp) != NULL;\n}}\n")


def generate(protos, table, wanted, renames=None):
    parts = ["/* Generated by test/tspice/cspice/instrument_c.py: do not edit. */", "",
             "#include <stdlib.h>", "#include <string.h>", '#include "f2c.h"', '#include "tspice_log.h"', ""]
    renames = renames or {}
    wrapped, symbols, skips, marks = [], [], [], []
    for name in sorted(set(wanted) | T.WRITE_MARKS):
        sym = renames.get(name, name + "_")
        res = instrument_routine(name, protos, table, sym) if name in wanted else "not in the routines file"
        if res.startswith("extern"):
            parts.append(res)
            wrapped.append(name)
            symbols.append(sym)
            continue
        if name in T.WRITE_MARKS:
            marks.append((name, res))
            res = mark_routine(name, protos, table, sym)
            if res.startswith("extern"):
                parts.append(res)
                wrapped.append(name)
                symbols.append(sym)
                continue
        skips.append((name, res))
    for h, (signature, body) in HOOKS.items():
        p = protos[h]
        f = renames.get(h, h + "_")
        if h not in wrapped:                 # else the wrapper declares it
            parts.append(f"extern {'int' if p.ret.startswith('int') else p.ret} {f}({decl(p.params)});")
        parts.append(f"{signature}\n{{\n    {body.format(f=f)}\n}}\n")
    parts.append("static int tlog_cmp(const void *a, const void *b)\n{\n"
                 "    return strcmp(*(const char *const *)a, *(const char *const *)b);\n}\n")
    parts.append(name_list("tlog_is_pure", sorted(n.upper() for n in T.PURE)))
    parts.append(name_list("tlog_is_query", sorted(n.upper() for n in T.QUERY)))
    return Generated("\n".join(parts), wrapped, symbols, skips, marks)


FAMILY = re.compile(r"^\s*(f_\w+__)\s*\(\s*&plain_ok\s*\)\s*;\s*$")
CFAMILY = re.compile(r"^\s*(f_\w+_c)\s*\(\s*&spice_ok\s*\)\s*;\s*$")
MAIN = re.compile(r"^(\s*)int\s+main\s*\(\s*int\s+argc\s*,\s*char\s*\*\*\s*argv\s*\)\s*$")


def driver(pgm):
    """The main program of tspice_c (tspice_c.pgm) made to run one f2c family by name. Its exit
    status is 1 if a test case failed (T_ANYBAD), 0 otherwise, 2 for an unknown family."""
    out, fams, replaced = ["#include <stdio.h>", "#include <string.h>"], [], False
    for line in pgm.splitlines():
        m = MAIN.match(line)
        if m:
            out.append(f"{m.group(1)}static void tlog_dispatch ( logical * plain_ok );")
            out.append(f"{m.group(1)}static int naif_main( int argc, char **argv )")
            continue
        m = FAMILY.match(line)
        if m:
            fams.append(m.group(1))
            if not replaced:
                # the flag of the C-API families, which this driver does not run
                out.append("   spice_ok = SPICETRUE;")
                out.append("   tlog_dispatch ( &plain_ok );")
                replaced = True
            continue
        if CFAMILY.match(line):
            continue
        out.append(line)
    assert replaced, "no family call in the main program"
    names = ["F_" + f[2:-2].upper() for f in fams]
    out.append("")
    out += [f"int {f} ( logical * plain_ok );" for f in fams]
    out.append("int t_anybad__ ( logical * ok );")
    out.append("static const char *const FAMILIES[] = {\n    " +
               ",\n    ".join(f'"{n}"' for n in names) + "\n};")
    out.append("static int chosen = -1;")
    out.append("static void tlog_dispatch ( logical * plain_ok )\n{\n   switch ( chosen )\n   {")
    out += [f"   case {k}: {f} ( plain_ok ); break;" for k, f in enumerate(fams)]
    out.append("   }\n}")
    out.append("""int main(int argc, char **argv)
{
   int k, n = (int)(sizeof FAMILIES / sizeof FAMILIES[0]);
   if (argc == 2 && strcmp(argv[1], "--list") == 0)
   {
      for (k = 0; k < n; k++) puts(FAMILIES[k]);
      return 0;
   }
   if (argc != 2)
   {
      fprintf(stderr, "usage: %s FAMILY | --list\\n", argv[0]);
      return 2;
   }
   for (k = 0; k < n; k++) if (strcmp(argv[1], FAMILIES[k]) == 0) chosen = k;
   if (chosen < 0)
   {
      fprintf(stderr, "unknown family %s\\n", argv[1]);
      return 2;
   }
   naif_main(1, argv);
   {
      /* failure if any test case failed (the flag of a family is that of its last
         check) */
      logical bad = 0;
      t_anybad__(&bad);
      return bad ? 1 : 0;
   }
}
""")
    return "\n".join(out)


def main(tspice_c, routines, outdir):
    protos = parse_prototypes(open(os.path.join(tspice_c, "include", "SpiceZfc.h"), encoding="latin-1").read())
    table = signatures.load(os.path.join(HERE, "signatures.tsv"))
    wanted = set()
    for line in open(routines):
        line = line.split("#", 1)[0].strip().lower()
        if line:
            wanted.add(line)
    os.makedirs(outdir, exist_ok=True)
    renames = parse_renames(open(os.path.join(tspice_c, "include", "SpiceZrnm.h"), encoding="latin-1").read())
    g = generate(protos, table, wanted, renames)
    open(os.path.join(outdir, "tlog_wrap.c"), "w").write(g.code)
    open(os.path.join(outdir, "wrap.syms"), "w").write("".join(f"{n} __wrap_{n}\n" for n in g.symbols))
    pgm = open(os.path.join(tspice_c, "src", "tspice_c", "tspice_c.pgm"), encoding="latin-1").read()
    open(os.path.join(outdir, "driver.c"), "w").write(driver(pgm))
    for n, why in g.marks:
        print(f"mark {n}: {why}", file=sys.stderr)
    for n, why in g.skips:
        print(f"skip {n}: {why}", file=sys.stderr)
    done = [n for n in g.wrapped if n in wanted]
    missing = sorted(wanted - set(done))
    print(f"instrumented {len(done)} routines; not instrumented: {len(missing)}", file=sys.stderr)
    if missing:
        print("  " + " ".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
