"""The signature table of the recorded routines, from NAIF's Fortran SPICELIB.

    python3 signatures.py SPICELIB_DIR ROUTINES_FILE OUT_TSV

SPICELIB_DIR holds the Fortran sources of SPICELIB N0067 (`*.f` and `*.inc`), for
example `tspice/src/spicelib` of NAIF's
https://naif.jpl.nasa.gov/pub/naif/misc/tspice/N0067/PC_Linux_64bit/tspice.tar.
For each routine of ROUTINES_FILE and of WRITE_MARKS (tables.py), OUT_TSV gets one
row per argument: its kind (d, i, l, c, or p for a procedure), its dimensions as declared
(named constants replaced by their values, except in the bounds of a cell) and its
direction from the routine's Brief_I/O section. instrument_c.py reads the table to write
the wrappers of CSPICE's f2c routines, whose C prototypes give no lengths. Run it by hand
when routines.txt changes; the table is committed, the sources are not.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COLUMNS = ["routine", "kind", "ftype", "arg", "akind", "dims", "intent"]

# Brief_I/O sections that state a direction the code contradicts: routine -> {arg: intent}.
# Each of these routines writes the argument that its Brief_I/O lists as an input (CKW04A
# packs the counts of each packet into PKTDAT and shortens PKTSIZ; the EK record routines
# assign SEGNO or RECNO), as the code of each routine shows. A cell
# listed as an output is also read (its size), which `cell_intent` handles.
OVERRIDES = {
    "ckw04a": {"pktsiz": "inout", "pktdat": "inout"},
    "dafra": {"iorder": "inout"},
    "ekdelr": {"recno": "inout"},
    "ekffld": {"rcptrs": "inout"},
    **{r: {"segno": "inout"} for r in ("ekaced", "ekacei", "ekrcec", "ekrced", "ekrcei",
                                       "ekucec", "ekuced", "ekucei")},
}


def cell_intent(dims, intent):
    """A cell (dimensions LBCELL:*) is read even as an output: its size is an input."""
    return "inout" if intent == "out" and dims.upper().startswith("LBCELL") else intent

TYPE = r"(DOUBLE\s*PRECISION|INTEGER|LOGICAL|CHARACTER\s*\*\s*(?:\(\s*[^)]*\)|\d+))"
UNIT = re.compile(r"^(?:" + TYPE + r"\s+)?(SUBROUTINE|FUNCTION)\s+([A-Z0-9_]+)\s*(?:\((.*)\))?\s*$")
ENTRY = re.compile(r"^ENTRY\s+([A-Z0-9_]+)\s*(?:\((.*)\))?\s*$")
DECL = re.compile(r"^" + TYPE + r"\s+(.*)$")
PARAM = re.compile(r"^PARAMETER\s*\((.*)\)\s*$")
INCLUDE = re.compile(r"^INCLUDE\s+'([^']+)'")
IO_ROW = re.compile(r"^ {5}([A-Za-z]\w*(?:\s*,\s*[A-Za-z]\w*)*)\s*(,?)\s*(I/O|I-O|I|O|P)?(?=\s|$)")
INTENT = {"I": "in", "O": "out", "I/O": "inout", "I-O": "inout"}


def kind_of(ftype):
    t = ftype.upper().replace(" ", "")
    if t.startswith("DOUBLE"):
        return "d"
    if t.startswith("INTEGER"):
        return "i"
    if t.startswith("LOGICAL"):
        return "l"
    return "c"


def split_top(s):
    """Split at the commas outside parentheses."""
    out, depth, cur = [], 0, ""
    for ch in s:
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
    return [x.strip() for x in out if x.strip()]


def statements(text):
    """(kind, text) items in order: ('stmt', upper-case statement with continuations
    joined) or ('comment', the comment line as written)."""
    items = []
    for line in text.splitlines():
        if line[:1] in ("C", "c", "*", "!"):
            items.append(("comment", line))
            continue
        if not line.strip():
            continue
        if len(line) > 5 and line[5] not in (" ", "0") and items and items[-1][0] == "stmt":
            # a continuation line, possibly after comment lines
            for k in range(len(items) - 1, -1, -1):
                if items[k][0] == "stmt":
                    items[k] = ("stmt", items[k][1] + " " + line[6:].strip().upper())
                    break
            continue
        items.append(("stmt", line[6:].strip().upper() if len(line) > 6 else ""))
    return items


def parameters(stmts, includes, params):
    """Add the integer named constants of the statements (and their includes) to params."""
    for kind, s in stmts:
        if kind != "stmt":
            continue
        m = INCLUDE.match(s)
        if m and m.group(1).lower() in includes:     # statements are upper case
            parameters(statements(includes[m.group(1).lower()]), includes, params)
            continue
        m = PARAM.match(s)
        if m:
            for item in split_top(m.group(1)):
                if "=" not in item:
                    continue
                name, expr = (x.strip() for x in item.split("=", 1))
                v = int_value(expr, params)
                if v is not None:
                    params[name] = v


def int_value(expr, params):
    e = re.sub(r"[A-Z][A-Z0-9_]*", lambda m: str(params[m.group(0)]) if m.group(0) in params else "?",
               expr.upper().replace(" ", ""))
    if "?" in e or not re.fullmatch(r"[-+*/()0-9]+", e):
        return None
    try:
        return int(eval(e.replace("/", "//")))
    except Exception:
        return None


def dims_text(dims, params):
    if dims is None:
        return ""
    d = dims.upper().replace(" ", "")
    if ":" in d:
        return d
    out = []
    for part in split_top(d):
        v = None if part == "*" else int_value(part, params)
        out.append(part if v is None else str(v))
    return ",".join(out)


def declared(s):
    """{name: (kind, dims)} of a type declaration, or {} for another statement."""
    m = DECL.match(s)
    if not m or re.match(r"^\S+\s+FUNCTION\b", s) or re.search(r"\bFUNCTION\s", s):
        return {}
    k = kind_of(m.group(1))
    out = {}
    for item in split_top(m.group(2)):
        n = re.match(r"([A-Z0-9_]+)\s*(?:\((.*)\))?$", item)
        if n:
            out[n.group(1).lower()] = (k, n.group(2))
    return out


def arg_names(text):
    return [a.strip().lower() for a in (text or "").split(",") if a.strip()]


def parse_file(text, includes):
    """{routine: {"kind", "ftype", "args": [(name, akind, dims, intent)]}} of one source."""
    stmts = statements(text)
    params = {}
    parameters(stmts, includes, params)
    units = []          # [name, kind, ftype, args, unit index]
    decls, procs, io = [], [], {}
    cur, section, pending = None, False, []
    for kind, s in stmts:
        if kind == "comment":
            body = s[1:]
            if body.startswith("$ ") or body.startswith("$P") or body.startswith("-&"):
                section = body.strip() == "$ Brief_I/O"
                pending = []
                continue
            if not section or cur is None:
                continue
            m = IO_ROW.match(body)
            if not m:
                if not body.strip():
                    pending = []
                continue
            names = [x.strip().lower() for x in m.group(1).split(",")]
            if m.group(3):
                for n in pending + names:
                    if m.group(3) in INTENT:
                        io.setdefault(cur, {}).setdefault(n, INTENT[m.group(3)])
                pending = []
            elif m.group(2):
                pending += names
            continue
        m = UNIT.match(s)
        if m and not s.startswith("END"):
            decls.append({})
            procs.append(set())
            ftype = kind_of(m.group(1)) if m.group(2) == "FUNCTION" and m.group(1) else (
                "d" if m.group(2) == "FUNCTION" else "")
            cur = m.group(3).lower()
            units.append([cur, m.group(2).lower(), ftype, arg_names(m.group(4)), len(decls) - 1])
            continue
        m = ENTRY.match(s)
        if m and units:
            cur = m.group(1).lower()
            units.append([cur, "entry", None, arg_names(m.group(2)), len(decls) - 1])
            continue
        if not decls:
            continue
        d = declared(s)
        if d:
            decls[-1].update(d)
            continue
        if s.startswith("EXTERNAL"):
            procs[-1].update(x.strip().lower() for x in s[len("EXTERNAL"):].split(","))
    table = {}
    for name, kind, ftype, args, u in units:
        if kind == "entry":
            parent = next(x for x in units if x[4] == u and x[1] != "entry")
            ftype = parent[2]
            if ftype and name in decls[u]:
                ftype = decls[u][name][0]
        rows = []
        for a in args:
            if a in procs[u]:          # a procedure, typed or not (LOGICAL UDBAIL)
                rows.append((a, "p", "", "in"))
                continue
            k, dims = decls[u].get(a, ("?", None))
            dims = dims_text(dims, params)
            intent = OVERRIDES.get(name, {}).get(a) or cell_intent(dims, io.get(name, {}).get(a, "?"))
            rows.append((a, k, dims, intent))
        table[name] = {"kind": kind, "ftype": ftype or "", "args": rows}
    return table


def write(table, path):
    with open(path, "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for name in sorted(table):
            t = table[name]
            rows = t["args"] or [("", "", "", "")]
            for a in rows:
                f.write("\t".join([name, t["kind"], t["ftype"], *a]) + "\n")


def load(path):
    table = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        assert head == COLUMNS, head
        for line in f:
            name, kind, ftype, arg, akind, dims, intent = line.rstrip("\n").split("\t")
            t = table.setdefault(name, {"kind": kind, "ftype": ftype, "args": []})
            if arg:
                t["args"].append((arg, akind, dims, intent))
    return table


def main(spicelib, routines, out):
    sys.path.insert(0, HERE)
    from tables import WRITE_MARKS
    wanted = set(WRITE_MARKS)
    for line in open(routines):
        line = line.split("#", 1)[0].strip().lower()
        if line:
            wanted.add(line)
    includes = {fn: open(os.path.join(spicelib, fn), encoding="latin-1").read()
                for fn in os.listdir(spicelib) if fn.endswith(".inc")}
    table = {}
    for fn in sorted(os.listdir(spicelib)):
        if fn.endswith(".f"):
            for name, t in parse_file(open(os.path.join(spicelib, fn), encoding="latin-1").read(),
                                      includes).items():
                if name in wanted:
                    table[name] = t
    write(table, out)
    missing = sorted(wanted - set(table))
    unknown = sorted(f"{n}.{a[0]}" for n, t in table.items() for a in t["args"] if "?" in (a[1], a[3]))
    print(f"{len(table)} routines; not found: {len(missing)}", file=sys.stderr)
    if missing:
        print("  " + " ".join(missing), file=sys.stderr)
    if unknown:
        print(f"arguments without a kind or a direction: {' '.join(unknown)}", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
