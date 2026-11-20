"""Tests of instrument_c.py on prototypes, signatures and a driver written for the tests."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import instrument_c as G  # noqa: E402
import instrument  # noqa: E402  (test/tspice, put on the path by instrument_c)

PROTOS = """\
extern int spkx_(char *targ, doublereal *et, doublereal *starg, ftnlen targ_len);
extern int twos_(char *name__, char *other, integer *n, ftnlen name_len, ftnlen other_len);
extern int vsum_(integer *n, doublereal *x, doublereal *s);
extern int wnvald_(integer *size, integer *n, doublereal *window);
extern int wnfltd_(doublereal *smlint, doublereal *window);
extern doublereal vnorm_(doublereal *v1);
extern int dafwdr_(integer *handle, integer *recno, doublereal *drec);
extern int dafus_(doublereal *sum, integer *nd, integer *ni, doublereal *dc, integer *ic);
extern logical failed_(void);
extern int getmsg_(char *option, char *msg, ftnlen option_len, ftnlen msg_len);
extern int reset_(void);
extern int zzddhnfo_(integer *handle, char *fname, integer *intarc, integer *intbff, integer *intamh, logical *found, ftnlen fname_len);
extern int dafhsf_(integer *handle, integer *nd, integer *ni);
extern int dafgh_(integer *handle);
extern /* Character */ VOID chfun_(char *ret_val, ftnlen ret_val_len, integer *n);
extern int unk_(integer *k, doublereal *x);
"""


def sig(kind, ftype, *args):
    return {"kind": kind, "ftype": ftype, "args": list(args)}


TABLE = {
    "spkx": sig("subroutine", "", ("targ", "c", "", "in"), ("et", "d", "", "in"), ("starg", "d", "6", "out")),
    "twos": sig("subroutine", "", ("name", "c", "", "in"), ("other", "c", "*", "in"), ("n", "i", "", "in")),
    "vsum": sig("subroutine", "", ("n", "i", "", "in"), ("x", "d", "*", "in"), ("s", "d", "", "out")),
    "wnvald": sig("subroutine", "", ("size", "i", "", "in"), ("n", "i", "", "in"),
                  ("window", "d", "LBCELL:*", "inout")),
    "wnfltd": sig("subroutine", "", ("smlint", "d", "", "in"), ("window", "d", "LBCELL:*", "inout")),
    "vnorm": sig("function", "d", ("v1", "d", "3", "in")),
    "dafwdr": sig("subroutine", "", ("handle", "i", "", "in"), ("recno", "i", "", "in"),
                  ("drec", "d", "128", "in")),
    "dafus": sig("subroutine", "", ("sum", "d", "*", "in"), ("nd", "i", "", "in"), ("ni", "i", "", "in"),
                 ("dc", "d", "*", "out"), ("ic", "i", "*", "out")),
    "chfun": sig("function", "c", ("n", "i", "", "in")),
    "unk": sig("subroutine", "", ("k", "i", "", "in"), ("x", "d", "*", "in")),
}


def wrapper_of(text, name):
    start = text.index(f" __wrap_{name}_(")
    return text[start:text.index("\n}\n", start)]


def in_order(test, text, parts):
    pos = 0
    for p in parts:
        k = text.find(p, pos)
        test.assertNotEqual(k, -1, f"{p!r} not found after position {pos} in:\n{text}")
        pos = k + len(p)


class GeneratorTests(unittest.TestCase):
    def setUp(self):
        self.protos = G.parse_prototypes(PROTOS)
        wanted = {"spkx", "twos", "vsum", "wnvald", "wnfltd", "vnorm", "dafus", "chfun", "unk"}
        self.out = G.generate(self.protos, TABLE, wanted)

    def test_scalar_and_array(self):
        w = wrapper_of(self.out.code, "spkx")
        in_order(self, w, ['tlog_begin("SPKX")', 'tlog_c("<targ", targ, targ_len, 1)',
                           'tlog_d("<et", et, 1)', "spkx_(targ, et, starg, targ_len)",
                           "tlog_mark()", 'tlog_d(">starg", starg, 6)', "tlog_end()"])

    def test_string_args(self):
        w = wrapper_of(self.out.code, "twos")
        in_order(self, w, ['tlog_c("<name", name__, name_len, 1)',
                           'tlog_c("<other", other, other_len, tlog_max0((1)*(long)(*n)))'])

    def test_assumed_size(self):
        w = wrapper_of(self.out.code, "vsum")
        self.assertIn('tlog_d("<x", x, tlog_max0((1)*(long)(*n)))', w)

    def test_fixed_lengths(self):
        self.assertEqual({k: set(v) for k, v in G.FIXED_LENGTHS_C.items()},
                         {k: set(v) for k, v in instrument.FIXED_LENGTHS.items()})

    def test_cell_arg(self):
        w = wrapper_of(self.out.code, "wnfltd")
        in_order(self, w, ['tlog_d("<window", window + 4, 2 + tlog_card_d(window))',
                           'tlog_d(">window", window + 4, 2 + tlog_card_d(window))'])
        # WNVALD's input window holds N endpoints that are not yet a window
        w = wrapper_of(self.out.code, "wnvald")
        self.assertIn('tlog_d("<window", window + 4, 2 + tlog_clamp(*n))', w)

    def test_function_result(self):
        w = wrapper_of(self.out.code, "vnorm")
        in_order(self, w, ["doublereal r;", "r = vnorm_(v1);", 'tlog_d(">result", &r, 1)',
                           "return r;"])

    def test_write_mark(self):
        w = wrapper_of(self.out.code, "dafwdr")
        in_order(self, w, ['tlog_begin("DAFWDR")', "tlog_writer()", 'tlog_i("<handle", handle, 1)',
                           "dafwdr_(handle, recno, drec)", "tlog_mark()", "tlog_end()"])
        self.assertNotIn("drec\", ", w)

    def test_hooks_call_cspice(self):
        # only the test program's calls are redirected to the wrappers (objcopy), so the
        # hooks call CSPICE's routines by their own names, wrapped or not
        hook = self.out.code[self.out.code.index("void tlog_hook_dafus"):]
        self.assertIn("dafus_((doublereal *)sum, &nd, &ni, dc, ic);", hook[:hook.index("}")])
        hook = self.out.code[self.out.code.index("logical tlog_hook_failed"):]
        self.assertIn("return failed_();", hook[:hook.index("}")])
        self.assertNotIn("__real_", self.out.code)

    def test_nested_passthrough(self):
        w = wrapper_of(self.out.code, "spkx")
        in_order(self, w, ["tlog_enter();", "if (tlog_depth() > 1 || !tlog_on()) {",
                           "r = spkx_(targ, et, starg, targ_len);", "tlog_leave();", "return r;",
                           'tlog_begin("SPKX")'])

    def test_wrapped_list(self):
        self.assertIn("spkx", self.out.wrapped)
        self.assertIn("dafwdr", self.out.wrapped)
        self.assertNotIn("failed", self.out.wrapped)

    def test_skips_with_reason(self):
        skips = dict(self.out.skips)
        self.assertIn("chfun", skips)
        self.assertIn("CHARACTER", skips["chfun"])
        self.assertIn("unk", skips)
        self.assertIn("length", skips["unk"])
        self.assertNotIn(" __wrap_chfun_(", self.out.code)
        self.assertNotIn(" __wrap_unk_(", self.out.code)

    def test_renamed_symbol(self):
        # SpiceZrnm.h renames a few f2c symbols (vnorm_ is zz_vnorm_): the wrapper takes the
        # symbol, the log keeps the routine's name
        g = G.generate(self.protos, TABLE, {"vnorm"}, G.parse_renames("   #define   vnorm_    zz_vnorm_\n"))
        w = wrapper_of(g.code, "zz_vnorm")
        in_order(self, w, ['tlog_begin("VNORM")', "r = zz_vnorm_(v1);"])
        self.assertIn("zz_vnorm_", g.symbols)

    def test_lists(self):
        in_order(self, self.out.code, ["int tlog_is_pure(const char *name)", '"VNORM"'])
        in_order(self, self.out.code, ["int tlog_is_query(const char *name)", '"SPKEZR"'])


PGM = """\
   #include "tutils_c.h"
   int main( int argc, char **argv )
   {
   int f_ab__ ( logical * plain_ok );
   int f_spk08__ ( logical * plain_ok );
   void f_xyz_c ( SpiceBoolean * spice_ok );
   logical plain_ok;
   putcml_c( argc, argv );
   tsetup_c ( cmline, "cspice{0-9}{0-9}.log", "1.0.0" );
   /*
   Invoke test families generated from tspice.
   */
   f_ab__ ( &plain_ok );
   f_spk08__ ( &plain_ok );
   f_xyz_c ( &spice_ok );
   tclose_c();
   return ( 0 );
   }
"""


class DriverTests(unittest.TestCase):
    def test_driver(self):
        d = G.driver(PGM)
        in_order(self, d, ["static int naif_main( int argc, char **argv )", "putcml_c( argc, argv );",
                           "tlog_dispatch ( &plain_ok );", "tclose_c();"])
        self.assertNotIn("f_xyz_c ( &spice_ok )", d)
        in_order(self, d, ['"F_AB"', '"F_SPK08"'])
        in_order(self, d, ["case 0: f_ab__ ( plain_ok ); break;", "case 1: f_spk08__ ( plain_ok ); break;"])
        in_order(self, d, ["int main(int argc, char **argv)", "--list", "naif_main(1, argv)"])
        # the status: any failed test case (T_ANYBAD), as PEPPER's driver; the flag of the
        # C-API families, which the driver does not run, is set
        in_order(self, d, ["spice_ok = SPICETRUE;", "tlog_dispatch ( &plain_ok );"])
        in_order(self, d, ["naif_main(1, argv);", "t_anybad__(&bad);", "return bad ? 1 : 0;"])


if __name__ == "__main__":
    unittest.main()
