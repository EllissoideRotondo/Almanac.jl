#define _DEFAULT_SOURCE
/* Test of tspice_log.c with stub hooks: no CSPICE. test_log.sh builds and runs it, then
   compares the log it writes with the expected text. */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "tspice_log.h"

/* the hooks that the generated wrappers define: nothing failed, no file is known */
logical tlog_hook_failed(void) { return 0; }
void tlog_hook_getmsg_short(char *out, ftnlen len) { memset(out, ' ', len); }
void tlog_hook_zzddhnfo(integer h, char *fname, ftnlen len, integer *arc, integer *bff,
                        integer *amh, logical *found) { *found = 0; }
void tlog_hook_dafhsf(integer h, integer *nd, integer *ni) { *nd = 0; *ni = 0; }
void tlog_hook_dafus(const doublereal *sum, integer nd, integer ni, doublereal *dc, integer *ic) {}
void tlog_hook_dafgh(integer *h) { *h = 0; }
void tlog_hook_reset(void) {}
int tlog_is_pure(const char *name) { return strcmp(name, "VNORM") == 0; }
int tlog_is_query(const char *name) { return strcmp(name, "BODN2C") == 0; }

static int fails = 0;

static void check_fmt(double x, const char *want) {
    char out[32];
    tlog_fmt_double(x, out);
    if (strcmp(out, want) == 0) {
        printf("ok   fmt %s\n", want);
    } else {
        printf("FAIL fmt [%s], want [%s]\n", out, want);
        fails++;
    }
}

/* one call as a wrapper makes it: begin, inputs, mark, outputs, end */
static void pure_call(double a) {
    doublereal v[3] = {a, 2.0, 3.0}, r = 1.0;
    tlog_enter();
    tlog_begin("VNORM");
    tlog_d("<v1", v, 3);
    tlog_mark();
    tlog_d(">result", &r, 1);
    tlog_end();
    tlog_leave();
}

static void query_call(void) {
    integer code = 399;
    logical found = 1;
    tlog_enter();
    tlog_begin("BODN2C");
    tlog_c("<name", "EARTH   ", 8, 1);
    tlog_mark();
    tlog_i(">code", &code, 1);
    tlog_l(">found", &found, 1);
    tlog_end();
    tlog_leave();
}

static void state_call(void) {
    tlog_enter();
    tlog_begin("CLPOOL");
    tlog_mark();
    tlog_end();
    tlog_leave();
}

int main(void) {
    static doublereal big[2001];
    int k;
    check_fmt(1.0, " 1.00000000000000000E+000");
    check_fmt(-2.5e-300, "-2.49999999999999998E-300");
    check_fmt(0.0, " 0.00000000000000000E+000");
    check_fmt(-0.0, "-0.00000000000000000E+000");
    check_fmt(NAN, "                      NaN");
    check_fmt(INFINITY, "                 Infinity");
    check_fmt(-INFINITY, "                -Infinity");
    check_fmt(1.0e300, " 1.00000000000000005E+300");
    check_fmt(4.9406564584124654e-324, " 4.94065645841246544E-324");
    check_fmt(0.1, " 1.00000000000000006E-001");

    pure_call(1.0);      /* written */
    pure_call(1.0);      /* a repeat: not written */
    pure_call(5.0);      /* a second distinct call, past TSPICE_MAXPURE=1: not written */
    query_call();        /* written */
    query_call();        /* the same state: not written */
    state_call();        /* written; changes the state */
    query_call();        /* written again */

    tlog_enter();        /* a string with a line feed, and a long array */
    tlog_begin("TKVRSN");
    tlog_c("<item", "A\nB ", 4, 1);
    for (k = 0; k < 2001; k++) big[k] = k + 1;
    tlog_d("<x", big, 2001);
    tlog_mark();
    tlog_end();
    tlog_leave();

    tlog_enter();        /* a nested call is the wrapper's business: depth 2 */
    tlog_enter();
    if (tlog_depth() != 2) { printf("FAIL depth %d\n", tlog_depth()); fails++; }
    tlog_leave();
    tlog_leave();

    printf("%s\n", fails ? "FAILED" : "PASSED");
    return fails != 0;
}
