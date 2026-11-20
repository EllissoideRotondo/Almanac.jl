/* The call log of the tspice_c recording (test/tspice/cspice): the C port of the
   Fortran module tspice_log that instrument.py writes for the PEPPER recording, with the
   same output (the format is described in LOGMOD of instrument.py). The wrappers that
   instrument_c.py generates call these functions. */

#ifndef TSPICE_LOG_H
#define TSPICE_LOG_H

#include "f2c.h"

#define TLOG_HUGE 2147483647L   /* huge(0) of a 32-bit integer: an array the log does not hold */

/* the state of the logger */
int tlog_on(void);
int tlog_depth(void);
void tlog_enter(void);
void tlog_leave(void);

/* a call: begin, inputs, mark, outputs, end */
void tlog_begin(const char *name);
void tlog_mark(void);
void tlog_end(void);

/* arguments: label is "<arg" for an input, ">arg" for an output */
void tlog_d(const char *label, const doublereal *x, long n);
void tlog_i(const char *label, const integer *x, long n);
void tlog_l(const char *label, const logical *x, long n);
void tlog_c(const char *label, const char *x, ftnlen len, long n);
void tlog_fmt_double(double x, char out[32]);

/* files */
void tlog_snapshot(void);
void tlog_snapshot_file(const char *name, ftnlen len);
void tlog_close_name(integer handle);
void tlog_close_snapshot(void);
void tlog_handle_file(integer handle);
void tlog_writer(void);

/* custom arguments and lengths */
void tlog_dafsum(integer handle, const doublereal *sum, char dir);
void tlog_dafsum_cur(const doublereal *sum, char dir);
void tlog_ek_note(integer handle, integer segno, integer ncols, integer nrows,
                  const char *cnames, ftnlen cnames_len, const char *decls, ftnlen decls_len);
long tlog_ek_nvals(integer handle, integer segno, const char *column, ftnlen column_len,
                   const integer *entszs);
long tlog_ek_nent(integer handle, integer segno, const char *column, ftnlen column_len);
long tlog_ek_nnul(integer handle, integer segno, const char *column, ftnlen column_len);
void tlog_spk14_note(integer handle, integer chbdeg);
long tlog_spk14_len(integer handle, integer ncsets);
long tlog_psum(const integer *npkts, const integer *subtps, long n, int s0, int s1, int s2, int s3);
long tlog_okn(long n);
long tlog_nsub(long nchars, long bpos, long epos, long slen);
long tlog_sum_i(const integer *x, long n);
long tlog_nint(double x);
long tlog_max0(long x);
long tlog_clamp(long n);                          /* max(0, min(n, 1000000)) */
long tlog_card_d(const doublereal *cell);         /* the cardinality of a cell, within 0:10**6 */
long tlog_card_i(const integer *cell);

/* defined by the generated code (tlog_wrap.c): CSPICE's routines (the wrappers take only
   the test program's calls), and the lists PURE and QUERY of instrument.py */
logical tlog_hook_failed(void);
void tlog_hook_getmsg_short(char *out, ftnlen len);
void tlog_hook_zzddhnfo(integer h, char *fname, ftnlen len, integer *arc, integer *bff,
                        integer *amh, logical *found);
void tlog_hook_dafhsf(integer h, integer *nd, integer *ni);
void tlog_hook_dafus(const doublereal *sum, integer nd, integer ni, doublereal *dc, integer *ic);
void tlog_hook_dafgh(integer *h);
void tlog_hook_reset(void);
int tlog_is_pure(const char *name);
int tlog_is_query(const char *name);

#endif
