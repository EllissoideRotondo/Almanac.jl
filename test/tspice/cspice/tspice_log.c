/* The call log of the tspice_c recording: a port of the Fortran module tspice_log of
   instrument.py (LOGMOD), routine by routine, with the same output. See tspice_log.h.

   The log goes to the file named by TSPICE_LOG (a FIFO that record.sh compresses);
   without it nothing is written. The depth of the calls is kept by the wrappers
   (tlog_enter, tlog_leave): only a call at depth 1, from the test program, is logged. */

#define _DEFAULT_SOURCE
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "tspice_log.h"

#define MAXLOG 2000                 /* longer arrays are logged as N = -1, not replayed */
#define MAXCOPY "268435456"          /* bytes of file copies per family (256 MB) */
#define NEK 4000
#define NS14 64
#define NCAP 400
#define FNV_OFFSET 0xcbf29ce484222325ULL
#define FNV_PRIME 0x100000001b3ULL

typedef unsigned long long u64;

static FILE *logfile = NULL;
static int checked = 0, opened = 0, keepall = 0, depth = 0, nsnap = 0, bincheck = 0;
static long nwritten = 0, maxpure = 20000;
static char logpath[1024] = "";
static char closename[1024] = "";
static long long binbytes = 0, maxbin = 268435456LL;
static char callname[33] = "";
static char *buf = NULL;
static size_t bcap = 0, blen = 0, inend = 0;
static u64 *seen = NULL;
static size_t nseen = 0, mseen = 0;
static u64 state = 0;               /* calls of routines that can change the state */
static char capname[NCAP][33];
static long capcount[NCAP];
static int ncapped = 0;
static integer hfh[4096];           /* tlog_handle_file: handle, hash of its file */
static u64 hff[4096];
static integer ekh[NEK], eks[NEK], eksz[NEK], ekrows[NEK];
static int eknul[NEK], nek = 0;
static char ekcol[NEK][33];
static integer s14h[NS14], s14d[NS14];

/* ---- the state ---- */

int tlog_on(void) {
    const char *v;
    if (!checked) {
        checked = 1;
        v = getenv("TSPICE_LOG");
        if (v != NULL && v[0] != '\0') {
            strncpy(logpath, v, sizeof logpath - 1);
            logfile = fopen(logpath, "wb");
            if (logfile != NULL) {
                opened = 1;
                bcap = 1048576;
                buf = malloc(bcap);
                mseen = 65536;
                seen = calloc(mseen, sizeof(u64));
                v = getenv("TSPICE_ALL");
                keepall = v != NULL && v[0] == '1';
                v = getenv("TSPICE_MAXPURE");
                if (v != NULL && v[0] != '\0') maxpure = atol(v);
            }
        }
    }
    return opened;
}

int tlog_depth(void) { return depth; }
void tlog_enter(void) { depth++; }
void tlog_leave(void) { depth--; }

static void emit(const char *s) {
    size_t n = strlen(s), need = blen + n + 1;
    if (need > bcap) {
        bcap = 2 * bcap > need ? 2 * bcap : need;
        buf = realloc(buf, bcap);
    }
    memcpy(buf + blen, s, n);
    buf[blen + n] = '\n';
    blen = need;
}

/* the length of a Fortran string without its trailing blanks */
static size_t len_trim(const char *x, size_t n) {
    while (n > 0 && x[n - 1] == ' ') n--;
    return n;
}

/* a Fortran string, trimmed, as a C string (truncated to size - 1) */
static void ctrim(char *out, size_t size, const char *x, size_t n) {
    n = len_trim(x, n);
    if (n > size - 1) n = size - 1;
    memcpy(out, x, n);
    out[n] = '\0';
}

static u64 fnv_bytes(u64 h, const unsigned char *s, size_t n) {
    size_t k;
    for (k = 0; k < n; k++) h = (h ^ s[k]) * FNV_PRIME;
    return h;
}

static u64 fnv(const char *s, size_t n) {
    u64 h = fnv_bytes(FNV_OFFSET, (const unsigned char *)s, n);
    return h == 0 ? 1 : h;
}

/* Add the key to the set of the inputs seen; false if it was there. */
static int insert(u64 key) {
    size_t k, m;
    if (2 * (nseen + 1) > mseen) {
        u64 *old = seen;
        size_t mold = mseen;
        mseen = 2 * mold;
        seen = calloc(mseen, sizeof(u64));
        nseen = 0;
        for (k = 0; k < mold; k++)
            if (old[k] != 0) insert(old[k]);
        free(old);
    }
    m = mseen;
    k = (size_t)(key % m);
    while (seen[k] != 0) {
        if (seen[k] == key) return 0;
        k = (k + 1) % m;
    }
    seen[k] = key;
    nseen++;
    return 1;
}

/* Count a new call of the routine without state NAME; false once the family has made
   maxpure distinct ones (TSPICE_MAXPURE, default 20000; 0: no limit). */
static int below_cap(const char *name) {
    int k;
    for (k = 0; k < ncapped; k++)
        if (strcmp(capname[k], name) == 0) break;
    if (k >= ncapped) {
        ncapped = ncapped + 1 < NCAP ? ncapped + 1 : NCAP;
        k = ncapped - 1;
        strcpy(capname[k], name);
        capcount[k] = 0;
    }
    capcount[k]++;
    return maxpure <= 0 || capcount[k] <= maxpure;
}

/* ---- calls ---- */

void tlog_begin(const char *name) {
    char ln[80];
    if (!tlog_on()) return;
    if (depth == 1) {
        blen = 0;
        inend = 0;
        strncpy(callname, name, 32);
        callname[32] = '\0';
    }
    sprintf(ln, "@ %s %d %c", name, depth, tlog_hook_failed() ? 'T' : 'F');
    emit(ln);
}

/* The inputs of the call are logged: they end here. */
void tlog_mark(void) {
    if (!tlog_on()) return;
    if (depth == 1) inend = blen;
}

void tlog_end(void) {
    char msg[65], ln[80];
    int failed, pure, query;
    if (!tlog_on()) return;
    failed = tlog_hook_failed() != 0;
    msg[0] = '\0';
    if (failed) {
        char raw[64];
        tlog_hook_getmsg_short(raw, 64);
        ctrim(msg, sizeof msg, raw, 64);
    }
    sprintf(ln, "! %c %s", failed ? 'T' : 'F', msg);
    emit(ln);
    if (depth > 1) return;
    pure = tlog_is_pure(callname);
    query = tlog_is_query(callname);
    if (!keepall && (pure || query)) {
        if (inend == 0) inend = blen;
        /* a query repeats only in the same state: since the last call of another routine */
        if (!insert(fnv(buf, inend) ^ (query ? state : 0))) return;     /* a repeat */
        /* a family that has made enough distinct calls of a routine without state */
        if (pure && !below_cap(callname)) return;
    } else {
        state++;
    }
    fwrite(buf, 1, blen, logfile);
    nwritten++;
    if (nwritten % 256 == 0) fflush(logfile);
}

/* ---- files ---- */

static void shell(const char *cmd) {
    if (system(cmd) == -1) perror("tspice_log");
}

/* Link the files of the working directory into a new snapshot directory. If no file
   changed since the last snapshot (no newer file, the same count), the new one is a
   symbolic link to it, relative so that the recording can move. */
void tlog_snapshot(void) {
    char dir[1100], prev[1100], cmd[4000];
    if (!tlog_on()) return;
    nsnap++;
    sprintf(dir, "%s.files/%d", logpath, nsnap);
    sprintf(prev, "%s.files/%d", logpath, nsnap - 1);
    sprintf(cmd, "d=\"%s\"; p=\"%s\"; "
            "n=$(ls -A . | wc -l); if [ -f \"$p/.stamp\" ] && [ \"$(cat \"$p/.count\")\" = \"$n\" ] "
            "&& [ -z \"$(find . -maxdepth 1 -type f -newer \"$p/.stamp\" -print -quit)\" ]; then "
            "ln -s \"$(basename \"$(readlink -f \"$p\")\")\" \"$d\"; else mkdir -p \"$d\" && ln -f ./* \"$d/\" 2>/dev/null; "
            "echo \"$n\" > \"$d/.count\"; touch \"$d/.stamp\"; fi", dir, prev);
    shell(cmd);
    sprintf(cmd, "S %s", dir);
    emit(cmd);
}

/* Copy the file NAME into a new snapshot directory, at a call from the test program. A
   file equal to its last copy is linked to that copy; past MAXCOPY bytes of copies in a
   family the file is linked, not copied. */
void tlog_snapshot_file(const char *name, ftnlen len) {
    char dir[1100], fname[1024], cmd[5000];
    if (!tlog_on() || depth != 1) return;
    nsnap++;
    sprintf(dir, "%s.files/%d", logpath, nsnap);
    ctrim(fname, sizeof fname, name, len);
    sprintf(cmd, "d=\"%s\"; f=\"%s\"; "
            "l=\"%s.files/latest\"; mkdir -p \"$d\" \"$l\"; b=$(basename \"$f\"); "
            "if [ -f \"$f\" ]; then if cmp -s \"$f\" \"$l/$b\"; then ln \"$l/$b\" \"$d/\"; "
            "else t=$(cat \"$l/.bytes\" 2>/dev/null || echo 0); if [ \"$t\" -gt " MAXCOPY " ]; then "
            "ln -f \"$f\" \"$d/\"; else cp \"$f\" \"$d/\" && ln -f \"$d/$b\" \"$l/$b\" && "
            "echo $((t + $(stat -c %%s \"$f\"))) > \"$l/.bytes\"; fi; fi; fi 2>/dev/null; true",
            dir, fname, logpath);
    shell(cmd);
    sprintf(cmd, "S %s", dir);
    emit(cmd);
}

/* Write `H HANDLE PATH` into the log, at once and outside the call, when a logged call
   takes a handle of a file open for reading whose file the log has not named. */
void tlog_handle_file(integer handle) {
    char raw[1024], fname[1025];
    integer arc = 0, bff = 0, amh = 0;
    logical found = 0;
    int k;
    u64 h;
    if (!tlog_on()) return;
    if (tlog_hook_failed()) return;
    tlog_hook_zzddhnfo(handle, raw, 1024, &arc, &bff, &amh, &found);
    if (tlog_hook_failed() || !found || amh != 1) return;
    ctrim(fname, sizeof fname, raw, 1024);
    h = fnv(fname, strlen(fname));
    k = ((handle % 4096) + 4096) % 4096;
    if (hfh[k] == handle && hff[k] == h) return;
    hfh[k] = handle;
    hff[k] = h;
    fprintf(logfile, "H %d %s\n", (int)handle, fname);
}

/* Note the file of HANDLE before a close from the test program, if the file is open for
   write (or new): tlog_close_snapshot copies it after the close. */
void tlog_close_name(integer handle) {
    char raw[1024];
    integer arc = 0, bff = 0, amh = 0;
    logical found = 0;
    closename[0] = '\0';
    if (!tlog_on() || depth != 1) return;
    if (tlog_hook_failed()) return;
    tlog_hook_zzddhnfo(handle, raw, 1024, &arc, &bff, &amh, &found);
    if (found && (amh == 2 || amh == 4)) ctrim(closename, sizeof closename, raw, 1024);
}

/* Copy the file of the close into a new snapshot directory. Past MAXCOPY bytes of copies
   the copy of the last close of the same file name is removed, so that only the last
   state of each file is kept. */
void tlog_close_snapshot(void) {
    char dir[1100], cmd[5000];
    if (!tlog_on() || depth != 1) return;
    if (closename[0] == '\0') return;
    nsnap++;
    sprintf(dir, "%s.files/%d", logpath, nsnap);
    sprintf(cmd, "d=\"%s\"; f=\"%s\"; "
            "l=\"%s.files/latest\"; mkdir -p \"$d\" \"$l\"; b=$(basename \"$f\"); "
            "if [ -f \"$f\" ]; then if cmp -s \"$f\" \"$l/$b\"; then ln \"$l/$b\" \"$d/\"; "
            "else t=$(cat \"$l/.bytes\" 2>/dev/null || echo 0); if [ \"$t\" -le " MAXCOPY " ]; then "
            "cp \"$f\" \"$d/\" && ln -f \"$d/$b\" \"$l/$b\" && "
            "echo $((t + $(stat -c %%s \"$f\"))) > \"$l/.bytes\"; "
            "else p=$(cat \"$l/.last.$b\" 2>/dev/null); [ -n \"$p\" ] && rm -f \"$p/$b\"; "
            "cp \"$f\" \"$d/\" && echo \"$d\" > \"$l/.last.$b\"; fi; fi; fi 2>/dev/null; true",
            dir, closename, logpath);
    shell(cmd);
    sprintf(cmd, "S %s", dir);
    emit(cmd);
    closename[0] = '\0';
}

/* The mark of a writer that is not recorded (WRITE_MARKS of instrument.py). */
void tlog_writer(void) {
    if (!tlog_on()) return;
    emit("W");
}

/* ---- arguments ---- */

static void head(const char *label, char kind, long n) {
    char ln[80];
    sprintf(ln, "%s %c %ld", label, kind, n);
    emit(ln);
}

/* gfortran's es25.17e3, trailing blanks removed (none: the field is full) */
void tlog_fmt_double(double x, char out[32]) {
    char m[40], *e;
    int ex;
    if (isnan(x)) {
        sprintf(out, "%25s", "NaN");
        return;
    }
    if (isinf(x)) {
        sprintf(out, "%25s", x > 0 ? "Infinity" : "-Infinity");
        return;
    }
    sprintf(m, "%.17E", x);
    e = strchr(m, 'E');
    ex = atoi(e + 1);
    *e = '\0';
    sprintf(out, "%s%sE%c%03d", m[0] == '-' ? "" : " ", m, ex < 0 ? '-' : '+', ex < 0 ? -ex : ex);
}

/* True if the array of N values of SIZE bytes, too long for the log, goes to a binary
   file: an argument of a call from the test program, of at most 50 million values and
   2 GB, while the family has written at most TSPICE_MAXBIN bytes of such files. */
static int bin_ok(long n, long size) {
    if (!bincheck) {
        char cmd[1200];
        const char *v = getenv("TSPICE_MAXBIN");
        bincheck = 1;
        if (v != NULL && v[0] != '\0') maxbin = atoll(v);
        sprintf(cmd, "mkdir -p \"%s.files/bin\"", logpath);
        shell(cmd);
    }
    return depth == 1 && n <= 50000000L && (long long)size * n <= 2000000000LL && binbytes <= maxbin;
}

/* The path of the binary file of values whose hash is H, of kind K; true if it exists. */
static int bin_path(u64 h, const char *k, char *path) {
    char gz[1300];
    struct stat st;
    sprintf(path, "%s.files/bin/%016llX.%s", logpath, h, k);
    sprintf(gz, "%s.gz", path);
    return stat(gz, &st) == 0;
}

/* Write the bytes to the binary file of hash H (gzip), unless it exists, and log its
   path in a `B` line. */
static void bin_write(u64 h, const char *kind, const void *x, size_t nbytes) {
    char path[1200], cmd[1300];
    struct stat st;
    if (!bin_path(h, kind, path)) {
        FILE *u = fopen(path, "wb");
        if (u != NULL) {
            fwrite(x, 1, nbytes, u);
            fclose(u);
        }
        sprintf(cmd, "gzip -1 -f \"%s\"", path);
        shell(cmd);
        sprintf(cmd, "%s.gz", path);
        if (stat(cmd, &st) == 0 && st.st_size > 0) binbytes += st.st_size;
    }
    sprintf(cmd, "B %s.gz", path);
    emit(cmd);
}

void tlog_d(const char *label, const doublereal *x, long n) {
    long k;
    char ln[32];
    if (!tlog_on()) return;
    if (n > MAXLOG) {
        if (bin_ok(n, 8)) {
            /* FNV-1a over the 8 bytes of each value, low byte first */
            u64 h = FNV_OFFSET, w;
            int b;
            for (k = 0; k < n; k++) {
                memcpy(&w, &x[k], 8);
                for (b = 0; b < 8; b++) h = (h ^ ((w >> (8 * b)) & 255)) * FNV_PRIME;
            }
            bin_write(h, "d", x, (size_t)n * 8);
        }
        head(label, 'd', -1);
        return;
    }
    head(label, 'd', n);
    for (k = 0; k < n; k++) {
        tlog_fmt_double(x[k], ln);
        emit(ln);
    }
}

void tlog_i(const char *label, const integer *x, long n) {
    long k;
    char ln[24];
    if (!tlog_on()) return;
    if (n > MAXLOG) {
        if (bin_ok(n, 4)) {
            u64 h = FNV_OFFSET, w;
            int b;
            for (k = 0; k < n; k++) {
                w = (u64)(long long)x[k];
                for (b = 0; b < 4; b++) h = (h ^ ((w >> (8 * b)) & 255)) * FNV_PRIME;
            }
            bin_write(h, "i4", x, (size_t)n * 4);
        }
        head(label, 'i', -1);
        return;
    }
    head(label, 'i', n);
    for (k = 0; k < n; k++) {
        sprintf(ln, "%d", (int)x[k]);
        emit(ln);
    }
}

void tlog_l(const char *label, const logical *x, long n) {
    long k;
    if (!tlog_on()) return;
    if (n > MAXLOG) {
        head(label, 'l', -1);
        return;
    }
    head(label, 'l', n);
    for (k = 0; k < n; k++) emit(x[k] ? "T" : "F");
}

/* A string as one line of the log: `|` and the string without its trailing blanks, or,
   if it holds a character outside printable ASCII, `%` and the hexadecimal codes. */
static void emit_str(const char *x, ftnlen len) {
    size_t n = len_trim(x, (size_t)len), k;
    char *line = malloc(2 * n + 2);
    int plain = 1;
    for (k = 0; k < n; k++)
        if ((unsigned char)x[k] < 32 || (unsigned char)x[k] > 126) plain = 0;
    if (plain) {
        line[0] = '|';
        memcpy(line + 1, x, n);
        line[n + 1] = '\0';
    } else {
        line[0] = '%';
        for (k = 0; k < n; k++) sprintf(line + 1 + 2 * k, "%02x", (unsigned char)x[k]);
        line[2 * n + 1] = '\0';
    }
    emit(line);
    free(line);
}

void tlog_c(const char *label, const char *x, ftnlen len, long n) {
    long k;
    if (!tlog_on()) return;
    if (n > MAXLOG) {
        if (bin_ok(n, len)) {
            char kind[16];
            sprintf(kind, "c%d", (int)len);
            bin_write(fnv_bytes(FNV_OFFSET, (const unsigned char *)x, (size_t)n * len), kind, x,
                      (size_t)n * len);
        }
        head(label, 'c', -1);
        return;
    }
    head(label, 'c', n);
    for (k = 0; k < n; k++) emit_str(x + k * len, len);
}

/* The summary SUM of a DAF array as its double and integer components, `<dc` and `<ic`
   (DIR '>': `>dc` and `>ic`): the decimal log does not keep the bits of the integers
   packed into doubles. No components for a handle that is not of a DAF. */
void tlog_dafsum(integer handle, const doublereal *sum, char dir) {
    char raw[1024], label[8];
    integer arc = 0, bff = 0, amh = 0, nd = 0, ni = 0;
    logical found = 0;
    doublereal dc[125];
    integer ic[250];
    if (!tlog_on()) return;
    if (!tlog_hook_failed()) tlog_hook_zzddhnfo(handle, raw, 1024, &arc, &bff, &amh, &found);
    if (found && arc == 1) tlog_hook_dafhsf(handle, &nd, &ni);
    if (tlog_hook_failed() || nd < 0 || nd > 124 || ni < 2 || ni > 250) {
        nd = 0;
        ni = 0;
    }
    if (nd + ni > 0) tlog_hook_dafus(sum, nd, ni, dc, ic);
    sprintf(label, "%cdc", dir);
    tlog_d(label, dc, nd);
    sprintf(label, "%cic", dir);
    tlog_i(label, ic, ni);
}

/* The summary of DAFWS or DAFRS, or with DIR '>' the one that DAFGS returned, for the DAF
   of the current search (DAFGH). */
void tlog_dafsum_cur(const doublereal *sum, char dir) {
    integer h = 0;
    if (!tlog_on()) return;
    if (!tlog_hook_failed()) {
        tlog_hook_dafgh(&h);
        /* no search: the call itself signals; this lookup leaves no error behind */
        if (tlog_hook_failed()) {
            tlog_hook_reset();
            h = 0;
        }
    }
    tlog_dafsum(h, sum, dir);
}

/* ---- lengths ---- */

/* N, or 0 after a failed call: the count of an output array that a failed call left unset. */
long tlog_okn(long n) {
    if (tlog_hook_failed()) return 0;
    return n > 0 ? n : 0;
}

/* Sum over the K = 1:N mini-segments of NPKTS(K) times the size of the subtype SUBTPS(K);
   huge for a subtype outside or a count past a million. */
long tlog_psum(const integer *npkts, const integer *subtps, long n, int s0, int s1, int s2, int s3) {
    long k;
    long long t = 0;
    int sz;
    if (n > 1000000) return TLOG_HUGE;
    for (k = 0; k < n; k++) {
        switch (subtps[k]) {
        case 0: sz = s0; break;
        case 1: sz = s1; break;
        case 2: sz = s2; break;
        case 3: sz = s3; break;
        default: sz = 0;
        }
        if (sz == 0) return TLOG_HUGE;
        t += (long long)(npkts[k] > 0 ? npkts[k] : 0) * sz;
    }
    return t < TLOG_HUGE ? (long)t : TLOG_HUGE;
}

/* Number of strings that NCHARS characters fill at the positions BPOS:EPOS of strings
   of length SLEN; huge for positions outside. */
long tlog_nsub(long nchars, long bpos, long epos, long slen) {
    if (bpos < 1 || epos > slen) return TLOG_HUGE;
    if (epos < bpos || nchars < 1) return 0;
    return (nchars + epos - bpos) / (epos - bpos + 1);
}

long tlog_sum_i(const integer *x, long n) {
    long k, t = 0;
    for (k = 0; k < n; k++) t += x[k];
    return t;
}

/* Fortran's NINT: half away from zero */
long tlog_nint(double x) { return x >= 0 ? (long)floor(x + 0.5) : -(long)floor(-x + 0.5); }

/* ---- EK fast loads and SPK type 14 ---- */

/* A column name as Fortran's upper32: left-adjusted, 32 characters, upper case
   (trailing blanks removed in the C string). */
static void upper32(char out[33], const char *x, size_t len) {
    size_t k = 0, n = 0;
    while (k < len && x[k] == ' ') k++;
    while (k < len && n < 32) {
        out[n++] = (char)toupper((unsigned char)x[k]);
        k++;
    }
    out[n] = '\0';
    n = len_trim(out, n);
    out[n] = '\0';
}

/* The columns of the EK segments that EKIFLD began: handle, segment, column name, entry
   size (-1: variable), whether nulls are allowed, and the rows of the segment. */
void tlog_ek_note(integer handle, integer segno, integer ncols, integer nrows,
                  const char *cnames, ftnlen cnames_len, const char *decls, ftnlen decls_len) {
    char u[513], *p;
    integer k;
    int j, l, sz;
    if (!tlog_on()) return;
    if (tlog_hook_failed()) return;
    for (k = 0; k < ncols; k++) {
        const char *d = decls + (size_t)k * decls_len;
        int n = (int)len_trim(d, (size_t)decls_len);
        l = 0;
        for (j = 0; j < n; j++) {
            if (d[j] == ' ' || l >= 512) continue;
            u[l++] = (char)toupper((unsigned char)d[j]);
        }
        u[l] = '\0';
        sz = 1;
        p = strstr(u, "SIZE=");
        if (p != NULL) {
            char item[513];
            char *c = strchr(p + 5, ',');
            size_t m = c == NULL ? strlen(p + 5) : (size_t)(c - (p + 5));
            memcpy(item, p + 5, m);
            item[m] = '\0';
            if (strcmp(item, "VARIABLE") == 0) {
                sz = -1;
            } else {
                char *end;
                long v = strtol(item, &end, 10);
                sz = (end != item && *end == '\0') ? (int)v : 1;
            }
        }
        nek = nek % NEK + 1;
        ekh[nek - 1] = handle;
        eks[nek - 1] = segno;
        upper32(ekcol[nek - 1], cnames + (size_t)k * cnames_len, (size_t)cnames_len);
        eksz[nek - 1] = sz;
        eknul[nek - 1] = strstr(u, "NULLS_OK=TRUE") != NULL;
        ekrows[nek - 1] = nrows;
    }
}

/* The note of a column of an EK fast load (1-based), or 0; the latest note first. */
static int ek_find(integer handle, integer segno, const char *column, ftnlen len) {
    char c[33];
    int j, k;
    upper32(c, column, (size_t)len);
    for (j = 0; j < NEK; j++) {
        k = (((nek - 1 - j) % NEK) + NEK) % NEK + 1;
        if (ekh[k - 1] == handle && eks[k - 1] == segno && strcmp(ekcol[k - 1], c) == 0 && ekrows[k - 1] > 0)
            return k;
    }
    return 0;
}

long tlog_ek_nvals(integer handle, integer segno, const char *column, ftnlen column_len,
                   const integer *entszs) {
    int k = ek_find(handle, segno, column, column_len);
    long long t = 0;
    integer r;
    if (k == 0) return TLOG_HUGE;
    if (eksz[k - 1] >= 0) {
        t = (long long)eksz[k - 1] * ekrows[k - 1];
    } else {
        for (r = 0; r < ekrows[k - 1]; r++) t += entszs[r] > 0 ? entszs[r] : 0;
    }
    return t < TLOG_HUGE ? (long)t : TLOG_HUGE;
}

long tlog_ek_nent(integer handle, integer segno, const char *column, ftnlen column_len) {
    int k = ek_find(handle, segno, column, column_len);
    if (k == 0) return TLOG_HUGE;
    return eksz[k - 1] < 0 ? ekrows[k - 1] : 0;
}

long tlog_ek_nnul(integer handle, integer segno, const char *column, ftnlen column_len) {
    int k = ek_find(handle, segno, column, column_len);
    if (k == 0) return TLOG_HUGE;
    return eknul[k - 1] ? ekrows[k - 1] : 0;
}

/* The degree of the type 14 segment that SPK14B began on HANDLE, for SPK14A. */
void tlog_spk14_note(integer handle, integer chbdeg) {
    int k;
    if (!tlog_on()) return;
    if (tlog_hook_failed()) return;
    for (k = 0; k < NS14; k++)
        if (s14h[k] == handle || s14h[k] == 0) break;
    if (k >= NS14) k = 0;
    s14h[k] = handle;
    s14d[k] = chbdeg;
}

long tlog_spk14_len(integer handle, integer ncsets) {
    int k;
    for (k = 0; k < NS14; k++) {
        if (s14h[k] == handle && handle != 0) {
            if (s14d[k] >= 0 && ncsets < 1000000)
                return (long)(ncsets > 0 ? ncsets : 0) * (2 + 6 * (s14d[k] + 1));
            return TLOG_HUGE;
        }
    }
    return TLOG_HUGE;
}
