"""The lists of routines and arguments that decide how a call is recorded: which calls a
family may log once, which file arguments are copied, which writers are only marked, and
the lengths that the arguments do not give. instrument_c.py and signatures.py read them;
the format of the log is described in tspice_log.h.
"""

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
wncard wncomd wncond wndifd wnelmd wnexpd wnextd wnfetd wnfild wnfltd wnincd wninsd wnintd
wnreld wnsumd wnunid wnvald
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

# Routines whose file arguments are copied when the call is recorded (the files of a
# family change and disappear while it runs): the snapshot of the working directory is
# taken before the call.
FILE_ROUTINES = {"furnsh", "ldpool", "spkcov", "spkobj", "ckcov", "ckobj", "pckcov", "pckfrm",
                 "dskobj", "dsksrf", "getfat", "spklef", "pcklof", "cklpf", "eklef",
                 "zzdsklsf"}
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
# Integer arguments of a marked writer that the replay maps to its files.
MARK_HANDLES = {"handle", "inhan", "han", "handl", "newh"}

# Assumed-size arrays counted by an argument other than the first count of the routine:
# the plates by NP (DSKRB2 and DSKMI2 also take NV), and the coordinate parameters,
# which have the fixed size NSYPAR.
NAMED_COUNTS = {"plates": "np", "vrtces": "nv", "qpnams": "qnpars", "qcpars": "qnpars", "qdpars": "qnpars",
                "qipars": "qnpars", "qlpars": "qnpars"}
FIXED_SIZES = {"corpar": "10"}

# The routines whose cells are logged: the GF searches, the window routines, and the
# routines that give the coverage or the objects of a file or the frames of a class.
# GFUDB, GFUDS, UDDC and UDDF take functions of the test program, which the replay
# cannot call.
CELL_ROUTINES = set("""
gfdist gfilum gfoclt gfpa gfposc gfrfov gfrr gfsep gfsntc gfsubc gftfov gfevnt gffove gfocce
wncard wncomd wncond wndifd wnelmd wnexpd wnextd wnfetd wnfild wnfltd wnincd wninsd wnintd
wnreld wnsumd wnunid wnvald
spkcov spkobj ckcov ckobj pckcov pckfrm dskobj dsksrf bltfrm kplfrm
""".split())
# Routines that take procedures: GFEVNT, GFFOVE and GFOCCE take the step, the
# refinement, the progress report and the interrupt check. Every call of tspice passes
# SPICELIB's GFSTEP and GFREFN (and report and interrupt routines, which do not change
# the result); the replay passes the step that GFSSTP set and Almanac's gfrefn. The
# procedures are not logged.
PROC_ROUTINES = {"gfevnt", "gffove", "gfocce"}
