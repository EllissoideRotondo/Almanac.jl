# DAF files and SPK segments

SPK, CK and binary PCK files are DAF (double precision array file) files. A file handle
is a [`DAFHandle`](@ref) object, not an integer. The calls that continue an array
search or an array being added ([`daffna`](@ref), [`daffpa`](@ref), [`dafgs`](@ref),
[`dafgn`](@ref), [`dafws`](@ref), [`dafrs`](@ref), [`dafrn`](@ref), [`dafada`](@ref),
[`dafena`](@ref)) take the handle, which the C calls do not, because the handle holds
the state of the search.

The identifier of an SPK or CK segment (`segid`, `ident`) is the name of its DAF array,
which [`dafgn`](@ref) returns.

## DAF files

```@docs
DAFHandle
dafopr
dafopw
dafonw
dafops
save_kernel
dafcls
dafhsf
dafrfr
dafbfs
dafbbs
daffna
daffpa
dafgs
dafgn
dafps
dafus
dafgda
dafgsr
dafbna
dafada
dafena
dafws
dafrs
dafrn
dafac
dafec
dafdc
```

## SPK segments

A loaded file has no handle, so [`spksfs`](@ref) returns the path of the file, which
[`spkpvn`](@ref) takes with the same context.

```@docs
spksfs
spkpvn
spkuds
spksub
```
