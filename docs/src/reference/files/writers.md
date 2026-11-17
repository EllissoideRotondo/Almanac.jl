# Kernel writers

The writers take the arguments of CSPICE in the order of CSPICE, but leave out most
counts and flags that follow from the input. For example, in [`ckw01`](@ref) and
[`ckw03`](@ref) `avvs` is `nothing` for a segment without angular velocity, in place of
the flag `avflag`. An input of the wrong shape throws `SpiceError(:badarg)`.
[Writing kernels](@ref) shows a writer in use.

## SPK writers

```@docs
spkopn
spkopa
spkcls
spkpds
spkw01
spkw02
spkw03
spkw05
spkw08
spkw09
spkw12
spkw13
spkw15
spkw17
spkw18
spkw19
spkw20
spkw21
spk14b
spk14a
spk14e
SPK14Writer
spkw10
```

## Two-line element sets

[`spkw10`](@ref) writes element sets to an SPK segment.

```@docs
getelm
evsgp4
```

## CK writers

```@docs
ckopn
ckcls
ckw01
ckw02
ckw03
ckw05
ckw06
```

## Binary PCK writer

```@docs
pckopn
pckw02
pckw20
pckcls
```
