# States and aberration

## Positions and states

```@docs
spkezr
spkpos
spkgeo
spkssb
spkez
spkezp
spkgps
```

## Derivatives of any order

```@docs
spkder
pxder
```

## Aberration corrections

```@docs
Aberration
spkltc
spkaps
spkacs
spkapo
ltime
stelab
stlabx
```

## Observers and targets with a constant state

```@docs
spkcpo
spkcpt
spkcvo
spkcvt
```

## Prepared state queries

```@docs
state_query
StateQuery
```

## Bodies defined by Julia functions

```@docs
add_ephemeris!
load_samples!
load_function!
```

## Tables for inner loops

```@docs
EphemerisTable
ephemeris_table
table_position
table_state
table_positions!
table_states!
table_jet
table_span
RotationTable
rotation_table
table_rotation
table_rotation_jet
table_to_spk
```
