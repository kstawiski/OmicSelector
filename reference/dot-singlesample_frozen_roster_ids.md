# Frozen Amendment \#4 roster identifiers

Independent enumeration of the 74 pre-specified method identifiers,
transcribed from the prespecification (SAP Amendment \#4 sections
4.1-4.4 and the round-2 survey, section B; plus Amendment \#5's
\`lrt-tcopula\`). Kept separate from the shipped manifest so a
completeness test can assert the manifest implements exactly this frozen
set (\`setequal\`): editing the CSV without editing this vector (or vice
versa) fails the gate.

## Usage

``` r
.singlesample_frozen_roster_ids()
```

## Value

Character vector of 74 method identifiers.
