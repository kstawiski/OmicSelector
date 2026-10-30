# Fingerprint frozen fitted state, including Python parameters and buffers

Computes SHA-256 of R state and deterministic content snapshots of
nested Python objects. Python snapshots include estimator preprocessing,
training context, neural parameters, buffers and evaluation flags.
Snapshots use numerical content rather than shared-object identities or
pickle bytes. Autograd ownership, TabICL checkpoint paths and temporary
inference workspaces are excluded; frozen context caches and estimator
inference settings remain included. Fingerprints from earlier package
versions must be regenerated. Unknown opaque state is rejected, rather
than hashed as a pointer placeholder.

## Usage

``` r
singlesample_state_fingerprint(model)
```

## Arguments

- model:

  A fitted model or single-sample deployment object.

## Value

One SHA-256 string. Equality tests detect state changes, not validity of
the fitted model or clinical performance.
