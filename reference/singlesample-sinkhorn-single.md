# Forced-single-sample entropic-OT (Sinkhorn) negative control

Forced \\n = 1\\, row-by-row entropic optimal-transport (Sinkhorn)
scorer (family NC), a deliberately weak **negative control** on the
single-sample discriminator roster. It scores each specimen
INDEPENDENTLY (as its own one-atom / Dirac distribution) against two
FROZEN training reference point clouds – the case cloud and the control
cloud – in per-sample robust CLR space, and reports the difference of
the two entropic-kernel transport energies oriented so that larger =
more case-like.

**Why a separate method, and what out-of-N coupling it removes.** The
package already ships a MULTI-ROW Sinkhorn-OT scorer
([`fit_sinkhorn_ot_scorer`](https://kstawiski.github.io/OmicSelector/reference/fit_sinkhorn_ot_scorer.md)
/
[`score_sinkhorn_ot_scorer`](https://kstawiski.github.io/OmicSelector/reference/score_sinkhorn_ot_scorer.md)).
That primitive is transductive / out-of-N: its batch projection
(`.sot_project_scaled`, the `nrow(z_scaled) > 1` branch) builds ONE
joint Sinkhorn transport plan whose source marginal is
`a <- rep(1 / nrow(z_scaled), nrow(z_scaled))` – it couples ALL scored
rows to the barycenter atoms and then renormalises each row by its plan
mass `row_mass <- rowSums(plan)`, so a specimen's projection (hence its
score) depends on which other specimens are scored alongside it. Its
cost normalisation (`.sot_cost`) also divides by the MEDIAN pairwise
cost of the scored batch – a second cross-row statistic. Both make the
primitive fail the single-sample deployability gate
([`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)).
This method removes BOTH couplings: every specimen is transported on its
own as a one-atom measure (the row-independent entropic-kernel softmin
the primitive itself uses only in its `nrow == 1` branch), and the cost
scale is FROZEN at fit from the training clouds, never recomputed from
the scored set.

**Entropic-kernel transport energy.** For one specimen with rCLR vector
\\z\\ (a Dirac \\\delta_z\\) and a frozen reference cloud
\\\\a_j\\\_{j=1}^m\\ with uniform weights, the entropic-kernel (Gibbs /
Sinkhorn) transport energy is the softmin of squared-Euclidean transport
costs \$\$E\_\varepsilon(z, A) = -\varepsilon \log\\\Big( \tfrac{1}{m}
\sum\_{j=1}^m \exp(-\\C(z, a_j) / \varepsilon) \Big),\qquad C(z, a_j) =
\lVert z - a_j \rVert_2^2 / \tau,\$\$ where \\\varepsilon\\ is the
entropic regularisation and \\\tau\\ the frozen cost scale (median
training pairwise cost). This is exactly the single-source entropic-OT /
Sinkhorn-potential value (the log-sum-exp the multi-row primitive
evaluates per atom in its single-sample branch); as \\\varepsilon \to
0\\ it tends to the nearest-atom transport cost. The score is
\\\mathrm{orientation} \cdot (E\_\varepsilon(z, \mathrm{control}) -
E\_\varepsilon(z, \mathrm{case}))\\: a specimen closer (in entropic-OT
energy) to the case cloud than to the control cloud scores higher.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`); this is exactly invariant
to per-specimen scaling and uses no cross-row statistic, exactly as the
other roster scorers do. Universe features absent from a specimen are
dropped from BOTH the specimen and the frozen reference atoms for that
one transport (the cost is computed on the present-feature subspace), so
a missing feature contributes no transport cost.

**Negative-control honesty.** This is EXPECTED to discriminate weakly on
real data and is not engineered to look strong; it nonetheless produces
a finite, deterministic, exactly row-equivariant score and separates a
strongly planted case-vs-control shift (so the wiring is testable).
Degenerate inputs (no present features, fewer than `hp$min_features`
present, an all-zero specimen, or an empty reference cloud) return the
neutral score `0`.

## References

Cuturi M. (2013) Sinkhorn Distances: Lightspeed Computation of Optimal
Transport. *Advances in Neural Information Processing Systems (NeurIPS)*
26:2292-2300.

Peyre G, Cuturi M. (2019) Computational Optimal Transport. *Foundations
and Trends in Machine Learning* 11(5-6):355-607.
