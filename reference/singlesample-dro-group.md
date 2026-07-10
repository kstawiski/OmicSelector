# Group-DRO kit x biofluid robust discriminator (canonical single-sample wrapper)

Canonical `fit_*`/`score_*` wrapper around the existing Group-DRO kit x
biofluid scorer primitive
([`fit_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/fit_group_dro_scorer.md)
/
[`score_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/score_group_dro_scorer.md),
family N, transfer estimand). The primitive fits a train-only robust-CLR
logistic classifier whose exponentiated-gradient sample weights are
pooled over kit x biofluid groups (Sagawa et al., ICLR 2020), so the
trained scorer is robust to the worst-performing provenance group rather
than only the pooled average. The SCORE path of the primitive is already
single-sample / row-equivariant: a specimen is mapped to the frozen
train-only CLR with a per-row pseudocount and per-row pivot centring,
scaled by the frozen median/MAD constants, and pushed through the frozen
logistic head; it carries no cross-row statistic and depends only on the
specimen and the frozen model. This wrapper does NOT reimplement any
Group-DRO math; it only adapts the primitive to the canonical contract
(`fit_dro_group(X_train, y_train, meta_train, hp)` /
`score_dro_group(model, X, meta)`) so the method is dispatchable via
[`singlesample_score_call`](https://kstawiski.github.io/OmicSelector/reference/singlesample_score_call.md).

**Auditable group engagement.** The Group-DRO objective only differs
from plain ERM when there is more than one kit x biofluid group. When
the canonical `meta_train` is `NULL` or lacks the configured kit /
biofluid columns, this wrapper DEGENERATES GRACEFULLY to a single-group
fit by synthesizing a constant single-group metadata frame, so the
primitive runs without crashing and still yields a valid single-sample
model – but the run is then mathematically ERM, not group-robust. To
never silently mislabel the estimand, the model records `n_groups` (the
number of distinct kit x biofluid groups actually engaged) and the
logical `groups_engaged` (`n_groups > 1`); these mirror the `n_cohorts`
/ `n_environments` engagement signals exposed by `reo-metaktsp` and
`sel-stablemate`.

**Orientation.** `score_group_dro_scorer` returns the logistic
probability of the case class, so larger already means more case-like;
the wrapper preserves this orientation unchanged.

**Metadata at score time.** `meta` passed to `score_dro_group` is
forwarded to the primitive only as a diagnostic (it attaches test-group
audit attributes); it never enters the numeric score, which is stripped
to a plain finite numeric vector before return.

## References

Sagawa S, Koh PW, Hashimoto TB, Liang P. (2020) Distributionally Robust
Neural Networks for Group Shifts. *ICLR*. arXiv:1911.08731.
