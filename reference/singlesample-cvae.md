# Counterfactual class-conditional VAE single-sample discriminator

Single-sample within-cohort discriminator (family N) implementing a
counterfactual CLASS-conditional variational autoencoder used as a
generative classifier. A class-conditional VAE (encoder \\q(z\|x,y)\\ +
decoder \\p(x\|z,y)\\, both plain BatchNorm-free MLPs conditioned on a
one-hot of the class label \\y\\) is trained on the TRAINING set by the
reparameterized ELBO via reticulate-python torch, then **FROZEN**; a new
specimen is scored by the COUNTERFACTUAL evidence log-likelihood-RATIO
\$\$s(x) = \mathrm{ELBO}(x, y=\mathrm{case}) - \mathrm{ELBO}(x,
y=\mathrm{control}),\$\$ computed in PURE base-R with the latent MEAN
\\z = \mu\\ (DETERMINISTIC, no sampling, no RNG at score). Larger score
= more case-like.

**PRESPEC DEPARTURE (kit-conditional -\> class-conditional).** The
present `fit_kit_conditional_vae`
(R/singlesample-02-learned-kit-aware-vae.R) conditions a VAE on KIT
metadata – an external cross-sample covariate that VIOLATES
single-sample deployability (the kit of a held-out specimen need not be
known/seen at score). The single-sample `cvae` instead conditions on the
CLASS label \\y\\ and is used as a COUNTERFACTUAL generative classifier:
it models \\p(x\|y)\\ for each class and scores a query by the
counterfactual evidence ratio \\\mathrm{ELBO}(x,\mathrm{case}) -
\mathrm{ELBO}(x,\mathrm{control})\\. This needs NO kit/cohort input at
score (only the query \\x\\; the class conditioning is internal – both
classes are evaluated and differenced), so it is single-sample-faithful.
This reinterpretation is flagged for the gate.

**Torch is confined to FIT.** The class-conditional VAE is trained with
torch on the rCLR'd, universe-aligned training matrix. After training,
the encoder, mu-head, logvar-head, decoder, and xhat-head weight
matrices and biases, plus the scalar decoder variance \\\sigma^2\\, are
EXPORTED to R as plain numeric matrices/vectors; the python module is
DISCARDED. The fitted model therefore holds NO external pointer and NO
python dependency at score time: the encoder forward, the decoder
forward, and the per-class ELBO are all pure base-R. This makes the
default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch forward by \\\sim 10^{-6}\\. The score uses ONLY the exported
weights through the SAME pure-R forwards (`.cvae_encode` /
`.cvae_decode`) at both fit-time verification and score, so there is no
second code path to drift (the deepmaha / proto-net / ssl-vicreg
fit/score-consistency lesson).

**The forwards.** The encoder is \\\mathrm{concat}\[x, e_y\](p+2) \to
\[\mathrm{Linear} \to \phi\]\_{1..H_e} \to (\mathrm{Linear}\to\mu,\\
\mathrm{Linear}\to\log\sigma^2_z)\\; the activation (`"relu"` default,
or `"tanh"`) follows every hidden Linear, and the two output heads (mu,
logvar) are bare Linear (no activation). The decoder is
\\\mathrm{concat}\[z, e_y\](L+2) \to \[\mathrm{Linear} \to
\phi\]\_{1..H_d} \to \mathrm{Linear} \to \hat{x}(p)\\; activation
follows every hidden Linear, the xhat-head is bare Linear. \\e_y =
\mathrm{onehot}(y)\\ (length 2).

**The counterfactual score (pure R, deterministic).** For a query rCLR
vector \\x\\, for EACH class \\c \in \\0,1\\\\: \$\$(\mu_c,
\log\sigma^2\_{z,c}) = \mathrm{encode}(x, e_c);\quad z_c = \mu_c\\
\mathrm{(latent\\ mean)};\quad \hat{x}\_c = \mathrm{decode}(\mu_c,
e_c);\$\$ \$\$\log N_c = -\tfrac12\big(p\\\log(2\pi\sigma^2) +
\tfrac{1}{\sigma^2}\sum_j (x_j - \hat{x}\_{c,j})^2\big);\$\$
\$\$\mathrm{KL}\_c = \tfrac12 \sum_l \big(\mu\_{c,l}^2 +
e^{\log\sigma^2\_{z,c,l}} - \log\sigma^2\_{z,c,l} - 1\big);\$\$
\$\$\mathrm{ELBO}\_c = \log N_c - \beta\\\mathrm{KL}\_c;\qquad s(x) =
\mathrm{ELBO}\_{\mathrm{case}} - \mathrm{ELBO}\_{\mathrm{ctrl}}.\$\$
Using \\z = \mu\\ (NOT a sample) makes the score DETERMINISTIC and
RNG-free, and it depends only on \\x\\ and the frozen model -\>
single-sample. Per-row -\> single-row == batch EXACTLY 0. \\\beta =
\mathrm{kl\\beta}\\ (default 1.0).

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score, so a specimen's score
depends only on that specimen and the frozen model – it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe tested on the (pre-rCLR) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin but is a VALID specimen and is
scored normally (NOT floored) – the empty-support test is on the
ORIGINAL abundances, NOT on `all(z == 0)`.

## References

Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes.
*International Conference on Learning Representations (ICLR)*.
arXiv:1312.6114.

Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation
using Deep Conditional Generative Models. *Advances in Neural
Information Processing Systems (NeurIPS)* 28.
