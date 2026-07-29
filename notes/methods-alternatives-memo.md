# Methods memo: alternatives to mean Hellinger enrichment

**Status:** research note (not shipped API)  
**Scope:** subject-level categorical compositions → group contrasts  
**Baseline:** Paul Edlefsen’s Hellinger between/within enrichment as implemented in `HellingerDistanceEnrichment`

This memo maps alternatives to the package’s **mean** pairwise Hellinger aggregation onto the same design, with primary emphasis on **Wilcox / WRS2** robust methods and a secondary survey of **Bayesian nonparametric** options (including where MCMC is required). No package API changes are proposed here—only a shortlist for later work.

Math outside tables uses GitHub’s `$` / `$$` delimiters; table cells use Unicode so they render reliably on GitHub.

---

## 1. Current procedure (baseline)

The unit of analysis is a **subject-level composition**: multinomial counts over categories (clusters, cell types, etc.), one row per subject, with a group label.

1. **Jeffreys softening.** Add pseudocounts (default $1/2$ per category) to each subject’s count vector.
2. **Pairwise Hellinger distance** on the simplex:

$$
H(p,q) = \frac{1}{\sqrt{2}}\,\bigl\|\sqrt{p}-\sqrt{q}\bigr\|_2 \in [0,1].
$$

3. **Group summary matrix.** For each pair of groups, within-block entries are the **mean** of pairwise Hellinger distances among subjects in that group; off-diagonal entries are the **mean** of cross-group distances.
4. **Effect size**

$$
R = \frac{\overline{d}_{\mathrm{between}}}{\overline{d}_{\mathrm{within}}}.
$$

5. **Inference**
   - **Permutation:** fix the softened distance matrix; shuffle labels; recompute $R$; observed-inclusive $p$-value.
   - **Bayes:** draw subject compositions from conjugate Dirichlet posteriors; rebuild Hellinger distances; nest the same permutation null; report posterior mean $R$ / $p$ and credible intervals.

Implemented in `R/hellinger-internals.R` and `R/compare-group-compositions.R`. The Dirichlet step is conjugate (`rgamma` normalization)—**no MCMC** in the current Bayes path.

---

## 2. Where “mean” enters

Two aggregation layers are easy to confuse:

| Layer | Current choice | What it controls |
|-------|----------------|------------------|
| **A. Pairwise distance summary** | Arithmetic mean of Hellinger distances within a set (or between two sets) | Sensitivity to outlier *subjects* (a single atypical composition inflates d̄) |
| **B. Contrast functional** | Ratio R = d̄_between / d̄_within | Scale of the effect size and behavior when within-group distance is near zero |

Alternatives can replace **A**, **B**, or both, while keeping Hellinger geometry and the label-permutation null. Wilcox/WRS2 ideas act mostly on **A** (robust location of the distance cloud) or bypass distances entirely (projection WMW on the composition vectors). Energy / probability-of-superiority functionals replace **B**.

```mermaid
flowchart LR
  counts[Subject counts]
  soft[Jeffreys soften]
  dist[Hellinger distance matrix]
  aggA[Layer A: mean or robust summary]
  aggB[Layer B: ratio energy or PS]
  null[Permutation or nested Bayes]
  counts --> soft --> dist --> aggA --> aggB --> null
```

---

## 3. Wilcox / WRS2 map

Primary references: Wilcox, *Introduction to Robust Estimation and Hypothesis Testing* (esp. Ch. 5 comparing two groups; Ch. 6 multivariate methods, including the projection-type Wilcoxon–Mann–Whitney extension); Mair & Wilcox, *Robust Statistical Methods Using WRS2* (CRAN package `WRS2` and vignette). We cite WRS2 as a convenient implementation of Wilcox’s toolkit—this package need not depend on it.

### 3.1 Robust aggregation of Hellinger distances (Layer A)

Replace the arithmetic mean of pairwise distances with a robust location:

- **Median** of within- or between-pair Hellinger distances.
- **20% trimmed mean** (Wilcox’s default trim level in many WRS2 routines).
- Optional **MOM / one-step M-estimator** of location on the distance list (WRS2 `est = "mom"` / `"onestep"` / `"median"` pattern used in `pb2gen`, `sppba`, etc.).

**Why it fits.** Outlier subjects produce large Hellinger distances that dominate $\overline{d}$ and can inflate or deflate $R$ even when most subjects are exchangeable. Robust Layer A keeps the Edlefsen geometry and the same permutation machinery: only `mean_pairwise_distance` / `mean_pairwise_distance_between` change.

**Inference.** Still label permutation on a fixed (softened) distance matrix, or nested Dirichlet draws that rebuild the matrix. Bootstrap CIs for the robust summary itself can use Winsorized variance ideas from Wilcox Ch. 3–5 if interval estimation (beyond permutation $p$) is desired.

### 3.2 Energy / DISCO-style contrast (Layer B)

Instead of the ratio $R$, use an energy-distance-style functional of the same Hellinger matrix. For two groups with within means $\overline{d}_{11}$, $\overline{d}_{22}$ and between mean $\overline{d}_{12}$,

$$
E = 2\,\overline{d}_{12} - \overline{d}_{11} - \overline{d}_{22}
$$

(with multi-group omnibus formed by averaging pairwise $E$ or via distance-components ANOVA / DISCO). Related energy statistics appear in Székely & Rizzo’s energy framework; the algebra is compatible with any metric, including Hellinger.

**Pros.** $E = 0$ under identical distributions of compositions (in population, for the underlying metric); no division by a near-zero within mean.  
**Cons.** Scale depends on the metric; less immediately “enrichment ratio” interpretable than $R$.  
**Null.** Same label permutation.

### 3.3 Probability of superiority on distances (Layer B)

Wilcox emphasizes nonparametric effect sizes built from stochastic ordering—classically

$$
\pi = P(X < Y)
$$

for univariate responses (Wilcox–Mann–Whitney / probability of superiority; WRS2 `dep.effect` reports a related **SIGN** $P(X < Y)$ for paired settings). Applied to the distance design:

$$
\mathrm{PS}_{\mathrm{dist}} = P\bigl(d_{\mathrm{between}} > d_{\mathrm{within}}\bigr)
$$

estimated by comparing the empirical cloud of between-pair Hellinger distances to the cloud of within-pair distances (ties handled by mid-ranks or half-credit, as in Vargha–Delaney $A$ / Ruscio’s $A$).

**Pros.** Bounded, interpretable (“how often is a random between-pair more separated than a random within-pair?”); robust to a few extreme distances without needing trimmed means.  
**Cons.** Depends on how within and between pairs are sampled when group sizes differ; needs a clear sampling convention for multi-group omnibus.  
**Null.** Permutation of labels (or bootstrap for CIs, as in Ruscio & Mullen on $A$).

### 3.4 Projection Wilcoxon–Mann–Whitney (bypass Layer A/B distances)

Wilcox Ch. 6 §6.10: *A two-sample, projection-type extension of the Wilcoxon–Mann–Whitney test* (also Wilcox, 2005, *British Journal of Mathematical and Statistical Psychology*). Rough procedure:

1. Soften subject compositions (or work on $\sqrt{p}$ / clr coordinates).
2. Estimate a **robust multivariate center** for each group (e.g. spatial/OP median; trimmed mean of composition vectors).
3. Project all subjects onto the line connecting the two centers.
4. Compute univariate WMW / $\pi$ (and optionally plot projected scores).

**Pros.** Directly multivariate; effect size is the familiar $\pi$; does not reduce the design to a single mean distance first; graphical projection is interpretable.  
**Cons.** Needs care with simplex geometry (closure constraint); centers on raw proportions vs Hellinger $\sqrt{p}$ can differ; multi-group omnibus needs pairwise projections or a multi-sample depth/rank extension (Wilcox Ch. 6 depth-based WMW analogs).  
**Null.** Label permutation is natural and avoids dependence issues Wilcox notes for naive tests on projected points.

### 3.5 Discrete multinomial toolkit (`discANOVA`, `discmcp`, `binband`) — scope mismatch

WRS2 provides:

- **`binband`:** compare two discrete distributions at each support point (Storer–Kim; Kulinskaya–Morgenthaler–Staudte options).
- **`discANOVA`:** global test that $J$ independent groups share **identical multinomial distributions** (Storer–Kim generalization; bootstrap under pooled multinomial).
- **`discmcp` / `discstep`:** multiple comparisons / step-down for discrete responses.

These methods assume each *observation* is a draw from a discrete sample space with **small cardinality** (e.g. Likert item, single category label). That is **not** the same design as this package:

| | WRS2 `discANOVA` | This package |
|--|------------------|--------------|
| Observation | One discrete outcome per subject | One **composition vector** (counts over many categories) per subject |
| Null | Identical multinomials for the scalar discrete response | Between-group Hellinger enrichment of subject compositions |
| Typical use | Likert / small support | Cluster frequencies, cell-type mixes |

**Legitimate adjacent uses** (not drop-in replacements for $R$):

- Dominant-category label per subject → `discANOVA` / `discmcp`.
- Category-wise binary presence/absence or per-category margins → `binband`-style contrasts (multiple-testing burden).

The memo should not claim `discANOVA` replaces Hellinger enrichment; it is a Wilcox tool for a different but related discrete problem.

---

## 4. Bayesian nonparametric options

### 4.1 Already exact: conjugate Dirichlet per subject

Current `method = "bayes"` places a Dirichlet–multinomial conjugate posterior on each subject’s composition, draws frequencies, and nests Edlefsen’s permutation null. This is **Bayesian**, but **parametric conjugate**, not a nonparametric prior on the *law of compositions* within a group. No MCMC.

### 4.2 Group-level BNP: are group composition distributions equal?

Target: random probability measures $F_g$ on the simplex (or on $\mathbb{R}^{K-1}$ after a log-ratio map), with subjects $p_{ig} \sim F_g$.

| Approach | What you get | Exact vs MCMC |
|----------|--------------|---------------|
| **Pólya tree two-sample test** (Holmes, Caron, Griffin & Stephens, 2015, *Bayesian Analysis*) | Bayes factor / posterior probability H₀: F₁ ≡ F₂ | Often **analytic** marginal likelihoods for Pólya-tree partitions (no MCMC for the BF itself) |
| **Dependent / hierarchical Dirichlet process** on group-level measures | Flexible F_g with shared atoms; posterior functionals of Hellinger / total variation between F_g | Typically **MCMC** (Gibbs / slice / stick-breaking) |
| **Bernstein / stick-breaking models on the simplex** (e.g. compositional BNP regression literature) | Density estimation and group comparison on compositions | **MCMC** |

**Fit to this package.** These ask a different primary question—“are the subject-composition *distributions* equal?”—than Edlefsen’s enrichment ratio. They can complement $R$ (especially for two-group contrasts) but do not automatically yield collapse/subset contrast tables. Dimension of the simplex ($K$ categories) makes generic continuous BNP harder as $K$ grows; Dirichlet process mixtures of Dirichlets or finite stick-breaking on the simplex are more natural than kernel densities in ambient Euclidean space.

### 4.3 Hellinger disparity + BNP density (MHB / BHM)

Literature combining **minimum Hellinger distance** estimation with Bayesian nonparametric density estimates (e.g. Hooker & Vidyashankar disparity posteriors; asymptotic MHB/BHM results such as Wu & Hooker-type lines of work) targets robust parametric inference when data are continuous and a kernel or histogram BNP estimate of $g$ is available.

**Applicability here is weaker:** subject compositions already live on a finite simplex with multinomial likelihoods; a continuous KDE/Hellinger disparity layer is usually overkill relative to Dirichlet or DP-Dirichlet models. If pursued, expect **MCMC or heavy numerical integration**, not a conjugate one-liner.

---

## 5. Comparison table

| Method | Effect size | Null | Exact / MCMC | Fit to subject-composition design | Impl. cost | Robust to outlier subjects |
|--------|-------------|------|--------------|-----------------------------------|------------|----------------------------|
| Current mean Hellinger R | d̄_B / d̄_W | Label perm (+ nested Dir.) | Exact (conjugate draws) | Native | Already shipped | Low (mean-sensitive) |
| Median / trimmed Hellinger R | Same R, robust Layer A | Same | Exact | Native plug-in | Low | High |
| Energy E on Hellinger | 2d̄_B − d̄_W1 − d̄_W2 | Label perm | Exact | Native plug-in | Low | Medium (still mean-based unless robustified) |
| PS on distances | P(d_B > d_W) | Label perm / bootstrap CI | Exact | Native plug-in | Low–medium | High |
| Projection WMW / π | π on projected scores | Label perm | Exact | Strong (multivariate) | Medium | High (robust centers) |
| WRS2 `discANOVA` etc. | Multinomial equality tests | Bootstrap under pooled mult. | Exact (bootstrap) | **Poor** as drop-in; OK for collapsed labels | Low (external pkg) | Designed for discrete support |
| Conjugate Dirichlet Bayes (current) | Posterior of R | Nested perm | Exact | Native | Shipped | Same as chosen Layer A/B |
| Pólya-tree two-sample BF | P(H₀ \| data) | Model comparison | Often analytic | Two-group distributional equality | Medium–high | Prior-dependent |
| HDP / DDP on simplex | Posterior distances between F_g | Posterior / decision | **MCMC** | Strong if K manageable | High | Prior-dependent |
| Hellinger disparity + BNP density | Disparity posterior for parametric θ | Model-based | Usually MCMC | Weak for finite simplex multinomials | High | Aimed at continuous contamination |

---

## 6. Recommended shortlist for later code

### Tier 1 — plug into the existing Hellinger pipeline (prefer first)

Keep Hellinger distances and Edlefsen contrasts; add arguments conceptually like:

- `summary = c("mean", "median", "trimmed")` for Layer A pairwise aggregation;
- `statistic = c("ratio", "energy", "ps")` for Layer B.

Permutation and nested Dirichlet Bayes reuse the same scaffolding. Highest leverage per line of code; closest to the current scientific claim (“Hellinger enrichment”).

### Tier 2 — Wilcox projection WMW / $\pi$

On softened compositions (or $\sqrt{p}$ coordinates), estimate robust group centers, project, report $\pi$ with label permutation. Complements Tier 1 when the scientific question is stochastic dominance of compositions along the primary between-group direction rather than a distance ratio. Multi-group: all pairwise projections + existing contrast/collapse machinery.

### Tier 3 — research / MCMC BNP

Only if conjugate subject-level Dirichlet nesting is judged insufficient for uncertainty in the *group-level laws* $F_g$:

- Start with **Pólya-tree two-sample** Bayes factors for two-group contrasts (often no MCMC).
- Escalate to **hierarchical / dependent DP** (or Dirichlet mixtures of Dirichlets) on the simplex when $K$ is moderate and sharing of mixture atoms across groups matters.

Do **not** prioritize continuous Hellinger-disparity BNP hybrids for this multinomial design unless a specific robustness-to-contamination story appears.

### Explicit non-goals for v1 alternatives

- Do not treat WRS2 `discANOVA` as a replacement for Hellinger enrichment.
- Do not add a WRS2 hard dependency; reimplement the small Tier 1–2 pieces needed, citing Wilcox/WRS2.
- Do not change vignette pedagogy until Tier 1 lands and is validated on example cohorts.

---

## References (selective)

**Baseline / package**

- Edlefsen, P. Original Hellinger between/within enrichment procedure (Jeffreys softening, ratio statistic, label permutation, collapse/subset contrasts); nested Dirichlet Bayes as planned composition-uncertainty extension. Implemented in `HellingerDistanceEnrichment`.

**Wilcox / WRS2**

- Wilcox, R. R. *Introduction to Robust Estimation and Hypothesis Testing* (4th/5th ed.). Elsevier. Especially Ch. 5 (two groups, effect sizes related to $\pi$); Ch. 6 (multivariate methods; projection-type WMW §6.10; depth-based WMW analogs).
- Wilcox, R. R. (2005). A multivariate projection-type analogue of the Wilcoxon–Mann–Whitney test. *British Journal of Mathematical and Statistical Psychology*, 58, 69–81.
- Mair, P., & Wilcox, R. (2019/ongoing). *Robust Statistical Methods Using WRS2*. CRAN package `WRS2` and vignette (incl. `binband`, `discANOVA`, `discmcp`, `discstep`, trimmed/MOM estimators, `dep.effect` SIGN/$\pi$-style summaries).
- Field, A. P., & Wilcox, R. R. (2017). Robust statistical methods: A primer for clinical psychology. *Clinical Psychology Review* (accessible overview of Wilcox’s program).

**Energy / distance functionals**

- Székely, G. J., & Rizzo, M. L. Energy statistics / distance components (DISCO) and related $k$-sample tests (R package `energy`).

**Probability of superiority**

- Ruscio, J. (2008). A probability-based measure of effect size. *Psychological Methods*.
- Ruscio, J., & Mullen, T. (2012). Confidence intervals for the probability of superiority. *Multivariate Behavioral Research*.
- Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistic. *Journal of Educational and Behavioral Statistics*.

**Bayesian nonparametric two-sample / compositional**

- Holmes, C. C., Caron, F., Griffin, J. E., & Stephens, D. A. (2015). Two-sample Bayesian nonparametric hypothesis testing. *Bayesian Analysis*, 10(2), 297–320. (Pólya tree; often analytic BF.)
- Teh, Y. W., Jordan, M. I., Beal, M. J., & Blei, D. M. (2006). Hierarchical Dirichlet processes. *JASA* (canonical HDP; MCMC).
- Dependent / graphical DP and simplex Bernstein-stick-breaking models for compositional responses (e.g. Barrientos, Jara, Quintana line of work)—MCMC for flexible $F_g$ on the simplex.

**Hellinger disparity + BNP (weaker fit here)**

- Hooker, G., & Vidyashankar, A. N. Bayesian disparity / Hellinger posterior constructions.
- Asymptotic MHB/BHM-type results combining minimum Hellinger distance with Bayesian nonparametric density estimates (continuous-data setting).

---

## Bottom line

The shipped method’s distinctive geometry is **Hellinger on subject compositions**; its fragile step is mostly the **arithmetic mean** in Layer A (and the ratio in Layer B when within-group distance is tiny). Wilcox/WRS2 give immediate, permutation-compatible upgrades: **robust Layer A**, **energy or PS Layer B**, and a **projection WMW** path that never averages distances. Discrete WRS2 multinomial tests are valuable for collapsed labels but are **not** substitutes for enrichment. Conjugate Dirichlet Bayes already covers subject-level uncertainty without MCMC; **group-level BNP** (Pólya tree first, HDP/DDP if needed) is the MCMC frontier if the scientific target shifts from “enrichment of pairwise Hellinger” to “equality of composition laws.”
