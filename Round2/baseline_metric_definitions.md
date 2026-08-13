# Evaluating OLS, PDP, and Forward Marginal Effects under the Direction and Strength Metrics

**Companion note for MS-INS-2024-07487.R2 — Theme 2 (benchmark contrast)**

Purpose: specify exactly how the baseline estimators are instantiated inside the alignment framework of §4.1 of
the round-2 manuscript, so that their direction and strength scores are numerically comparable with those already
reported for SHAP and LIME.

---

## Contents

- [0. Notation](#0-notation)
- [1. The two metrics, as defined in the round-2 manuscript](#1-the-two-metrics-as-defined-in-the-round-2-manuscript)
  - [1.1 Direction alignment](#11-direction-alignment)
  - [1.2 Strength alignment](#12-strength-alignment)
- [2. What each baseline must supply](#2-what-each-baseline-must-supply)
- [3. Baseline instantiations](#3-baseline-instantiations)
  - [3.1 OLS (linear probability model)](#31-ols-linear-probability-model)
  - [3.2 Partial dependence](#32-partial-dependence)
  - [3.3 Forward marginal effects](#33-forward-marginal-effects)
  - [3.4 Summary](#34-summary)
- [4. Why a second ground-truth target is required](#4-why-a-second-ground-truth-target-is-required)
  - [4.1 The rule, and what it implies for each method](#41-the-rule-and-what-it-implies-for-each-method)
  - [4.2 Why this is not merely inelegant](#42-why-this-is-not-merely-inelegant)
  - [4.3 The addition](#43-the-addition)
- [5. Magnitude accuracy: MAPE, and the SHAP estimand problem](#5-magnitude-accuracy-mape-and-the-shap-estimand-problem)
  - [5.1 What MAPE is](#51-what-mape-is)
  - [5.2 Why SHAP has no magnitude estimand](#52-why-shap-has-no-magnitude-estimand)
  - [5.3 Three reporting options](#53-three-reporting-options)
  - [5.4 A constructive resolution on the linear datasets](#54-a-constructive-resolution-on-the-linear-datasets)
  - [5.5 Guard rails](#55-guard-rails)
  - [5.6 A bonus that closes a repeated request](#56-a-bonus-that-closes-a-repeated-request)
- [6. Estimator multiplicity: bootstrapping the baselines](#6-estimator-multiplicity-bootstrapping-the-baselines)
- [7. Alternative baselines](#7-alternative-baselines)
- [8. Producing the Figure 3-style alignment PDFs for the baselines](#8-producing-the-figure-3-style-alignment-pdfs-for-the-baselines)
  - [8.1 Decisions to settle before running](#81-decisions-to-settle-before-running)
  - [8.2 Procedure](#82-procedure)
  - [8.3 Two exclusions that must be handled, not ignored](#83-two-exclusions-that-must-be-handled-not-ignored)
  - [8.4 Stratification is required, not optional](#84-stratification-is-required-not-optional)
  - [8.5 Plotting](#85-plotting)
  - [8.6 Expected patterns, and what each would mean](#86-expected-patterns-and-what-each-would-mean)
  - [8.7 On restricting to a linear subset](#87-on-restricting-to-a-linear-subset)
- [9. Practical settings](#9-practical-settings)
- [10. Corrections to the current plan](#10-corrections-to-the-current-plan)
- [Appendix A. Issues in the 81-dataset design, independent of the baselines](#appendix-a-issues-in-the-81-dataset-design-independent-of-the-baselines)
- [References](#references)

---
## 0. Notation

$G$ is the ground-truth data-generating process, $\mathcal{M}$ the fitted predictive model, $d$ the number of
features, $S$ the evaluation set, $\mathbf{x}^{(i)}\in\mathbb{R}^d$ an instance. Features are standardised to
z-scores. Perturbed instances and the perturbation grid are as in the manuscript:

$$\mathbf{x}^{(i)}_{j,\delta}=\bigl(x^{(i)}_1,\dots,x^{(i)}_j+\delta,\dots,x^{(i)}_d\bigr),
\qquad
\Delta=\{m\cdot \mathrm{Std}\;\mid\; m\in\{-M,\dots,-1,1,\dots,M\}\},\quad M=3 .$$

**Estimand.** All quantities below are *associational and ceteris paribus*: the change in $\mathbb{E}[Y\mid X]$
when coordinate $j$ moves and all other coordinates are held fixed.

---

## 1. The two metrics, as defined in the round-2 manuscript

### 1.1 Direction alignment

For a method $\mathcal{E}$ with attribution $e_\mathcal{E}(\mathbf{x}^{(i)},j)$,

$$\Delta e^{(i)}_j(\delta)=e_\mathcal{E}\bigl(\mathbf{x}^{(i)}_{j,\delta},j\bigr)-e_\mathcal{E}\bigl(\mathbf{x}^{(i)},j\bigr),
\qquad
\Delta Y^{(i)}_j(\delta)=G\bigl(\mathbf{x}^{(i)}_{j,\delta}\bigr)-G\bigl(\mathbf{x}^{(i)}\bigr),$$

$$I^{\mathrm{dir}}_{ij}(\delta)=
\begin{cases}
\mathbb{1}\!\left[\Delta e^{(i)}_j(\delta)\,\Delta Y^{(i)}_j(\delta)>0\right] & \text{if } \bigl|\Delta Y^{(i)}_j(\delta)\bigr|>\varepsilon_1\\[4pt]
\mathbb{1}\!\left[\bigl|\Delta e^{(i)}_j(\delta)\bigr|<\varepsilon_2\right] & \text{otherwise,}
\end{cases}
\qquad
\rho_{\mathrm{dir},j}=\frac{1}{|S|\,|\Delta|}\sum_{i=1}^{|S|}\sum_{\delta\in\Delta} I^{\mathrm{dir}}_{ij}(\delta).$$

**The reference quantity $\Delta Y$ is method-independent.** Direction alignment is therefore *already*
cross-method comparable; a baseline is added simply by supplying its $\Delta e^{(i)}_j(\delta)$.

### 1.2 Strength alignment

With $I_\mathcal{E}(j)$ the global importance under $\mathcal{E}$, $I_G(j)$ the ground-truth importance, and
$R_\mathcal{E},R_G$ the induced rankings,

$$\rho_{\mathrm{strength}}(R_\mathcal{E},R_G)=1-\frac{6\sum_{j=1}^{d}\bigl(R_\mathcal{E}(j)-R_G(j)\bigr)^2}{d(d^2-1)} .$$

Strength requires *two* supplied objects: $I_\mathcal{E}(j)$ and a target $I_G(j)$. The manuscript defines $I_G$
**method-specifically**. §4 explains why that becomes a problem once baselines are added; §5 adds the magnitude
metric R1 asked for, having observed that *"rank correlation alone is [not] enough to capture alignment. A method
can recover rankings well while misestimating magnitudes, or vice versa."*

---

## 2. What each baseline must supply

| Slot | Symbol | Used by |
|---|---|---|
| Perturbation response | $\Delta e^{(i)}_j(\delta)$ | direction alignment |
| Global importance | $I_\mathcal{E}(j)$ | strength alignment |
| Ground-truth target | $I_G(j)$ | strength alignment |

---

## 3. Baseline instantiations

### 3.1 OLS (linear probability model)

Fit $\hat{Y}=\hat\beta_0+\sum_j\hat\beta_j X_j$ on standardised features.

$$\boxed{\;\Delta e^{(i)}_j(\delta)=\hat\beta_j\,\delta\;}
\qquad
\boxed{\;I_{\mathrm{OLS}}(j)=\bigl|\hat\beta_j\bigr|\;}$$

**Interpretation.** $\hat\beta_j\delta$ is the predicted change in the outcome from moving feature $j$ by $\delta$
standard deviations with all other features held fixed. It is constant across instances, so OLS asserts *one sign
per feature for the entire population*; the direction metric consequently asks of OLS only whether that single
global sign is right, a weaker demand than the instance-level claim it imposes on SHAP. This asymmetry must be
stated when the table is reported. $|\hat\beta_j|$ ranks features by expected outcome change per standard
deviation — an effect-size ranking, not a predictive-contribution ranking. That distinction is precisely the one
R1 pressed in round 1: *"Feature importance typically refers to a feature's ability to predict the dependent
variable, whereas marginal effects quantify how the dependent variable changes when the feature itself changes.
These are fundamentally different concepts."* Having a baseline whose importance *is* an effect size makes the
contrast legible rather than merely asserted.

Two technical notes. First, $\hat\beta_j\delta$ is an identity rather than an approximation: for a linear model
$\hat{Y}(\mathbf{x}_{j,\delta})-\hat{Y}(\mathbf{x})=\hat\beta_j\delta$ exactly, so the plug-in prediction
difference and the coefficient rule coincide. This mirrors the manuscript's LIME instantiation
($\Delta e = w^{(i)}_j\delta$), letting OLS enter as the *global* counterpart to LIME's local surrogate under an
identical rule. Second, standardisation is load-bearing: without z-scoring, $|\hat\beta_j|$ ranks by units of
measurement. If features are left on raw scales, use $|\hat\beta_j|\hat\sigma_j$.

### 3.2 Partial dependence

With $\mathrm{PD}_j(v)=\mathbb{E}_{X_{-j}}\bigl[\mathcal{M}(v,X_{-j})\bigr]$, estimated as
$\widehat{\mathrm{PD}}_j(v)=\frac1n\sum_{i=1}^{n}\mathcal{M}\bigl(v,\mathbf{x}^{(i)}_{-j}\bigr)$
([Friedman 2001](https://doi.org/10.1214/aos/1013203451)):

$$\boxed{\;\Delta e^{(i)}_j(\delta)=\widehat{\mathrm{PD}}_j\bigl(x^{(i)}_j+\delta\bigr)-\widehat{\mathrm{PD}}_j\bigl(x^{(i)}_j\bigr)\;}$$

$$\boxed{\;I_{\mathrm{PD}}(j)=\sqrt{\frac{1}{K-1}\sum_{k=1}^{K}\Bigl(\widehat{\mathrm{PD}}_j(v_k)-\overline{\mathrm{PD}}_j\Bigr)^{2}},
\qquad \overline{\mathrm{PD}}_j=\frac{1}{K}\sum_{k=1}^{K}\widehat{\mathrm{PD}}_j(v_k)\;}$$

the **variation (dispersion) measure** of
[Greenwell, Boehmke & McCarthy (2018)](https://arxiv.org/abs/1805.04755): the standard deviation of the partial
dependence function of $X_j$, evaluated at the $K$ unique observed values $v_1,\dots,v_K$ of $X_j$.

**Interpretation.** $\Delta e$ is how the *population-average* prediction moves when feature $j$ is set to a new
value — it answers "what happens to the average account," not "what happens to this account," because $X_{-j}$ has
been integrated out. $I_{\mathrm{PD}}(j)$ is how much that average prediction swings across the observed range of
$X_j$: a feature whose PD curve is nearly constant moves the average prediction very little and is ranked
unimportant. Note what this cannot express — because $\widehat{\mathrm{PD}}_j$ depends on $i$ only through
$x^{(i)}_j$, instances sharing a value of $X_j$ receive identical $\Delta e$, so PD direction is a property of
position on one marginalised curve and is structurally blind to heterogeneity driven by $X_{-j}$. Evaluate the
curve once on a grid and interpolate; cost is then independent of $|\Delta|$.

*Terminology.* Greenwell et al. motivate the measure by the heuristic that a **flat** PD curve indicates an
unimportant feature — "it is this notion of 'flatness' which we will use as a basis to define our variable
importance measure" — but the quantity computed is a standard deviation, which increases with *departure* from
flatness. It is a variation measure, correctly oriented as importance (larger = more important) and requiring no
inversion. Describing it as "a flatness measure" reverses the sign of the description; the `vip` package's
phrasing of relative "flatness" propagates the same looseness and is best not quoted.

*Evaluation points.* The strength measure uses the empirical distribution of $X_j$ so that dispersion is weighted
by where data actually lie; a uniform grid over the range would overweight sparse regions and inflate
$I_{\mathrm{PD}}$ for skewed features. This differs from the grid used for the **direction** metric, which must
extend beyond the observed range by $\max|\delta|$ so perturbed points remain interpolable. Use the extended grid
for $\Delta e$ and the observed values for $I_{\mathrm{PD}}$, and say so.

### 3.3 Forward marginal effects

For step $h$ ([Scholbeck et al. 2022](https://arxiv.org/abs/2201.08837);
[2024](https://doi.org/10.1007/s10618-023-00993-x)),
$\mathrm{FME}^{(i)}_j(h)=\mathcal{M}\bigl(\mathbf{x}^{(i)}_{j,h}\bigr)-\mathcal{M}\bigl(\mathbf{x}^{(i)}\bigr)$:

$$\boxed{\;\Delta e^{(i)}_j(\delta)=\mathcal{M}\bigl(\mathbf{x}^{(i)}_{j,\delta}\bigr)-\mathcal{M}\bigl(\mathbf{x}^{(i)}\bigr)\;}$$

$$\mathrm{AME}_j=\frac{1}{n}\sum_{i=1}^{n}\frac{\mathrm{FME}^{(i)}_j(h)}{h},
\qquad
\boxed{\;I_{\mathrm{FME}}(j)=\frac{1}{n}\sum_{i=1}^{n}\frac{\bigl|\mathrm{FME}^{(i)}_j(h)\bigr|}{h}\;}$$

**Why the mean absolute FME rather than the AME.** Three reasons, all decisive:

1. *The manuscript's own definition of strength is undirected.* §4.1.3 states that strength alignment "evaluates
   whether an explainer correctly captures the relative and **undirected** importance of features."
   $\mathrm{AME}_j$ is a signed average and therefore answers a different question.
2. *Signed averaging cancels.* A feature whose effect is positive for some instances and negative for others
   averages toward zero, so $\mathrm{AME}_j\approx 0$ ranks it least important even when its effects are the
   largest anywhere in the sample. This is not hypothetical for the simulation designs in the paper: quadratic
   terms make $\partial G/\partial x_j$ change sign across the support, and interaction terms make the sign depend
   on a second feature. Under exactly the conditions the paper is built to study, AME systematically understates
   the features of greatest interest.
3. *Like-for-like comparability.* The competing quantities — $\mathbb{E}|\phi_j|$ for SHAP, $\mathbb{E}|w_j|$ for
   LIME, $|\hat\beta_j|$ for OLS, the PD standard deviation — are all non-cancelling. Ranking FME by a cancelling
   statistic while ranking every rival by a non-cancelling one would make the comparison an artefact of the
   aggregation rule rather than of the estimators.

$\mathrm{AME}_j$ is retained, but for a different slot: it is the natural signed estimate of the true average
partial effect, and so belongs in the **magnitude** comparison (MAPE against $\mathbb{E}[\partial G/\partial x_j]$)
that R1 requested — *"Including another metric that captures magnitude accuracy—such as mean absolute percentage
error (MAPE)—would provide a more complete and interpretable assessment."* Reporting $\mathrm{AME}_j$ there and $I_{\mathrm{FME}}(j)$ here uses each where its sign
convention matches the target. The gap between them is itself informative: $I_{\mathrm{FME}}(j)\gg|\mathrm{AME}_j|$
flags a feature with strongly heterogeneous effects, which is a warning that *any* single-number global summary —
including a SHAP bar chart — will mislead.

**Interpretation.** $\Delta e$ is the actual change this model predicts for this instance: no surrogate, no
averaging, no attribution scheme. $I_{\mathrm{FME}}(j)$ is the typical absolute size of the model's local response
per unit step — the model's own answer to "how much does moving feature $j$ move the prediction." This makes FME
the pivotal row of the comparison: it interrogates $\mathcal{M}$ with no explainer in the loop, so Oracle→FME is
model approximation error and FME→SHAP/LIME is what the explainer layer adds. That decomposition answers R3's
first priority item, which asks the paper to distinguish a property of the model and data from a property of the
attribution method: *"Separate identification from explainer behavior. Your central result like in 5.2 is that
predictive accuracy does not identify G. This is a model/data property. Clarify what is failing in SHAP/LIME
beyond this."*

### 3.4 Summary

| Method | $\Delta e^{(i)}_j(\delta)$ | $I_\mathcal{E}(j)$ | Locality | Reads as |
|---|---|---|---|---|
| OLS | $\hat\beta_j\delta$ | $\lvert\hat\beta_j\rvert$ | global, instance-invariant | one sign, one effect size, for everyone |
| PDP | $\widehat{\mathrm{PD}}_j(x^{(i)}_j+\delta)-\widehat{\mathrm{PD}}_j(x^{(i)}_j)$ | sd of PD curve (variation) | global curve, varies with $x_j$ only | effect on the average account |
| FME | $\mathcal{M}(\mathbf{x}^{(i)}_{j,\delta})-\mathcal{M}(\mathbf{x}^{(i)})$ | mean $\lvert\mathrm{FME}\rvert/h$ | fully local | effect on *this* account, per the model |
| SHAP | $\phi_j(\mathbf{x}^{(i)}_{j,\delta})-\phi_j(\mathbf{x}^{(i)})$ | mean $\lvert\phi_j\rvert$ | local, aggregated | credit allocated to feature $j$ |
| LIME | $w^{(i)}_j\delta$ | mean $\lvert w_j\rvert$ | local surrogate | local linear slope |

---

## 4. Why a second ground-truth target is required

### 4.1 The rule, and what it implies for each method

The manuscript sets $I_G$ **method-specifically**: the ground-truth importance is obtained by applying the *same*
method to $G$. So each method is asked its own question:

| Method | Question $\rho_{\mathrm{strength}}$ actually answers under the current rule |
|---|---|
| SHAP | Does SHAP-of-$\mathcal{M}$ rank features as SHAP-of-$G$ does? |
| LIME | Does LIME-of-$\mathcal{M}$ rank features as LIME-of-$G$ does? |
| OLS | Does OLS-on-data rank features as the population projection of $G$ does? |
| PDP | Does PD-of-$\mathcal{M}$ rank features as PD-of-$G$ does? |
| FME | Does FME-of-$\mathcal{M}$ rank features as FME-of-$G$ does? |

Every row is a sensible question. But they are **five different questions with five different answer keys**, so
the numbers cannot be placed in one column and compared. If SHAP scores $0.95$ and OLS scores $0.90$, nothing
follows about which better recovers the data-generating process: they were graded on different exams.

### 4.2 Why this is not merely inelegant

The deeper problem is that a self-referential target is **blind to the exact error the paper exists to document**.

Suppose — as the manuscript argues — that mean absolute Shapley value is a poor proxy for the strength of feature
effects in a data-generating process. Then $I_G^{\text{self}}$ for SHAP, being *itself* a mean absolute Shapley
value computed on $G$, inherits that distortion. SHAP-of-$\mathcal{M}$ can then reproduce SHAP-of-$G$ perfectly,
scoring $\rho_{\mathrm{strength}}=1.00$, while both rank $G$'s actual feature effects badly. The metric would
certify the pipeline as flawless in precisely the case the paper calls misinterpretation.

The analogy is grading a translation by comparing it against another translation by the same translator. Perfect
agreement demonstrates *consistency*, not *fidelity to the original*. The self-referential target measures
**pipeline fidelity** — does substituting $\mathcal{M}$ for $G$ corrupt this method's output? — but is silent on
**estimand validity** — does this method's notion of importance track the DGP's feature effects at all?

This also explains why the issue surfaces only now. For a paper studying explainers alone, the self-referential
target was a defensible choice: it isolates the model-substitution step, which was the object of study. Adding
baselines is what breaks it, because the baselines' whole purpose is to answer AE point 2 — *"how much of the
observed misalignment is specific to the post hoc explanation pipeline and how much reflects the broader
difference between prediction and inference"* — and that question cannot be answered by a target that moves with
the method being scored.

Note also the asymmetry with §1.1: direction alignment already uses a common reference, $\Delta Y$ computed from
$G$ itself. Strength is the only metric with this defect, and fixing it makes the two metrics structurally
consistent.

### 4.3 The addition

Report both targets:

$$I_G^{\text{self}}(j)=\text{method applied to } G,
\qquad
I_G^{\text{common}}(j)=\mathbb{E}\!\left[\left|\frac{\partial G(\mathbf{X})}{\partial x_j}\right|\right].$$

$I_G^{\text{common}}$ is the mean absolute partial effect of the true process: unsigned (matching every
$I_\mathcal{E}$ in §3), method-free, and available analytically for a simulated $G$. Define

$$\text{estimand-mismatch gap}_\mathcal{E}
=\rho^{\text{self}}_{\mathrm{strength}}(\mathcal{E})-\rho^{\text{common}}_{\mathrm{strength}}(\mathcal{E}).$$

The gap should be near zero for OLS, PDP, and FME — all of which estimate partial effects, so their own yardstick
nearly coincides with the common one — and materially positive for SHAP, which does not. That contrast *is* the
paper's central claim expressed as a measured quantity rather than an argument, and it converts R1's round-1
objection into a result rather than a concession: *"The most critical issue is the paper's conflation of feature
importance with marginal effects… Critically, the most important features might not correspond to the features
with the largest marginal effects."*

For $I_G^{\text{self}}$ with the baselines: OLS on $G$ is the population least-squares projection
$\boldsymbol\beta_G=\arg\min_{a,\mathbf b}\mathbb{E}\bigl[(G(\mathbf{X})-a-\mathbf{b}^\top\mathbf{X})^2\bigr]$;
PDP and FME on $G$ substitute $G$ for $\mathcal{M}$ in §3.2–3.3.

---

## 5. Magnitude accuracy: MAPE, and the SHAP estimand problem

R1's request in round 2: *"A method can recover rankings well while misestimating magnitudes, or vice versa… Including another metric that captures magnitude accuracy—such as mean absolute percentage error (MAPE)—would provide a more complete and interpretable assessment."* This is Theme 2, and it is coupled to Theme 1, because a percentage error is undefined until the estimand is declared.

### 5.1 What MAPE is

For true values $\theta_j$ and estimates $\hat\theta_j$ across $d$ features,

$$\mathrm{MAPE}=\frac{1}{d}\sum_{j=1}^{d}\frac{\bigl|\hat\theta_j-\theta_j\bigr|}{\bigl|\theta_j\bigr|}\;(\times 100).$$

Each error is expressed relative to the size of the quantity being estimated, so the measure is scale-free and answers a question rank correlation cannot: *by what percentage are the estimated effects wrong?* Its known pathologies — undefined at $\theta_j=0$, explosive for small $|\theta_j|$, unbounded above but bounded below by zero, hence asymmetric between over- and under-estimation — are documented in [Hyndman & Koehler (2006)](https://doi.org/10.1016/j.ijforecast.2006.03.001) and will be raised by a referee if not pre-empted.

The natural target here is the true average partial effect, $\mathrm{APE}_j=\mathbb{E}\!\left[\partial G/\partial x_j\right]$.

### 5.2 Why SHAP has no magnitude estimand

An estimand is a population quantity such that "$\hat\theta_j$ estimates $\theta_j$" is a meaningful claim. The asymmetry across methods:

| Method | What it estimates | APE estimand? |
|---|---|---|
| OLS $\hat\beta_j$ | $\partial\mathbb{E}[Y\mid X]/\partial x_j$ | yes |
| PD slope | derivative of the marginalised curve | yes (exactly, under independence) |
| $\mathrm{AME}_j$ | average forward difference of $\mathcal{M}$ | yes |
| LIME $w^{(i)}_j$ | local slope of $\mathcal{M}$ | yes, locally |
| SHAP $\phi_j$ | a **share** of $f(\mathbf{x})-\mathbb{E}[f]$ | **no** |

Shapley values are defined by an allocation axiom: they divide a prediction's departure from baseline among features subject to efficiency, $\sum_j\phi_j=f(\mathbf{x})-\mathbb{E}[f]$. A share of a total is not a rate of change. Two symptoms make the distinction concrete:

1. **Baseline collapse.** $\phi_j\approx 0$ whenever $x_j$ sits at its baseline value, *however large that feature's effect is*. No derivative behaves this way. R1 noted the same phenomenon in round 1 and judged it reasonable **for an importance measure**: *"the authors
   suggest it is undesirable for a feature to have a SHAP value of 0 if its value equals the average. However, from
   a feature importance perspective, this is reasonable because an average feature value is unlikely to indicate a
   larger or smaller dependent variable value."* That is precisely the point — reasonable for an importance
   measure, unreasonable for an effect estimate.
2. **Unit sensitivity.** Rescale $x_j\mapsto 2x_j$. For a linear model $\beta_j\mapsto\beta_j/2$ while $(x_j-\mu_j)\mapsto 2(x_j-\mu_j)$, so $\phi_j$ is **unchanged** while the APE **halves**. Two quantities that respond differently to a change of units cannot be estimates of the same thing.

So $\mathbb{E}|\phi_j|$ is a legitimate ranking device and not an effect-magnitude estimator: MAPE against APE is not *difficult* for SHAP, it is *undefined*. The practical difficulty is presentational — reporting MAPE for four baselines and a blank cell for SHAP answers R1 for every method except the one under critique, which reads as evasion.

### 5.3 Three reporting options

1. **MAPE in APE units.** Defined for OLS, logit, PD slope, AME, and LIME; blank for SHAP. Honest, and the blank cell *is* the paper's argument — but it must be argued, not merely displayed.
2. **MAPE on normalised importance shares.** Divide each $I_\mathcal{E}$ by its sum and compare against the normalised $I_G^{\text{common}}$. Defined for **every** method including SHAP, at the cost of measuring profile shape rather than effect size.
3. **MAPE after optimal rescaling.** Choose $c^\star=\arg\min_c\sum_j\bigl(c\,I_\mathcal{E}(j)-I_G(j)\bigr)^2$ and score $c^\star I_\mathcal{E}$. This hands each method its most favourable scale and tests shape alone; if SHAP still fails after being given the best possible scaling, the result is far stronger than failure on units it never claimed.

Report 1 and 2 together, with 3 as the generous robustness check.

### 5.4 A constructive resolution on the linear datasets

For linear $G$ with independent features, Shapley values admit the closed form $\phi_j(\mathbf{x})=\beta_j\bigl(x_j-\mathbb{E}[X_j]\bigr)$, so

$$\mathbb{E}\bigl|\phi_j\bigr|=\bigl|\beta_j\bigr|\cdot\mathbb{E}\bigl|X_j-\mu_j\bigr|
\;=\;\bigl|\beta_j\bigr|\,\sigma_j\sqrt{2/\pi}\quad\text{(Gaussian marginals).}$$

On the linear subset, therefore, $\mathbb{E}|\phi_j|$ **is** proportional to $|\mathrm{APE}_j|$ with an analytically known constant ($\sqrt{2/\pi}\approx 0.798$ for standardised features). Divide it out, place SHAP on the APE scale, and compute a genuine MAPE; residual error is then attributable to $\mathcal{M}\neq G$ rather than to the estimand.

This yields the paper's story in one exhibit: *SHAP magnitudes are interpretable as effect sizes under linearity, here is the exact conversion, and here is how it degrades as quadratic and interaction terms enter.* R1 is answered fully where an answer exists, and the impossibility is demonstrated rather than asserted where it does not.

Two caveats. The constant assumes Gaussian marginals, which the balancing step in Algorithm 1 violates (Appendix A.2) — compute $\mathbb{E}|X_j-\mu_j|$ empirically instead. And the closed form requires linear $G$, so it applies only to the nonlinear-free, interaction-free cells.

### 5.5 Guard rails

- **Exclude near-zero denominators.** Drop features with $|\mathrm{APE}_j|<\tau$ for a stated $\tau$ and report how many were dropped per dataset. This is not optional here: quadratic terms contribute $2\gamma_j x_j$ to $\partial G/\partial x_j$, which averages to $\approx 0$ under symmetric features, so quadratically-driven features have $|\mathrm{APE}_j|\approx 0$ by construction (Appendix A.8).
- **Report a bounded companion** — sMAPE, or the total-variation distance between normalised profiles from §4 — so no single feature can dominate.
- **Beware the $G$ versus $\mathbb{E}[Y\mid X]$ scale gap.** With $G=\sigma(z)$ but $\mathbb{E}[Y\mid X]=\Phi(z/s)$, the baselines and $\mathcal{M}$ estimate derivatives of $\Phi(z/s)$ while the target is a derivative of $\sigma(z)$; at $s=0.1$ these differ by a large multiplicative factor. Direction and rank are unaffected (both are monotone in $z$), but raw-units MAPE is not. Either declare $\mathbb{E}[Y\mid X]$ the magnitude target, or use options 2–3 above, which are scale-free by construction. See Appendix A.3.

### 5.6 A bonus that closes a repeated request

R1's underlying complaint — misranking matters more when the magnitude gaps are large — was made in round 1 as well, about the concordance and relevance metrics: *"weighting by the size of the discrepancy would better reflect the severity of the mistake."* A **weighted Kendall's $\tau$**, with rank discrepancies weighted by the magnitude gaps between the features involved, implements exactly that and costs one line of code. Adopting it visibly closes a suggestion R1 has now made in two consecutive rounds.

---

## 6. Estimator multiplicity: bootstrapping the baselines

Refitting OLS on $B$ resamples is nearly free and buys considerably more than error bars. Use the **pairs
bootstrap** (resample $(\mathbf{x}^{(i)},y^{(i)})$ jointly), since $X$ is random here; the residual bootstrap
would impose a homoskedastic linear model that the DGP violates. Let $\hat{\boldsymbol\beta}^{(b)}$,
$b=1,\dots,B$ denote the replicate estimates.

**6.1 Alignment distributions rather than point estimates.** Compute $\rho_{\mathrm{dir}}$ and
$\rho_{\mathrm{strength}}$ separately for each replicate. This yields a *distribution* of alignment scores for the
baseline, directly comparable to the distributions the manuscript already plots for SHAP and LIME. Since the
paper's Finding 2 is about **long left tails** rather than means, showing whether the OLS alignment distribution
has a comparably heavy tail is the sharpest possible form of the benchmark contrast — stated in the paper's own
currency.

**6.2 Sign stability as a direction diagnostic.** For each feature define
$\pi_j=\frac1B\sum_b \mathbb{1}\bigl[\mathrm{sign}(\hat\beta^{(b)}_j)=\mathrm{sign}(\bar\beta_j)\bigr]$.
Under strong collinearity $\pi_j$ falls well below 1, quantifying the sign instability that the notebook already
detects indirectly (OLS strength alignment collapsing to $0.539$ at $\rho=0.9$). $\pi_j$ makes it a reported
statistic rather than an inference from a degraded score.

**6.3 Bagged versus single-fit estimates.** Scoring the averaged $\bar{\boldsymbol\beta}$ against scoring a single
fit reproduces, for a parametric baseline, the manuscript's Robustness Check 1 (averaging explanations over a
Rashomon set). Identical experimental logic, different model class — which lets you say whether averaging helps
*generally* or only for post hoc explainers.

**6.4 The payoff for Theme 5: a second multiplicity diagnostic.** Define bootstrap agreement

$$A_{\mathrm{boot}}=\frac{2}{B(B-1)}\sum_{b<b'}\mathrm{Spearman}\Bigl(R_{|\hat\beta^{(b)}|},\,R_{|\hat\beta^{(b')}|}\Bigr),$$

the exact analogue of the explanation agreement in §7 of the manuscript, and test whether it predicts
$\rho_{\mathrm{strength}}$ as Table 1 tests Rashomon agreement. Either outcome is publishable. If it does,
agreement-as-reliability-signal is a **general property of estimator multiplicity**, which broadens the
contribution and speaks to the AE's demand that the diagnostic be positioned against prior work: *"The paper also
needs to distinguish this idea carefully from existing work on model multiplicity, explanation instability, and
disagreement across near-optimal models."* If it does not,
the Rashomon version is doing something specific to model multiplicity, which is a stronger novelty claim. Given
that the AE called the Rashomon diagnosis *"the most promising part of the paper"* and *"the clearest basis for a
distinctive contribution,"* resolving this is high value for very little compute.

**6.5 Keep the two multiplicities distinct.** Bootstrap multiplicity is *sampling* multiplicity — one estimator,
many samples. A Rashomon set is *specification* multiplicity — one sample, many estimators. They are different
sources and must not be conflated in the write-up; the honest framing is that both are instances of "the data do
not pin down a unique answer," which is R3's point. A crossed design (bootstrap $\times$ model class) separates
them if you want the full decomposition, and connects naturally to model class reliance
([Fisher, Rudin & Dominici 2019](https://jmlr.org/papers/v20/18-760.html)) and variable importance clouds
([Dong & Rudin 2020](https://arxiv.org/abs/1901.03209)).

**6.6 Extending beyond OLS.** The same scheme applies to PDP and FME by refitting $\mathcal{M}$ on each resample,
but the cost is $B$ model fits plus $B$ explanation passes. With XGBoost this is affordable at $B\approx 25$–$50$
for PD and FME and prohibitive for SHAP; report the reduced $B$ honestly rather than silently dropping the check.
$B=500$ is trivial for OLS alone.

---

## 7. Alternative baselines

If the team drops one of the three, these are the defensible substitutes. The middle columns record what each can
supply, which determines whether it enters both metrics or only one.

| Candidate | Supplies $\Delta e$? | Supplies $I_\mathcal{E}$? | Why choose it |
|---|---|---|---|
| **ALE** ([Apley & Zhu 2020](https://doi.org/10.1111/rssb.12377)) | yes: $\mathrm{ALE}_j(x_j+\delta)-\mathrm{ALE}_j(x_j)$ | yes: sd of the ALE curve | Designed to stay on-manifold under correlated features — the one condition the paper identifies as the dominant driver of misalignment. The strongest substitute for PDP, and arguably worth including *alongside* it. |
| **GAM / spline additive model** | yes: $s_j(x_j+\delta)-s_j(x_j)$ | yes: sd of $s_j$ over observed $x_j$ | A *correctly specified* nonparametric benchmark for additive DGPs. Recovers nonlinearity OLS misses while retaining a unique fit, sharpening R3's contrast between direct estimation and the explanation pipeline. EBMs ([Nori et al. 2019](https://arxiv.org/abs/1909.09223)) are the ML-flavoured version. |
| **Permutation feature importance** ([Breiman 2001](https://doi.org/10.1023/A:1010933404324)) | no | yes | The default practitioner baseline, so its absence would be conspicuous. But it measures *predictive* contribution, not effect size — making it a useful foil for R1's importance-vs-effect distinction rather than a competitor on strength. Strength only. |
| **Model class reliance** ([Fisher et al. 2019](https://jmlr.org/papers/v20/18-760.html)) | no | yes, as an interval | Importance bounds across all near-optimal models rather than a point estimate; the natural bridge to Theme 5's positioning requirement. |
| **LOCO** ([Lei et al. 2018](https://doi.org/10.1080/01621459.2017.1307116)) | no | yes | Refit-without-$j$ importance with distribution-free inference; expensive ($d$ refits) but yields valid confidence intervals. |
| **SAGE** ([Covert et al. 2020](https://arxiv.org/abs/2004.00668)) | no | yes | Shapley values of the *loss* rather than of a prediction — a global, model-level importance that is the fair Shapley-family comparator to mean $\lvert\phi_j\rvert$. |
| **Shapley effects / Sobol indices** ([Owen 2014](https://doi.org/10.1137/130936233); [Song et al. 2016](https://doi.org/10.1137/15M1048070)) | no | yes | Variance-based sensitivity indices defined on the *function*, computable directly on $G$. A principled fallback for $I_G^{\text{common}}$ if a referee objects to the mean absolute derivative. |
| **Double/debiased ML** ([Chernozhukov et al. 2018](https://doi.org/10.1111/ectj.12097)) | yes | yes | Modern partial-effect estimation, but standardly motivated causally — it reopens the causal flank Theme 1 closes. Include only if framed strictly as an associational estimator. |

Practical guidance. If only one baseline survives, keep **FME**: it is the only row that isolates the explainer
layer. That is the request most likely to hold the paper at Major Revision if left unanswered, because R3
conditions the significance of the whole contribution on it: *"The contribution is significant if the authors can
isolate what is genuinely specific to post hoc explainers beyond this generic identification problem; as written,
a reader can come away feeling the core result is an XAI-flavored restatement of omitted-variable/identification
concerns."* R3 then names the experiment required: *"Show a case where the model faithfully recovers the data's
associational structure yet the explainer still misleads."*

If two baselines survive, add **OLS**. It is the one R1 named — *"At a minimum, including a simple baseline—such
as linear regression with normalized features—would make the results more interpretable"* — and the one R3
expects to strengthen rather than weaken the argument: *"I would also encourage at least one benchmark contrast
(e.g., what a correctly specified regression recovers on the same DGP) so readers can see that the failure is
shared by any method that does not impose identifying assumptions, which sharpens rather than weakens the paper's
message."*

If PDP is dropped over the correlation concerns
in §3.2, replace it with **ALE** rather than deleting the global-curve row, since losing it would leave no global
effect-shape estimator in the comparison.

---

## 8. Producing the Figure 3-style alignment PDFs for the baselines

Figure 3 plots densities of direction and strength alignment across the 81 dataset–model pairs. Extending it to baselines requires **no new datasets, no new model fits, and no new SHAP computation**: each of OLS, PD, and FME is cheap, and re-scoring existing $\mathbb{E}|\phi_j|$ vectors against a new target costs nothing. This is the highest value-per-hour item in Theme 2. The design decisions below must be settled first, because several are not neutral.

### 8.1 Decisions to settle before running

**(a) Which features are scored.** §4.3 averages direction scores "over the five features with the largest mean explanation values" — i.e. each method is graded on the features *it* nominated. That mirrors practice, but it means OLS and SHAP may be scored on different feature sets, and a method that wrongly elevates a feature is graded on its own error. Make the **ground-truth top-5** the primary rule and retain method-selected as a practice-mirroring variant. This is the same logical move as §4: comparability requires a fixed yardstick.

**(b) Whether strength uses all $d$ features or only 5.** The averaging sentence in §4.3 is well defined for direction, which is per-feature ($\rho_{\mathrm{dir},j}$), but strength alignment is already a single global Spearman per dataset. As written it is ambiguous which set it covers. Resolve it explicitly: AE point 5 asks for exactly this kind of
pass — *"The authors should conduct a full consistency check of the main paper, the appendix, the analyses, the
figures, the captions, and all reported results."*

**(c) The magnitude target.** Per §5.5, declare whether the target is $G=\sigma(z)$ or $\mathbb{E}[Y\mid X]=\Phi(z/s)$. Direction and strength are insensitive to the choice; MAPE is not.

**(d) Ground truth is computed on the delivered data.** The balancing step conditions on $Y$, so the feature distribution in $X_{\mathrm{final}}$ is not the generating distribution (Appendix A.2). Compute $I_G^{\text{common}}$, APEs, and $\Delta Y$ from the balanced sample actually used, not from the nominal $\mathcal{N}(0,1)$ design.

### 8.2 Procedure

For each of the 81 dataset–model pairs $(D_t,\mathcal{M}_t)$:

1. **Recover the design record** $(m,\rho,p,q)$ and the sampled coefficients, so ground truth is exact rather than re-estimated.
2. **Compute ground truth** on the evaluation split: $\Delta Y^{(i)}_j(\delta)$ for $\delta\in\Delta$; $\mathrm{APE}_j$; and $I_G^{\text{common}}(j)=\mathbb{E}\lvert\partial G/\partial x_j\rvert$ analytically from the index.
3. **Fit the baselines.** OLS on the training split (standardised); PD curves on a grid extended by $\max|\delta|$ with the observed-value grid retained for $I_{\mathrm{PD}}$; FME at step $h$. Total cost is seconds per dataset.
4. **Assemble** $\Delta e^{(i)}_j(\delta)$ and $I_\mathcal{E}(j)$ per §3, plus the cached $\mathbb{E}|\phi_j|$ and LIME coefficients.
5. **Score** $\rho_{\mathrm{dir},j}$ per feature, average over the selected five per 8.1(a); compute $\rho_{\mathrm{strength}}$ against **both** $I_G^{\text{self}}$ and $I_G^{\text{common}}$ (§4.3).
6. **Record** dataset identifiers $(m,\rho,p,q)$ alongside each score so the densities can be stratified.

This yields, per method, 81 direction scores and 81 strength scores — the same object Figure 3 plots.

### 8.3 Two exclusions that must be handled, not ignored

**$\rho=1$ produces an exactly duplicated feature.** Algorithm 1 line 5 sets $x_m\leftarrow\rho x_k+\sqrt{1-\rho^2}\,\epsilon$; at $\rho=1$ the noise term vanishes and $x_m=x_k$ identically. $X^\top X$ is then singular, OLS is not identified, and $\hat\beta_k,\hat\beta_m$ are individually meaningless — only their sum is estimable. Statsmodels will silently drop a column or return unstable values. This is **27 of the 81 datasets**. Options, in order of preference:

- Report the **54 datasets with $\rho\in\{0,0.5\}$** as the primary baseline panel, with the $\rho=1$ stratum shown separately and labelled a non-identifiability limiting case.
- Or move that factor level to $\rho=0.9$, which preserves the "strong correlation" contrast without singularity. This requires regenerating 27 datasets and re-running SHAP on them, so it is only worth doing if the schedule allows.
- Or use a pseudoinverse/ridge fit and state plainly that the reported OLS is regularised at that level. Least preferred: it changes the estimator mid-table.

**Do not silently drop the stratum.** The correlation-strength coefficient in Figure 6 is the paper's largest reported effect on strength misalignment, and part of it is driven by this level, where two identical columns carry independently drawn true coefficients that no estimator can distinguish. Adding OLS makes that visible, because OLS fails loudly where SHAP failed quietly. Better to own it than to have a referee find it.

### 8.4 Stratification is required, not optional

Under a random ranking, Spearman's $\rho$ has standard deviation $\approx 1/\sqrt{d-1}$: $0.33$ at $d=10$ versus $0.19$ at $d=30$. Pooling all 81 datasets into one density therefore mixes three different null distributions, and part of the left tail is a dimensionality artefact rather than a property of any method. R2 made this exact critique in round 1 about the relevance metric — *"it is inherently true that relevance is higher with fewer features… the authors should represent this criterion in relative terms"* — so it will be recognised if repeated. Either facet the strength panel by $d\in\{10,20,30\}$, or report a normalised score (e.g. $\rho$ divided by its null standard deviation at that $d$).

### 8.5 Plotting

Use **ECDFs (or survival curves) as the primary display**, with KDEs retained for visual continuity with Figure 3. The headline claim is Finding 2, the long left tails, and KDEs are the wrong instrument for tails: alignment scores are bounded above at 1.0, so a Gaussian kernel leaks probability mass past the boundary, and the apparent tail shape depends on a bandwidth a referee can contest. ECDFs have no bandwidth and read tails directly.

Accompany the curves with a small table of **mean, median, P05, P10, and minimum** per method. That table, not the curves, is what actually substantiates a claim about tail risk.

On layout, R2's round-2 complaints apply directly — *"the y-axis ranges differ across panels even though the
underlying quantities appear comparable… The absence of grid lines further reduces readability. More generally,
the figure is not sufficiently self-contained"* — so use consistent y-axis ranges, gridlines, and captions that
define the plotted quantity without reference to the main text. Two panels, direction and strength, with methods
overlaid will read better than a multi-panel grid with varying axes.

### 8.6 Expected patterns, and what each would mean

State these as predictions before running, so the write-up is not retrofitted.

- **Every method's density shifts left as $\rho$ rises.** This is the pattern R3 predicts and asks the paper to
  make explicit — that the failure is *"shared by any method that does not impose identifying assumptions"* — so
  it is an identification result rather than an explainer defect, and reporting it as such is what R3 says
  *"sharpens rather than weakens the paper's message."*
- **OLS and logit degrade most in *magnitude* while holding up in *rank*.** Collinearity inflates coefficient variance without necessarily reordering features. If so, the paper must not claim regression is uniformly superior.
- **PD understates features whose contribution is dominated by interaction rather than by their own linear coefficient.** For an interaction $\lambda x_j x_k$, the PD curve retains only $\lambda v\,\mathbb{E}[x_k]$, which vanishes when $\mathbb{E}[x_k]=0$, whereas the ground-truth $\partial z/\partial x_j=\beta_j+2\gamma_j x_j+\lambda x_k$ retains it. Expect PD's strength density to have a heavier left tail on the $q>0$ cells. Note that balancing can push $\mathbb{E}[x_k]$ off zero and partially leak the interaction back into the curve, so the effect will be attenuated rather than exact.
- **Where SHAP or LIME fall below FME**, the deficit is attributable to the explainer layer, since FME interrogates the same model with no explainer in the loop. Those cells are the empirical content of R3's request for *"an example where the model recovers the data's
  associational structure but the explainer still misleads."* If no such cells appear at any $\rho$, that is itself
  a finding, and the honest response is to lean the contribution on the Rashomon diagnostic rather than on
  explainer-specific critique.

### 8.7 On restricting to a linear subset

A tempting simplification is to run baselines only on datasets with no nonlinear and no interaction terms. Note the count: with a $3\times3\times3\times3$ design, cells with **both** $p=0$ and $q=0$ number $3\times3\times1\times1=\mathbf{9}$, not 27 — 27 is the count with $p=0$ *or* with $q=0$ taken separately. Of those 9, three sit at $\rho=1$, leaving **6** that are both linear and non-singular: far too few for a density plot.

More importantly, the linear subset discards exactly the conditions that motivate black-box models, and R1 asked
in round 1 for nonlinear cases specifically: *"there should also be an example with non-linear relationships that
are better captured by a machine learning model, illustrating that while linear regression may perform poorly,
SHAP and LIME still fail to accurately capture marginal effects in these more complex scenarios."* Since the baselines are computationally trivial, run all 81 and stratify. If a reduced first pass is wanted, use the 54 datasets with $\rho\in\{0,0.5\}$, which sidesteps the singularity while preserving all nonlinearity and interaction variation. Reserve the 9 linear cells for the closed-form SHAP conversion in §5.4, where linearity is a requirement rather than a convenience.

---

## 9. Practical settings

**Thresholds.** $\varepsilon_1$ is in outcome units; a change below it is substantively negligible (e.g. $10^{-3}$
in probability). $\Delta e$ for OLS, PDP, FME, and LIME is *also* in outcome units, so $\varepsilon_2=\varepsilon_1$
is coherent for those four. SHAP attributions are not, so use a scale-matched
$\varepsilon_2=\varepsilon_1\cdot\mathrm{med}|\Delta e|/\mathrm{med}|\Delta Y|$. Validate the harness by scoring
$G$ against itself: it must return $\rho_{\mathrm{dir}}=\rho_{\mathrm{strength}}=1$ exactly.

**Ties.** The closed-form Spearman expression in §1.2 assumes no ties; PD variation and $|\hat\beta|$ can tie in
principle. Use a tie-corrected implementation and report which.

**Step sizes.** $\delta$ enters the metric; $h$ enters $I_{\mathrm{FME}}$. Keep them distinct and report
sensitivity to both — tree ensembles are piecewise constant, so small $\delta$ yields exactly zero predicted change
for many instances, depressing instance-level direction scores for FME and SHAP for reasons attributable to the
model class rather than the explainer.

**Correlation.** $\widehat{\mathrm{PD}}_j$ averages over the empirical marginal of $X_{-j}$ and so evaluates
$\mathcal{M}$ off the data manifold when features are dependent
([Apley & Zhu 2020](https://doi.org/10.1111/rssb.12377);
[Hooker et al. 2021](https://arxiv.org/abs/1905.03151)). Since feature correlation is the manuscript's headline
driver of misalignment, PD's degradation across a $\rho$-sweep is expected and should be reported as a finding,
with ALE as the robustness check.

---

## 10. Corrections to the current plan

**C1 — PDP and FME directionality as drafted are the same estimator.** The plan defines the PDP version via ICE
curves ([Goldstein et al. 2015](https://doi.org/10.1080/10618600.2014.907095)):
$f^{(i)}_j(x^{(i)}_j+\delta)-f^{(i)}_j(x^{(i)}_j)$. But the ICE curve is
$f^{(i)}_j(v)=\mathcal{M}(v,\mathbf{x}^{(i)}_{-j})$, so that difference equals
$\mathcal{M}(\mathbf{x}^{(i)}_{j,\delta})-\mathcal{M}(\mathbf{x}^{(i)})$ — exactly the forward marginal effect. The
two baselines would return identical direction scores at every instance and every $\delta$, and a referee would
notice. Fix: PDP direction uses the **marginalised PD curve** (§3.2); the ICE-difference form is FME and should be
labelled as such. This also restores a meaningful global-vs-local contrast between the rows.

**C2 — the OLS correlation worry largely dissolves; the causal-ancestor worry is a design choice.** The metric's
reference quantity $\Delta Y^{(i)}_j(\delta)$ is *itself* a ceteris paribus perturbation — only coordinate $j$
moves. Under correct specification $\hat\beta_j$ targets that same object, so correlation among features creates no
estimand mismatch; it inflates $\mathrm{Var}(\hat\beta_j)$ and produces finite-sample sign instability, which is a
*result to measure* (see §6.2) rather than a problem to solve. On causal ancestry: in a simulation that draws
$\mathbf{X}$ from a joint distribution and then sets $Y=G(\mathbf{X})+\epsilon$, features are correlated but none
is a causal ancestor of another — state this explicitly as a scope condition. Were feature-to-feature edges
introduced, OLS, PD, FME *and* $\Delta Y$ would all continue to target direct ceteris paribus effects and so remain
mutually comparable; what none would deliver is a total effect. Given Theme 1's associational commitment, that is
a footnote, not a redesign.

**C3 — rank FME by the mean absolute effect, not the AME.** Full argument in §3.3; in brief, AME is a signed
average that cancels under sign-heterogeneous effects, while every competing importance measure is non-cancelling.
Retain AME for the magnitude (MAPE) comparison.

**C4 — add a common strength target.** Per §4, a self-referential $I_G$ cannot support the cross-method benchmark
AE point 2 requires: *"Without a clear comparison point, it is difficult to know how much of the observed
misalignment is specific to the post hoc explanation pipeline and how much reflects the broader difference between
prediction and inference."* Report both targets and the gap.

**C5 — z-score before ranking on $|\hat\beta_j|$**, or rank on $|\hat\beta_j|\hat\sigma_j$.

---

## Appendix A. Issues in the 81-dataset design, independent of the baselines

These concern Algorithm 1 (Appendix OA 3) itself and would apply even if no baseline were added. Several bear on results already reported in the manuscript, so they are worth resolving before the third submission rather than after a referee raises them.

**A.1 — $\rho=1$ yields an exactly duplicated feature.** Line 5 gives $x_m=\rho x_k+\sqrt{1-\rho^2}\,\epsilon$, so at $\rho=1$, $x_m=x_k$ identically while the two carry independently drawn coefficients from $\mathrm{Uniform}(-1,1)$. The ground truth distinguishes them; no estimator can. This affects 27 datasets and inflates the correlation-strength coefficient in Figure 6 — the paper's largest reported driver of strength misalignment — by mixing a genuine correlation effect with a non-identifiability artefact. Recommend either relabelling that level as a deliberate limiting case or moving it to $\rho=0.9$.

**A.2 — Balancing induces dependence among all features.** The final step samples $N$ positives and $N$ negatives (Algorithm 1, section 8, lines 13–14). Since $Y$ is a common child of every feature, conditioning on it is collider conditioning: features generated independently become dependent in $X_{\mathrm{final}}$, and the marginals are no longer $\mathcal{N}(0,1)$. Three consequences. The $\rho=0$ stratum is not an uncorrelated control, which weakens its role as a baseline in Figure 6. Any ground-truth quantity must be computed under the delivered distribution. And a convenient property is lost: under symmetric $X$, omitted quadratic terms are orthogonal to linear ones, so a linear fit would be unbiased for the linear coefficients despite misspecification — balancing breaks that symmetry.

**A.3 — $G$ is not $\mathbb{E}[Y\mid X]$.** With $G=\sigma(z)$ and $Y=\mathbb{1}[\sigma(z+\epsilon)>0.5]=\mathbb{1}[z+\epsilon>0]$, the conditional expectation is $\mathbb{E}[Y\mid X]=\Phi(z/s)$. Any model trained on $Y$ targets $\Phi(z/s)$, not $\sigma(z)$. The two are monotone transforms of the same index, so **direction alignment and rank-based strength are unaffected** — but their derivatives differ by a large factor at $s=0.1$, so any magnitude comparison in raw units is affected. State which object is the target in each metric.

**A.4 — The labelling rule is nearly deterministic, which bears on §5.1 and on R3's request.** Because $Y=\mathbb{1}[z+\epsilon>0]$ with $s=0.1$ while $\mathrm{SD}(z)$ is of order 1–2, the Bayes-optimal classifier achieves accuracy well above the $[0.85,0.90]$ top bracket in Figure 4 — the irreducible error is small. The reported accuracy range is therefore produced almost entirely by deliberate model degradation (varying training size and hyperparameters) rather than by noise in the DGP. Two implications. First, §5.1's claim that accuracy is necessary for alignment is being tested in a regime where even the best models sit well below the achievable ceiling; a referee may reasonably ask whether misalignment persists at near-Bayes accuracy — which is close to R3's
request for *"an example where the model recovers the data's associational structure but the explainer still
misleads."* Second, [Semenova et al. (2022)](https://doi.org/10.1145/3531146.3533232) tie Rashomon set size to noise, so a near-deterministic labelling rule implies a *small* Rashomon set at the optimum; the multiplicity observed in §5.2 is then attributable to finite-sample and model-class effects rather than to label noise. That is a defensible position but should be argued explicitly, since §5.2 currently leans on the noise intuition.

**A.5 — Exactly one correlated pair regardless of dimension.** Line 4 correlates a single feature with a single other, so the share of correlated features is $1/m$: $1/10$ at $m=10$ but $1/30$ at $m=30$. The correlation factor's bite therefore shrinks mechanically as dimension grows, confounding the two factors. This is the most likely explanation for the "counterintuitive" positive dimension coefficient noted in §5.3, and it is a cleaner account than the one currently given. Consider scaling the number of correlated pairs with $m$, or reporting correlation effects within each $m$ stratum.

**A.6 — Spearman's null variance depends on $d$.** Under random ranking, $\mathrm{SD}(\rho)\approx 1/\sqrt{d-1}$, so
$d=10$ and $d=30$ have materially different null distributions. Pooling across $m\in\{10,20,30\}$ mixes them, and
part of the left tail in the strength panel of Figure 3 is a dimensionality artefact. R2 raised the identical
objection in round 1 about the relevance metric — *"it is inherently true that relevance is higher with fewer
features (e.g., at k=3) as there are fewer 'competitors' for the top three. This holds true even if feature
importances are drawn at random… I think the authors should represent this criterion in relative terms"* — so a
repeat in a different metric will be recognised. Stratify or normalise.

**A.7 — Algorithm 1 uses two interleaved numbering schemes.** Bold section headers run 1–8 while algorithm lines run 1–15, so "step 8" is ambiguous between "Balance Dataset" and "Select $q$ pairs, add interaction terms." Renumber, or reference lines only. This is small but exactly the class of issue AE point 5 asks to be swept out.

**A.8 — Quadratic terms drive true APEs toward zero.** $\partial z/\partial x_j$ contains $2\gamma_j x_j$, which averages to $\approx 0$ under symmetric features. Features whose contribution is chiefly quadratic therefore have $\mathrm{APE}_j\approx 0$ while their mean absolute influence is large. This is harmless for MAI-based strength but destroys MAPE denominators (§5.5), and it is worth a sentence in the text since it also means a signed average effect is a poor summary for those features — the same cancellation argument that governs the AME choice in §3.3.

**A.9 — Unstated: whether generated terms appear as columns.** OA 3 lines 7–8 say squared and interaction terms are "added," without stating whether they enter only the index $z$ or also appear as columns of $X_{\mathrm{final}}$. The distinction is consequential: if they are columns, a linear model is correctly specified and would win the benchmark comparison trivially. Specify explicitly.

**A.10 — Model-selection asymmetry between sections.** Section 5 uses 12 XGBoost models per dataset spanning accuracy bins, while Section 7 uses 10 high-accuracy models drawn from nine algorithm families. Figures built on these two sets are not directly comparable, which is the likely source of R2's confusion between Figures 3 and 4: *"In Figure 3, SHAP appears clearly superior
to LIME with respect to both direction and strength alignment… However, in the next section, where model accuracy
is varied, Figure 4 seems to suggest that in the highest-performance bracket SHAP and LIME perform rather
similarly on both dimensions. I found this somewhat confusing."* Label each figure with the model set it uses.

---

## References

- Apley, D. W., & Zhu, J. (2020). Visualizing the effects of predictor variables in black box supervised learning models. *JRSS-B*, 82(4), 1059–1086. [doi.org/10.1111/rssb.12377](https://doi.org/10.1111/rssb.12377)
- Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. [doi.org/10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324)
- Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68. [doi.org/10.1111/ectj.12097](https://doi.org/10.1111/ectj.12097)
- Covert, I., Lundberg, S., & Lee, S.-I. (2020). Understanding global feature contributions with additive importance measures. *NeurIPS 33*. [arxiv.org/abs/2004.00668](https://arxiv.org/abs/2004.00668)
- Dong, J., & Rudin, C. (2020). Exploring the cloud of variable importance for the set of all good models. *Nature Machine Intelligence*, 2, 810–824. [arxiv.org/abs/1901.03209](https://arxiv.org/abs/1901.03209)
- Efron, B. (1979). Bootstrap methods: Another look at the jackknife. *Annals of Statistics*, 7(1), 1–26. [doi.org/10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552)
- Fisher, A., Rudin, C., & Dominici, F. (2019). All models are wrong, but many are useful: Learning a variable's importance by studying an entire class of prediction models simultaneously. *JMLR*, 20(177), 1–81. [jmlr.org/papers/v20/18-760.html](https://jmlr.org/papers/v20/18-760.html)
- Friedman, J. H. (2001). Greedy function approximation: A gradient boosting machine. *Annals of Statistics*, 29(5), 1189–1232. [doi.org/10.1214/aos/1013203451](https://doi.org/10.1214/aos/1013203451)
- Goldstein, A., Kapelner, A., Bleich, J., & Pitkin, E. (2015). Peeking inside the black box: Visualizing statistical learning with plots of individual conditional expectation. *JCGS*, 24(1), 44–65. [doi.org/10.1080/10618600.2014.907095](https://doi.org/10.1080/10618600.2014.907095)
- Greenwell, B. M., Boehmke, B. C., & McCarthy, A. J. (2018). A simple and effective model-based variable importance measure. [arxiv.org/abs/1805.04755](https://arxiv.org/abs/1805.04755)
- Hooker, G., Mentch, L., & Zhou, S. (2021). Unrestricted permutation forces extrapolation: Variable importance requires at least one more model, or there is no free variable importance. *Statistics and Computing*, 31, 82. [arxiv.org/abs/1905.03151](https://arxiv.org/abs/1905.03151)
- Hyndman, R. J., & Koehler, A. B. (2006). Another look at measures of forecast accuracy. *International Journal of Forecasting*, 22(4), 679–688. [doi.org/10.1016/j.ijforecast.2006.03.001](https://doi.org/10.1016/j.ijforecast.2006.03.001)
- Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R. J., & Wasserman, L. (2018). Distribution-free predictive inference for regression. *JASA*, 113(523), 1094–1111. [doi.org/10.1080/01621459.2017.1307116](https://doi.org/10.1080/01621459.2017.1307116)
- Nori, H., Jenkins, S., Koch, P., & Caruana, R. (2019). InterpretML: A unified framework for machine learning interpretability. [arxiv.org/abs/1909.09223](https://arxiv.org/abs/1909.09223)
- Owen, A. B. (2014). Sobol' indices and Shapley value. *SIAM/ASA Journal on Uncertainty Quantification*, 2(1), 245–251. [doi.org/10.1137/130936233](https://doi.org/10.1137/130936233)
- Scholbeck, C. A., Casalicchio, G., Molnar, C., Bischl, B., & Heumann, C. (2022/2024). Marginal effects for non-linear prediction functions. *Data Mining and Knowledge Discovery*, 38, 2997–3042. [arxiv.org/abs/2201.08837](https://arxiv.org/abs/2201.08837) · [doi.org/10.1007/s10618-023-00993-x](https://doi.org/10.1007/s10618-023-00993-x)
- Semenova, L., Rudin, C., & Parr, R. (2022). On the existence of simpler machine learning models. *ACM FAccT 2022*, 1827–1858. [doi.org/10.1145/3531146.3533232](https://doi.org/10.1145/3531146.3533232)
- Song, E., Nelson, B. L., & Staum, J. (2016). Shapley effects for global sensitivity analysis: Theory and computation. *SIAM/ASA Journal on Uncertainty Quantification*, 4(1), 1060–1083. [doi.org/10.1137/15M1048070](https://doi.org/10.1137/15M1048070)
- Wooldridge, J. M. (2010). *Econometric Analysis of Cross Section and Panel Data* (2nd ed.). MIT Press. [mitpress.mit.edu/9780262232586](https://mitpress.mit.edu/9780262232586/)
