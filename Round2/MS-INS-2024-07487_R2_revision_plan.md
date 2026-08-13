# Revision Plan: MS-INS-2024-07487.R2 (Third Submission)

**Paper:** From Model Explanation to Data Misinterpretation: A Cautionary Analysis of Post Hoc Explainers in Business Research
**Decision on R1:** Major Revision (18-Jul-2026), DE Hemant Bhargava
**Goal:** A revision decisive enough that the review team's only remaining requests are cosmetic — i.e., a Minor Revision decision.
**Status:** Consolidated plan (merges both working drafts). Prepared for coauthor discussion, 20-Jul-2026.

---

## 1. Reading the decision letter strategically

The tone of this round is materially warmer than round 1. All three reviewers acknowledge substantial improvement; R1 says they "enjoyed reading the paper and learned from it"; R2 explicitly drops their novelty objection; R3 calls the Rashomon-agreement diagnostic "genuinely useful." The DE's hesitation is that he "would have liked a more definitive movement towards Accept/Reject," which tells us the bar for the next round is *decisiveness*: every AE point must be closed so completely that no reviewer can reasonably ask for another experiment.

The AE report is the roadmap. Its five points absorb nearly every reviewer comment, so we organize the plan around them. Two of the five (AE-1, AE-5) are described by the DE himself as "direct to address." Two (AE-2, AE-3) are conceptual and will decide the paper's fate. One (AE-4) is our opportunity: the AE is inviting us to elevate the Rashomon diagnostic into the paper's headline contribution.

A useful north star for the reframing that AE-1/AE-2/R3 are all circling: the paper's deepest point is not "SHAP and LIME are broken" but "predictive accuracy does not identify attribution." That is precisely the explanation/prediction distinction of Shmueli (2010), surfaced through a modern XAI pipeline, with Rashomon multiplicity (Breiman 2001; Fisher et al. 2019; Semenova et al. 2022) as the mechanism and Rashomon agreement as the practical diagnostic. If the revision commits to that framing consistently, AE points 1–4 largely resolve into one coherent story.

---

## 2. Theme-by-theme response options

For each theme we give the reviewer basis, two or three response options ordered by ambition, the methodological grounding for each, and a recommendation. Effort estimates assume the existing simulation infrastructure and cached SHAP explanations are reusable (per Appendix OA 6, regenerating SHAP from scratch is the binding compute constraint — roughly 27 server-days — so every option below is designed to avoid full regeneration).

### Theme 1 — Define the estimand and name what is failing (AE-1; R1 writing comments; R3 asks 1–2)

All three parties converge on the same demand: state precisely whether X→Y is associational or causal, and be precise about *which* component fails — the model, the data's identifiability, or the explainer.

**Option 1A (minimal — definitional surgery).** Add an early "Estimands and Scope" subsection declaring X→Y strictly associational: the objects of interest are properties of the conditional expectation E[Y|X] under the DGP (direction of local/average partial change; ranking of average partial effect magnitudes). Sweep the manuscript for causal verbs ("impact," "drives," "effect of X on Y") and replace or qualify them. Add the notation table R3 requested. Ground the associational/causal boundary in Shmueli (2010, *Statistical Science*), Zhao & Hastie (2021, *JBES*, "Causal Interpretations of Black-Box Models" — directly on point: it delineates when PDP/ICE-style objects admit causal readings), and Hernán & Robins (2020, *Causal Inference: What If*) for what we explicitly do *not* claim.

**Option 1B (recommended — failure-source taxonomy).** Everything in 1A, plus a short conceptual section decomposing misalignment into three separable sources, and then *tagging every empirical section by which source it isolates*:

1. **Approximation failure:** M ≉ G (Section 5.1 isolates this by varying accuracy).
2. **Identification failure:** accuracy does not pin down attribution — the Rashomon mechanism (Sections 5.2, 7 isolate this).
3. **Explainer-induced distortion:** the attribution method itself distorts even a faithful model — off-manifold perturbations, conditional-vs-interventional value functions, surrogate misspecification (the E1 experiment under Theme 3 isolates this).

This taxonomy is well supported in published work: Chen, Janizek, Lundberg & Lee (2020, "True to the Model or True to the Data?") for the value-function distinction; Kumar et al. (2020, ICML, "Problems with Shapley-value-based explanations as feature importance measures") for structural limits of Shapley attribution; Molnar et al. (2022, "General Pitfalls of Model-Agnostic Interpretation Methods," Springer LNAI) as an organizing citation. The taxonomy directly answers AE-1's "these are related issues, but they are not the same," and it converts a criticism into an expositional asset.

**Option 1C (ambitious).** 1B plus a formal proposition decomposing an alignment-error bound into the three terms. Attractive but risky: a loose bound invites new technical review. Recommend only if the derivation falls out cleanly in week 1; otherwise defer to future work.

**Recommendation: 1B.** Low compute, high leverage, and it is the connective tissue for Themes 3 and 5.

### Theme 2 — Benchmarks and practical significance (AE-2; R1 method comments; R3 minor)

R1 and R3 both ask: misaligned *relative to what*? And the AE asks: what decisions go wrong?

**Option 2A (minimal — regression baseline).** Run OLS/logit with standardized covariates on all 81 DGPs and score it with the identical direction/strength metrics. This is nearly free computationally and gives the reader the comparison R3 predicts: "a correctly specified regression returns a unique estimate and recovers the population relationship, whereas the post hoc explanation workflow need not." Report it as sharpening, not weakening, the message (R3's own words). Average-partial-effect machinery per Wooldridge (2010, ch. 2) grounds the estimand. Run this in two forms — a *naive* standardized regression (generically misspecified; what a researcher who skipped the ML pipeline would report) and a *correctly specified oracle* regression per DGP (the unique-estimate upper bound R3 describes) — so the reader sees both the floor and the ceiling.

**Option 2B (recommended — full benchmark suite + magnitude metric + consequences vignette).** 2A plus:

- **ML-native benchmarks:** partial dependence (Friedman 2001 — R1 explicitly requests it) and Accumulated Local Effects (Apley & Zhu 2020, *JRSS-B*). ALE is the right foil because it was designed specifically to remain valid under correlated features — the very condition we identify as the dominant driver of misalignment. If ALE also degrades under high correlation, that strengthens the identification story; if it degrades less than SHAP/LIME, that quantifies explainer-specific loss. Either result is publishable and both close AE-2.
- **A magnitude-accuracy metric** (R1): alongside Spearman, report normalized magnitude error of attributions versus true average partial effects (MAPE or normalized RMSE). R1 proposed MAPE by name; adopting it verbatim earns goodwill cheaply.
- **An aggregate direction variant** (R1): compute direction alignment at the population level (sign of the averaged explainer-implied change vs. sign of the averaged true change) and report alongside the instance-level metric. Frame the instance-level metric as the stricter notion and the aggregate one as matching common research practice. This defuses R1's "may overstate misalignment" concern without abandoning our metric.
- **A consequences vignette** (AE-2, R2): one worked example translating misalignment into a decision error — e.g., a manager allocates retention budget to the top-3 SHAP-ranked features; compute the realized lift versus targeting by the true top-3, and report the regret. One figure, one paragraph. This is the "why it matters" the AE asked for. Pair it with a plain-language translation of alignment scores (e.g., "0.8 direction alignment ⇒ roughly one sign reading in five is wrong").

**Option 2C (ambitious).** 2B plus a double/debiased ML partial-effects benchmark (Chernozhukov et al. 2018, *Econometrics Journal*) as the modern gold standard. Impressive but opens a causal-inference flank we just closed in Theme 1 (DML is usually motivated causally). Recommend against unless framed purely as an associational estimator.

**Recommendation: 2B.** The compute is light — all benchmarks are closed-form or cheap, run on existing datasets, and the metrics recompute from cached attributions.

### Theme 3 — Isolate explainer-specific failure from identification failure (R3's central ask; AE-1)

R3's sharpest request: "Show a case where the model faithfully recovers the data's associational structure yet the explainer still misleads." Unanswered, this comment keeps the paper at Major Revision, because R3 has said that without it the contribution "reads as a restatement of the well-known prediction vs inference gap."

**Option 3A (argumentative only).** Concede that most observed misalignment is identification-driven and reposition the contribution entirely around diagnosing it. Honest, but declines R3's explicit ask; not recommended alone.

**Option 3B (recommended — the E1 "explainer-on-G" experiment).** Apply SHAP and LIME *directly to the ground-truth generator G* — the limiting case of a perfectly estimated model (M = G), where approximation error is zero and the Rashomon set is a single point by construction. Score their direction/strength alignment against G's *analytic* partial responses and importance ordering, sweeping feature correlation ρ ∈ {0, .5, .9}. Any misalignment here is explainer-intrinsic with no confound: not the model (it is exact), not identification (the optimum is unique). This is the cleanest possible isolation of row 3 of the Theme-1 taxonomy, and it is cheap — we already apply explainers to G when constructing the strength ground truth I_G (§4.1), so the machinery exists; the only change is scoring E(G) against G's analytic structure rather than using it as the reference.

The *mechanism* is what makes this bite, and we should name it precisely — this is where the technical grounding lives. Under feature dependence, interventional SHAP evaluates the value function off the data manifold while conditional SHAP splits credit among correlated features, so even an exact model yields mis-ranked strength and, at high ρ, flipped local directions (Chen, Janizek, Lundberg & Lee 2020, "True to the Model or True to the Data?"; Aas, Jullum & Løland 2021, *Artificial Intelligence*; Janzing, Minorics & Blöbaum 2020, AISTATS). Sundararajan & Najmi (2020, ICML, "The Many Shapley Values") shows the choice among Shapley variants alone changes attributions — i.e., the explainer, not the model, injects the degree of freedom. For LIME, show local-coefficient sensitivity to kernel bandwidth and perturbation distribution on G (Hooker, Mentch & Zhou 2021, *Statistics and Computing*; Slack et al. 2020, already cited).

**Option 3C (recommended complement — faithful *fitted* model, matches R3's literal wording).** Because R3 asked specifically for a "model" that recovers structure yet is misexplained, pair E1 with a fitted-model version: on a DGP whose class we know, fit a *correctly specified* model (OLS or a GAM matching the DGP class), verify it recovers the structure (coefficients ≈ truth, unique optimum, trivially small Rashomon set), then apply SHAP/LIME to *that* model and reproduce the E1 failure. E1 gives the clean theoretical limit; 3C gives the concrete object R3 named. Together they leave no opening for "but that's not the pipeline." Small compute (one DGP family, one or two faithful models, existing explainer code).

**Recommendation: 3B (E1) as the headline result, with 3C as its fitted-model companion.** Both populate row 3 of the Theme-1 taxonomy, and they let us honestly grant R3 that rows 1–2 dominate in the wild — which we can now *say* because we have separated the rows empirically.

### Theme 4 — Novelty framing, prevalence candor, and the literature audit (AE-3)

**Option 4A (recommended core — candor + audit).**

- **Prevalence reframing:** lead with the leading-journal figures (14–17%; 16.1% across the union of lists) *before* the 42.5% full-sample figure, and motivate the paper as (i) prevention — the practice is diffusing from adjacent applied literatures toward top journals, and (ii) protection of the much larger practitioner audience the AE himself flags. R3 asked for exactly this candor; giving it prominently costs little because even ~1 in 6 papers in elite outlets is a defensible motivation.
- **Positioning paragraphs:** two tight paragraphs distinguishing our claim from Slack et al. (2020) (adversarial manipulability of explanations ≠ non-recoverability of DGP structure under good-faith use) and Fernández-Loría et al. (2022) (decision-level counterfactual critique ≠ data-level alignment evaluation). We committed to this in round 1; R2 has stood down, so the task is preservation, not persuasion. Anchor it with the "two questions" contrast: CS critiques test whether E(M) is faithful to M; we test whether E(M) is valid evidence about G — different question, different ground truth.
- **Full audit of the 181-paper review:** re-verify every classification, quotation, and severity code, per the AE's explicit warning that criticized authors will scrutinize us. Method: two coders independently re-code the full sample against the written codebook; report inter-rater reliability (Cohen's κ; interpret per Landis & Koch 1977, or Krippendorff's α per Krippendorff 2004); adjudicate disagreements with the documented LLM-assisted protocol already in OA 1; log every quotation against the source PDF. This is standard systematic-review practice (PRISMA-style transparency; Page et al. 2021, *BMJ*) and makes the review bulletproof. Note: the AE asked for a *recheck*, not a submitted evidence sheet — the coding log is internal insurance; the response letter simply states the recheck was performed.

**Option 4B (defensive supplement).** Additionally soften the in-text presentation: keep the anonymized "boxes" of prototypical claims in the main paper, and move the citable, named examples to an appendix table with exact quotations and page numbers so the evidence is verifiable but the main text reads as pattern-documentation rather than indictment. Recommended — it lowers litigational temperature while *increasing* verifiability, which is exactly the AE's dual demand ("fair" and "defensible").

**Recommendation: 4A + 4B.** Pure labor, no compute; the audit is the schedule's critical path for human-hours, so it starts week 1.

### Theme 5 — Elevate the Rashomon-agreement diagnostic into a usable tool (AE-4; R3; R1 endorsement)

The AE calls this "the most promising part of the paper" and asks two things: make it *usable*, and position it against the multiplicity literature.

**Option 5A (minimal — positioning only).** Add a focused related-work passage distinguishing our diagnostic from: model class reliance and Rashomon-set variable importance (Fisher, Rudin & Dominici 2019, *JMLR*, "All Models Are Wrong, but Many Are Useful"); variable importance clouds (Dong & Rudin 2020, *Nature Machine Intelligence*); predictive multiplicity (Marx, Calmon & Ustun 2020, ICML; Watson-Daniels, Parkes & Ustun 2023, ICML; Hsu & Calmon 2022, already cited); and the Rashomon Importance Distribution (Donnelly et al. 2023, NeurIPS, already cited). The delta to state exactly: prior work *characterizes or aggregates* cross-model disagreement; we *validate agreement as an empirical predictor of ground-truth alignment* (Table 1's 0.79 correlation) and derive an ex-ante reliability signal from it. That validation-against-known-DGP step is, to our knowledge, new — say so in one sentence and no more.

**Option 5B (recommended — the protocol).** 5A plus operationalize the diagnostic as a boxed, numbered procedure a researcher can follow — a "Rashomon Agreement Protocol":

1. Fit K diverse models within an ε-accuracy band of the best model (we used K = 10, ε = 3%).
2. Compute each model's global attribution profile; take pairwise Spearman correlations; average → agreement score A.
3. Compare A to calibrated reference points from our 81-DGP experiments — e.g., report the empirical probability that strength alignment exceeded 0.9 conditional on A falling in each band (a simple calibration curve or ROC-style figure derived from data we already have).
4. Report A alongside any explanation-based claims; treat low-A settings as hypothesis-generation-only.

Everything needed for the calibration figure already exists in the Section 7 experiment outputs — this is analysis, not new simulation. The protocol box plus one calibration figure converts "interesting empirical result" into "something a researcher can apply," which is the AE's phrasing verbatim.

**Option 5C (ambitious — real-data case study).** 5B plus a short end-to-end application on one public business dataset (e.g., a churn or marketing-response dataset): run the pipeline, compute A, show what a researcher would report and how the conclusion changes in a low-A versus high-A regime. This is the single strongest possible answer to AE-4 and gives the paper the applied anchor R1's "so what" comment gestures at. Moderate compute (one dataset, K = 10 models, SHAP on a subsample — days, not weeks, on the existing server).

**Recommendation: 5B firmly; adopt 5C if the compute finishes by week 7.** If 5C threatens the timeline, hold it in reserve for the response letter as "additionally, we applied the protocol to..." — but note the AE has signaled this is where the accept-level contribution lives, so 5C deserves real consideration.

### Theme 6 — Consistency, structure, figures, prose (AE-5; R2 throughout; R1 structural suggestions)

Not optional and not option-based — an execution checklist. Every item below is directly traceable to a reviewer sentence.

- **Reconcile Figure 3 vs. Figure 4** (R2): explain — or eliminate — the apparent contradiction that SHAP dominates LIME overall (Fig 3, best models) while the top accuracy bracket shows them similar (Fig 4). If the resolution is that Fig 3 uses best-model-per-dataset while Fig 4 bins across deliberately degraded models, say so in both captions; if the numbers genuinely conflict, rerun and fix.
- **Figure overhaul** (R2): common y-axis ranges within each figure family (Figs 4, 6, 7, 8), gridlines, and self-contained captions that define plotted quantities and annotate that boxed numbers are univariate OLS slopes with significance stars.
- **Fix the "Feature 0" reference** (R2): the text on p. 25 references an index absent from Figure 5a; align text and figure labels.
- **Linearize the empirical narrative** (R2) and adopt R1's structural suggestion: move Section 5.3 (data properties → alignment) into Section 7.3's discussion or an appendix, and compress Section 6 (robustness) to a summary paragraph with details in OA 4. This simultaneously satisfies R1 ("remove or move Sections 5.3 and 6") and R2 ("more linear structure"). Resulting flow: metrics → alignment results (best models) → why (accuracy → multiplicity → diagnostic) → robustness summary → guidance.
- **Notation table** (R3, twice-requested across rounds — do not ship without it).
- **Full consistency audit:** main-vs-appendix terminology (e.g., "directionality"/"direction alignment," ν vs. 𝜈 usage, Δ definition typo "{−, 3, 3}" in OA 4.3), figure/caption/number cross-checks, and the specific typos R2 lists (p. 3 "made substantive SHAP or LIME," Fig 3c "alighment").
- **Professional proofreading pass** at the end (R2 requested it in both rounds).
- **Causal-language sweep** merges with Theme 1's estimand work (R1).

---

## 3. Recommended package (for discussion)

Adopt: **1B + 2B + 3B/3C (E1 explainer-on-G + faithful-model companion) + 4A/4B + 5B (5C strongly considered) + Theme 6 in full.**

This package answers every AE point with new evidence rather than argument alone, adds four bounded new analyses (regression/PDP/ALE benchmarks; magnitude + aggregate-direction metrics; the E1 explainer-on-G experiment plus its fitted-model companion; agreement calibration), none of which requires regenerating the expensive SHAP corpus, and repositions the paper around its strongest asset per the AE's own signal. The main scope risk is 5C; the main schedule risk is the 181-paper audit, which is human-time bound and must start immediately.

---

## 4. Master change table (priority-ranked)

Priority key — **P1 (critical):** DE-flagged make-or-break items (AE-2, AE-3) plus anything the decision plausibly turns on. **P2 (high):** explicitly requested by a reviewer; skipping needs a written justification in the response. **P3 (medium):** polish and goodwill; cheap, do them all. "Theme" ties each row to §2.

| # | Change | Theme | Type | Comments | Priority |
|---|---|---|---|---|---|
| 1 | Formal estimand definition (X→Y strictly associational) | T1 | Framing | AE1, R3.2, R1.8 | **P1** |
| 2 | Three-way failure taxonomy (approximation / identification / explainer) | T1 | Framing | AE1, R3.1, R3.7 | **P1** |
| 3 | **E1: explainer-on-G experiment** (M=G limit; score vs G's analytic structure; sweep ρ) | T3 | Experiment | R3.6, R3.1, AE1 | **P1** |
| 4 | Faithful fitted-model companion (correctly specified model, then explain it) | T3 | Experiment | R3.6, R3.2 | **P2** |
| 5 | Naive standardized-regression benchmark (the floor) | T2 | Experiment | AE2, R1.1, R1.3, R1.6 | **P1** |
| 6 | Oracle correctly-specified regression benchmark (the ceiling) | T2 | Experiment | R3.5, R3.7, AE2 | **P1** |
| 7 | Partial dependence benchmark | T2 | Experiment | R1.3 | **P2** |
| 8 | ALE benchmark (correlated-feature foil) | T2 | Experiment | AE2, R3.5 | **P2** |
| 9 | Two-questions positioning: E(M)↔M faithfulness vs E(M)↔G validity | T4 | Writing | AE3, R2.1 | **P1** |
| 10 | Multiplicity-literature positioning (MCR, VIC, predictive multiplicity) | T5 | Writing | AE4, R3.3, R3.8 | **P1** |
| 11 | 181-paper audit (dual coding, κ/α, quotation verification) | T4 | Process | AE3 | **P1** |
| 12 | Candid prevalence framing (lead with ~16% top-journal rate) | T4 | Writing | AE3, R3.4 | **P2** |
| 13 | Rashomon Agreement Protocol (boxed procedure) | T5 | Analysis+writing | AE4, R1.2, R1.6 | **P1** |
| 14 | Agreement→alignment calibration curve (from §7 outputs) | T5 | Analysis | AE4, R1.1 | **P2** |
| 15 | Real-data diagnostic case study | T5 | Experiment | AE4, R2.2, R1.2 | **P2** |
| 16 | Code release (protocol package/notebook) | T5 | Process | AE4 | **P3** |
| 17 | Aggregate-direction metric (population-level sign) | T2 | Experiment | R1.4 | **P2** |
| 18 | Magnitude strength metric (MAPE / normalized RMSE) | T2 | Experiment | R1.5 | **P2** |
| 19 | Managerial consequences vignette (retention-budget regret) | T2 | Writing+analysis | AE2, R2.2, R1.2 | **P2** |
| 20 | Decision-relevant translation of alignment scores | T2 | Writing | AE2, R2.2 | **P2** |
| 21 | "What to do instead" decision tree | T2 | Writing | R1.2, AE2 | **P2** |
| 22 | Fig 3 vs Fig 4 reconciliation | T6 | Analysis | R2.3, AE5 | **P2** |
| 23 | Figure overhaul (uniform y-axes, gridlines, self-contained captions) | T6 | Presentation | R2.4, R2.5, AE5 | **P2** |
| 24 | Linear restructure (§5.3 → §7.3/appendix; §6 → summary) | T6 | Writing | R2.5, R1.10, R1.7, AE5 | **P2** |
| 25 | Consistency audit (main vs appendix; numbers/terms/figures) | T6 | Process | AE5, R2.5 | **P2** |
| 26 | Notation table | T6 | Writing | R3.5 | **P3** |
| 27 | Shmueli (2010) explain-vs-predict framing | T1 | Writing | R1.9 | **P3** |
| 28 | De-duplication + professional proofread | T6 | Writing | R1.7, R2.6, AE5 | **P3** |

---

## 5. Timeline

Target: submission in **12 weeks (week of 12-Oct-2026)**, comfortably inside the journal's one-year window and fast enough to retain the same review team's momentum. Assumes three coauthors plus the existing compute server; audit assumes two coders working in parallel.

| Weeks | Workstream | Owner-type | Deliverable |
|---|---|---|---|
| 1 | Kickoff: lock option package; draft estimand section + failure taxonomy (1B); notation table | Lead author | New Section 2.x draft |
| 1–4 | 181-paper re-audit: dual coding, κ/α reliability, quotation verification, appendix evidence table (4A/4B) | 2 coders | Audited OA 1 + evidence appendix |
| 2–4 | Benchmark suite: naive + oracle regression, PDP, ALE on 81 DGPs; MAPE + aggregate-direction metrics recomputed from cached attributions (2B) | Empirical lead | Benchmark figures + tables |
| 2–5 | E1 explainer-on-G experiment: explainers applied to G, scored vs analytic structure; interventional/conditional SHAP mechanism + LIME bandwidth sensitivity; faithful fitted-model companion (3C) | Empirical lead | New subsection + 1–2 figures |
| 4–7 | Rashomon protocol: positioning passage, protocol box, calibration curve from Section 7 outputs (5B); launch 5C case study in background by week 5 if adopted | Lead + empirical | Revised Section 7 |
| 5–6 | Consequences vignette (regret example) + score translation (2B) | Any | One figure + paragraph |
| 6–8 | Rewrite Intro, Related Work (positioning + multiplicity literature), Conclusion; prevalence reframing; restructure per Theme 6 (move 5.3, compress 6) | Lead author | Full-text draft v3 |
| 8–10 | Figure overhaul; Fig 3/Fig 4 reconciliation; full main-appendix consistency audit | Empirical + one fresh-eyes coauthor | Camera-consistent draft |
| 9–11 | Point-by-point response letter (every comment quoted, response, manuscript location); summary of changes; professional proofread | All | Response package |
| 12 | Buffer: cold read-through by the coauthor least involved in drafting; final checks; submit | All | Submission |

**Aggressive variant (9 weeks):** drop 3C and 5C, run the audit with three coders, overlap the rewrite with the audit's tail. Feasible but leaves 5C — the AE's clearest accept-lever — on the table; we'd recommend the 12-week plan.

**Compute note:** nothing above requires regenerating the 27-day SHAP corpus. E1 reuses the existing explainer-on-G machinery from §4.1's I_G construction; new SHAP work is confined to the small 3C fitted-model DGP and (if adopted) the single 5C dataset — both fit inside normal week-long windows on the existing 80-core server.

---

## 6. Response-letter strategy

- Quote every comment verbatim; respond point-by-point with manuscript page/section pointers (the MS instructions require this and the AE will check coverage of non-highlighted reviewer details too).
- Open the AE response by mapping our five biggest changes onto his five points, one-to-one.
- Concede early and cleanly where R3 is right: yes, the dominant failure mode is identification, and the revision now says so explicitly — then show the taxonomy and the E1 experiment as what the paper adds *beyond* that concession. Reviewers reward candor here far more than defense (and R3 has essentially drafted our framing for us: "sharpens rather than weakens the paper's message").
- Flag R1's adopted suggestions by name (MAPE, aggregate direction, PDP benchmark, Shmueli framing, restructuring) — R1 handed us a checklist; visibly completing it invites a "my concerns are addressed" review.
- Where we decline a suggestion (e.g., DML benchmark, formal bound), give a one-paragraph reasoned decline with a citation, never silence.

---

## 7. Comment-to-change traceability map (for drafting R2R.tex)

Each reply slot in `Round2/R2R.tex` maps to change-table rows in §4.

| Reply | Change rows |
|---|---|
| SE (both points) | Overview of package; highlights of rows 1–16 |
| AE 1 | 1, 2, 3, 4 |
| AE 2 | 5, 6, 7, 8, 19, 20 |
| AE 3 | 9, 11, 12 |
| AE 4 | 10, 13, 14, 15, 16 |
| AE 5 | 22, 23, 24, 25, 26, 28 |
| R1.1 | 5, 6, 14 |
| R1.2 | 13, 15, 19, 21 |
| R1.3 | 5, 7 |
| R1.4 | 17 |
| R1.5 | 18 |
| R1.6 | 5, 13 |
| R1.7 | 24, 28 |
| R1.8 | 1 |
| R1.9 | 27 |
| R1.10 | 24 |
| R2.1 | 9 |
| R2.2 | 15, 19, 20 |
| R2.3 | 22 |
| R2.4 | 23 |
| R2.5 | 23, 24, 25 |
| R2.6 | 28 |
| R3.1 | 2, 3 |
| R3.2 | 1, 4 |
| R3.3 | 10, 13 |
| R3.4 | 12 |
| R3.5 | 6, 8, 26 |
| R3.6 | 3, 4 |
| R3.7 | 2, 6 |
| R3.8 | 10 |

---

## 8. Open questions for this morning

1. Do we commit to the 5C real-data case study now, or hold as reserve? (Recommendation: commit — it is the accept-lever.)
2. Who are the two audit coders, and can they start this week?
3. Do we run both E1 (explainer-on-G) and its fitted-model companion (3C), or is E1 alone sufficient? (Recommendation: both — 3C matches R3's literal wording, E1 is the clean limit.)
4. Any coauthor concerns with leading the prevalence framing with the top-journal (16.1%) figure?
