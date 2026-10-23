# Rung 0 — measure the error rate before trying to change it

*A one-page research protocol for the clinical team (MVA Hackathon 2026, Track 2). It gives no treatment and changes
nothing in the child's care; it needs the clinical team's judgement, ethics approval and the family's consent. Code
and assumptions: `scripts/t2/nof1.py`; results: `reports/t2/nof1.md`.*

**Why.** Every candidate in the Track 2 report needs a readout that moves in weeks, and the usual ones do not: a
lymphocyte karyotype reports aneuploidy accumulated over years. The frequency of **micronucleated reticulocytes
(%MN-RET)** counts micronuclei — the direct product of chromosome missegregation — in red cells that circulate for
1–2 days, so it reports a *recent* error rate [1, 2]. It has a paediatric reference (upper bound of the 90 %
tolerance interval in controls, **0.42 %**) [3], it rose within 3–5 days in 9 of 10 adults starting chemo- or
radiotherapy [4], and it was 5.2-fold raised in untreated high-grade glioma [7] — but it has **never been measured in
MVA**. Two questions
first: is it raised in this child, and how much does it vary between visits (σ)? Every n-of-1 design depends on σ.

| Item | Specification | Basis |
|---|---|---|
| Samples | 6 over 8 weeks (weeks 0, 1, 2, 4, 6, 8) plus 1 spare | about 10 % of paediatric samples failed QC [3] |
| Material | 240 µL of whole blood from a K2-EDTA tube, ideally drawn at a routine blood test (no extra venepuncture) | [3] |
| Processing | fixed in ultracold methanol per kit instructions, stored at −75 to −85 °C, **all samples analysed in one batch** | [3]; batching removes run-to-run variation |
| Events | ≥ 20,000 CD71⁺ reticulocytes per sample (minimum 10,000); **≥ 100,000 if the baseline is below 1 %** | [3]; below 1 % counting error starts to cost power |
| Off treatment | *suggested:* no cytotoxic chemotherapy or radiotherapy in the previous 4 weeks (to agree with the clinical team) | the rise peaks 3–5 days after starting therapy, up to 26-fold [4]; time back to baseline is not established |
| Record each visit | infection (fever, CRP), nutritional status, folate and B12, reticulocyte count, transfusions | infection and malnutrition raise %MN-RET in children [5]; folate depletion raises micronuclei [6] |

**Pre-declared analysis.** Report the geometric-mean %MN-RET with its 90 % interval and σ̂, the SD of log %MN-RET
between visits. **Plan with the upper 80 % bound of σ̂, not the point estimate**: with six samples σ̂ is
imprecise. The trial itself re-estimates σ (t posterior), but how many visits it needs depends on the true σ.

| Result | Reading | Next step |
|---|---|---|
| **A** · mean ≤ 0.42 % | the endpoint does not see this child's instability in blood | do not use it for an n-of-1; report it as the first MVA measurement — it is also evidence against the disease model the measurement plan rests on |
| **B** · raised, σ̂ ≤ 0.5 | an n-of-1 on this endpoint is feasible | size it from the table below, reading the column at or above the upper bound of σ̂ |
| **C** · raised, σ̂ > 0.5 | too noisy for a decision within weeks | sample more often, or choose another endpoint |

If an mTOR inhibitor is ever started (tier 3), the same assay becomes the safety monitor, because rapamycin was
aneugenic in human lymphocytes in one in-vitro study [8]: pre-declare a stop rule, e.g. P(increase ≥ 25 %) ≥ 90 %.

**Visits per condition (on/off) for 80 % power**, 20,000 reticulocytes per sample; each cell gives the visits for a
true 35 % / 50 % reduction, «—» = more than 30.

| Responder rule (false-positive rate) | σ = 0.20 | σ = 0.35 | σ = 0.50 |
|---|---|---|---|
| reduction ≥ 25 %, 95 % sure (< 1 %) | 28 / 5 | — / 11 | — / 20 |
| any reduction, 95 % sure (4 %) | 4 / 3 | 9 / 5 | 18 / 8 |
| any reduction, 90 % sure (10 %) | 3 / 2 | 7 / 4 | 13 / 6 |

The rule is a choice for the clinical team and the family, **made before the first sample**: a strict rule needs a
large effect or many visits, because its threshold sits close to the effect one hopes for; a loose rule accepts more
false positives.

**What this is not.** Not a validated diagnostic test (a first disease-activity observation); not validated in MVA; it reads missegregation in the erythroid lineage,
not in muscle or other tissues; a fall in %MN-RET is a pharmacodynamic signal, not proof of clinical benefit.

**References (PMID).** [1] Dertinger 2003, 14644356 · [2] Dertinger 2011, 21164195 · [3] Baig 2020, 32761646
(paediatric) · [4] Dertinger 2007, 17059891 · [5] Cervantes-Ríos 2012, 22119781 · [6] Everson 1988, 3367390 ·
[7] Barrera 2018, 30335243 · [8] Bonatti 1998, 9914383.
