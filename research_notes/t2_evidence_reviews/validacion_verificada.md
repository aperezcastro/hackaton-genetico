# How a repurposed drug is actually validated for an ultra-rare disease — verified citations

Coordinator's verification pass (2 Oct 2026). Every PMID/NCT below was retrieved from PubMed or ClinicalTrials.gov
in this session. Items the subagent cited from memory and that did **not** verify are listed at the end.

## 1. Verified precedents: the realistic ladder

| Precedent | Citation (verified) | Evidence that justified use | Why it is relevant to us |
|---|---|---|---|
| **Sirolimus in idiopathic multicentric Castleman disease (iMCD-TAFRO)** | Fajgenbaum DC et al., *J Clin Invest* 2019 — **PMID 31408438**; trial **NCT03933904** | Proteomic/mechanistic signal (PI3K/AKT/mTOR activation) measured **in the patients themselves** → an already-approved drug → sustained remission → then a prospective trial | **The closest model to our project**: measure a pathway in the patient, repurpose an approved drug against it, then formalise it. Same drug class we propose for H3 |
| **Milasen (individualised ASO)** | Kim J et al., *N Engl J Med* 2019;381:1644 — **PMID 31597037** | Patient-specific molecular diagnosis → drug designed and tested **in the patient's own cells** → single-patient IND → treatment within ~1 year | Established that patient-derived-cell evidence can support first-in-human use in one child |
| **Framework for individualised splice-switching ASOs** | Kim J et al., *Nature* 2023 — **PMID 37438524** | Generalises the milasen route into a reusable framework | Shows regulators/field accept a *platform* logic, which is what our gene-agnostic pipeline offers |
| **Lonafarnib in Hutchinson-Gilford progeria** | Gordon LB et al., *JAMA* 2018;319:1687 — **PMID 29710166**; trials NCT00425607, NCT00916747; everolimus+lonafarnib NCT02579044 | Single-arm treated cohort vs **natural-history external controls**; mortality 3.7 % vs 33.3 % | **The key design lesson**: with no possible placebo arm, a rigorously assembled natural-history control can support approval (Zokinvy, 2020). Also a progeroid disease, like the BubR1 phenotype |
| **Patient-specific in vivo base editing (CPS1 deficiency)** | Musunuru K et al., *N Engl J Med* 2025 — **PMID 40373211** | Bespoke therapy for a single infant | Current outer edge of n-of-1 individualised therapy |
| **Baricitinib in paediatric autoinflammatory interferonopathies** | Sanchez GAM et al., *J Clin Invest* 2018 — **PMID 29649002**; trials NCT03921554 (AGS), NCT04517253 (CANDLE/NNS, paediatric) | Open-label, mechanism-driven JAK inhibition in children with a chronic interferon signature | Direct precedent for our data-driven H6 line: JAK inhibition **is** done in children when interferon drives the disease |
| **NR in Li-Fraumeni syndrome** | **NCT03789175**, COMPLETED, **enrolment = 1** (true n-of-1), NCI | NR effect on skeletal-muscle mitochondrial function in one patient with a cancer-predisposition syndrome | **Highly relevant precedent**: a NAD⁺ precursor has already been given, under a protocol, to a patient with a cancer-predisposition syndrome — and as a formal n-of-1 |

## 2. Verified picture of the therapeutic area

- **No clinical trial has ever been registered for mosaic variegated aneuploidy** (ClinicalTrials.gov search for
  "mosaic variegated aneuploidy OR premature chromatid separation": **0 studies**). Everything here is first-in-disease.
- Progeria, by contrast, has **12 registered studies**, including one of everolimus in children (NCT02579044). The
  progeria community is the operational template: patient registry → natural history → single-arm trial → approval.
- **BUB1B has no drugs in DGIdb and no known drugs in Open Targets.** Open Targets tractability for BUB1B: small
  molecule "High-Quality Ligand" and "Druggable Family", but every compound annotated against the checkpoint genes
  (BUB1, TTK/MPS1) is an **inhibitor**, i.e. the wrong direction for a patient who already has too little BUBR1.
  This is a strong, citable justification for working through mechanism and consequences instead of "gene → drug".

## 3. Verified human senolytic evidence (tempering our H4)

| Trial | Citation (verified) | Result |
|---|---|---|
| D+Q in postmenopausal women, phase 2 RCT, n=60 | Farr JN et al., *Nat Med* 2024 — **PMID 38956196** | **Primary endpoint not met**: change in CTx did not differ (D+Q −4.1 % vs control −7.7 %, P = 0.611). Bone-formation marker P1NP rose transiently at 2 and 4 weeks (+16 %), not sustained |
| D+Q in diabetic kidney disease, n=9, open-label | Hickson LJ et al., *EBioMedicine* 2019 — **PMID 31542391** | Senescent-cell markers fell in adipose tissue and skin; no placebo group |
| D+Q in idiopathic pulmonary fibrosis, first-in-human pilot | Justice JN et al., *EBioMedicine* 2019 — **PMID 30616998** | Physical function improved; **lung function unchanged** |
| D+Q in aged mice | Xu M et al., *Nat Med* 2018 — **PMID 29988130** | Improved physical function and increased lifespan |
| ABT-263/navitoclax in aged mice | Chang J et al., *Nat Med* 2016 — **PMID 26657143** | Cleared senescent cells, rejuvenated haematopoietic stem cells |
| Genetic clearance in **BubR1-hypomorphic** mice | Baker DJ et al., *Nature* 2011 — **PMID 22048312** | Delayed sarcopenia, cataract, lipodystrophy — but via the **INK-ATTAC transgene, not a drug** |

**Conclusion for the report:** the best human senolytic RCT to date **missed its primary endpoint**, and no senolytic
drug has ever been tested in a BubR1 or aneuploidy model. Our H4 therefore stays where we put it: a preclinical
experiment to run, not a treatment to propose. This strengthens the honesty of the report rather than weakening it.

## 4. Context verified for MVA itself

- **Malumbres M, *Nat Rev Genet* 2024 — PMID 39169218**, "Mosaic variegated aneuploidy in development, ageing and
  cancer": MVA involves "cellular and systemic stress, including **chronic immune response** and the associated
  clinical implications". This is independent, human-level support for the interferon signature we found in LINCS,
  and therefore for H6 being a real feature of the disease rather than a cell-line artefact.
- **TRIM37 / Mulibrey nanism** can also present with mosaic variegated aneuploidy (*J Med Genet* 2026, PMID 41702694).
  Checked against our data: no rare functional variant in TRIM37, Exomiser score 0. Differential excluded.
- BubR1 insufficiency damages the **heart** (*Aging Cell* 2025, PMID 40607964) and **axon myelination/motor function**
  (*Aging* 2016, PMID 27922816) in mice → two further arguments for the safety card, since anthracyclines are
  cardiotoxic and vincristine is neurotoxic.

## 5. Not verified — do not cite without checking

- Hawksworth 2024 *Trials* review of 74 randomised n-of-1 trials (DOI given from memory; did not resolve in PubMed).
- EMA compassionate-use guideline date/reference, EMA small-populations guideline code.
- FDA Rare Disease Evidence Principles (Sep 2025) and the Feb 2026 "plausible mechanism" draft guidance: read only
  through secondary summaries. **If cited in the final report, cite them as secondary or fetch the FDA page first.**
- Final (vs draft) status of the 2021 individualised-ASO guidances.
- The subagent attributed the 2023 *Nature* ASO framework to "Yu"; the verified first author is **Kim J**.
