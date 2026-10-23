# Read-through of the BUB1B nonsense allele: allele-specific context analysis + verified evidence

Coordinator's analysis and verification pass (2 Oct 2026). The sequence work below uses only the **public reference
CDS** of BUB1B (Ensembl ENST00000287598, 3,153 nt / 1,050 aa + stop) and the single known coding change; no genotype
table or patient file was used.

## 1. The patient allele's stop-codon context (computed here)

```
reference codons 735-740 : CCA GAG TTA AGT GCC TCT
mutant    codons 735-740 : CCA GAG TGA AGT GCC TCT
                                   ^^^ premature stop (Leu737Ter)
stop + 4 nucleotide      : TGA | A      →  "UGA-A"
```

| Feature | Value | Read-through outlook |
|---|---|---|
| Stop codon identity | **UGA** | **Favourable** — UGA is the most read-through-permissive stop codon |
| +4 nucleotide (first base after the stop) | **A** | **Unfavourable** — the +4 preference is **C > U > A/G**, so A is in the weakest group |
| Position in transcript | exon 17 of 23, >50 nt upstream of the last exon–exon junction | **NMD predicted**, so most of the template is destroyed before translation |

**Verified source for the context rule:** Mangkalaphiban K et al., "Extended stop codon context predicts nonsense
codon readthrough efficiency in human cells", *Nat Commun* 2024 — **PMID 38509072**. Machine-learning analysis of
HEK293T ribosome-profiling data; identity of the stop codon, its context and 3'-UTR length are conserved
determinants.

**Net conclusion:** this allele sits in a **mixed-to-poor** read-through context — the right stop codon, the wrong
neighbouring base, and NMD on top. This is a concrete, allele-level argument that H2 should be **downgraded**, and it
is exactly the kind of analysis that distinguishes a specific proposal from a generic one.

## 2. Verified evidence against our previous H2 drugs

| Drug | Verified finding | Consequence for our report |
|---|---|---|
| **Ataluren (PTC124)** | McElroy SP et al., *PLoS Biol* 2013 — **PMID 23824517**: "A lack of premature termination codon read-through efficacy of PTC124 (Ataluren) in a diverse array of reporter assays". Its original activity is attributed to an **off-target effect on the firefly luciferase reporter** used to discover it. Phase 3 in cystic fibrosis failed (PMID 24836205); Cochrane found evidence insufficient (PMID 36866921); the ACT DMD trial missed significance (13.0 m, p = 0.213, PMID 28728956); EMA's CHMP did not renew its conditional authorisation. | **Remove ataluren from the candidate list.** Keeping it would be a rigour failure. |
| **Amlexanox** | In the SMG1-inhibitor study (McHugh DR et al., *Int J Mol Sci* 2021, DOI 10.3390/ijms22010344) amlexanox showed **no synergy** with the SMG1 inhibitor and **did not increase G542X CFTR mRNA**; in a factor V nonsense study (PMID 37866515) amlexanox produced **no measurable effect** while DAP gave 10.8-fold and G418 7.0-fold increases. Its action is mainly NMD-related, not read-through. | **Downgrade amlexanox**; present it as a weak, NMD-directed option, not as "the drug for the null allele". |
| **Gentamicin / aminoglycosides** | Real read-through activity (G418 is the laboratory benchmark), but our own pipeline already excluded gentamicin for this patient (boxed nephrotoxicity vs nephrocalcinosis). | Keep excluded, now with a second reason: efficacy in this context is unlikely to justify the renal risk. |
| **2,6-diaminopurine (DAP)** | Trzaska C et al., *Nat Commun* 2020 — **PMID 32198346**: a highly potent corrector **specifically of UGA** nonsense mutations, acting by inhibiting FTSJ1. Strongest effect in the factor V comparison. | **Not an approved drug** → outside a repurposing-of-approved-drugs submission, but the single most mechanistically apt molecule for a UGA allele. Belongs in "future directions", not in the candidate table. |
| **ELX-02** | Investigational; results available only through company/SEC releases, no peer-reviewed phase 2; monotherapy sweat-chloride change 5.4 mmol/L; no significant benefit added to ivacaftor. | Mention as investigational only; do not score. |

## 3. Why read-through attempts generally fail (verified themes)

- **Assay artefacts** (the ataluren/luciferase story) — the reason to demand a non-luciferase readout in any go/no-go.
- **NMD destroys the template**, so read-through alone has little substrate; combination with NMD inhibition (SMG1
  inhibitors) is what works in models, and SMG1 inhibitors are **toxic at high doses** and not approved.
- **Mutation- and patient-specific responses**: only ~23 % of NF1 patient lines responded to ataluren.
- **Narrow therapeutic window** and cytotoxicity for the effective compounds (G418, DAP).
- Most published effects are reported as **fold-change**, not percent of wild-type protein, which hides whether the
  restored amount could ever cross a biological threshold.

## 4. What this means for our recommendation

H2 moves from "Tier 4, exploratory candidate" to **"mechanistically attractive but currently undruggable with
approved agents"**:

1. The target allele has the right stop codon but a poor +4 context and predicted NMD.
2. The only approved-ish candidates are discredited (ataluren) or weak and off-mechanism (amlexanox).
3. The compounds that actually work on UGA (DAP, G418, SMG1 inhibitors) are **not approved** and some are toxic.
4. Therefore the honest recommendation is: **do not propose a read-through drug now**; instead, define the
   experiment that would make it decidable — a non-luciferase reporter carrying the native BUB1B UGA-A context,
   ± NMD inhibition, in patient-derived cells — and watch DAP-class compounds as they develop.

This is a line we are closing on evidence, and saying so is a strength of the submission.

## 5. Gaps

- No study of read-through of any BUB1B nonsense variant exists (searched; none found).
- Whether this specific exon-17 stop triggers NMD is inferred from the standard >50 nt rule, not measured.
- Some of the subagent's secondary citations (PMID 22938201 for amlexanox, the ataluren regulatory history) were
  taken from search snippets; the primary claims above were verified directly.
