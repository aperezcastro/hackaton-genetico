# Pharmacogenomic safety note — biallelic *BUB1B* (mosaic variegated aneuploidy type 1)

*For the clinical team. Hypothesis-based guidance from the MVA Hackathon 2026 (Track 2); it does not replace
specialist judgement and should be reviewed by a paediatric oncologist and a clinical pharmacologist before use.
Evidence: **C** = case reports in MVA (number of cases given), **G** = guideline/consensus, **M** = mechanism or
cell/mouse data, **T** = transcriptomic connectivity (this work), **L** = FDA label.*

**The risk of this note is under-treatment.** It flags drugs for discussion and monitoring. It never recommends
replacing or de-intensifying a curative regimen on its own: that decision belongs to the treating oncologist.

**Why this matters.** The child has a weak spindle assembly checkpoint (low functional BUBR1). For drugs that act
through mitotic arrest the prediction runs **both ways**: checkpoint-weak cells may slip out of arrest and survive
with more missegregation (less tumour kill, more chromosomal damage in normal dividing tissue), or fail to sustain
arrest and die faster (more sensitivity and more toxicity). Which dominates in this child is unknown; both argue for
closer monitoring, not for avoidance. Mouse data add one organ-specific caution: BubR1 insufficiency causes cardiac
hypertrophy, fibrosis and senescence (PMID 40607964). BubR1-hypomorphic mice also show central-nervous-system
hypomyelination with motor deficits (PMID 27922816); peripheral nerves were not studied, so this does not by itself
predict vincristine neuropathy.

**A · Drugs this child may receive in oncology or general care**

| Drug / class | Concern | Suggested action | Evidence |
|---|---|---|---|
| **Curative chemotherapy as a whole** (VAC/VI and others) | one MVA child with rhabdomyosarcoma was treated with reduced-intensity chemotherapy; one biallelic *BUB1B* child had a severe reaction to dactinomycin | do not change a curative regimen on this basis alone; agree intensity and monitoring with the treating oncologist | C (2) |
| **Vinca alkaloids** (vincristine, vinblastine, vinorelbine) | checkpoint-dependent mechanism; direction of effect unresolved (see above) | standard dosing unless the oncologist decides otherwise; close monitoring of neuropathy and counts | M |
| **Taxanes, eribulin** | same mechanism | flag for discussion; no MVA-specific data | M |
| **Dactinomycin** | severe adverse reaction in one biallelic *BUB1B* child | close monitoring; dose decisions with the oncologist | C (1) |
| **Anthracyclines, topoisomerase inhibitors** (doxorubicin, irinotecan, etoposide) | reproduce the *BUB1B*-loss cellular state in LINCS; BubR1-insufficient mouse hearts are hypertrophic and fibrotic | cardiac monitoring per protocol, with a low threshold for echocardiography | M, T |
| **Ionising radiation** | exquisite sensitivity reported in several genomic-instability disorders (mainly DNA-repair defects); unknown for MVA | specialist review of fields and dose | G (indirect) |
| **Haematopoietic stem cell transplant** (if marrow failure or MDS arises) | one MVA1 report: reduced-intensity conditioning tolerated, then graft rejection | reduced-intensity conditioning; specialist centre | C (1) |
| **HSP90 inhibitors** (investigational) | HSP90 folds kinase-domain BUBR1 mutants; inhibiting it severely lowers mutant protein in cells | avoid unless the missense allele is shown not to be low-abundance | M |
| **Nephrotoxic drugs** (aminoglycosides, etc.) | the child has nephrocalcinosis | prefer alternatives; levels and renal monitoring | clinical |
| **Growth hormone** | poor or absent growth response in four MVA patients (three *CEP57*, one unspecified; **none *BUB1B***); consensus: generally contraindicated in cancer-predisposition syndromes | only with oncology/endocrinology agreement and a pre-defined trial-of-response endpoint | G, C (4) |
| **NAD⁺ precursors**, including over-the-counter supplements | rhabdomyosarcoma depends on NAMPT/NAD⁺ salvage; a NAD⁺ supplement would antagonise a NAMPT inhibitor should one ever be used | ask about supplements; avoid while disease is active | M |

**B · If a candidate from the Track 2 report is ever considered** (label items verified in current FDA labels)

| Candidate | Items to plan for | Suggested action | Evidence |
|---|---|---|---|
| **Dasatinib** | pleural effusion; pulmonary arterial hypertension; paediatric effects on bone growth and epiphyses — relevant to a child with short stature | baseline and periodic echocardiography; growth and bone monitoring | L |
| **Sirolimus, everolimus** | rapamycin caused chromosome malsegregation in human lymphocytes in one in-vitro study; proteinuria (child has nephrocalcinosis); hyperlipidaemia; impaired wound healing; trough-level monitoring; CYP3A4 interactions; avoid live vaccines | %MN-RET as a safety readout; renal, lipid and trough monitoring; plan around surgery | L, M |
| **JAK inhibitors** | serious infections including herpes zoster and tuberculosis reactivation; thrombosis; avoid live vaccines (tofacitinib, baricitinib labels); fedratinib: boxed warning for encephalopathy, including Wernicke's (thiamine) | TB screening and VZV immunity before starting; vaccinate beforehand; thiamine levels before and during fedratinib | L |
| **NAD⁺ precursors** (niacin, nicotinamide) | see section A | only in documented remission, with intensified surveillance | M |

**Surveillance:** follow the AACR 2024 recommendations for MVA (Nakano et al., *Clin Cancer Res* 2024, PMID 39264246).

**Key references (PMIDs):** 31081598 and 31184400 (MVA cases), 31053147 (HSCT), 39264246, 35319491 (growth hormone
consensus), 35434947, 41700350, 36588761, 12116237 (growth hormone in MVA), 20516114 (HSP90), 37616468 (NAD⁺/RMS),
9914383 (rapamycin), 40607964 (heart), 27922816 (CNS myelin). Connectivity analysis: `scripts/t2/cmap.py`.
