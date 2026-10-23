# Translational readthrough of PTCs and NMD inhibition as a strategy in genetic disease

Scope note: compiled 2 Oct 2026 from PubMed E-utilities abstracts, Europe PMC, EMA and company/SEC releases. Where only an abstract was read, numbers are abstract-level. Items marked (unverified) came from a search snippet or memory and were not confirmed on a primary page.

## 1. How efficient is readthrough by stop codon identity and +4 nucleotide context?

### Takeaway
UGA is the most readthrough-permissive stop, UAA the least, UAG intermediate. The +4 base matters, and cytosine is the strongest promoter (UGAC best, UAAG worst). Context rules are probabilistic, and basal readthrough is typically well under a few percent, so context only sets the ceiling for any drug.

### Cited Findings
- In human cells, UGA is the most readthrough-permissive stop and UAA the least. With the stop held constant, +4 order is C > U > A/G. UGAC is the most permissive tetranucleotide and UAAG the most inhibitory. A model trained on ribosome profiling predicted readthrough of CFTR nonsense variants under G418 (Spearman 0.88 after excluding outliers). Longer 3'UTRs raise readthrough (PMID 38509072, Mangkalaphiban et al., Nat Commun 2024; DOI 10.1038/s41467-024-46703-z) — [PMC10954755](https://pmc.ncbi.nlm.nih.gov/articles/PMC10954755)
- Ribosome profiling after aminoglycoside (G418) showed stop identity, flanking nucleotides and wider sequence context all shape drug-stimulated readthrough. This is a genome-wide effect on native stops, not only on disease PTCs (PMID 31971508, Wangen & Green, eLife 2020; DOI 10.7554/eLife.52611) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/31971508/)
- General ranking in the literature: UAA most efficient termination, then UAG, then UGA. Cytosine at +4 is the strongest inducer for all three stops. One summarised ordering for UGA is C > A > G > U, which differs from Mangkalaphiban's C > U > A/G. eRF1 is reported to sense +4 C poorly (search-engine summary of older work, primary papers not opened) — [PMC1084031](https://pmc.ncbi.nlm.nih.gov/articles/PMC1084031), [PMC4748822](https://pmc.ncbi.nlm.nih.gov/articles/PMC4748822)
- Position-specific effects extend beyond +4. In mammalian cells +8 is also influential, and +5/+6 act in a more complex way (search summary of [PMC10954755 refs](https://pmc.ncbi.nlm.nih.gov/articles/PMC10954755)).
- The inserted amino acid depends on the stop (near-cognate tRNAs: Trp/Cys/Arg for UGA, Gln/Tyr/Lys for UAA/UAG). Readthrough therefore yields a missense-like full-length product whose function depends on the residue. "Readthrough effectiveness depends on genetic context and amino acid incorporation type" (PMID 40474765, Vitale, IUBMB Life 2025) — abstract via PubMed efetch.
- Response is mutation-specific even within one gene. COL17A1 nonsense variants responded variably to readthrough drugs, with a personalised approach urged (PMID 42824587, Tan, Mol Ther Nucleic Acids 2026). Ataluren restored function in only about 23% of NF1 nonsense fibroblast lines (PMID 41256734, Kim, MedComm 2025; DOI 10.1002/mco2.70485).
- NMD itself is position-dependent and not binary. Saturation genome editing of 722 PTCs showed complex position-dependent NMD (bioRxiv preprint, PMID 41509423, Cortázar 2025). EJC-enhanced NMD lowered RNA 2-9 fold and protein 7-12 fold versus EJC-independent NMD (PMID 38375914, Kolakada, NAR 2024; DOI 10.1093/nar/gkae121).

### Inferences
- A UGA stop followed by C is the best case for any readthrough drug. A UAA or UAG followed by A/G is the worst. The case's BUB1B stop codon identity and +4 base should be read off the reference sequence (gene-level, no genotype table needed) before proposing readthrough. This was not done here.
- Mixed amino-acid insertion means that even successful readthrough may give a hypomorphic BUBR1. This is consistent with MVA1 already being a null/hypomorph disease, but no direct evidence for BUB1B readthrough was found.

### Gaps
- No primary +4 percentage tables were extracted. Absolute readthrough figures per context (e.g. Manuvakhova 2000, Dabrowski 2015) were not fetched.
- No study of readthrough of BUB1B PTCs in any system was found.

## 2. Compound classes: aminoglycosides, ataluren, amlexanox, SMG1 inhibitors, DAP, ELX-02

### Takeaway
Aminoglycosides (G418, gentamicin) have genuine readthrough activity but are toxic and weak in vivo. Ataluren's mechanism is disputed and two phase 3 programs failed (EMA non-renewal). Amlexanox acts mainly via NMD, DAP is the most potent in vitro for UGA, and ELX-02 is a toxicity-reduced aminoglycoside with modest clinical effect. SMG1 inhibition synergises with aminoglycosides but not with non-aminoglycoside agents, and is toxic.

### Cited Findings
**Aminoglycosides**
- Review stresses toxicity, bioavailability and cellular interactions as key obstacles (PMID 30134808, Dabrowski, Mol Med 2018) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/30134808/)
- Designer aminoglycosides aim to keep readthrough while cutting toxicity (derivative with enhanced eukaryotic selectivity; PMID 37312846, Guchhait, ACS Med Chem Lett 2023).
- TRPC channel inhibition (AC1903) reduced G418 uptake and readthrough, showing uptake matters (PMID 34999117, Baradaran-Heravi, J Biol Chem 2022).
- Minipig NF1 model: gentamicin alone gave unreliable neurofibromin rescue, and gentamicin plus NMD inhibition restored it in primary cells. Gentamicin reached the optic nerve (PMID 37520682, Osum, Mol Ther Nucleic Acids 2023; DOI 10.1016/j.omtn.2023.06.018).

**Ataluren (PTC124)**
- Original mechanism was challenged: ataluren stabilises firefly luciferase, a false positive in the FLuc reporter. In non-luciferase assays it showed "no evidence of activity" and did not match G418 (PMID 23824517, McElroy et al., PLoS Biol 2013; DOI 10.1371/journal.pbio.1001593). Context: [Nature blog](https://blogs.nature.com/spoonful/2013/06/researchers-question-read-through-mechanism-of-muscular-dystrophy-drug-ataluren.html). The follow-up literature is partial, and some later papers still report ataluren effects (below).
- Ataluren has a bell-shaped dose-response with a narrow effective window (MPS I-H mouse, PMID 35857082, Wang, J Mol Med 2022).
- Restored SBDS in Shwachman-Diamond cells in vitro (PMID 37876306, Cipolli, Br J Haematol 2024), partial BMPR1A restoration with G418 and PTC124 (PMID 38192234, Wang, Ann Hum Genet 2024). Both are cell-level.
- Factor V: PTC-124 and amlexanox were ineffective (PMID 37866515, Todaro, J Thromb Haemost 2024; DOI 10.1016/j.jtha.2023.10.007).

**Amlexanox**
- Increased nonsense mRNA and full-length protein; first described in Gonzalez-Hilarion et al., Orphanet J Rare Dis 2012 (PMID 22938201; abstract not opened, from search snippet). It raised COL7A1 transcript and UPF1 phosphorylation in RDEB cells, suggesting NMD inhibition.
- At 100 uM for 20 h it stabilised GDAP1 UGA-PTC mRNA and restored mitochondrial morphology in CMT2H patient neurons (PMID 37513945, Benslimane, Pharmaceuticals 2023).
- Did not increase G542X-Cftr quantity and gave no synergy with SMG1i (McHugh et al., below).

**SMG1 inhibitors / NMD inhibition**
- SMG1i synergised with G418, gentamicin and paromomycin, but not amikacin, tobramycin, PTC124, escin or amlexanox. SMG1i plus 100 uM G418 raised G542X-Cftr mRNA 7.6-fold (72.9% of WT mRNA); with VX-661 an estimated 8.8% of WT CFTR function. At 200 uM G418 plus 1 uM SMG1i, toxicity reduced the response. SMG-1 acts outside NMD too (McHugh, Cotton, Hodges, Int J Mol Sci 2021; DOI 10.3390/ijms22010344; PMCID PMC7794695, I did not confirm the PMID) — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7794695/)
- CFTR PTC mRNAs are decayed through SMG6, and NMD inhibition improves readthrough therapy (PMID 35487895, Sanderlin, Nat Commun 2022; DOI 10.1038/s41467-022-29935-9).
- DMD: readthrough drugs plus NMD inhibitors act synergistically on dystrophin mRNA (PMID 37145099, Amar-Schwartz, Hum Mol Genet 2023; DOI 10.1093/hmg/ddad072).
- G418 plus SMG1i was additive for EJC-independent and synergistic for EJC-enhanced NMD reporters (PMID 38375914).
- Reviews: PMID 38543100 (Benslimane, Pharmaceuticals 2024), PMID 36979640 (Temaj, Biomedicines 2023), PMID 40420818 (Ricci, IUBMB Life 2025), methods PMID 40699434 (Nedjma 2025).
- Quintuple combination in CF organoids (readthrough + NMD inhibitor + correctors + potentiator) achieved clinically relevant restoration (PMID 34666947, de Poel, J Cyst Fibros 2022).

**2,6-Diaminopurine (DAP)**
- Corrects UGA (and UAA) nonsense by inhibiting FTSJ1 (tRNA-Trp C34 2'-O-methylation). Raised p53 in UGA-TP53 cancer cells and reduced tumour growth in mice (PMID 32198346, Trzaska, Nat Commun 2020) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/32198346/). Note the effect is UGA-specific, i.e. depends on Trp-tRNA near-cognate decoding.
- CF mouse in utero/breastfeeding: "clinically relevant" correction (PMID 36641622, Leroy, Mol Ther 2023; DOI 10.1016/j.ymthe.2023.01.014).
- Factor V: DAP raised activity 10.8-fold, versus G418 7.0-fold and ELX-02 3.1-fold, from 0.7% baseline; G418 and DAP showed cytotoxicity (PMID 37866515).
- Restored ACTN3 577X alpha-actinin-3 at micromolar levels (PMID 41931225, Harada, Biochem Genet 2026).

**ELX-02**
- G542X CFTR organoids: dose-dependent function, 5-fold CFTR mRNA in homozygotes (PMID 33558100, Crawford, J Cyst Fibros 2021) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/33558100/)
- Phase 2 monotherapy (G542X): sweat chloride fell by a statistically significant 5.4 mmol/L at 1.5 mg/kg/day (company release, [SEC 8-K 2021](https://www.sec.gov/Archives/edgar/data/1035354/000110465921140321/tm2133239d1_8k.htm)).
- Phase 2 ELX-02 plus ivacaftor missed statistical significance on sweat chloride and ppFEV1 (Sept 2022 company release, [Benzinga](https://www.benzinga.com/general/biotech/22/09/28873111/eloxx-shares-plummet-as-mid-stage-cystic-fibrosis-study-results-miss-efficacy-endpoint), [Global Genes](https://globalgenes.org/?p=621814)). Sources are press/SEC rather than peer review.
- Triple combination ELX-02 + VX-809 + CC-90009 gave 3-fold CFTR rescue in cells, increased to 15-fold with IL-4 (PMID 41381228, Borrelli, Eur Respir J 2026).

**Other mechanisms**
- SRI-41315 reduces eRF1 and synergises with aminoglycosides; identified from 771,345 compounds (PMID 34272367, Sharma, Nat Commun 2021).
- Heldin 2023: C47/C61 plus G418 or eRF3-targeting restored p53 and PTEN (PMID 37339263, RNA Biol).
- Suppressor tRNAs (non-drug alternative): Lueck 2019 (PMID 30778053), Wang 2022 AAV (PMID 35322228), Albers 2023 (PMID 37258671).
- Targeted pseudouridylation suppressed NMD and promoted readthrough (PMID 36764303, Adachi, Mol Cell 2023).

### Inferences
- A rational readthrough regimen for a UGA PTC would be DAP and/or an aminoglycoside (G418-like) with NMD inhibition. Ataluren and amlexanox have the weakest and least reproducible evidence. No agent here is approved for this purpose, and DAP and ELX-02 are investigational, not repurposable approved drugs.
- Of approved drugs, only gentamicin/paromomycin-type aminoglycosides plus amlexanox (approved for aphthous ulcers in the US) and ataluren (not FDA approved; EU authorisation lapsed) are candidates, with the evidence caveats above.

### Gaps
- No peer-reviewed publication of ELX-02 phase 2 data was found; only company releases.
- DAP human safety data were not found. DAP is not an approved drug in the sources reviewed.

## 3. Reported protein restoration levels

### Takeaway
Reported rescue is usually a few percent to a few tens of percent of wild type in vitro, rarely quantified in vivo, and often reported as fold-change over a very low PTC baseline, which exaggerates the therapeutic impact. Clinical benefit is small and mostly not statistically convincing.

### Cited Findings
- G418 + SMG1i: 72.9% of WT Cftr mRNA (mRNA, not protein), about 8.8% of WT CFTR function with VX-661 in murine G542X organoids (McHugh 2021, DOI 10.3390/ijms22010344).
- Factor V: 3.1-10.8-fold over a 0.7% baseline, so absolute activity stays in the single-digit percent range (PMID 37866515; the percent range is my arithmetic).
- ELX-02: 5-fold CFTR mRNA in G542X homozygous organoids; function comparable to F508del-modulator treatments (PMID 33558100).
- Clinical ELX-02: sweat chloride -5.4 mmol/L (SEC 8-K 2021). This is small relative to CFTR modulators.
- Ataluren NF1: functional restoration (RAS-GTP, pERK) in about 23-24% of patient fibroblast lines (PMID 41256734).
- DAP in CF mice: "clinically relevant" without a percent in the abstract (PMID 36641622).
- Ataluren DMD study 041: 6MWD decline reduced 21% (p=0.0248) in ITT (n=359), 30% in the 300-400 m subgroup (PMID 40843507, Vlodavets, J Comp Eff Res 2025). The earlier ACT DMD trial had a non-significant 13.0 m difference (p=0.213) in ITT, significant only in the prespecified subgroup (42.9 m, p=0.007) (PMID 28728956, McDonald, Lancet 2017).

### Inferences
- Because BUBR1 is dosage-sensitive (hypomorphs), even 5-15% restoration could matter biologically, but no evidence of the required threshold was found in the sources here.

### Gaps
- Most abstracts give fold-changes, not percent of WT protein. Absolute percentages would need full texts.

## 4. Why many attempts failed, including negative and contrary results

### Takeaway
Failures stem from weak in vivo efficacy, toxicity (nephro/ototoxicity, SMG1 pleiotropy), artefactual assays (luciferase), heterogeneous responses by PTC context and patient, NMD degrading the target mRNA, and low-dose or narrow-window pharmacology. Regulators did not accept ataluren.

### Cited Findings
- Luciferase artefact and non-reproducibility of ataluren readthrough (PMID 23824517; see above).
- Ataluren CF phase 3: no overall benefit; post hoc benefit only in patients not on inhaled tobramycin; creatinine rises in 15% versus <1% of placebo (PMID 24836205, Kerem, Lancet Respir Med 2014; DOI 10.1016/S2213-2600(14)70100-6). The Cochrane review (2 trials, n=517) found no difference in quality of life or lung function, higher renal events (RR 12.81), and the tobramycin subgroup was not reproduced in a later prospective trial; evidence insufficient (PMID 36866921, Aslam, 2023).
- Ataluren DMD: phase 2b missed the primary endpoint; the Genzyme deal ended; FDA declined filing in 2016 (summary: [Wikipedia PTC124 via search](https://www.wikipedia.com/wiki/PTC124), secondary source). EMA CHMP confirmed non-renewal of the conditional authorisation: confirmatory studies 020 and 041 "failed to confirm" phase 2 study 007 — [EMA](https://www.ema.europa.eu/en/news/ema-confirms-recommendation-non-renewal-authorisation-duchenne-muscular-dystrophy-medicine-translarna). Contrary view: Study 041 authors report significant 21% slowing of decline and describe benefit (PMID 40843507; sponsor-affiliated authors).
- ELX-02: monotherapy modest, combination with ivacaftor non-significant (company releases above).
- Aminoglycosides: nephro/ototoxicity, poor tissue penetration, need for high doses. Cellular uptake limits effect (PMID 34999117). Designer analogues remain preclinical (PMID 37312846).
- NMD works against readthrough: PTC mRNA is degraded (SMG6 for CFTR, PMID 35487895), so readthrough alone yields little protein, but NMD inhibitors (SMG1i) are cytotoxic at high dose (McHugh 2021) and SMG1 is needed for other processes; effect depends on EJC context (PMID 38375914).
- No synergy of SMG1i with non-aminoglycoside agents (ataluren, amlexanox, escin) (McHugh 2021). The two proposed mechanisms are therefore not interchangeable.
- Heterogeneity across variants and patients (COL17A1 PMID 42824587; NF1 PMID 41256734). Some host genes modulate readthrough (e.g. TRPC uptake PMID 34999117).
- Readthrough may have immunogenic or off-target consequences. Hemophilia A: natural readthrough linked to inhibitor development (PMID 42374266, Testa, Mol Med 2026). Native stop codons are also read through (Wangen 2020).
- Readthrough products may be mutated proteins with a wrong residue at the PTC position.
- Alternative approaches have emerged because of these limits: suppressor tRNAs and AAV delivery (PMIDs 30778053, 35322228, 37258671), ASO exon skipping for CFTR W1282X (PMID 35017301).

### Inferences
- For a transient or vascular-limited indication a safety review would be needed, and readthrough in a child with a rhabdomyosarcoma-prone, aneuploidy-prone syndrome has unclear safety (genome-instability, not directly addressed by any source found).
- The strongest evidence for benefit is aminoglycoside plus NMD inhibitor synergy in vitro/ex vivo; none of it is for BUB1B.

### Gaps
- No evidence on NMD of BUB1B nonsense transcripts in MVA patient cells beyond a search snippet reporting NMD-reduced BUB1B transcripts in Japanese PCS families (monoallelic BUB1B mutations; primary paper not opened; [PMC3508359](https://pubannotation.org/docs/sourcedb/PMC/sourceid/3508359)). Exon 17/23 location and downstream EJC rule: a stop in exon 17 of 23 is expected to trigger NMD (general rule, not source-confirmed here).
- Whether complete loss of BUBR1 is lethal is relevant: MVA1 requires a hypomorphic second allele (from project memory, not re-sourced). Restoring PTC-allele protein could therefore be harmful or helpful, and no source was found.
- Wangen/Mangkalaphiban were not tested for BUB1B stops specifically.
