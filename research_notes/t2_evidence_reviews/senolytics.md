# Senolytics and senomorphics in progeroid / chromosomal-instability models and their translation to humans

Method note: ~19 tool calls (WebSearch + ClinicalTrials.gov API v2). Europe PMC REST returned HTTP 503 and PubMed pages did not render in WebFetch, so most claims rest on search-result summaries, not full-text reading. DOIs marked "(DOI from memory, not verified in this session)" should be checked before citing. No patient data used.

## 1. Which compounds were tested in which progeroid / CIN models, and with what outcome?

### Takeaway
The only direct test in a chromosomal-instability (CIN) model is genetic (INK-ATTAC) clearance of p16Ink4a+ cells in BubR1-hypomorphic mice (Baker 2011), which delayed tissue-specific ageing phenotypes but did not extend lifespan (to be confirmed, see Gaps). Pharmacological senolytics (fisetin, CK2 inhibition, ABT-263) have been tested in Ercc1 and Zmpste24 progeroid mice and in aged wild-type mice; I found NO study of D+Q, fisetin or navitoclax in BubR1-hypomorphic or other aneuploidy mouse models.

### Cited Findings
- BubR1 H/H hypomorphic mice (short lifespan, infertility, lordokyphosis, sarcopenia, cataracts, fat loss, arterial stiffening, impaired wound healing). The INK-ATTAC transgene (p16Ink4a promoter driving FKBP-caspase-8, activated by AP20187) was crossed in; life-long removal of p16Ink4a+ cells delayed onset of phenotypes in adipose tissue, skeletal muscle and eye; late-life clearance attenuated progression of established disorders. Baker et al., Nature 2011, DOI 10.1038/nature10600 — [Springer](https://link.springer.com/article/10.1038/nature10600); [PMC3468323](https://pmc.ncbi.nlm.nih.gov/articles/PMC3468323)
- Selectivity limits in that paper (from my background knowledge, not re-verified here): benefit was tissue-specific (tissues where p16Ink4a drives pathology); clearance did not rescue all tissues (e.g., effect was absent in organs where p16 was not the driver). Source above should be checked for exact wording.
- Fisetin in Ercc1-/Δ progeroid mice (carrying p16-luciferase reporter) and aged WT mice: most potent of 10 flavonoids screened; acute or intermittent ("hit-and-run") treatment reduced senescence markers in multiple tissues and age-related pathology; late-life fisetin extended median and maximum lifespan in WT mice. Yousefzadeh et al., EBioMedicine 2018 — [PMC6197652](https://pmc.ncbi.nlm.nih.gov/articles/PMC6197652) (DOI 10.1016/j.ebiom.2018.09.015, DOI from memory, not verified)
- CONTRARY: NIA Interventions Testing Program found fisetin did NOT significantly extend lifespan in either sex at the doses/schedules used — [Fight Aging summary of the ITP report (secondary source)](https://www.fightaging.org/archives/2023/12/the-nia-interventions-testing-program-shows-that-fisetin-does-not-extend-life-in-mice/). The primary ITP paper was not retrieved.
- D+Q in naturally aged mice: intermittent oral dosing alleviated physical dysfunction and increased post-treatment survival by 36%; same paper used Ercc1-/Δ mice (from memory). Xu et al., Nat Med 2018, PMID 29988130 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/29988130/) (DOI 10.1038/s41591-018-0092-9, from memory)
- Zmpste24-/- (HGPS-like) mice: fisetin in vivo reduced senescent fibro-adipogenic progenitors, restored muscle stem cell number and improved muscle pathology; CK2 inhibition eliminates senescent cells and prolongs lifespan in the same model — [PMC11139886](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11139886/) (the fisetin muscle result is cited within this paper's search summary; original fisetin-Zmpste24 paper not identified here)
- Navitoclax (ABT-263; BCL-2/BCL-xL inhibitor): oral dosing depleted senescent cells in irradiated and aged mice and rejuvenated aged hematopoietic stem cells and muscle stem cells — Chang et al., Nat Med 2016, PMID 26657143 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/26657143/) (DOI 10.1038/nm.4010, [Nature](https://www.nature.com/articles/nm.4010)). Search found no navitoclax study in a progeroid or CIN mouse model.
- Senomorphic/mTOR arm in HGPS: genetic mTOR reduction extends lifespan in LmnaG609G mice — [PMC8441492](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8441492/); rapamycin increased lifespan >50% in Lmna-/- mice (as summarized in the same search results; those are cardiomyopathy-type laminopathy mice, not true HGPS). Ruxolitinib/metformin in progeroid models: search returned nothing relevant (see Gaps).
- Aneuploidy-adjacent: aneuploid cells with complex karyotypes arrest, show senescence features, secrete pro-inflammatory cytokines and are eliminated by NK cells — Santaguida et al., Dev Cell 2017;41:638-651, [PMC5536848](https://pmc.ncbi.nlm.nih.gov/articles/PMC5536848) (DOI 10.1016/j.devcel.2017.05.022, from memory)

### Inferences
- BubR1 hypomorphy is the best mechanistic precedent for a BUB1B hypomorph, but the evidence is genetic ablation in a mouse that is progeroid via p16-driven senescence; BUB1B-mutant MVA1 is not obviously the same (human MVA patients show cancer predisposition, not clear progeroid ageing); translation of senolytic benefit is hypothetical.
- A BubR1 senolytic gap exists: no one (that I found) has tested D+Q/fisetin/navitoclax in BubR1 H/H mice.

### Gaps
- Baker 2011 lifespan outcome and exact tissue list not verified from full text.
- Original Zmpste24/fisetin paper and ITP primary paper not retrieved.
- Other senolytic studies in BubR1 mice (e.g., later work by Baker/van Deursen on BubR1 + NAD+ or senolytics) not searched properly; Europe PMC down.
- Ruxolitinib/JAK, metformin, rapalog senomorphic data in progeroid/CIN models not retrieved.

## 2. Human trial results to date

### Takeaway
Human senolytic trials are small, mostly open-label and short; the clearest biological signal is the 9-patient DKD biopsy study (p16+ cells down 35% adipose, 20% epidermis). The one RCT in bone (D+Q, n=60) missed its primary endpoint. No trial in children or in a genetic chromosomal-instability disorder was found.

### Cited Findings
- Diabetic kidney disease, open-label phase 1 (n=9, D 100 mg + Q 1000 mg x3 days): p16INK4a+ cells reduced 35% in adipose and 20% in epidermis, reduced macrophage accumulation and circulating SASP factors; no placebo. Hickson et al., EBioMedicine 2019, PMID 31542391, NCT02848131 (NCT from memory, not verified) — [PMC6796530](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6796530/)
- IPF, open-label pilot (n=14; D 100 mg + Q 1250 mg, 3 days/week x 3 weeks): 100% retention, improved 6-min walk, gait speed, chair stands; pulmonary function, FI-LAB and clinical chemistry UNCHANGED. Justice et al., EBioMedicine 2019, PMID 30616998 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/30616998/), [PMC6413337](https://pmc.ncbi.nlm.nih.gov/articles/PMC6413337); NCT02874989 (listed as completed phase 1, n=26, in ClinicalTrials.gov API query)
- Postmenopausal women, phase 2 RCT (n=60, D+Q intermittent): primary endpoint (20-week % change in CTx) NOT met (median D+Q -4.1% vs control -7.7%, P=0.611); post hoc subgroup with high baseline senescent burden (e.g., high-CTx/ older) showed increased bone formation, decreased resorption, wrist BMD gain. Farr et al., Nat Med 2024, PMID 38956196 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/38956196/); NCT04313634 (completed, has results, n=74 listed) — [ClinicalTrials.gov](https://classic.clinicaltrials.gov/show/NCT04313634)
- Alzheimer's SToMP-AD, open-label phase 1 pilot, 12 weeks intermittent D+Q; primary aim CNS penetrance (CSF levels). Gonzales et al., Nat Med 2023, [PMC10168460](https://pmc.ncbi.nlm.nih.gov/articles/PMC10168460); NCT04063124 — [CT.gov](https://classic.clinicaltrials.gov/ct2/show/NCT04063124). (Detail from memory, needs verification: dasatinib detected in CSF at low level, quercetin not; small n=5.)
- Systemic sclerosis, prolonged dasatinib reduced SASP/senescence markers in skin biopsies (as reported in search summary; trial not identified).
- Fisetin trials (ClinicalTrials.gov API v2 query, status as of Oct 2026): STOP-Sepsis NCT05758246 (phase 2, recruiting, n=220); vascular function in older adults NCT06133634 (active, primary completion ~Mar 2026); frailty in older women NCT03430037; COVID-19 fisetin NCT04771611 (completed, results posted), NCT04537299 (terminated, results posted), NCT04476953; carpal tunnel NCT05416515 (completed, results); PAD NCT06399809; fisetin + exercise in breast cancer survivors NCT06113016; mild AD NCT07279714 (completed, n=5). Published results for most not retrieved. A 2025 Aging Cell paper reports intermittent fisetin improved physical function and reduced senescence in skeletal muscle (preclinical mouse study, PMID 40437670 if the PubMed hit is that paper — not verified).
- Other D+Q trials from CT.gov: knee OA NCT04210986 (completed, results), liver fibrosis NCT05506488 (phase 1/2, completed), HSCT survivors NCT02652048? [NCT02652052 listed as completed phase 1 n=9; adult HSCT survivors; not pediatric as far as retrievable], D+Q with CAR-T in myeloma NCT06940297, ILD in CVID NCT05593588, ALSENLITE NCT04785300, Alzheimer's-risk adults NCT05422885 (completed, results), osteoporosis NCT06018467 (SENIOR). ROPE osteoarthritis NCT04770064 was WITHDRAWN (n=0); GlaucoCetin NCT04784234 TERMINATED; VIAging NCT04994561 withdrawn. — [CT.gov API query](https://clinicaltrials.gov/api/v2/studies?query.term=senolytic%20OR%20dasatinib%20quercetin%20OR%20fisetin)
- Navitoclax: dose-limiting on-target thrombocytopenia (BCL-xL in platelets) and neutropenia; D: anaemia, thrombocytopenia, neutropenia, pleural effusion, GI effects — [PMC10952275 table](https://pmc.ncbi.nlm.nih.gov/articles/PMC10952275/table/febs16351-tbl-0001); platelet-sparing PROTAC and galacto-prodrug versions exist in preclinical work — [Fight Aging](https://www.fightaging.org/archives/2020/05/turning-navitoclax-into-a-protac-senolytic-with-a-better-safety-profile/).

### Inferences
- Most positive human signals are open-label and surrogate-based; the only RCT here was negative on its primary endpoint.
- No senolytic has been tested in pediatric or syndromic-CIN populations (none found).

### Gaps
- Results for SToMP-AD detail, STOP-Sepsis, COVID-FIS results not read.
- NCT02848131 and NCT02652048 identifiers are from memory/unconfirmed; pediatric trial listings not specifically searched.
- Dasatinib pediatric safety (it is approved for pediatric Ph+ ALL/CML) not searched; a reasonable inference is its pediatric safety profile is characterized in oncology, but this needs a source.

## 3. Risks of removing senescent cells in a growing organism

### Takeaway
Senescence has essential roles in mammalian embryogenesis (limb AER, mesonephros, inner ear), in wound healing and in liver endothelium; senescent-cell clearance can worsen injury in some cell types. No study of senolytic effects on childhood growth was found.

### Cited Findings
- Developmental senescence: p21-dependent programmed senescence in mesonephros and endolymphatic sac (Muñoz-Espín et al., Cell 2013) and in the apical ectodermal ridge and neural roof plate (Storer et al., Cell 2013); p21-deficient mice show patterning defects; senescent cells are cleared by macrophages — [Nat Rev Mol Cell Biol commentary](https://www.nature.com/articles/nrm3731), [UPF repository entry](https://repositori.upf.edu/handle/10230/23507). (Note: those are embryonic stages; relevance to postnatal growth is an extrapolation.)
- Liver: clearing p16Ink4a+ macrophages mitigates hepatocellular damage, whereas eliminating p16Ink4a+ sinusoidal endothelial cells (~70% of p16+ cells in 12-month livers) aggravates liver injury and fibrosis; Grosse et al., Cell Metab 2020 — [Aging-US/related review PMC12867169](https://pmc.ncbi.nlm.nih.gov/articles/PMC12867169) and [PubMed 39368477](https://pubmed.ncbi.nlm.nih.gov/39368477/) (the latter being a review on functional roles of senescence across cell types)
- Transient senescence after injury promotes repair (SASP factors); only persistent senescent cells are harmful — see search summary in the same reviews. Demaria 2014 (PDGF-AA from senescent fibroblasts aids wound closure; ablation delayed it) is from my knowledge and was not retrieved: DOI 10.1016/j.devcel.2014.11.012 (from memory, unverified).
- Senolytic-specific toxicities: thrombocytopenia (navitoclax), cytopenias / pleural effusion (dasatinib) as above.
- INK-ATTAC BubR1 data: BubR1 H/H mice themselves have impaired wound healing and dermal thinning — [Baker 2011 PMC3468323](https://pmc.ncbi.nlm.nih.gov/articles/PMC3468323).

### Inferences
- In a child with a BUB1B hypomorph, the main extra risks vs adults are (i) interference with growth/tissue modelling, (ii) cytopenias with marrow already stressed by aneuploidy, and (iii) loss of tumour-suppressive senescence of aneuploid cells that otherwise arrest (Santaguida above); the third is a hypothesis inferred from findings, not demonstrated.

### Gaps
- No direct data on senolytic dosing during childhood growth, puberty or bone development.
- No source retrieved on senolytics and tumour risk in the context of CIN.

## 4. Has senescence burden been measured in CIN / aneuploidy disorders?

### Takeaway
Cell-culture and constitutive-trisomy data show senescence/inflammatory signatures (cGAS-STING) in aneuploid cells, and MVA is described as having chronic immune response and cellular stress, but I found no study quantifying senescent-cell burden (e.g., p16, SA-beta-gal in tissues) in MVA patients.

### Cited Findings
- Constitutive trisomies (primary trisomic fibroblasts, Down syndrome): genotoxic stress activates cGAS-STING, autophagy and innate immune response independent of which chromosome is extra — [PMC8253785](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8253785/)
- MVA review: MVA accompanied by tumour susceptibility and cellular and systemic stress including chronic immune response. "Mosaic variegated aneuploidy in development, ageing and cancer", PMID 39169218 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/39169218/) (abstract not retrieved; content from search snippet only)
- Experimentally induced missegregation produces arrested, senescence-featured, SASP-secreting cells cleared by NK cells (Santaguida 2017, above).
- Senescent fibroblasts activate inflammatory programs through cGAS/STING; non-canonical pathway proposed to drive ageing — bioRxiv preprint 2025, [biorxiv 10.1101/2025.04.03.645994](https://www.biorxiv.org/content/10.1101/2025.04.03.645994.full.pdf) (preprint, not peer reviewed)

### Inferences
- If MVA cells are cleared by immune surveillance via senescence-like arrest, senescent burden may actually be low in patients (cells removed) — speculative; no data.

### Gaps
- No measurement of senescence markers in MVA (BUB1B, CEP57, TRIP13) patient tissues or patient-derived fibroblasts was found; Europe PMC 503 prevented an exhaustive search. A rerun with targeted terms (e.g., "BUB1B mutant fibroblasts senescence SA-β-gal") is advisable.
