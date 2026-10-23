"""formulario_metodos.py — rellena la hoja «Track 1 methods» del formulario oficial (docs/methods_description_form.xlsx)
y guarda la copia en results/envios/t1_final/aperezcastro_track1_methods_form.xlsx. Solo texto de métodos.
"""
import pathlib

import openpyxl

ROOT = pathlib.Path(__file__).resolve().parents[1]
LLM = ("Anthropic Claude via the Claude Code CLI (Claude Opus models), Claude Max (consumer tier), model-training "
       "opt-out enabled ('Help improve Claude' off), no feedback provided on any output, no content flagged for safety "
       "review - per the organisers' confirmed guidance of 28 Sep 2026 (discussion #2). Applicable terms: opted-out "
       "conversations are not used for training and are retained 30 days. Additional data minimisation: the model "
       "only saw code, aggregate counts, HPO terms and short candidate lists (<= 50 annotated variants); raw data and "
       "genome-scale derivatives were blocked from the agent by configuration. No patient variant was sent to any third-party API.")
T1 = {
    "Team name": "aperezcastro (Agustín Pérez Castro, individual)",
    "Model number": "1 (BUB1B compound heterozygous; aperezcastro_bub1b_comphet.csv)",
    "Please describe your model/approach in detail.": (
        "Three-arm pipeline on the provided VCF, run as time-capped jobs on Modal (scripts in the GitHub repo). "
        "(1) bcftools norm (split multiallelics, left-align vs Ensembl GRCh38) -> Ensembl VEP 116 offline (MANE, HGVS, "
        "gnomAD v4.1 exomes+genomes, ClinVar with review status, SIFT, PolyPhen-2, AlphaMissense). "
        "(2) Unbiased genome-wide funnel: PASS -> max gnomAD AF < 1% -> HIGH/MODERATE/splice -> gene-level recessive "
        "model (>=2 het or 1 hom) = 339 variants in 99 genes, each reviewed. "
        "(3) Hypothesis-driven arm: MVA/spindle-assembly-checkpoint panel (BUB1B, CEP57, TRIP13, CENATAC, BUB1, BUB3, "
        "MAD1L1, MAD2L1BP, SMC5/6 complex...) at any consequence and AF < 5%; of the 42 rare PASS variants, the 32 "
        "inside annotated gene bodies were scored with SpliceAI 1.3.1 (50 and 2000 nt windows; all max delta <= 0.03). "
        "(4) Exomiser 15.1 (data 2512) with the 8 HPO terms, with and without a phenotype pre-filter: BUB1B ranks #1 "
        "genome-wide in both. "
        "(5) ACMG/AMP classification with ClinGen SVI refinements and Bayesian points (Tavtigian 2020) for each "
        "submitted variant; phase feasibility tested directly (only 3 heterozygous sites in the 10.9 kb interval, gaps "
        "6.8 and 4.1 kb: not phasable with short reads) and argued from biology. "
        "(6) Genome-level QC from aggregate statistics: allelic balance per autosome and in 5-Mb bins (excess "
        "overdispersion confined to pericentromeric/segmental-duplication regions, i.e. technical; no detectable "
        "whole-chromosome imbalance, as expected for variegated low-level mosaicism), balanced allele fractions around "
        "BUB1B, and no run of homozygosity > 2 Mb (consistent with compound heterozygosity rather than autozygosity). Answer: BUB1B c.2210T>G p.(Leu737*) "
        "[Likely pathogenic, 9 points; Pathogenic with ClinVar corroboration] + c.3006T>G p.(Asn1002Lys) [VUS, 4-5 points], MVA type 1."),
    "(required)": LLM,
    "Please state whether the submission file is the automated output": (
        "Candidate generation is automated (funnel, panel, SpliceAI, Exomiser); the final 3-row file and EPCR values "
        "were curated manually from the automated evidence."),
    "Please describe any downstream manual review in detail.": (
        "Manual review of the 99 recessive candidate genes (artefact-prone regions, ClinVar-benign, low depth, "
        "phenotype fit), of Exomiser's top hits, and of dominant rhabdomyosarcoma predisposition genes. Row 1 = the pair "
        "(EPCR 0.95). Rows 2-3 are low-probability alternative pairings with a third rare deep-intronic BUB1B variant "
        "(SpliceAI 0.00; EPCR 0.03 and 0.01). EPCR reflects residual uncertainty about phase and unobserved SV/regulatory "
        "alleles. Validated with the organisers' evaluation.py."),
    "Please state whether your approach only used publicly available data": "Only publicly available data.",
    "Please describe any public data used in detail": (
        "GRCh38 (Ensembl 116 primary assembly), Ensembl VEP 116 cache, GENCODE (SpliceAI grch38), gnomAD v4.1, "
        "ClinVar (NCBI VCF, Oct 2026), AlphaMissense hg38, Exomiser data release 2512 (REVEL, MVP, AlphaMissense, "
        "SpliceAI, HPO/MGI/ZFIN phenotypes), Human Phenotype Ontology (hp.obo), UniProt O60566 and the AlphaFold model, "
        "PubMed literature (cited in the report)."),
    "Please describe any proprietary data used in detail.": "None.",
    "Is your approach able to output proposed pairs": (
        "Yes: pairs are produced by the gene-level recessive model (all qualifying pairs within a gene) and ranked; "
        "single variants are produced for dominant models."),
    "How did you handle secondary or incidental findings": (
        "No secondary findings were submitted. Incidental findings were handled following the organisers' IRB "
        "guidance (discussion #20) and are not reported."),
    "Please provide an estimate of run time and cost": (
        "About 1 hour wall-clock of cloud compute (annotation parallelised over 25 containers ~20 min, Exomiser ~7 min, "
        "SpliceAI ~4 min) plus reference downloads; total Modal cost ~1.4 USD (scripts/modal/gastos.csv). Reproducible from the repository."),
    "Provide a method abstract": (
        "We combined an unbiased genome-wide search with a hypothesis-driven analysis of genes linked to mosaic "
        "variegated aneuploidy (MVA). The single-sample WGS VCF was normalised and annotated offline with Ensembl VEP 116 "
        "(gnomAD v4.1, ClinVar, AlphaMissense) on a private cloud environment. A recessive gene-level funnel reduced "
        "4.96 M variants to 339 variants in 99 genes, of which only BUB1B belongs to the MVA/spindle-checkpoint panel. "
        "Because the second, hypomorphic allele in MVA1 is often subtle, the panel was searched at every consequence "
        "class and scored with SpliceAI; independently, Exomiser 15.1 with the patient's HPO terms ranked BUB1B first "
        "genome-wide with and without a phenotype filter (score 0.59 dominant / 0.56 recessive vs <= 0.22 for the next gene). The proband carries "
        "BUB1B c.2210T>G p.(Leu737*) (nonsense, NMD predicted; ClinVar P/LP; ACMG 9 points, Likely pathogenic, Pathogenic with ClinVar corroboration) and "
        "c.3006T>G p.(Asn1002Lys) (ultra-rare kinase-domain missense, AlphaMissense 0.92, buried high-confidence residue "
        "in the AlphaFold model 12 Å from the known MVA1 allele L1012P; 4-5 points, VUS; its class, low abundance or impaired function, is unmeasured). This "
        "null + kinase-domain hypomorph architecture is the classical MVA1 pattern and mirrors a published mouse model "
        "(BubR1 X753 + L1002P, mouse L1002 = human L1012). A third rare deep-intronic BUB1B variant was assessed and deprioritised (SpliceAI 0.00). "
        "Phase could not be observed: the 10.9 kb interval contains only three heterozygous sites (gaps 6.8 and 4.1 kb), "
        "so we tested and documented that short-read phasing is infeasible and argued trans configuration from "
        "biology (monoallelic loss does not cause MVA), absence of any other qualifying variant and balanced allele "
        "fractions, and the genome shows no run of homozygosity > 2 Mb. Strengths: three independent arms converging on the same pair, explicit ACMG evidence with Bayesian "
        "points, sensitivity to subtle alleles, full reproducibility with cost-capped scripts, and strict data "
        "minimisation (no third-party variant APIs). Limitations: no BAM/CRAM, so structural variants were not called "
        "genome-wide; phase is inferred rather than observed; in silico predictors disagree in magnitude for the "
        "missense (REVEL 0.47); single proband without segregation data."),
}

T2 = {
    "Team name": "aperezcastro (Agustín Pérez Castro, individual)",
    "Please describe your approach in detail": (
        "Allele-aware, patient-aware pipeline (scripts/t2/). (1) Mechanism model from Track 1: null allele "
        "p.(Leu737*) (UGA, NMD) + kinase-domain missense p.(Asn1002Lys) -> partial loss of function in a "
        "dosage-sensitive checkpoint. The CLASS of the missense allele (low abundance vs stable-but-dysfunctional) "
        "is unmeasured; every disease-associated substitution tested in or near the kinase domain was 5-10x less abundant and the "
        "known stable alleles lie outside that domain, so we carry both scenarios at a 0.75/0.25 prior and price "
        "the uncertainty into mechanistic fit, not evidence. AlphaFold structure: Asn1002 buried high-confidence residue near "
        "the MVA1 allele L1012P (destabilisation). (2) Knowledge base (hipotesis.yaml): 7 therapeutic lines + a safety "
        "card, each with PMID-cited evidence graded by level and direction (supporting/contrary): H1 stabilise BUBR1 "
        "(NAD+/SIRT2), H2 read-through/NMD inhibition of the nonsense allele, H3 mTORC1 (sarcopenia), H4 senolytics, "
        "H5 aneuploid-selective stress, H6 interferon/JAK-STAT (biomarker-gated), H7 proteostasis amplification/antioxidants (evaluated and declined). (3) Connectivity (cmap.py): consensus of "
        "32 BUB1B CRISPR-knockout L1000 signatures (20 cell lines) via SigCom LINCS, two-sided rank enrichment against "
        "chemical perturbation signatures, cross-cell-line consistency, Enrichr interpretation (interferon/NF-kB/p53 up, "
        "E2F/G2M down). Used to check H6 for class consistency (JAK inhibitors reverse the signature), to support the safety card "
        "(DNA-damaging drugs mimic it) and to challenge H3 (mTOR inhibitors mimic it). (4) candidatos.py: seeded "
        "approved drugs (ChEMBL target expansion is available but added none in this run), annotated with ChEMBL and openFDA labels (paediatric use, boxed "
        "warnings) and penalised by the patient's HPO-derived risk profile; transparent score 0.30 evidence + 0.25 "
        "mechanistic fit + 0.20 paediatric feasibility + 0.25 patient-specific safety, with Dirichlet weight "
        "sensitivity; only approved medicines are ranked (discussion #5). (5) Clinical tiering by benefit-risk: tier 0 "
        "safety note plus rung-0 measurement; tier 1 intermittent dasatinib (highest-scoring approved drug), to be tested first in the "
        "BubR1 H/L1002P mouse; tier 2 an NAD+ precursor - niacin (nicotinamide ties on score but inhibits sirtuins) - "
        "gated on the allele class AND documented remission; tier 3 low-dose sirolimus for muscle with %MN-RET as a "
        "safety readout; tier 4 JAK inhibitors only if an interferon signature is measured; read-through "
        "closed on current allele-level evidence; proteostasis declined because raising HSF1 is a documented "
        "aneuploidy-TOLERATING event (Donnelly 2014). (6) nof1.py: an n-of-1 power simulation on micronucleated "
        "reticulocytes (%MN-RET), the only endpoint for this disease that moves in weeks. Validation ladder: rung 0 "
        "measure %MN-RET and its variability; patient cells (decides the allele scenario); iPSC myotubes; the BubR1 "
        "allelic mouse; then an n-of-1 protocol with pre-declared responder criteria."),
    "(required)": LLM,
    "Was candidate identification automated": (
        "Hybrid: candidate expansion, annotation, safety flags, connectivity and scoring are automated; the hypotheses, "
        "evidence grading, mechanistic-fit values and final clinical tiering were curated manually (LLM-assisted "
        "literature search with every reference verified against PubMed)."),
    "Please describe any manual literature review or expert curation in detail.": (
        "PubMed searches at gene/pathway/drug level only (no patient data). Each claim in hipotesis.yaml cites a PMID "
        "whose title/abstract was retrieved and checked; contrary evidence is recorded explicitly (e.g., SIRT2 does not "
        "extend lifespan in wild-type mice; aneuploid-cell immune clearance argues against JAK inhibition). "
        "Mechanistic-fit values (0-1) are justified in the YAML. Automated label flags were reviewed manually."),
    "Please state whether you relied solely on publicly available": "Solely public databases and literature.",
    "Please describe the public data sources used in detail": (
        "PubMed/NCBI E-utilities; ChEMBL REST API (molecules, mechanisms, targets); openFDA drug labels; "
        "ClinicalTrials.gov API v2; LINCS L1000 via SigCom LINCS (l1000_xpr knockout and l1000_cp chemical "
        "signatures); Enrichr (MSigDB Hallmark 2020); UniProt O60566; AlphaFold DB AF-O60566-F1 v6; OMIM/HPO."),
    "Please describe any proprietary data sources used in detail.": "None.",
    "How did you characterize the variant's mechanism": (
        "Loss of function by reduced dosage: the nonsense allele is predicted to undergo NMD (no protein). For the "
        "kinase-domain missense we DO NOT assert a class: Suijkerbuijk 2010 splits MVA alleles into low-abundance and "
        "functionally-impaired classes, and although all disease-associated substitutions tested in or near the kinase domain were "
        "5-10x less abundant, this residue has never been characterised, so scenario A (0.75) and scenario B (0.25) "
        "are both carried and a BUBR1 immunoblot in patient cells is the pre-declared gate. Pathway: spindle assembly checkpoint / mitotic checkpoint "
        "complex (APC/C-CDC20 inhibition) and kinetochore-microtubule attachment. Downstream: chromosome "
        "missegregation -> micronuclei -> cGAS-STING interferon/NF-kB and p53 responses (confirmed in the knockout "
        "signature), senescence/SASP, mTORC1 hyperactivity in muscle (mouse allelic series), tumour predisposition."),
    "Please provide an estimate of time or effort": (
        "About 2 working days of analysis and writing; computation runs in minutes on a laptop (public APIs, cached)."),
    "Provide a method abstract": (
        "We present an allele-aware, patient-aware repurposing pipeline for a child with biallelic BUB1B mosaic "
        "variegated aneuploidy (null p.Leu737* + kinase-domain p.Asn1002Lys). The mechanism is partial loss of function "
        "in a dosage-sensitive checkpoint; the missense allele class is unmeasured and both scenarios are carried. Three evidence streams were "
        "combined: (i) a PMID-cited knowledge base of seven therapeutic lines plus a safety card, scoring supporting "
        "and contrary evidence explicitly and weighting target validation below demonstrated benefit; (ii) "
        "connectivity analysis of 32 BUB1B-knockout L1000 signatures across 20 cell lines; (iii) automated "
        "annotation from ChEMBL and FDA labels, penalised by the patient's HPO-derived risk profile, with Dirichlet "
        "weight-sensitivity analysis. The knockout signature is interferon/NF-kB/p53-up and E2F/G2M-down; "
        "DNA-damaging chemotherapeutics reproduce it and mTOR inhibitors mimic it; JAK inhibitors reverse it with "
        "class consistency but at ranks 15-154 of 539 reversers, a targeted search rather than a discovery. Two evidence "
        "reviews changed most lines: read-through was closed, on current evidence, by the stop-codon context (UGA, unfavourable +4 A, NMD) "
        "plus evidence that ataluren's activity is a luciferase artefact; "
        "nicotinamide fell because rhabdomyosarcoma depends on NAMPT/NAD+ and nicotinamide inhibits sirtuins; senolysis "
        "became a preclinical experiment after the best senolytic RCT missed its endpoint; and sirolimus fell from "
        "first to tenth when its mouse support proved to come from monoallelic carriers. "
        "THE MAIN RESULT is a decision rule: classifying every line by what it does to a chromosomally "
        "unstable cell (correct fewer errors / clear the abnormal cells / tolerate help them cope) shows that no axis "
        "wins on evidence and that each corrective line is blocked by a different nameable fact about this child. "
        "Today the best-supported option is intermittent dasatinib (clear axis, narrowly ahead of the JAK inhibitor "
        "ruxolitinib); the mechanism-matched candidate, an NAD+ precursor (niacin), ranks eighth, and rises to first if the clinical team establishes the "
        "allele class and documented remission - so the first recommendation is a measurement, not a prescription. "
        "In a cancer-predisposition syndrome a premalignant clone is an abnormal cell too. Deliverables: a two-page "
        "pharmacogenomic safety card, usable now; %MN-RET, a blood measurement never made in MVA, diagnostic "
        "first and then an endpoint with its power calculation; a ranked and "
        "conditioned candidate list, and a four-rung validation ladder (patient cells, iPSC myotubes, the published "
        "BubR1 allelic mouse model, then an n-of-1 protocol with pre-declared responder criteria) built on verified "
        "precedents (lonafarnib, sirolimus in Castleman disease). Strengths: "
        "per-allele strategies, explicit contrary evidence that demoted our own candidates, patient-specific safety "
        "beyond drug labels, and reproducible gene-agnostic code. Limitations: efficacy evidence is preclinical; "
        "connectivity signatures come from tumour cell lines after acute knockout; mechanistic-fit values are "
        "expert judgements and the ranking is sensitive to two of them, so the counterfactual is published rather than a "
        "robustness claim; every ranked drug was hand-seeded (automated target expansion added none); the %MN-RET "
        "baseline and variability used in the power calculation have never been measured in MVA; regulatory routes "
        "are orientation, not legal advice. An independent blind adversarial "
        "review of an earlier draft drove these corrections."),
}


def rellena(ws, respuestas):
    for fila in range(7, ws.max_row + 1):
        pregunta = str(ws.cell(fila, 1).value or "")
        for clave, texto in respuestas.items():
            if pregunta.startswith(clave):
                ws.cell(fila, 2).value = texto
                ws.cell(fila, 2).alignment = openpyxl.styles.Alignment(wrap_text=True, vertical="top")
    return [ws.cell(f, 1).value[:50] for f in range(7, ws.max_row + 1) if ws.cell(f, 1).value and not ws.cell(f, 2).value]


def anexo_md(ws, titulo):
    """Anexo A: preguntas del formulario oficial con nuestras respuestas, para adjuntarlo al informe (no hay campo
    para subir el .xlsx en el formulario de envío)."""
    l = [f"## Appendix A — Methods description form ({titulo})", "",
         "*Answers to the organisers' official methods description template (version with the required LLM question),"
         " reproduced here because the submission form has no separate upload field for it.*", ""]
    for fila in range(7, ws.max_row + 1):
        q, a = ws.cell(fila, 1).value, ws.cell(fila, 2).value
        if q and a:
            l += [f"**{str(q).strip()}**", "", str(a).strip(), ""]
    return "\n".join(l)


def esc(t):
    """Escapa texto plano para LaTeX."""
    rep = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
           "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "<=": r"$\leq$",
           ">=": r"$\geq$", "->": r"$\rightarrow$"}
    out = ""
    i = 0
    while i < len(t):
        for k in sorted(rep, key=len, reverse=True):
            if t.startswith(k, i):
                out += rep[k]; i += len(k); break
        else:
            out += t[i]; i += 1
    return out


def anexo_tex(ws, titulo):
    l = [r"\section{Methods description form (" + esc(titulo) + ")}",
         r"\textit{Answers to the organisers' official methods description template (version with the required LLM "
         r"question), reproduced here because the submission form has no separate upload field for it.}", "",
         r"\begin{description}[leftmargin=0pt,style=nextline,itemsep=4pt,font=\sffamily\bfseries\small\color{branddark}]"]
    for fila in range(7, ws.max_row + 1):
        q, a = ws.cell(fila, 1).value, ws.cell(fila, 2).value
        if q and a:
            l.append(r"\item[" + esc(str(q).strip()) + "] " + esc(str(a).strip()))
    l.append(r"\end{description}")
    return "\n".join(l) + "\n"


for hoja, resp, carpeta, nombre in (("Track 1 methods", T1, "t1_final", "aperezcastro_track1_methods_form.xlsx"),
                                    ("Track 2 methods", T2, "t2_final", "aperezcastro_track2_methods_form.xlsx")):
    wb = openpyxl.load_workbook(ROOT / "docs" / "methods_description_form.xlsx")
    vacias = rellena(wb[hoja], resp)
    dst = ROOT / "results" / "envios" / carpeta / nombre
    dst.parent.mkdir(parents=True, exist_ok=True)
    wb.save(dst)
    tr = "track1" if "1" in hoja else "track2"
    (ROOT / "reports" / f"methods_description_{tr}.md").write_text(anexo_md(wb[hoja], hoja) + "\n")
    (ROOT / "reports" / "latex" / f"methods_description_{tr}.tex").write_text(anexo_tex(wb[hoja], hoja))
    n_pal = len(resp["Provide a method abstract"].split())
    if n_pal > 500:                      # el formulario oficial pide «up to 500 words»
        raise SystemExit(f"✗ {dst.name}: el resumen tiene {n_pal} palabras (máximo 500)")
    print("✓", dst.relative_to(ROOT), "· sin respuesta:", vacias or "ninguna",
          f"· palabras del resumen: {n_pal}/500")
