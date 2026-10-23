# MVA Hackathon 2026 — biallelic *BUB1B*: variant interpretation and allele-aware drug repurposing

Submission by **Agustín Pérez Castro** (`aperezcastro`) to
[Rare Disease, Real Kid: The MVA Hackathon 2026](https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/).

| Track | Result | Report |
|---|---|---|
| 1 — Variant prediction | *BUB1B* compound heterozygous: c.2210T>G p.(Leu737\*) + c.3006T>G p.(Asn1002Lys) → MVA type 1; Exomiser rank 1 genome-wide | [`reports/latex/track1.pdf`](reports/latex/track1.pdf) |
| 2 — Drug repurposing | allele- and patient-aware pipeline. **Main result: a decision rule, not a single drug.** Every line is classified by what it does to a chromosomally unstable cell (correct / clear / tolerate). Today intermittent dasatinib leads; the mechanism-matched NAD⁺ line ranks 7th, gated on two measurable facts about this child (missense allele class and remission). Plus a pharmacogenomic safety card, a weeks-scale endpoint never measured in MVA (%MN-RET) with its n-of-1 power calculation, and a one-page rung-0 protocol | [`reports/latex/track2.pdf`](reports/latex/track2.pdf), [`reports/t2/safety_card.pdf`](reports/t2/safety_card.pdf), [`reports/t2/rung0_protocol.pdf`](reports/t2/rung0_protocol.pdf) · pitch: https://youtu.be/YzJo-CZt5Ok |

**No patient data are in this repository.** The controlled dataset (VCF, FASTQ) and every genome-scale derivative live
in `data/` and `results/`, which are git-ignored and are deleted ≤ 30 days after the hackathon, as the data-use terms
require.

## Repository layout

```
scripts/
  hpo_extraer.py          HPO term IDs from the clinical summary (IDs only; labels resolved against a local hp.obo)
  modal/                  Track 1 compute on Modal (each job: hard timeout, ≤ 2 retries, cost ceiling, cost log)
    comun.py              credentials, cost ceiling, watchdog, gastos.csv logging
    refs.py               public references: VEP 116 cache (gnomAD v4.1, ClinVar), Ensembl fasta, AlphaMissense, ClinVar
    exomiser_refs.py      Exomiser 15.1 + data release 2512
    anotar.py             bcftools norm + VEP 116 annotation, one container per chromosome
    filtra.py             genome-wide funnel, recessive gene model, MVA/SAC panel (all consequences)
    splice_ai.py          SpliceAI 1.3.1 on rare panel variants (50 and 2000 nt windows)
    exomiser.py           Exomiser with and without a phenotype pre-filter
    fase_viable.py        tests whether the pair is phasable with short reads (it is not)
    balance_alelico.py    genome-level QC from aggregates (allelic balance, runs of homozygosity)
    estado_vol.py         lists the Modal volumes (housekeeping)
    gastos.csv            cost log of every Modal run
  t2/                     Track 2 (public gene/drug-level data only; runs on a laptop)
    hipotesis.yaml        mechanism model + therapeutic lines, PMID-cited, evidence level and direction
    candidatos.py         ChEMBL + openFDA annotation, patient-specific safety flags, scoring + sensitivity
    cmap.py               LINCS L1000 connectivity of the BUB1B-knockout signature (SigCom LINCS, Enrichr)
    estructura.py         AlphaFold structural context of the BUBR1 alleles
    contexto_stop.py      sequence context of the nonsense codon (read-through assessment)
    revision_casos.py     review of published MVA cases that received treatment
    nof1.py               n-of-1 power simulation on micronucleated reticulocytes (%MN-RET), the endpoint
                          proposed for measuring effect in this child, and the design table used by the one-page
                          rung-0 protocol (reports/t2/rung0_protocol.md)
  formulario_metodos.py   fills the official methods description form
  md_a_pdf.py             Markdown report → PDF (headless Chrome)
reports/latex/            Track 1 and Track 2 reports (LuaLaTeX, IBM Plex fonts, figures in fig/)
reports/t2/               Track 2 tables, safety card, connectivity results (no genome-scale data)
scripts/figuras/          report figures (matplotlib, validated colour palette)
scripts/compila_informes.py   rebuilds figures, the methods-form annex and both reports
reports/methods_description_track{1,2}.{md,pdf}   answers to the official methods description form
reports/submission/       the Track 1 CSV as submitted (organisers' format)
research_notes/           working notes of the Track 2 evidence reviews (superseded by the report where they differ)
docs/                     official templates (submission CSV, methods form) and the LLM-use statement
.claude/settings.json     agent configuration: data/ blocked, feedback commands and error reporting disabled
```

## Reproduce

**Requirements.** Python 3.12 with `uv`; a [Modal](https://modal.com) account (Track 1 only); access to the gated
dataset `SageBio/mva-hackathon-2026-data`; for the documents, LuaLaTeX with `latexmk` and Google Chrome (headless,
for the Markdown PDFs), plus `pdfinfo`/`pdftotext` (poppler).

**Local files that are not in the repository** (git-ignored on purpose): `data/` (the controlled dataset and its
genome-scale derivatives), `results/` (candidate lists and figure inputs), `refs/` (public reference downloads and
the Track 2 API cache) and `.env`, a plain `key=value` file with `modal-token-id`, `modal-token-secret` and
`hf-token`.

```bash
uv venv && source .venv/bin/activate
uv pip install modal pyyaml openpyxl markdown biopython numpy scipy matplotlib huggingface_hub
# 0. the case data, into data/raw/ (needs the dataset access and hf-token)
python -c "from huggingface_hub import hf_hub_download as d; [d('SageBio/mva-hackathon-2026-data', f, repo_type='dataset', local_dir='data/raw') for f in ('WGS_EX2312012_HGWCNDSX7.vcf.gz', 'WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi')]"
#    (the phenotype .docx is downloaded the same way; see the dataset README for its name)
# Track 1 (≈ 1.4 USD of cloud compute in total, logged in scripts/modal/gastos.csv)
python scripts/hpo_extraer.py
python scripts/modal/refs.py && python scripts/modal/exomiser_refs.py
python scripts/modal/anotar.py --subir          # VEP per chromosome, then the funnel and panel (filtra.py runs inside)
python scripts/modal/splice_ai.py && python scripts/modal/exomiser.py && python scripts/exomiser_top.py
python scripts/modal/fase_viable.py
python scripts/modal/balance_alelico.py && python scripts/modal/balance_alelico.py --diagnostico \
  && python scripts/modal/balance_alelico.py --roh
# Track 2 (minutes, laptop; public API responses are cached in refs/t2_cache/ on first run)
python scripts/t2/estructura.py && python scripts/t2/contexto_stop.py && python scripts/t2/revision_casos.py
python scripts/t2/cmap.py && python scripts/t2/cmap.py --comparar BUB1B,CEP57,TRIP13,BUB3,MAD2L1,BUB1,TTK
python scripts/t2/candidatos.py && python scripts/t2/nof1.py
# Reports: figures, the methods-form annex and both LaTeX reports
python scripts/compila_informes.py
python scripts/md_a_pdf.py reports/t2/safety_card.md reports/t2/safety_card.pdf --compacto
python scripts/md_a_pdf.py reports/t2/rung0_protocol.md reports/t2/rung0_protocol.pdf --compacto
python scripts/md_a_pdf.py reports/methods_description_track1.md reports/methods_description_track1.pdf
python scripts/md_a_pdf.py reports/methods_description_track2.md reports/methods_description_track2.pdf
```

**Reproducibility limits.** The Track 2 inputs come from live public services (ChEMBL, openFDA, SigCom LINCS,
Enrichr). The cache in `refs/t2_cache/` is not published, and Enrichr calls are not cached, so a fresh run can differ
slightly as those databases change; the committed `reports/t2/*.csv` and `*.md` are the outputs used in the report.
The submission CSV (`reports/submission/`) was validated with the organisers' public scorer `evaluation.py` from the
[hackathon Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).

## Data protection and LLM use

- Analysis assisted by Anthropic Claude (Claude Code, Claude Max consumer plan, "Help improve Claude" disabled, no
  feedback submitted). Data minimisation: the model saw code, aggregate counts, HPO terms and short candidate lists
  only. The raw data and genome-scale derivatives in `data/` were blocked from the agent by configuration
  (`.claude/settings.json`); intermediate tables in `results/` were processed by scripts, never loaded whole. See
  [`docs/llm_use.md`](docs/llm_use.md).
- Patient variants were never sent to third-party APIs: annotation used downloaded public resources on a private
  Modal volume. Track 2 queried public services at gene, pathway and drug level only.

## Licence and acknowledgement

Results and code: CC BY 4.0, as required by the hackathon rules.

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA
Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize
sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed
their data and their story to advance research into this rare disease. We acknowledge their trust in making this
Hackathon possible.
