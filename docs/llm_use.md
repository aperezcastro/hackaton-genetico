# LLM use and data handling

**Declaration** (organisers' format, Space discussion #2, 28 September 2026): Anthropic Claude via the Claude Code
CLI, Claude Max (consumer tier), model-training opt-out enabled ("Help improve Claude" off), no feedback provided on
any output, no content flagged for safety review.

## Why this is permitted

On 28 September 2026 the organisers confirmed in discussion #2 that consumer-tier accounts may be used with
patient-specific or derived information provided that (1) model training is opted out in the account settings and
(2) no feedback is given on outputs and no content is flagged for safety review. Both conditions held for the
whole project. Under the applicable terms (discussion #21), opted-out conversations are not used for training and
are retained for 30 days; the two exceptions (feedback and safety flags) were avoided.

## Safeguards in this repository

- `.claude/settings.json` disables the feedback, bug-report and transcript-survey commands and error reporting, and
  denies the agent read access to `data/`, where the controlled dataset and its largest derivatives live.
- `results/` holds candidate lists and intermediate tables (for example the ~1,000 rare functional variants left
  after filtering, and Exomiser's full output). From 1 October 2026 the agent could read this folder, under the
  rule that no more than ~50 annotated variants enter the model context at once; the larger tables were processed by
  scripts and never loaded whole.
- **Data minimisation**, as a precaution beyond the organisers' requirement: the model saw code, aggregate counts,
  HPO terms, the VCF header and short candidate lists. Genome-scale files were processed only by scripts.
- No patient variant was sent to a third-party API. Annotation (VEP, gnomAD, ClinVar, AlphaMissense, SpliceAI,
  Exomiser) ran on downloaded public resources in a private Modal volume (Modal acts as a data processor and does
  not train on customer data). Track 2 queried public services at gene, pathway and drug level only.
- `data/` and `results/` are git-ignored; the dataset, its derivatives and the Modal volume holding the case VCF are
  deleted within 30 days of the hackathon close, as the data-use terms require.

## Other models

No other LLM saw any patient-derived information. The voice of the Track 2 pitch video is a text-to-speech model
reading a script written from the published findings, with no genotype data.
