"""revision_casos.py — revisión sistemática de lo que la literatura publicada dice sobre TRATAMIENTOS recibidos por
pacientes con MVA (literatura pública; nada del paciente).

1. PubMed: ("mosaic variegated aneuploidy"[tiab] OR "premature chromatid separation"[tiab]) → todos los registros.
2. Cribado automático de título + resumen con términos de tratamiento (quimioterapia y fármacos concretos,
   radioterapia, trasplante, hormona de crecimiento, toxicidad, intensidad reducida).
3. Los registros con algún término se leen a mano; el resultado curado está en el informe del T2 (reports/latex/track2.pdf, revisión de casos MVA).

    .venv/bin/python scripts/t2/revision_casos.py

Escribe reports/t2/revision_casos.md (recuento por término y PMIDs). Caché en refs/t2_cache/.
"""
import collections
import json
import pathlib
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parents[2]
CACHE = ROOT / "refs" / "t2_cache"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = '("mosaic variegated aneuploidy"[tiab] OR "premature chromatid separation"[tiab])'
TERMINOS = {
    "chemotherapy (any)": r"chemotherap",
    "vincristine": r"vincristin",
    "dactinomycin": r"actinomycin|dactinomycin",
    "cyclophosphamide": r"cyclophosphamid",
    "anthracyclines": r"doxorubicin|anthracyclin|daunorubicin|epirubicin",
    "etoposide / irinotecan": r"etoposid|irinotecan",
    "radiotherapy": r"radiotherap|irradiat|radiation therap",
    "haematopoietic stem cell transplant": r"stem cell transplant|bone marrow transplant|\bhsct\b",
    "growth hormone": r"growth hormone|somatropin",
    "toxicity / adverse reaction": r"toxicit|adverse|severe reaction|intoleran",
    "reduced intensity / dose reduction": r"reduced[- ]intensity|dose reduc|reduced dose",
}


def get(url, nombre, data=None):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / nombre
    if not f.exists():
        req = urllib.request.Request(url, data=data.encode() if data else None)
        f.write_bytes(urllib.request.urlopen(req, timeout=120).read())
    return f.read_bytes()


def main():
    ids = json.loads(get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=1000&term={urllib.parse.quote(QUERY)}",
                         "mva_esearch.json"))["esearchresult"]["idlist"]
    xml = get(f"{EUTILS}/efetch.fcgi", "mva_efetch.xml", data=f"db=pubmed&rettype=abstract&retmode=xml&id={','.join(ids)}")
    arts = []
    for a in ET.fromstring(xml).iter("PubmedArticle"):
        ti = "".join(a.find(".//ArticleTitle").itertext())
        ab = " ".join("".join(x.itertext()) for x in a.iter("AbstractText"))
        arts.append((a.find(".//PMID").text, ti, ab))
    hits = collections.defaultdict(list)
    for pmid, ti, ab in arts:
        for k, rx in TERMINOS.items():
            if re.search(rx, (ti + " " + ab).lower()):
                hits[k].append(pmid)
    con_resumen = sum(1 for _, _, ab in arts if ab.strip())
    l = ["# Treatments reported in the published MVA literature — automated screen",
         "", f"PubMed query: `{QUERY}` → **{len(arts)} records** ({con_resumen} with an abstract). Title + abstract "
         "screened for treatment terms; every hit was then read manually (curated synthesis: Track 2 report, reports/latex/track2.pdf, MVA case review).",
         "", "| Term | Records | PMIDs |", "|---|---|---|"]
    l += [f"| {k} | {len(hits.get(k, []))} | {', '.join(hits.get(k, [])) or '—'} |" for k in TERMINOS]
    tot = sorted({p for v in hits.values() for p in v})
    l += ["", f"Records with at least one treatment term: **{len(tot)}** of {len(arts)}.",
          "", "Limitation: abstracts only; treatment details in full texts of case reports may be missed."]
    (ROOT / "reports" / "t2" / "revision_casos.md").write_text("\n".join(l) + "\n")
    print("\n".join(l))


if __name__ == "__main__":
    main()
