#!/usr/bin/env python3
"""hpo_extraer.py — extrae del documento de fenotipo del reto (data/raw/*.docx) SOLO los términos HPO
(HP:nnnnnnn + etiqueta si la hay) a results/t1/hpo.txt, para Exomiser/LIRICAL. No copia el resto del texto clínico.
Términos HPO permitidos en el contexto del agente (minimización de datos, docs/llm_use.md)."""
import re, sys, zipfile, pathlib
R = pathlib.Path(__file__).resolve().parents[1]
doc = next((R / "data" / "raw").glob("*Phenotype*.docx"))
xml = zipfile.ZipFile(doc).read("word/document.xml").decode("utf8")
txt = re.sub(r"<w:p[ >].*?</w:p>", lambda m: re.sub(r"<[^>]+>", "", m.group(0)) + "\n", xml, flags=re.S)
txt = re.sub(r"<[^>]+>", "", txt)
terms = {}
for line in txt.splitlines():
    for m in re.finditer(r"HP[:_ ]?(\d{7})", line):
        hp = f"HP:{m.group(1)}"
        lab = re.sub(r"HP[:_ ]?\d{7}", "", line).strip(" \t:-–—()[]|,;")
        terms.setdefault(hp, lab[:80])
out = R / "results" / "t1" / "hpo.txt"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("".join(f"{k}\t{v}\n" for k, v in terms.items()))
print(f"✓ {len(terms)} términos HPO → {out}")
