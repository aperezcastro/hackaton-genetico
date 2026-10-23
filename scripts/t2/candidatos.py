"""candidatos.py — Track 2: de un modelo de mecanismo (hipotesis.yaml) a una lista puntuada de fármacos aprobados.

Solo usa información a nivel gen/vía/fármaco (nunca genotipos del paciente): ChEMBL (estado de aprobación, avisos,
mecanismos y expansión por diana) y las fichas técnicas de la FDA vía openFDA (uso pediátrico, recuadro de
advertencia, advertencias). Las respuestas se guardan en refs/t2_cache/ para que la ejecución sea reproducible.

Puntuación (0–1), transparente y documentada en el informe:
  E  evidencia: peso del mejor estudio a favor (niveles `target_*` = la diana existe, no que el fármaco funcione) (+0,1 por cada estudio extra, máx. +0,2; −0,3 × Σ nivel de los estudios en contra, máx. −0,3)
  M  ajuste mecanístico al genotipo de este paciente (curado en hipotesis.yaml, justificado en `razon`); los
     fármacos añadidos por expansión de diana (no sembrados) llevan M × 0,9, porque nadie los ha revisado a mano
  P  viabilidad pediátrica: ficha FDA con uso pediátrico establecido 1 · ficha sin uso pediátrico establecido 0,6 ·
     aprobado fuera de EE. UU. / sin ficha 0,5 · suplemento o no aprobado 0,35
  S  seguridad para ESTE paciente: 1 − penalización por avisos de la ficha que chocan con su perfil de riesgo (HPO):
     0,25 por bandera en el recuadro de advertencia, 0,1 en advertencias, más `riesgo_curado` de hipotesis.yaml
     (riesgos documentados en la literatura que la ficha no recoge), hasta un máximo de 0,7
  score = 0,30·E + 0,25·M + 0,20·P + 0,25·S  (+ análisis de sensibilidad de los pesos, Dirichlet)

    .venv/bin/python scripts/t2/candidatos.py [--sin-red]

Salidas: reports/t2/candidatos.csv, reports/t2/candidatos.md (tabla de candidatos + tarjeta de seguridad).
"""
import csv
import hashlib
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[2]
CACHE = ROOT / "refs" / "t2_cache"
OUT = ROOT / "reports" / "t2"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
OPENFDA = "https://api.fda.gov/drug/label.json"
SIN_RED = "--sin-red" in sys.argv

# Pesos de evidencia. Distinguimos evidencia de BENEFICIO de la intervención (la que puntúa alto) de evidencia de
# VALIDACIÓN DE DIANA (la vía está alterada en pacientes): saber que una vía está activa no demuestra que
# modularla ayude, así que pesa como un estudio celular.
NIVEL = {"human_same_gene": 1.0, "human_other": 0.7, "mouse_same_gene": 0.8, "mouse_other": 0.5,
         "animal_other": 0.4,            # modelo animal no mamífero (p. ej. Drosophila) de la misma vía
         "feasibility_other": 0.3,       # tolerabilidad/viabilidad en otra enfermedad, sin beneficio en este mecanismo
         "cell_same_gene": 0.5, "cell_other": 0.3,
         "target_same_gene": 0.5, "target_other": 0.3}
BANDERAS = {
    "malignancy": r"malignan|lymphoma|carcinogen|second(ary)? (primary )?(cancer|tumou?r|neoplasm)",
    "nephrotoxicity": r"nephrotox|renal (failure|toxicity|impairment)|kidney injury",
    "myopathy": r"(?<!cardio)myopath|rhabdomyolysis|muscle weakness",
    "myelosuppression": r"myelosuppress|neutropenia|bone marrow suppression|pancytopenia",
    "growth": r"growth (retardation|suppression|delay)|epiphys|linear growth",
}


def get_json(url):
    """GET con caché en disco (las APIs públicas cambian; la caché fija la ejecución del informe)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / (hashlib.sha1(url.encode()).hexdigest()[:16] + ".json")
    if f.exists():
        return json.loads(f.read_text())
    if SIN_RED:
        return None
    for intento in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}),
                                        timeout=60) as r:
                d = json.loads(r.read())
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                d = None
                break
            if intento == 3:
                raise
        except urllib.error.URLError:
            if intento == 3:
                raise
        time.sleep(2 ** intento)
    time.sleep(0.25)
    f.write_text(json.dumps(d))
    return d


def chembl_molecula(nombre):
    q = urllib.parse.quote(nombre)
    for filtro in ("pref_name__iexact", "molecule_synonyms__molecule_synonym__iexact"):
        d = get_json(f"{CHEMBL}/molecule.json?{filtro}={q}&limit=5")
        mols = (d or {}).get("molecules", [])
        if mols:
            return max(mols, key=lambda m: float(m.get("max_phase") or 0))
    return None


def chembl_mecanismos(chembl_id):
    """Mecanismos de la molécula o, si están anotados en otra forma (sal/parental), de sus formas hermanas."""
    formas = [chembl_id] + [f["molecule_chembl_id"] for f in
                            (get_json(f"{CHEMBL}/molecule_form/{chembl_id}.json") or {}).get("molecule_forms", [])]
    for cid in dict.fromkeys(formas):
        d = get_json(f"{CHEMBL}/mechanism.json?molecule_chembl_id={cid}&limit=20")
        mec = [m["mechanism_of_action"] for m in (d or {}).get("mechanisms", []) if m.get("mechanism_of_action")]
        if mec:
            return list(dict.fromkeys(mec))
    return []


def chembl_aprobados_por_diana(simbolo):
    """Fármacos aprobados (max_phase 4, no retirados) cuyo mecanismo en ChEMBL actúa sobre la diana humana `simbolo`."""
    d = get_json(f"{CHEMBL}/target/search.json?q={urllib.parse.quote(simbolo)}&limit=20")
    dianas = [t["target_chembl_id"] for t in (d or {}).get("targets", [])
              if t.get("organism") == "Homo sapiens" and t.get("target_type") == "SINGLE PROTEIN"
              and any(s.get("component_synonym", "").upper() == simbolo
                      for c in t.get("target_components", []) for s in c.get("target_component_synonyms", []))]
    nombres = set()
    for tid in dianas:
        mec = get_json(f"{CHEMBL}/mechanism.json?target_chembl_id={tid}&limit=100")
        for m in (mec or {}).get("mechanisms", []):
            mol = get_json(f"{CHEMBL}/molecule/{m['molecule_chembl_id']}.json")
            if mol and float(mol.get("max_phase") or 0) >= 4 and not mol.get("withdrawn_flag") and mol.get("pref_name"):
                nombres.add(mol["pref_name"])
    return sorted(nombres)


VIAS_SISTEMICAS = {"ORAL", "INTRAVENOUS", "SUBCUTANEOUS", "INTRAMUSCULAR"}


def fda_fichas(nombre):
    """Hasta 10 fichas (marca y genéricos): los genéricos a veces omiten («carve-out») la indicación pediátrica.
    Si el fármaco tiene fichas por vía sistémica, solo se usan esas: la búsqueda por nombre genérico mezcla
    presentaciones locales (p. ej. la crema de ruxolitinib, con un recuadro de advertencia que el comprimido no lleva)."""
    q = urllib.parse.quote(f'openfda.generic_name:"{nombre}"')
    d = get_json(f"{OPENFDA}?search={q}&limit=10")
    fichas = (d or {}).get("results", [])
    sistemicas = [f for f in fichas if {r.upper() for r in f.get("openfda", {}).get("route", [])} & VIAS_SISTEMICAS]
    return sistemicas or fichas


PED_SI = (r"(safety and (effectiveness|efficacy)|effectiveness and safety)[^.]{0,120}(has|have) been established"
          r"|pediatric use information is approved for|indicated (for|in) (the treatment of )?pediatric"
          r"|pediatric patients (\d|aged|from|weighing|1 year)")


def pediatria(fichas, mol):
    if fichas:
        for f in fichas:
            txt = " ".join(f.get("pediatric_use", []) + f.get("indications_and_usage", [])
                           + f.get("dosage_and_administration", [])).lower()
            txt = re.sub(r"(have|has) not been established", "", txt)
            if re.search(PED_SI, txt):
                return 1.0, "FDA label: paediatric use/dosing established"
        return 0.6, "FDA label: paediatric use not established"
    if mol and float(mol.get("max_phase") or 0) >= 4:
        return 0.5, "approved (ChEMBL phase 4), no current US label found"
    return 0.35, "supplement / not an approved drug"


def estado_regulatorio(mol):
    """Track 2 se limita a medicamentos aprobados y comercializados (organizadores, discusión #5 del Space)."""
    fase = float((mol or {}).get("max_phase") or 0)
    if fase >= 4 and not (mol or {}).get("withdrawn_flag"):
        return "approved medicine"
    if fase >= 4:
        return "approved, withdrawn in some markets"
    return "NOT an approved medicine (supplement / investigational)"


FORMAS_PED = ("SOLUTION", "SUSPENSION", "FOR SUSPENSION", "SYRUP", "ELIXIR", "DISPERSIBLE", "FOR SOLUTION",
              "GRANULE", "PELLET", "POWDER, FOR")


def formulacion_pediatrica(nombre):
    """¿Existe en EE. UU. una forma oral líquida o dispersible (dosificable por peso en niños)? Directorio NDC de la FDA."""
    q = urllib.parse.quote(f'generic_name:"{nombre}"')
    res = (get_json(f"https://api.fda.gov/drug/ndc.json?search={q}&limit=100") or {}).get("results", [])
    formas = {r.get("dosage_form", "").upper() for r in res if r.get("dosage_form")}
    liquidas = sorted(f for f in formas if any(k in f for k in FORMAS_PED))
    return ", ".join(liquidas).lower() if liquidas else ("solid oral/other forms only" if formas else "not in US NDC directory")


def banderas(fichas, riesgos):
    if not fichas:
        return 0.0, []
    caja = " ".join(x for f in fichas for x in f.get("boxed_warning", [])).lower()
    aviso = " ".join(x for f in fichas for x in f.get("warnings_and_cautions", []) + f.get("warnings", [])).lower()
    pen, hits = 0.0, []
    for b, rx in BANDERAS.items():
        if b not in riesgos:
            continue
        if re.search(rx, caja):
            pen, hits = pen + 0.25, hits + [f"{b} (boxed)"]
        elif re.search(rx, aviso):
            pen, hits = pen + 0.1, hits + [b]
    return min(pen, 0.6), hits


def evidencia(h):
    """E = mejor estudio a favor + 0,1 por cada uno más (máx. +0,2) − Σ 0,3·nivel de los estudios en contra (máx. −0,3).

    Solo cuenta evidencia sobre si FUNCIONA (`support`/`contrary`). Lo que dice si DAÑA a este niño (`harm`) va a
    la seguridad S, no aquí: antes se restaba en los dos sitios (corregido el 4 oct tras la revisión independiente).
    `note` se documenta y no puntúa. Las contrarias se ponderan por su nivel y llevan tope, para que E no dependa
    de cuánto buscó el curador evidencia en contra de cada línea."""
    a_favor = sorted((NIVEL[e["nivel"]] for e in h["evidencia"] if e["dir"] == "support"), reverse=True)
    en_contra = [NIVEL[e["nivel"]] for e in h["evidencia"] if e["dir"] == "contrary"]
    e = (a_favor[0] if a_favor else 0) + min(0.1 * (len(a_favor) - 1), 0.2) - min(0.3 * sum(en_contra), 0.3)
    return max(0.0, min(1.0, e))


def conectividad():
    """Resumen CMap por fármaco (de scripts/t2/cmap.py): nº de líneas celulares en que invierte (R) o imita (M) la
    firma del knockout de BUB1B. Informativo: no entra en la puntuación (la firma es de líneas tumorales)."""
    res = {}
    for clase, letra in (("inversores", "R"), ("imitadores", "M")):
        f = OUT / f"cmap_{clase}.csv"
        for r in csv.DictReader(f.open()) if f.exists() else []:
            res.setdefault(r["compuesto"], {})[letra] = int(r["n_lineas"])
    alias = {"tg-101348": "fedratinib"}                 # nombres de código en LINCS → nombre genérico
    return {alias.get(k, k): f"reverses {v.get('R', 0)} / mimics {v.get('M', 0)} lines" for k, v in res.items()}


PESOS = (0.30, 0.25, 0.20, 0.25)   # E, M, P, S: la seguridad pesa tanto como el ajuste (niño con predisposición al cáncer)


def puntua(E, M, P, S, w=PESOS):
    return w[0] * E + w[1] * M + w[2] * P + w[3] * S


def rangos(cand, n=2000, semilla=7, ruido_m=0.0, ruido_s=0.0):
    """n sorteos: pesos ~ Dirichlet centrado en PESOS; si ruido_m > 0, M de cada LÍNEA ± ruido_m (uniforme); si
    ruido_s > 0, S de cada LÍNEA ± ruido_s. Se perturba por línea porque M y los riesgos curados son juicios sobre la
    hipótesis, no sobre el fármaco. Los empates se rompen al azar en cada sorteo (antes los decidía el orden de la
    lista y dos fármacos idénticos salían con robustez distinta). Devuelve fármaco → lista de puestos."""
    import random
    rnd = random.Random(semilla)
    rk = {f["farmaco"]: [] for f in cand}
    lineas = sorted({f["hipotesis"] for f in cand})
    for _ in range(n):
        g = [rnd.gammavariate(20 * w, 1) for w in PESOS]
        w = [x / sum(g) for x in g]
        dm = {h: rnd.uniform(-ruido_m, ruido_m) for h in lineas}
        ds = {h: rnd.uniform(-ruido_s, ruido_s) for h in lineas}
        orden = sorted(cand, key=lambda f: (-puntua(float(f["E"]), min(1, max(0, float(f["M"]) + dm[f["hipotesis"]])),
                                                    float(f["P"]), min(1, max(0, float(f["S"]) + ds[f["hipotesis"]])), w),
                                            rnd.random()))
        for i, f in enumerate(orden, 1):
            rk[f["farmaco"]].append(i)
    return rk


def sensibilidad(cand, n=2000, semilla=7, ruido_m=0.0, ruido_s=0.0):
    """Por fármaco: (puesto medio, % de sorteos en el top 3, % de sorteos en el puesto 1)."""
    rk = rangos(cand, n, semilla, ruido_m, ruido_s)
    return {k: (sum(v) / n, 100 * sum(r <= 3 for r in v) / n, 100 * sum(r == 1 for r in v) / n) for k, v in rk.items()}


# Escenarios contrafactuales: ¿qué cambiaría el orden si el equipo clínico estableciera los dos hechos que hoy
# bloquean la línea correctora? Se recalcula con los mismos E y P; solo cambian M o S de H1.
ESCENARIOS = [
    ("as scored today", {}, {}),
    ("allele class A confirmed (rung 1 blot)", {"H1": 0.8}, {}),
    ("documented remission (NAMPT concern 0.35 → 0.10)", {}, {"H1": 0.10}),
    ("both", {"H1": 0.8}, {"H1": 0.10}),
]


def contrafactual(cfg, cand):
    lin = {h["id"]: h for h in cfg["hipotesis"]}
    salida = []
    for nombre, nuevo_m, nuevo_riesgo in ESCENARIOS:
        filas = []
        for f in cand:
            h = lin[f["hipotesis"]]
            if f["hipotesis"] not in nuevo_m and f["hipotesis"] not in nuevo_riesgo:
                # línea sin cambios: se reutiliza su puntuación, no se recalcula desde los valores redondeados
                filas.append((round(float(f["score"]), 3), f))
                continue
            M, S = float(f["M"]), float(f["S"])
            if f["hipotesis"] in nuevo_m:
                M = nuevo_m[f["hipotesis"]] + (h.get("ajuste_farmaco") or {}).get(f["farmaco"].upper(), 0.0)
            if f["hipotesis"] in nuevo_riesgo:
                S = S + (h.get("riesgo_curado") or {}).get("penalizacion", 0) - nuevo_riesgo[f["hipotesis"]]
            filas.append((round(puntua(float(f["E"]), M, float(f["P"]), S), 3), f))
        orden = {id(f): i for i, f in enumerate(cand)}       # empates: el mismo orden que la tabla principal
        filas.sort(key=lambda x: (-x[0], orden[id(x[1])]))
        salida.append((nombre, filas))
    return salida


def main():
    cfg = yaml.safe_load((pathlib.Path(__file__).parent / "hipotesis.yaml").read_text())
    riesgos = set(cfg["caso"]["riesgos"])
    cmap = conectividad()
    filas = []
    for h in cfg["hipotesis"]:
        E = evidencia(h)
        farmacos = [(n, "seed") for n in h["semillas"]]
        for diana in h.get("dianas_expandir", []):
            farmacos += [(n, f"ChEMBL target {diana}") for n in chembl_aprobados_por_diana(diana)
                         if n not in h["semillas"]]
        for nombre, origen in farmacos:
            mol = chembl_molecula(nombre)
            fichas = fda_fichas(nombre)
            P, p_txt = pediatria(fichas, mol)
            pen, hits = banderas(fichas, riesgos)
            M = (h["ajuste"] + (h.get("ajuste_farmaco") or {}).get(nombre, 0.0)) * (1.0 if origen == "seed" else 0.9)
            # riesgo documentado en literatura que la ficha técnica NO recoge (p. ej. una dependencia metabólica
            # del tumor del paciente): se declara en hipotesis.yaml y se suma a la penalización de seguridad
            rc = h.get("riesgo_curado") or {}
            if rc:
                pen, hits = min(pen + rc["penalizacion"], 0.7), hits + [f"{rc['motivo'][:60]}… (literature)"]
            S = 1 - pen
            filas.append({
                "hipotesis": h["id"], "modo": h["modo"], "eje": h.get("eje", ""),
                "eje_riesgo": h.get("eje_riesgo", ""), "farmaco": nombre.lower(), "origen": origen,
                "chembl_id": (mol or {}).get("molecule_chembl_id", ""),
                "max_phase": (mol or {}).get("max_phase") or "", "primera_aprobacion": (mol or {}).get("first_approval") or "",
                "mecanismo_chembl": "; ".join(chembl_mecanismos(mol["molecule_chembl_id"]))[:120] if mol else "",
                "estado_regulatorio": estado_regulatorio(mol),
                "formulacion_pediatrica": formulacion_pediatrica(nombre),
                "pediatria": p_txt, "banderas_paciente": ", ".join(hits),
                "condiciones": h.get("condiciones", ""),
                "E": round(E, 2), "M": round(M, 2), "P": P, "S": round(S, 2),
                "score": round(puntua(E, M, P, S), 3) if h["modo"] != "avoid" else "",
                "cmap": cmap.get(nombre.lower(), "not profiled / n.s."),
                "pmids": ",".join(str(e["pmid"]) for e in h["evidencia"] if "pmid" in e),
            })
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "candidatos.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    escribe_md(cfg, filas)
    print(f"✓ {len(filas)} fármacos ({sum(f['modo'] != 'avoid' for f in filas)} candidatos, "
          f"{sum(f['modo'] == 'avoid' for f in filas)} en la tarjeta de seguridad) → reports/t2/")


def escribe_md(cfg, filas):
    nombres = {h["id"]: h["nombre"] for h in cfg["hipotesis"]}
    todos = sorted((f for f in filas if f["modo"] != "avoid"), key=lambda f: -f["score"])
    cand = [f for f in todos if f["estado_regulatorio"].startswith("approved")]
    fuera = [f for f in todos if not f["estado_regulatorio"].startswith("approved")]
    sens = sensibilidad(cand)
    sens_t = sensibilidad(cand, ruido_m=0.30, ruido_s=0.15)
    l = ["# Track 2 — scored candidates (generated by `scripts/t2/candidatos.py`)", "",
         "Score = 0.30·Evidence + 0.25·Mechanistic fit + 0.20·Paediatric feasibility + 0.25·Patient-specific safety. "
         "Only **approved medicines** are ranked (Track 2 scope, organisers' discussion #5). E counts only evidence "
         "about whether a line works; evidence that it could harm this child lowers S instead (it used to be "
         "subtracted in both places). Robustness: 2,000 draws of the weights (Dirichlet centred on the chosen "
         "weights); the last three columns add, in the same draws, a uniform ±0.30 perturbation of each line's "
         "mechanistic fit M and ±0.15 of its safety S — the two most curated inputs — with ties broken at random.", "",
         "| # | Drug | Line | Axis | E | M | P | S | Score | % top-3 (weights) | Mean rank (weights + M ±0.30 + S ±0.15) | % top-3 (all) | % ranked #1 (all) | Liquid/dispersible form (FDA NDC) | Label/literature flags vs patient risks | CMap (BUB1B-KO signature) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, f in enumerate(cand, 1):
        _, t3, _ = sens[f["farmaco"]]
        mr, t3t, p1 = sens_t[f["farmaco"]]
        eje = f["eje"] + (f" (risk: {f['eje_riesgo']})" if f["eje_riesgo"] else "")
        l.append(f"| {i} | {f['farmaco']} | {f['hipotesis']} | {eje} | {f['E']} | {f['M']} | {f['P']} | {f['S']} | "
                 f"**{f['score']}** | {t3:.0f} | {mr:.1f} | {t3t:.0f} | {p1:.0f} | {f['formulacion_pediatrica']} | "
                 f"{f['banderas_paciente'] or '—'} | {f['cmap']} |")
    l += ["", "### Not ranked: products that are not approved medicines", "",
          "Listed for transparency (they appear in the mechanistic literature) but outside the Track 2 scope of "
          "market-approved medications.", "",
          "| Drug | Line | Regulatory status (ChEMBL max phase) | Score if it were eligible |", "|---|---|---|---|"]
    l += [f"| {f['farmaco']} | {f['hipotesis']} | {f['estado_regulatorio']} (phase {f['max_phase'] or '—'}) | {f['score']} |"
          for f in fuera]
    # transparencia de la curación: cuántos estudios hay en cada dirección por línea
    l += ["", "## Lines, and how much evidence was curated for each", "",
          "| Line | Name | Axis | E | for | against (scored) | harm → S | notes (not scored) |", "|---|---|---|---|---|---|---|---|"]
    for h in cfg["hipotesis"]:
        cuenta = {d: sum(e["dir"] == d for e in h["evidencia"]) for d in ("support", "contrary", "harm", "note")}
        l.append(f"| {h['id']} | {h['nombre']} | {h.get('eje', '')} | {round(evidencia(h), 2)} | {cuenta['support']} | "
                 f"{cuenta['contrary']} | {cuenta['harm']} | {cuenta['note']} |")
    # eje terapéutico POR LÍNEA (no media por fármaco: la media dependería de cuántos fármacos se sembraron)
    l += ["", "## Therapeutic axis, line by line", "",
          "`correct` fewer missegregation events · `clear` remove the abnormal cells · `tolerate` help abnormal "
          "cells cope · `avoid` do not add insult. Summarised per **line** (its best approved drug), because a mean "
          "over drugs would depend on how many drugs were seeded per line.", "",
          "| Line | Axis | Best approved drug | Score | Rank | Blocked by |", "|---|---|---|---|---|---|"]
    bloqueo = {h["id"]: h.get("bloqueo", "—") for h in cfg["hipotesis"]}
    for h in cfg["hipotesis"]:
        fs = [f for f in cand if f["hipotesis"] == h["id"]]
        if not fs:
            continue
        b = fs[0]
        l.append(f"| {h['id']} | {h.get('eje', '')} | {b['farmaco']} | {b['score']} | {cand.index(b) + 1} | "
                 f"{bloqueo[h['id']]} |")
    # escenarios contrafactuales
    l += ["", "## What would change the ranking: two facts the clinical team can establish", "",
          "Same E and P; only the inputs that encode today's two patient-specific uncertainties about the corrective "
          "NAD⁺ line change — the allele class (M of H1 = 0.8 × P(class A), now 0.6) and the active-tumour concern "
          "(curated S penalty of H1, now 0.35).", "",
          "| Scenario | #1 | #2 | #3 | Best corrective drug (rank) | Top 7 by axis (correct / clear / tolerate) |",
          "|---|---|---|---|---|---|"]
    for nombre, fs in contrafactual(cfg, cand):
        top = [f for _, f in fs]
        corr = next(((i + 1, f) for i, f in enumerate(top) if f["eje"] == "correct"), (None, None))
        c7 = [f["eje"] for f in top[:7]]
        l.append(f"| {nombre} | {top[0]['farmaco']} | {top[1]['farmaco']} | {top[2]['farmaco']} | "
                 f"{corr[1]['farmaco'] if corr[1] else '—'} ({corr[0]}) | "
                 f"{c7.count('correct')} / {c7.count('clear')} / {c7.count('tolerate')} |")
    l += ["", "## Safety card (S1): drugs to avoid or dose-adapt", "",
          "| Drug | ChEMBL mechanism | Label flags vs patient risks | CMap (BUB1B-KO signature) |", "|---|---|---|---|"]
    l += [f"| {f['farmaco']} | {f['mecanismo_chembl'] or '—'} | {f['banderas_paciente'] or '—'} | {f['cmap']} |"
          for f in filas if f["modo"] == "avoid"]
    (OUT / "candidatos.md").write_text("\n".join(l) + "\n")


if __name__ == "__main__":
    main()
