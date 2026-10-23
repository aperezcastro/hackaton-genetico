"""cmap.py — línea de conectividad (tipo CMap) para el Track 2 con SigCom LINCS (datos públicos L1000, nivel gen).

1. Busca todas las firmas de knockout CRISPR de BUB1B (base l1000_xpr) y descarga su ranking de genes.
   Convención comprobada: rango alto = gen subido (CDKN1A y MDM2 arriba, BUB1B abajo).
2. Consenso: rango medio normalizado por gen entre firmas → top N subidos / bajados.
3. Busca en l1000_cp (perturbaciones químicas) las firmas que se le parecen (imitadores) o que la invierten
   (inversores) con la prueba de rango a dos colas de SigCom (enrich/ranktwosided).
4. Agrega por compuesto (nº de firmas significativas, mejor z) y escribe reports/t2/cmap_{inversores,imitadores}.csv.

Control de validez: los IMITADORES deberían incluir fármacos antimitóticos/inhibidores de quinasas mitóticas.

    .venv/bin/python scripts/t2/cmap.py [--n 150]
    .venv/bin/python scripts/t2/cmap.py --comparar BUB1B,CEP57,TRIP13   # escalabilidad: firmas de otros genes MVA

El modo --comparar calcula la firma de consenso de cada gen (mismas reglas), su correlación de Spearman con la de
BUB1B sobre los ~12.000 genes del ranking, y las vías Hallmark principales de cada una; escribe
reports/t2/cmap_comparacion.md.
"""
import collections
import csv
import hashlib
import json
import pathlib
import sys
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(pathlib.Path(__file__).parent))
CACHE = ROOT / "refs" / "t2_cache"
OUT = ROOT / "reports" / "t2"
META = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"
N = int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 150


def post(url, cuerpo):
    """POST JSON con caché en disco y reintentos (las respuestas fijan la ejecución del informe)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    clave = hashlib.sha1((url + json.dumps(cuerpo, sort_keys=True)).encode()).hexdigest()[:16]
    f = CACHE / f"cmap_{clave}.json"
    if f.exists():
        return json.loads(f.read_text())
    for intento in range(4):
        try:
            req = urllib.request.Request(url, data=json.dumps(cuerpo).encode(),
                                         headers={"Content-Type": "application/json"}, method="POST")
            with urllib.request.urlopen(req, timeout=300) as r:
                d = json.loads(r.read())
            break
        except urllib.error.HTTPError as e:
            if e.code in (307, 308):          # urllib no reenvía POST en redirecciones: seguirla a mano
                url = e.headers["Location"]
                continue
            if intento == 3:
                raise
            time.sleep(2 ** intento)
        except urllib.error.URLError:
            if intento == 3:
                raise
            time.sleep(2 ** intento)
    f.write_text(json.dumps(d))
    return d


def firmas_knockout(gen):
    d = post(f"{META}/signatures/find", {"filter": {"where": {"meta.pert_name": gen}, "limit": 500}})
    return [x for x in d if x["meta"].get("pert_type", "").startswith("CRISPR")]


def main():
    firmas = firmas_knockout("BUB1B")
    print(f"{len(firmas)} firmas CRISPR KO de BUB1B en {len({f['meta'].get('cell_line') for f in firmas})} líneas")
    r = post(f"{DATA}/fetch/rank", {"database": "l1000_xpr", "signatures": [f["id"] for f in firmas], "entities": []})
    ent, maxr = r["entities"], r["maxrank"]
    media = [0.0] * len(ent)
    for s in r["signatures"]:
        for i, x in enumerate(s["ranks"]):
            media[i] += x / maxr / len(r["signatures"])
    orden = sorted(range(len(ent)), key=lambda i: media[i])
    abajo, arriba = [ent[i] for i in orden[:N]], [ent[i] for i in orden[-N:]]

    sim = {x["id"]: x["meta"]["symbol"] for x in post(f"{META}/entities/find", {
        "filter": {"where": {"id": {"inq": abajo + arriba}}, "fields": ["id", "meta.symbol"]}})}
    print("top subidos:", ", ".join(sim.get(e, "?") for e in arriba[-15:][::-1]))
    print("top bajados:", ", ".join(sim.get(e, "?") for e in abajo[:15]))

    res = post(f"{DATA}/enrich/ranktwosided", {"up_entities": arriba, "down_entities": abajo, "limit": 1000,
                                                "database": "l1000_cp"})
    hits = res.get("results", [])
    meta = {x["id"]: x["meta"] for x in post(f"{META}/signatures/find", {
        "filter": {"where": {"id": {"inq": [h["uuid"] for h in hits]}}, "limit": len(hits) + 10}})}
    grupos = {"imitadores": collections.defaultdict(list), "inversores": collections.defaultdict(list)}
    for h in hits:
        m = meta.get(h["uuid"], {})
        nombre = (m.get("pert_name") or "?").lower()
        z_up, z_dn = h.get("z-up", 0), h.get("z-down", 0)
        tipo = h.get("type", "")
        clase = "imitadores" if tipo == "mimickers" or (z_up > 0 and z_dn < 0) else "inversores"
        grupos[clase][nombre].append({"z_sum": h.get("z-sum", z_up - z_dn),
                                      "fdr": max(h.get("fdr-up", 1), h.get("fdr-down", 1)),
                                      "linea": m.get("cell_line", ""), "dosis": m.get("pert_dose", "")})
    OUT.mkdir(parents=True, exist_ok=True)
    import candidatos                                   # reutiliza la anotación ChEMBL (con caché)
    lineas = {c: {k: {x["linea"] for x in v if x["fdr"] < 0.05} for k, v in g.items()} for c, g in grupos.items()}
    resumen = {}
    for clase, g in grupos.items():
        otra = "inversores" if clase == "imitadores" else "imitadores"
        filas = []
        for k, v in g.items():
            sig = [x for x in v if x["fdr"] < 0.05]
            if not sig or k.startswith("brd-"):
                continue
            n_l, n_otra = len(lineas[clase][k]), len(lineas[otra].get(k, ()))
            filas.append({"compuesto": k, "n_firmas_fdr05": len(sig), "n_lineas": n_l, "n_lineas_contrario": n_otra,
                          "consistencia": round((n_l - n_otra) / (n_l + n_otra), 2),
                          "mejor_z": round(max(abs(x["z_sum"]) for x in sig), 2)})
        filas.sort(key=lambda f: (-(f["n_lineas"] - f["n_lineas_contrario"]), -f["mejor_z"]))
        for f in filas[:40]:
            mol = candidatos.chembl_molecula(f["compuesto"].upper())
            f["max_phase"] = (mol or {}).get("max_phase") or ""
            f["mecanismo_chembl"] = "; ".join(candidatos.chembl_mecanismos(mol["molecule_chembl_id"]))[:90] if mol else ""
        with (OUT / f"cmap_{clase}.csv").open("w", newline="") as fh:
            # dict.fromkeys: filas[0] ya trae max_phase y mecanismo_chembl (antes salían duplicadas en la cabecera)
            w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(list(filas[0]) + ["max_phase", "mecanismo_chembl"])),
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(filas)
        resumen[clase] = filas
        print(f"\n{clase}: {len(filas)} compuestos con FDR < 0,05 · top 25 por líneas netas (aprobados marcados ★):")
        for x in filas[:25]:
            ap = "★" if str(x.get("max_phase", "")).startswith("4") else " "
            print(f"  {ap} {x['compuesto']:<24} líneas={x['n_lineas']:<2} contrario={x['n_lineas_contrario']:<2} "
                  f"cons={x['consistencia']:<5} |z|={x['mejor_z']:<6} {x.get('mecanismo_chembl', '')[:60]}")
    enriquecimiento([sim.get(e) for e in arriba if sim.get(e)], "SUBIDOS")
    enriquecimiento([sim.get(e) for e in abajo if sim.get(e)], "BAJADOS")


def enriquecimiento(genes, etiqueta, libreria="MSigDB_Hallmark_2020", silencioso=False):
    """Enrichr (lista de genes de la firma pública; nada del paciente)."""
    import uuid
    caja = uuid.uuid4().hex
    cuerpo = (f"--{caja}\r\nContent-Disposition: form-data; name=\"list\"\r\n\r\n" + "\n".join(genes) +
              f"\r\n--{caja}\r\nContent-Disposition: form-data; name=\"description\"\r\n\r\nbub1b_ko\r\n--{caja}--\r\n")
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=cuerpo.encode(), method="POST",
                                 headers={"Content-Type": f"multipart/form-data; boundary={caja}"})
    lid = json.loads(urllib.request.urlopen(req, timeout=60).read())["userListId"]
    r = json.loads(urllib.request.urlopen(
        f"https://maayanlab.cloud/Enrichr/enrich?userListId={lid}&backgroundType={libreria}", timeout=60).read())
    if not silencioso:
        print(f"\nEnriquecimiento {libreria} · {etiqueta}:")
        for t in r[libreria][:6]:
            print(f"  {t[1]:<40} p_aj={t[6]:.1e}  genes={','.join(t[5][:8])}")
    return [t for t in r[libreria] if t[6] < 0.05]

def consenso(gen):
    """Firma de consenso (rango medio normalizado por gen del ranking) de todas las firmas CRISPR KO de `gen`."""
    firmas = firmas_knockout(gen)
    if not firmas:
        return None, 0, 0
    r = post(f"{DATA}/fetch/rank", {"database": "l1000_xpr", "signatures": [f["id"] for f in firmas], "entities": []})
    media = {}
    for s_ in r["signatures"]:
        for e, x in zip(r["entities"], s_["ranks"]):
            media[e] = media.get(e, 0) + x / r["maxrank"] / len(r["signatures"])
    return media, len(firmas), len({f["meta"].get("cell_line") for f in firmas})


def spearman(a, b):
    comunes = sorted(set(a) & set(b))
    def rangos(v):
        orden = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for k, i in enumerate(orden):
            r[i] = k
        return r
    ra, rb = rangos([a[e] for e in comunes]), rangos([b[e] for e in comunes])
    n = len(comunes)
    return 1 - 6 * sum((x - y) ** 2 for x, y in zip(ra, rb)) / (n * (n * n - 1)), n


def hallmark_top(media, k=3):
    orden = sorted(media, key=media.get)
    abajo, arriba = orden[:N], orden[-N:]
    sim = {x["id"]: x["meta"]["symbol"] for x in post(f"{META}/entities/find", {
        "filter": {"where": {"id": {"inq": abajo + arriba}}, "fields": ["id", "meta.symbol"]}})}
    out = {}
    for nombre, lista in (("up", arriba), ("down", abajo)):
        out[nombre] = enriquecimiento([sim[e] for e in lista if e in sim], nombre, silencioso=True)[:k]
    return out


def comparar(genes):
    base, *_ = consenso(genes[0])
    l = ["# Scalability check — knockout signatures of MVA genes in LINCS L1000", "",
         f"Consensus signature per gene (CRISPR knockout, l1000_xpr; top/bottom {N} genes), Spearman correlation with "
         f"the {genes[0]} consensus over all ranked genes, and top MSigDB Hallmark sets (Enrichr).", "",
         "| Gene | Signatures | Cell lines | Spearman ρ vs " + genes[0] + " | Top Hallmark up | Top Hallmark down |",
         "|---|---|---|---|---|---|"]
    for g in genes:
        m, nf, nl = consenso(g)
        if m is None:
            l.append(f"| {g} | 0 | 0 | — (not profiled in L1000 CRISPR) | — | — |")
            continue
        rho, n = spearman(base, m)
        h = hallmark_top(m)
        fmt = lambda xs: "; ".join(f"{t[1]} ({t[6]:.0e})" for t in xs) or "—"
        l.append(f"| {g} | {nf} | {nl} | {rho:.2f} (n = {n}) | {fmt(h['up'])} | {fmt(h['down'])} |")
        print(l[-1])
    (OUT / "cmap_comparacion.md").write_text("\n".join(l) + "\n")


if __name__ == "__main__":
    if "--comparar" in sys.argv:
        comparar(sys.argv[sys.argv.index("--comparar") + 1].split(","))
    else:
        main()
