"""exomiser_top.py — de las salidas completas de Exomiser (data/derivados/exomiser_{A,B}.genes.tsv, a escala
genómica) saca solo los 12 primeros puestos de cada configuración a results/t1/exomiser_top.json, que es lo que lee
la figura t1_exomiser (que a su vez solo nombra los 3 primeros). No imprime nombres de genes, y los que figuran
en results/sensibles.txt (fuera de git) se tapan con un marcador antes de escribir el JSON.

    .venv/bin/python scripts/exomiser_top.py [--comprobar]

  --comprobar  no escribe: dice si el resultado coincide con el results/t1/exomiser_top.json existente
"""
import csv
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "derivados"
DST = ROOT / "results" / "t1" / "exomiser_top.json"
N = 12
MASCARA = "(not named: privacy)"


def col(fila, *nombres):
    for n in nombres:
        if n in fila and fila[n] not in ("", None):
            return fila[n]
    raise KeyError(nombres)


def top(tsv):
    with tsv.open() as f:
        filas = list(csv.DictReader(f, delimiter="\t"))
    salida = []
    for fila in filas[:N]:
        fila = {k.lstrip("#").strip(): v for k, v in fila.items()}
        salida.append({"rank": int(col(fila, "RANK")), "gene": col(fila, "GENE_SYMBOL"), "moi": col(fila, "MOI"),
                       "combined": round(float(col(fila, "EXOMISER_GENE_COMBINED_SCORE")), 4),
                       "pheno": round(float(col(fila, "EXOMISER_GENE_PHENO_SCORE")), 4),
                       "variant": round(float(col(fila, "EXOMISER_GENE_VARIANT_SCORE")), 4)})
    return salida


def main():
    lista = ROOT / "results" / "sensibles.txt"
    sensibles = set(lista.read_text().split()) if lista.exists() else set()
    res = {c: top(SRC / f"exomiser_{c}.genes.tsv") for c in ("A", "B")}
    # los genes de results/sensibles.txt no llegan ni al JSON: se tapan con un marcador
    for filas in res.values():
        for t in filas:
            if t["gene"] in sensibles:
                t["gene"] = MASCARA
    if "--comprobar" in sys.argv:
        viejo = json.loads(DST.read_text())
        difs = []
        for c in res:
            if len(res[c]) != len(viejo[c]):
                difs.append(f"{c}: {len(res[c])} filas frente a {len(viejo[c])}")
            genes_viejos = [t["gene"] for t in viejo[c]]
            for i, (n, v) in enumerate(zip(res[c], viejo[c]), 1):
                # solo posiciones, campos y diferencias numéricas: nunca el nombre del gen
                campos = [f"{k} ({n[k] - v[k]:+.4f})" for k in ("combined", "pheno", "variant") if abs(n[k] - v[k]) > 1e-4]
                campos += [k for k in ("rank", "moi") if n[k] != v[k]]
                tapado_igual = n["gene"] == MASCARA and not re.fullmatch(r"[A-Z0-9-]{2,15}", v["gene"])
                if n["gene"] != v["gene"] and not tapado_igual:  # un marcador antiguo cuenta como tapado
                    otro = genes_viejos.index(n["gene"]) + 1 if n["gene"] in genes_viejos else None
                    campos.append(f"gene (el nuevo {'estaba en la fila ' + str(otro) if otro else 'no estaba en el top 12'}"
                                  f"{', tapado' if n['gene'] == MASCARA else ''}; el viejo "
                                  f"{'en sensibles.txt' if v['gene'] in sensibles else 'no está en sensibles.txt'})")
                if campos:
                    difs.append(f"{c} fila {i}: {', '.join(campos)}")
        print("✓ coincide con el JSON existente" if not difs else "✗ difiere del JSON existente:\n  " + "\n  ".join(difs))
        return
    DST.parent.mkdir(parents=True, exist_ok=True)
    DST.write_text(json.dumps(res, indent=1))
    print(f"✓ top {N} de A y B → {DST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
