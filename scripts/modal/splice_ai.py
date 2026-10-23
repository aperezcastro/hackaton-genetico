"""splice_ai.py (no llamarlo spliceai.py: taparía al paquete en el contenedor) — puntúa con SpliceAI (Illumina, pesos abiertos) las variantes raras del panel MVA/SAC en Modal (CPU).

Entrada: results/t1/panel_mva.tsv filtrado a PASS y AF máx. gnomAD < 0,1 % (≈ 40 variantes: una lista corta, no una
tabla a escala genómica). Se escribe un VCF mínimo, se sube al contenedor junto con nada más y se puntúa contra el
fasta Ensembl GRCh38 del volumen público `mva-refs`, con la anotación GENCODE grch38 incluida en SpliceAI.

Dos distancias: -D 50 (por defecto, sitios cercanos) y -D 2000 (activación de sitios crípticos lejanos, útil para la
variante intrónica profunda de BUB1B). Delta ≥ 0,2 = posible efecto, ≥ 0,5 = probable, ≥ 0,8 = alta confianza
(Jaganathan et al. 2019, Cell).

    .venv/bin/python scripts/modal/splice_ai.py [--max-usd 0.5] [--timeout 30]

Baja a results/t1/: spliceai.tsv (una fila por variante y gen, deltas DS_AG/DS_AL/DS_DG/DS_DL y posiciones).
"""
import csv
import io
import pathlib
import sys
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
arg = lambda n, d: type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
TIMEOUT_MIN = arg("--timeout", 30)
MAX_USD = arg("--max-usd", 0.5)
CPU, MEM = 4.0, 8
ROOT = pathlib.Path(__file__).resolve().parents[2] if modal.is_local() else pathlib.Path("/")
FASTA = "/refs/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

refs = modal.Volume.from_name("mva-refs")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("tensorflow-cpu==2.15.1", "spliceai==1.3.1", "pyfaidx", "pysam", "biopython"))
app = modal.App("mva-spliceai", image=image)


@app.function(volumes={"/refs": refs}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=0)
def puntua(vcf: str):
    import subprocess
    t0 = time.time()
    open("/tmp/in.vcf", "w").write(vcf)
    out = {}
    for d in (50, 2000):
        r = subprocess.run(["spliceai", "-I", "/tmp/in.vcf", "-O", f"/tmp/out{d}.vcf", "-R", FASTA, "-A", "grch38",
                            "-D", str(d)], capture_output=True, text=True)
        out[d] = (open(f"/tmp/out{d}.vcf").read() if r.returncode == 0 else "", r.stderr[-1500:])
    return {"out": out, "secs": time.time() - t0}


def vcf_minimo(rows):
    h = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c}>" for c in sorted({r["chrom"] for r in rows})]
    h.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    body = sorted(rows, key=lambda r: (r["chrom"], int(r["pos"])))
    return "\n".join(h + [f"{r['chrom']}\t{r['pos']}\t{r['gene'] or '.'}\t{r['ref']}\t{r['alt']}\t.\tPASS\t."
                          for r in body]) + "\n"


def parsea(texto, d):
    filas = []
    for l in texto.splitlines():
        if l.startswith("#"):
            continue
        c = l.split("\t")
        info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
        for s in info.get("SpliceAI", "").split(","):
            p = s.split("|")
            if len(p) < 10:
                continue
            ds = [float(x) for x in p[2:6]]
            filas.append({"chrom": c[0], "pos": c[1], "ref": c[3], "alt": c[4], "gen": p[1], "D": d,
                          "DS_AG": ds[0], "DS_AL": ds[1], "DS_DG": ds[2], "DS_DL": ds[3], "DS_max": max(ds),
                          "DP_AG": p[6], "DP_AL": p[7], "DP_DG": p[8], "DP_DL": p[9]})
    return filas


if __name__ == "__main__":
    rows = list(csv.DictReader(open(ROOT / "results/t1/panel_mva.tsv"), delimiter="\t"))
    sel = [r for r in rows if r["filter"] == "PASS" and (r["af_max"] == "" or float(r["af_max"]) < 0.001)]
    print(f"{len(sel)} variantes raras PASS del panel → SpliceAI")
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, MAX_USD, label="spliceai")
    comun.vigilante(TIMEOUT_MIN * 60 + 300, "spliceai")
    t0 = time.time()
    with modal.enable_output(), app.run():
        r = puntua.remote(vcf_minimo(sel))
    filas = []
    for d, (texto, err) in r["out"].items():
        if not texto:
            print(f"✗ D={d}:\n{err}")
        filas += parsea(texto, d)
    dst = ROOT / "results/t1/spliceai.tsv"
    if filas:
        with dst.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0]), delimiter="\t")
            w.writeheader()
            w.writerows(filas)
    usd = comun.coste(CPU, MEM, r["secs"])
    print(f"✓ {len(filas)} puntuaciones → {dst.relative_to(ROOT)} · {r['secs'] / 60:.1f} min · ~{usd:.3f} $")
    comun.apunta("spliceai panel", CPU, MEM, r["secs"], time.time() - t0, usd, bool(filas))
