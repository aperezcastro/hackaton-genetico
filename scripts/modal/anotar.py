"""anotar.py — paso 2 del pipeline T1: normaliza y anota el VCF del caso con VEP 116 en Modal, un contenedor por
cromosoma, y filtra allí mismo. Solo bajan a local tablas cortas (results/) y conteos del embudo.

    .venv/bin/python scripts/modal/anotar.py [--max-usd 2] [--timeout 60] [--solo 22] [--subir]

  --subir   sube (una vez) el VCF + índice de data/raw/ al volumen privado `mva-data` de Modal.
  --solo C  anota solo el cromosoma C (prueba barata).

Volúmenes: `mva-refs` (referencias públicas, ver refs.py) y `mva-data` (datos del caso; Modal actúa como
procesador: no entrena con ellos; se borra al terminar el hackatón, ver docs/llm_use.md).

Por cromosoma: bcftools norm (-m -any, left-align contra la referencia Ensembl, contigs sin «chr») → VEP 116
offline (--everything: consecuencias por transcrito, MANE, HGVS, gnomAD v4.1 exomas/genomas, ClinVar de la caché,
SIFT/PolyPhen, regulación) + plugin AlphaMissense + ClinVar NCBI (estrellas) → mva-data/anot/<chrom>.vcf.gz.
Después filtra.py (en Modal) produce las tablas de results/.
"""
import json
import pathlib
import sys
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
arg = lambda n, d: type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
TIMEOUT_MIN = arg("--timeout", 60)
MAX_USD = arg("--max-usd", 2.0)
CPU, MEM = 4.0, 8
ROOT = pathlib.Path(__file__).resolve().parents[2] if modal.is_local() else pathlib.Path("/")
VCF = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
CROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

refs = modal.Volume.from_name("mva-refs", create_if_missing=True)
datos = modal.Volume.from_name("mva-data", create_if_missing=True)
image = (modal.Image.from_registry("ensemblorg/ensembl-vep:release_116.2", add_python="3.12")
         .apt_install("bcftools", "tabix")
         .add_local_file(str(pathlib.Path(__file__).parent / "filtra.py"), "/work/filtra.py"))
app = modal.App("mva-anotar", image=image)

FASTA = "/refs/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"


@app.function(volumes={"/refs": refs, "/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60)
def anota(crom: str):
    import os
    import subprocess
    t0 = time.time()
    os.makedirs("/data/anot", exist_ok=True)
    out = f"/data/anot/{crom}.vcf.gz"
    vcf_crom = "M" if crom == "MT" else crom          # el VCF del caso llama «M» al mitocondrial
    plug = os.environ.get("VEP_DIR_PLUGINS", "/plugins")
    cmd = f"""
set -euo pipefail
cd /tmp
bcftools view -r {vcf_crom} /data/{VCF} -Oz -o in.vcf.gz
{"bcftools annotate --rename-chrs <(echo 'M MT') in.vcf.gz -Oz -o in2.vcf.gz && mv in2.vcf.gz in.vcf.gz" if crom == "MT" else ""}
bcftools norm -m -any -f {FASTA} --check-ref w in.vcf.gz -Oz -o norm.vcf.gz 2> norm.log || (cat norm.log; exit 1)
N=$(bcftools view -H norm.vcf.gz | wc -l)
echo "variantes_norm $N"
vep --offline --cache --dir_cache /refs/vep --assembly GRCh38 --cache_version 116 --fasta {FASTA} \
    --input_file norm.vcf.gz --format vcf --vcf --compress_output bgzip --output_file vep.vcf.gz --force_overwrite \
    --everything --fork {int(CPU)} --buffer_size 20000 --no_stats \
    --dir_plugins {plug} --plugin AlphaMissense,file=/refs/alphamissense/AlphaMissense_hg38.tsv.gz \
    --custom file=/refs/clinvar/clinvar.vcf.gz,short_name=ClinVarNCBI,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
tabix -p vcf vep.vcf.gz
cp vep.vcf.gz {out}; cp vep.vcf.gz.tbi {out}.tbi
"""
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    datos.commit()
    log = (r.stdout + r.stderr)[-2500:]
    n = next((int(l.split()[1]) for l in r.stdout.splitlines() if l.startswith("variantes_norm")), None)
    return {"ok": r.returncode == 0, "crom": crom, "n": n, "secs": time.time() - t0, "log": log}


@app.function(volumes={"/refs": refs, "/data": datos}, cpu=4.0, memory=16 * 1024, timeout=40 * 60)
def filtra(croms: list):
    import subprocess
    t0 = time.time()
    r = subprocess.run(["python", "/work/filtra.py", "/data/anot", "/data/filtrado", ",".join(croms)],
                       capture_output=True, text=True)
    datos.commit()
    res = {}
    import os
    for f in sorted(os.listdir("/data/filtrado")) if os.path.isdir("/data/filtrado") else []:
        if f.endswith((".tsv", ".json")):
            res[f] = open(f"/data/filtrado/{f}", "rb").read()
    return {"ok": r.returncode == 0, "log": (r.stdout + r.stderr)[-4000:], "files": res, "secs": time.time() - t0}


@app.function(volumes={"/data": datos}, timeout=600)
def prueba_plugins():
    import os
    import subprocess
    plug = os.environ.get("VEP_DIR_PLUGINS", "/plugins")
    return subprocess.run(["bash", "-c", f"echo {plug}; ls {plug} | grep -i -E 'alpha|splice|revel|cadd' ; ls /data"],
                          capture_output=True, text=True).stdout


def main():
    if "--subir" in sys.argv:
        with datos.batch_upload(force=True) as b:
            for f in (VCF, VCF + ".tbi"):
                b.put_file(str(ROOT / "data" / "raw" / f), f"/{f}")
        print("✓ VCF subido al volumen privado mva-data")
    croms = [arg("--solo", "")] if "--solo" in sys.argv else CROMS
    print(prueba_plugins.remote())
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, MAX_USD, n=len(croms), label=f"anotar ×{len(croms)}")
    comun.vigilante(TIMEOUT_MIN * 60 + 45 * 60, "anotar")
    t0, busy, ok, tot = time.time(), 0.0, True, 0
    for r in anota.map(croms, order_outputs=False, return_exceptions=True):
        if isinstance(r, Exception):
            print("✗", r); ok = False; continue
        busy += r["secs"]
        print(f"[{'✓' if r['ok'] else '✗'}] chr{r['crom']}: {r['n']} variantes · {r['secs'] / 60:.1f} min", flush=True)
        if not r["ok"]:
            print(r["log"]); ok = False
        else:
            tot += r["n"] or 0
    usd = comun.coste(CPU, MEM, busy)
    print(f"anotación: {tot} variantes normalizadas · {busy / 3600:.2f} h·contenedor · ~{usd:.3f} $")
    comun.apunta(f"anotar ×{len(croms)}", CPU, MEM, busy, time.time() - t0, usd, ok)
    if ok:
        r = filtra.remote(croms)
        print(r["log"])
        dst = ROOT / "results" / "t1"
        dst.mkdir(parents=True, exist_ok=True)
        for f, b in r["files"].items():
            (dst / f).write_bytes(b)
        u2 = comun.coste(4, 16, r["secs"])
        comun.apunta("filtra", 4, 16, r["secs"], r["secs"], u2, r["ok"])
        print(f"✓ {len(r['files'])} ficheros en results/t1/ · filtrado ~{u2:.3f} $")


if __name__ == "__main__":
    with modal.enable_output(), app.run():
        main()
