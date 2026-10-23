"""refs.py — descarga UNA VEZ las referencias públicas del pipeline al volumen de Modal `mva-refs`.
Ningún dato del paciente pasa por aquí.

    .venv/bin/python scripts/modal/refs.py [--max-usd 1.5] [--timeout 180]

Contenido del volumen (/refs):
  vep/homo_sapiens/116_GRCh38/      caché indexada de VEP 116 (transcritos Ensembl/MANE, gnomAD v4.1, ClinVar, dbSNP)
  fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz(.fai,.gzi)   referencia Ensembl (contigs sin «chr», como el VCF)
  alphamissense/AlphaMissense_hg38.tsv.gz(.tbi)              plugin AlphaMissense
  clinvar/clinvar.vcf.gz(.tbi)                              ClinVar GRCh38 (NCBI), para estrellas de revisión

Robusto frente a cortes (1 oct 2026: Modal perdió el contenedor a mitad de una descarga en streaming de 27 GB):
  · primero lo pequeño; cada fichero se marca con <fichero>.ok al terminar y no se repite;
  · el tar de VEP se baja como FICHERO al volumen con `aria2c -c -x 8` (reanuda, 8 conexiones) y se hace commit del volumen cada
    5 min, así un contenedor nuevo continúa donde quedó; se descomprime después y se borra el tar;
  · la app se lanza con detach=True: no muere si se corta la conexión del portátil.
"""
import os
import subprocess
import sys
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
arg = lambda n, d: type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
TIMEOUT_MIN = arg("--timeout", 180)
MAX_USD = arg("--max-usd", 1.5)
CPU, MEM = 2.0, 4

vol = modal.Volume.from_name("mva-refs", create_if_missing=True)
image = (modal.Image.from_registry("ensemblorg/ensembl-vep:release_116.2", add_python="3.12")
         .apt_install("wget", "curl", "aria2"))
app = modal.App("mva-refs", image=image)

E = "https://ftp.ensembl.org/pub/release-116"
PEQUENOS = [
    ("fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz", f"{E}/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"),
    ("fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai", f"{E}/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai"),
    ("fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi", f"{E}/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi"),
    ("alphamissense/AlphaMissense_hg38.tsv.gz", "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"),
    ("clinvar/clinvar.vcf.gz", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"),
    ("clinvar/clinvar.vcf.gz.tbi", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"),
]
TAR = ("vep/homo_sapiens_vep_116_GRCh38.tar.gz", f"{E}/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz",
       27644657162)


def sh(c):
    r = subprocess.run(["bash", "-c", c], capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr)[-1500:]


@app.function(volumes={"/refs": vol}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=2)
def descarga():
    t0, log = time.time(), []
    os.chdir("/refs")
    for d in ("vep", "fasta", "alphamissense", "clinvar"):
        os.makedirs(d, exist_ok=True)
    for rel, url in PEQUENOS:
        if not os.path.exists(rel + ".ok"):
            rc, out = sh(f"curl -sfL --retry 5 -o {rel} {url} && touch {rel}.ok")
            log.append(f"{rel}: rc={rc} {out[-200:]}")
    if os.path.exists("alphamissense/AlphaMissense_hg38.tsv.gz.ok") and not os.path.exists("alphamissense/AlphaMissense_hg38.tsv.gz.tbi"):
        rc, out = sh("tabix -s 1 -b 2 -e 2 -f -S 1 alphamissense/AlphaMissense_hg38.tsv.gz")
        log.append(f"tabix AM: rc={rc} {out[-200:]}")
    vol.commit()
    rel, url, size = TAR
    if not os.path.exists("vep/homo_sapiens/116_GRCh38/.ok"):
        if not os.path.exists(rel + ".ok"):
            # aria2c: 8 conexiones en paralelo (Ensembl limita ~1,6 MB/s por conexión) y -c reanuda el fichero parcial
            p = subprocess.Popen(["bash", "-c", f"until aria2c -q -c -x 8 -s 8 -k 64M --file-allocation=none "
                                  f"-d {os.path.dirname(rel)} -o {os.path.basename(rel)} {url}; do sleep 10; done"])
            while p.poll() is None:            # commit periódico: un contenedor nuevo reanuda desde aquí
                time.sleep(300)
                vol.commit()
                log.append(f"[{(time.time() - t0) / 60:.0f} min] tar {os.path.getsize(rel) / 1e9:.1f}/{size / 1e9:.1f} GB")
            if os.path.getsize(rel) == size:
                open(rel + ".ok", "w").close()
            vol.commit()
        if os.path.exists(rel + ".ok"):
            sh("rm -rf vep/homo_sapiens/116_GRCh38")
            rc, out = sh(f"tar -xzf {rel} -C vep && touch vep/homo_sapiens/116_GRCh38/.ok && rm -f {rel} {rel}.ok")
            log.append(f"tar -x: rc={rc} {out[-300:]}")
            vol.commit()
    rc, out = sh("du -sh /refs/* ; ls /refs/vep/homo_sapiens/116_GRCh38 | head -40; ls /refs/*/")
    log.append(out)
    ok = os.path.exists("vep/homo_sapiens/116_GRCh38/.ok") and all(os.path.exists(r + ".ok") for r, _ in PEQUENOS)
    return {"ok": ok, "log": "\n".join(log)[-5000:], "secs": time.time() - t0}


if __name__ == "__main__":
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60 * 3, MAX_USD, label="refs (con 2 reintentos)")
    comun.vigilante(TIMEOUT_MIN * 60 * 3 + 600, "refs")
    t0 = time.time()
    with modal.enable_output(), app.run(detach=True):
        r = descarga.remote()
    print(r["log"])
    usd = comun.coste(CPU, MEM, r["secs"])
    print(f"{'✓' if r['ok'] else '✗'} refs · {r['secs'] / 60:.1f} min · ~{usd:.3f} $")
    comun.apunta("refs", CPU, MEM, r["secs"], time.time() - t0, usd, r["ok"])
