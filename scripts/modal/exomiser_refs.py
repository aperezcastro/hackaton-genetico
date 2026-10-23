"""exomiser_refs.py — descarga UNA VEZ Exomiser 15.1.0 y sus datos 2512 (hg38 + fenotipo) al volumen `mva-exomiser`.
Solo datos públicos.

    .venv/bin/python scripts/modal/exomiser_refs.py [--max-usd 1] [--timeout 150]
"""
import sys
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
arg = lambda n, d: type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
TIMEOUT_MIN = arg("--timeout", 150)
MAX_USD = arg("--max-usd", 1.0)
CPU, MEM = 4.0, 8
vol = modal.Volume.from_name("mva-exomiser", create_if_missing=True)
image = modal.Image.debian_slim(python_version="3.12").apt_install("curl", "unzip")
app = modal.App("mva-exomiser-refs", image=image)

PASOS = r"""
set -euo pipefail
cd /ex
if [ ! -d exomiser-cli-15.1.0 ]; then
  curl -sfL -o cli.zip https://github.com/exomiser/Exomiser/releases/download/15.1.0/exomiser-cli-15.1.0-distribution.zip
  unzip -q cli.zip && rm cli.zip
fi
mkdir -p data
for f in 2512_hg38 2512_phenotype; do
  if [ ! -f data/$f.ok ]; then
    echo "[exomiser] $f…"
    curl -sfL -o data/$f.zip https://data.monarchinitiative.org/exomiser/data/$f.zip
    unzip -q -o data/$f.zip -d data && rm data/$f.zip && touch data/$f.ok
  fi
done
ls -la data; du -sh /ex
"""


@app.function(volumes={"/ex": vol}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60)
def descarga():
    import subprocess
    t0 = time.time()
    r = subprocess.run(["bash", "-c", PASOS], capture_output=True, text=True)
    vol.commit()
    return {"ok": r.returncode == 0, "log": (r.stdout + r.stderr)[-3000:], "secs": time.time() - t0}


if __name__ == "__main__":
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, MAX_USD, label="exomiser-refs")
    comun.vigilante(TIMEOUT_MIN * 60 + 300, "exomiser-refs")
    t0 = time.time()
    with modal.enable_output(), app.run():
        r = descarga.remote()
    print(r["log"])
    usd = comun.coste(CPU, MEM, r["secs"])
    print(f"{'✓' if r['ok'] else '✗'} exomiser-refs · {r['secs'] / 60:.1f} min · ~{usd:.3f} $")
    comun.apunta("exomiser-refs", CPU, MEM, r["secs"], time.time() - t0, usd, r["ok"])
