"""estado_vol.py — tamaño de lo descargado en los volúmenes de referencias (diagnóstico rápido)."""
import os, pathlib
_env = pathlib.Path(__file__).resolve().parents[2] / ".env"   # opcional: si no existe, Modal usa su propia config
for l in (open(_env) if _env.exists() else []):
    k, _, v = l.strip().partition("=")
    if k == "modal-token-id": os.environ["MODAL_TOKEN_ID"] = v
    if k == "modal-token-secret": os.environ["MODAL_TOKEN_SECRET"] = v
import modal
r = modal.Volume.from_name("mva-refs")
print("vep:", len(r.listdir("vep/homo_sapiens/116_GRCh38")), "entradas;",
      {d: len(r.listdir(d)) for d in ["fasta", "alphamissense", "clinvar"]})
e = modal.Volume.from_name("mva-exomiser")
print("exomiser:", [(x.path, round(x.size / 1e6)) for x in e.listdir("data")])
