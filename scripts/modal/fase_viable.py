"""fase_viable.py — ¿se podría fasar el par de BUB1B realineando los FASTQ? Comprueba en el VCF (en Modal) si entre
las dos variantes hay una cadena de heterocigotas con huecos menores que el tamaño de inserto (~500 pb), que es lo
que necesita WhatsHap para enlazarlas con lecturas cortas. Solo devuelve recuentos y distancias, no genotipos.

    .venv/bin/python scripts/modal/fase_viable.py
"""
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
CPU, MEM, TIMEOUT_MIN = 1.0, 2, 10
VCF = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
INI, FIN = 40209701, 40220612          # alelo 1 (stop) y alelo 2 (missense), GRCh38

datos = modal.Volume.from_name("mva-data")
app = modal.App("mva-fase-viable", image=modal.Image.debian_slim().apt_install("bcftools"))


@app.function(volumes={"/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=0)
def mira():
    import subprocess
    t0 = time.time()
    out = subprocess.run(["bcftools", "query", "-r", f"15:{INI}-{FIN}", "-i", 'GT="het" && FILTER="PASS"',
                          "-f", "%POS\t[%PGT]\t[%PID]\n", f"/data/{VCF}"], capture_output=True, text=True, check=True)
    pos, pids = [], set()
    for l in out.stdout.splitlines():
        p, pgt, pid = l.split("\t")
        pos.append(int(p))
        if pid not in (".", ""):
            pids.add(pid)
    pos = sorted(set(pos) | {INI, FIN})
    huecos = [b - a for a, b in zip(pos, pos[1:])]
    return {"n_het": len(pos), "huecos": huecos, "bloques_pid": len(pids), "secs": time.time() - t0}


if __name__ == "__main__":
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, 0.05, label="fase viable")
    with app.run():
        r = mira.remote()
    h = r["huecos"]
    print(f"heterocigotas PASS en el intervalo (incluidas las 2): {r['n_het']} · bloques PID: {r['bloques_pid']}")
    print(f"hueco máximo entre heterocigotas consecutivas: {max(h)} pb · huecos > 500 pb: {sum(x > 500 for x in h)}")
    print("huecos ordenados:", sorted(h, reverse=True)[:10])
    comun.apunta("fase viable", CPU, MEM, r["secs"], r["secs"], comun.coste(CPU, MEM, r["secs"]), True)
