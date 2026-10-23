"""balance_alelico.py — control de calidad a escala de cromosoma (solo agregados, nada de variantes):

1. Equilibrio alélico por autosoma. En una célula euploide, una heterocigota tiene una fracción alélica (VAF)
   ~0,5. Si una fracción f de las células de la sangre tuviera una trisomía de ese cromosoma, las VAF se separarían
   en dos modos (~0,5 ± f/(2·(2+f))) y la VARIANZA observada superaría a la esperada por muestreo binomial. Se
   calcula por autosoma: n heterocigotas, mediana de VAF, sobredispersión = var_obs / var_binomial_esperada y
   fracción de heterocigotas con |VAF − 0,5| > 0,2. Una aneuploidía variegada de bajo nivel (cada cromosoma en
   pocas células) NO es detectable así: el resultado esperado es nulo y se informa como tal.
2. Regiones de homocigosidad (bcftools roh): longitud total y nº de ROH > 5 Mb. Muchas ROH largas sugerirían
   consanguinidad (y un modelo homocigoto); pocas apoyan un heterocigoto compuesto.
3. Comprobación local en 15q15 (BUB1B ± 5 Mb): equilibrio alélico y ausencia de pérdida de heterocigosidad.

Solo SNV bialélicas PASS, GT 0/1, DP ≥ 20. Cromosomas sexuales excluidos a propósito (no se informa el sexo).

    .venv/bin/python scripts/modal/balance_alelico.py

Baja a results/t1/balance_alelico.json (≈ 25 filas de estadísticos).
"""
import json
import pathlib
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
CPU, MEM, TIMEOUT_MIN = 2.0, 4, 20
VCF = "WGS_EX2312012_HGWCNDSX7.vcf.gz"
ROOT = pathlib.Path(__file__).resolve().parents[2] if modal.is_local() else pathlib.Path("/")

datos = modal.Volume.from_name("mva-data")
app = modal.App("mva-balance", image=modal.Image.debian_slim(python_version="3.12").apt_install("bcftools"))


@app.function(volumes={"/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=0)
def calcula():
    import statistics
    import subprocess
    t0 = time.time()
    res = {"autosomas": {}, "region_bub1b": {}, "roh": {}}
    for c in [str(i) for i in range(1, 23)]:
        out = subprocess.run(
            ["bcftools", "query", "-r", c, "-i", 'FILTER="PASS" && TYPE="snp" && N_ALT=1 && GT="het" && FMT/DP>=20',
             "-f", "%POS\t[%AD]\n", f"/data/{VCF}"], capture_output=True, text=True, check=True).stdout
        vaf, esp, cerca = [], [], []
        for l in out.splitlines():
            pos, ad = l.split("\t")
            r, a = (int(x) for x in ad.split(",")[:2])
            if r + a < 20:
                continue
            v = a / (r + a)
            vaf.append(v)
            esp.append(0.25 / (r + a))
            if c == "15" and 35_200_000 <= int(pos) <= 45_200_000:
                cerca.append(v)
        var_obs = statistics.pvariance(vaf)
        res["autosomas"][c] = {"n_het": len(vaf), "mediana_vaf": round(statistics.median(vaf), 4),
                               "sobredispersion": round(var_obs / statistics.mean(esp), 3),
                               "frac_desviadas_0.2": round(sum(abs(v - 0.5) > 0.2 for v in vaf) / len(vaf), 4)}
        if cerca:
            res["region_bub1b"] = {"intervalo": "15:35.2-45.2 Mb", "n_het": len(cerca),
                                   "mediana_vaf": round(statistics.median(cerca), 4)}
    # ROH en autosomas (frecuencias alélicas por defecto de bcftools: se estiman de la propia muestra → sólo
    # orientativo; se informa el total y las ROH largas, que son robustas a ese supuesto)
    roh = subprocess.run(["bash", "-c", f"bcftools roh -G30 --AF-dflt 0.4 -r {','.join(str(i) for i in range(1, 23))} "
                          f"-i 'FILTER=\"PASS\" && TYPE=\"snp\"' /data/{VCF} 2>/dev/null | awk '$1==\"RG\"'"],
                         capture_output=True, text=True).stdout
    largos = [int(l.split("\t")[5]) for l in roh.splitlines() if l.strip()]
    res["roh"] = {"n_total": len(largos), "total_mb": round(sum(largos) / 1e6, 1),
                  "n_mayor_5mb": sum(x > 5e6 for x in largos), "mayor_mb": round(max(largos, default=0) / 1e6, 1),
                  "nota": "bcftools roh -G30 --AF-dflt 0.4 (sin panel poblacional: orientativo)"}
    res["secs"] = time.time() - t0
    return res


@app.function(volumes={"/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=0)
def diagnostico(cromos_bins: list, bin_mb: int = 5):
    """¿Artefacto o aneuploidía en mosaico? (1) profundidad mediana por autosoma (todas las SNV PASS, het y hom):
    una ganancia en una fracción f de células sube la profundidad ~(1 + f/2); (2) bins de bin_mb Mb con
    sobredispersión, profundidad y mediana de VAF, para ver si el exceso es uniforme (cromosoma entero) o
    focal (zonas difíciles); (3) histograma de VAF de heterocigotas (20 cajas) por cromosoma pedido."""
    import statistics
    import subprocess
    t0 = time.time()
    res = {"dp_mediana": {}, "bins": {}, "hist": {}}
    for c in [str(i) for i in range(1, 23)]:
        out = subprocess.run(
            ["bcftools", "query", "-r", c, "-i", 'FILTER="PASS" && TYPE="snp" && N_ALT=1 && (GT="het" || GT="AA")',
             "-f", "%POS\t[%GT]\t[%AD]\n", f"/data/{VCF}"], capture_output=True, text=True, check=True).stdout
        dps, bins, hist = [], {}, [0] * 20
        for l in out.splitlines():
            pos, gt, ad = l.split("\t")
            try:
                r, a = (int(x) for x in ad.split(",")[:2])
            except ValueError:
                continue
            dps.append(r + a)
            if c in cromos_bins and gt in ("0/1", "0|1", "1|0") and r + a >= 20:
                v = a / (r + a)
                b = int(pos) // (bin_mb * 1_000_000)
                bins.setdefault(b, []).append((v, r + a))
                hist[min(int(v * 20), 19)] += 1
        res["dp_mediana"][c] = statistics.median(dps)
        if c in cromos_bins:
            res["hist"][c] = hist
            res["bins"][c] = {}
            for b, xs in sorted(bins.items()):
                if len(xs) < 200:
                    continue
                vs = [x[0] for x in xs]
                esp = statistics.mean(0.25 / x[1] for x in xs)
                res["bins"][c][f"{b * bin_mb}-{(b + 1) * bin_mb}"] = {
                    "n": len(xs), "sobredisp": round(statistics.pvariance(vs) / esp, 2),
                    "med_vaf": round(statistics.median(vs), 3), "med_dp": statistics.median(x[1] for x in xs)}
    res["secs"] = time.time() - t0
    return res


@app.function(volumes={"/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60, retries=0)
def homocigosidad(min_mb: float = 2.0, min_hom_kb: float = 0.05):
    """Regiones de homocigosidad (ROH) robustas: tramos entre heterocigotas consecutivas (PASS, DP ≥ 20) de más de
    min_mb Mb que SÍ contienen homocigotas alternativas llamadas (densidad > min_hom_kb por kb), para no confundir
    huecos del ensamblaje o centrómeros (sin llamadas) con homocigosidad. Además el stderr de bcftools roh, para
    entender por qué devolvió 0."""
    import subprocess
    t0 = time.time()
    rohs = []
    for c in [str(i) for i in range(1, 23)]:
        out = subprocess.run(
            ["bcftools", "query", "-r", c, "-i", 'FILTER="PASS" && TYPE="snp" && N_ALT=1 && FMT/DP>=20 && (GT="het" || GT="AA")',
             "-f", "%POS\t[%GT]\n", f"/data/{VCF}"], capture_output=True, text=True, check=True).stdout
        het_prev, homs = None, 0
        for l in out.splitlines():
            pos, gt = l.split("\t")
            pos = int(pos)
            if gt in ("1/1", "1|1"):
                homs += 1
                continue
            if het_prev is not None:
                largo = pos - het_prev
                if largo > min_mb * 1e6 and homs / (largo / 1e3) > min_hom_kb:
                    rohs.append((c, het_prev, pos, round(largo / 1e6, 2), homs))
            het_prev, homs = pos, 0
    err = subprocess.run(["bash", "-c", f"bcftools roh -G30 --AF-dflt 0.4 -r 22 /data/{VCF} 2>&1 >/dev/null | head -5"],
                         capture_output=True, text=True).stdout
    return {"rohs": rohs, "bcftools_roh_stderr": err, "secs": time.time() - t0}


if __name__ == "__main__" and "--roh" in __import__("sys").argv:
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, 0.1, label="homocigosidad")
    with app.run():
        r = homocigosidad.remote()
    largos = [x[3] for x in r["rohs"]]
    res = {"n_roh_mayor_2mb": len(largos), "total_mb": round(sum(largos), 1),
           "n_mayor_5mb": sum(x > 5 for x in largos), "mayor_mb": max(largos, default=0),
           "por_cromosoma": {c: sum(1 for x in r["rohs"] if x[0] == c) for c in {x[0] for x in r["rohs"]}}}
    j = ROOT / "results/t1/balance_alelico.json"
    d = json.loads(j.read_text()); d["roh_robusto"] = res; j.write_text(json.dumps(d, indent=1))
    print("ROH robustas (> 2 Mb, con homocigotas llamadas):", res)
    print("stderr de bcftools roh:", r["bcftools_roh_stderr"][:400])
    comun.apunta("homocigosidad", CPU, MEM, r["secs"], r["secs"], comun.coste(CPU, MEM, r["secs"]), True)
    raise SystemExit


if __name__ == "__main__" and "--diagnostico" in __import__("sys").argv:
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, 0.1, label="balance alélico · diagnóstico")
    comun.vigilante(TIMEOUT_MIN * 60 + 300, "diagnóstico")
    t0 = time.time()
    with app.run():
        r = diagnostico.remote(["1", "5", "9", "17", "20", "21", "22"])
    (ROOT / "results/t1/balance_alelico_diagnostico.json").write_text(json.dumps(r, indent=1))
    g = __import__("statistics").median(r["dp_mediana"].values())
    print("profundidad mediana por autosoma (relativa a la mediana genómica):")
    print("  " + "  ".join(f"{c}:{v / g:.3f}" for c, v in r["dp_mediana"].items()))
    for c, bins in r["bins"].items():
        print(f"\nchr{c}  bins de 5 Mb (n · sobredisp · med_VAF · med_DP):")
        print("  " + " | ".join(f"{k}: {v['n']}·{v['sobredisp']}·{v['med_vaf']}·{v['med_dp']:.0f}" for k, v in bins.items()))
    print("\nhistograma de VAF (20 cajas de 0,05):")
    for c, h in r["hist"].items():
        t = sum(h)
        print(f"  chr{c:>2}: " + " ".join(f"{100 * x / t:4.1f}" for x in h))
    comun.apunta("balance alélico diagnóstico", CPU, MEM, r["secs"], time.time() - t0, comun.coste(CPU, MEM, r["secs"]), True)
    raise SystemExit


if __name__ == "__main__":
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, 0.1, label="balance alélico")
    comun.vigilante(TIMEOUT_MIN * 60 + 300, "balance alélico")
    t0 = time.time()
    with app.run():
        r = calcula.remote()
    (ROOT / "results/t1/balance_alelico.json").write_text(json.dumps(r, indent=1))
    so = [v["sobredispersion"] for v in r["autosomas"].values()]
    print(f"{'crom':>4} {'n_het':>7} {'med_VAF':>8} {'sobredisp':>9} {'|d|>0.2':>8}")
    for c, v in r["autosomas"].items():
        print(f"{c:>4} {v['n_het']:>7} {v['mediana_vaf']:>8} {v['sobredispersion']:>9} {v['frac_desviadas_0.2']:>8}")
    print("15q15 (BUB1B ± 5 Mb):", r["region_bub1b"])
    print("ROH:", r["roh"])
    print(f"sobredispersión: mín {min(so)}, máx {max(so)}, mediana {sorted(so)[len(so) // 2]}")
    usd = comun.coste(CPU, MEM, r["secs"])
    comun.apunta("balance alélico", CPU, MEM, r["secs"], time.time() - t0, usd, True)
    print(f"✓ {r['secs'] / 60:.1f} min · ~{usd:.3f} $")
