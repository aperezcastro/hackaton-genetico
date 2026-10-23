"""figuras_t2.py — figuras del informe del Track 2 (PDF vectorial en reports/latex/fig/).

Solo datos públicos a nivel gen/fármaco: firmas LINCS (caché de scripts/t2/cmap.py), ranking de
scripts/t2/candidatos.py. Ningún dato del paciente.

    .venv/bin/python scripts/figuras/figuras_t2.py
"""
import csv
import math
import random
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parents[1] / "t2"))
from estilo import *  # noqa: E402,F403

T2 = ROOT / "reports" / "t2"
aplica()
LINEA = {"H1": ("H1 · stabilise BUBR1", AZUL), "H2": ("H2 · read-through", GRIS), "H3": ("H3 · mTORC1", NARANJA),
         "H4": ("H4 · senolysis", AQUA), "H5": ("H5 · aneuploid stress", GRIS), "H6": ("H6 · interferon / JAK", VIOLETA),
         "H7": ("H7 · proteostasis (declined)", GRIS)}


def firma():
    import cmap
    m, nf, nl = cmap.consenso("BUB1B")
    h = cmap.hallmark_top(m, k=6)
    lim = lambda n: n.replace("Signaling via ", "via ").replace(" Response", "").replace("Targets V1", "targets")
    up = [(lim(t[1]), -math.log10(t[6])) for t in h["up"]]
    dn = [(lim(t[1]), -math.log10(t[6])) for t in h["down"]]
    # conectividad por clase
    lin = {}
    for clase, f in (("R", "cmap_inversores.csv"), ("M", "cmap_imitadores.csv")):
        for r in csv.DictReader((T2 / f).open()):
            lin.setdefault(r["compuesto"], {})[clase] = int(r["n_lineas"])
    grupos = [("DNA-damaging chemotherapy", ROJO, ["idarubicin", "mitoxantrone", "etoposide", "irinotecan", "doxorubicin",
                                                  "gemcitabine", "cytarabine"]),
              ("mTOR inhibitors", NARANJA, ["sirolimus", "everolimus", "temsirolimus"]),
              ("JAK inhibitors", VIOLETA, [("tg-101348", "fedratinib"), "baricitinib", "tofacitinib", "ruxolitinib"]),
              ("our other candidates", AZUL, ["nicotinamide", "dasatinib"])]
    fig, (a, b) = plt.subplots(1, 2, figsize=(ANCHO, 3.1), gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.95))
    filas = [(n, v, ROJO) for n, v in up[::-1]] + [(n, -v, AZUL) for n, v in dn[::-1]]
    ys = np.arange(len(filas))
    for y, (n, v, c) in zip(ys, filas):
        a.barh(y, v, color=c, height=0.62, edgecolor="white")
        a.text(v + (0.3 if v > 0 else -0.3), y, n, va="center", ha="left" if v > 0 else "right", fontsize=6.9,
               color=TEXTO)
    a.axvline(0, color=TEXTO2, linewidth=0.6)
    a.set_yticks([])
    a.spines["left"].set_visible(False)
    a.set_xlim(-22, 22)
    a.set_xticks([-15, -10, -5, 0, 5, 10, 15], ["15", "10", "5", "0", "5", "10", "15"])
    a.set_xlabel("−log₁₀ adjusted p (MSigDB Hallmark)")
    a.text(-21, len(filas) - 0.2, "down in knockout", fontsize=7, color=AZUL, fontweight="semibold", va="bottom")
    a.text(21, len(filas) - 0.2, "up in knockout", fontsize=7, color=ROJO, fontweight="semibold", va="bottom",
           ha="right")
    a.set_title(f"a  BUB1B knockout: {nf} signatures, {nl} lines", fontsize=8.4, pad=14)
    y = 0
    etiquetas = []
    for nombre, color, drogas in grupos[::-1]:
        for d in drogas[::-1]:
            clave, et = (d if isinstance(d, tuple) else (d, d))
            v = lin.get(clave, {})
            neto = v.get("R", 0) - v.get("M", 0)
            b.barh(y, neto, color=color, height=0.62, edgecolor="white")
            etiquetas.append((y, et))
            y += 1
        b.annotate(nombre, xy=(-0.02, y + 0.05), xycoords=("axes fraction", "data"), ha="right", va="center",
                   fontsize=6.9, color=color, fontweight="semibold")
        y += 1.3
    b.set_yticks([e[0] for e in etiquetas], [e[1] for e in etiquetas], fontsize=7)
    b.tick_params(axis="y", length=0)
    b.axvline(0, color=TEXTO2, linewidth=0.6)
    b.set_xlim(-9, 9.5)
    b.set_xlabel("cell lines reversing − mimicking (FDR < 0.05)")
    b.text(-8.8, -1.4, "← mimics the knockout", fontsize=6.8, color=TEXTO2)
    b.text(8.8, -1.4, "reverses it →", fontsize=6.8, color=TEXTO2, ha="right")
    b.set_ylim(-1.8, y)
    b.spines["left"].set_visible(False)
    b.set_title("b  Connectivity by drug class", fontsize=8.4, pad=14)
    guarda(fig, "t2_firma")


def ranking():
    filas = [r for r in csv.DictReader((T2 / "candidatos.csv").open())
             if r["modo"] != "avoid" and r["estado_regulatorio"].startswith("approved")]
    filas.sort(key=lambda r: -float(r["score"]))
    sens = muestras(filas)
    fig, (a, b) = plt.subplots(1, 2, figsize=(ANCHO, 3.0), sharey=True,
                               gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.08))
    y = np.arange(len(filas))[::-1]
    for yy, r in zip(y, filas):
        et, c = LINEA[r["hipotesis"]]
        a.barh(yy, float(r["score"]), color=c, height=0.64, edgecolor="white")
        a.text(float(r["score"]) + 0.008, yy, f"{float(r['score']):.2f}", va="center", fontsize=6.9, color=TEXTO2)
        rk = sens[r["farmaco"]]
        b.boxplot([rk], positions=[yy], vert=False, widths=0.55, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=GRIS_CLARO, edgecolor=TEXTO2, linewidth=0.6),
                  medianprops=dict(color=TEXTO, linewidth=1.0), whiskerprops=dict(color=TEXTO2, linewidth=0.6),
                  capprops=dict(color=TEXTO2, linewidth=0.6))
        t3 = 100 * sum(x <= 3 for x in rk) / len(rk)
        b.text(len(filas) + 3.8, yy, f"{t3:.0f} %", va="center", ha="right", fontsize=6.9,
               color=TEXTO if t3 >= 50 else TENUE, fontweight="semibold" if t3 >= 50 else "normal")
    a.set_yticks(y, [r["farmaco"] for r in filas], fontsize=7.4)
    a.tick_params(axis="y", length=0)
    a.set_xlim(0.45, 0.80)
    a.set_xlabel("score (0.30·E + 0.25·M + 0.20·P + 0.25·S)")
    a.set_title("a  Approved medicines, ranked", fontsize=8.4)
    a.spines["left"].set_visible(False)
    vistos = []
    for k, (et, c) in LINEA.items():
        if any(r["hipotesis"] == k for r in filas) and c != GRIS:
            vistos.append(plt.Rectangle((0, 0), 1, 1, color=c, label=et))
    a.legend(handles=vistos, loc="lower right", fontsize=6.5, handlelength=0.9, borderaxespad=0.1)
    b.set_xlim(0.5, len(filas) + 4.0)
    b.set_xticks([1, 3, 6, 9, 12, 15])
    b.axvline(len(filas) + 0.7, color=GRIS_CLARO, linewidth=0.6)
    b.axvline(3.5, color=GRIS, linewidth=0.5, linestyle=(0, (3, 3)))
    b.set_xlabel("rank across 2,000 draws (weights, M ± 0.30, S ± 0.15)")
    b.set_title("b  Robustness · % of draws in top 3", fontsize=8.4)
    b.spines["left"].set_visible(False)
    b.tick_params(axis="y", length=0)
    guarda(fig, "t2_ranking")


def muestras(filas, n=2000, semilla=7):
    import candidatos
    return candidatos.rangos(filas, n=n, semilla=semilla, ruido_m=0.30, ruido_s=0.15)


def cambio():
    """Gráfico de pendiente: puesto antes (solo bibliografía + CMap) y después de las revisiones de evidencia."""
    v2 = {"dasatinib": 1, "nicotinamide": 2, "quercetin": 3, "sirolimus": 4, "niacin": 5, "nicotinamide riboside": 6,
          "fisetin": 7, "amlexanox": 8, "ataluren": 9, "everolimus": 10, "gentamicin": 11, "ruxolitinib": 12,
          "tofacitinib": 13, "fedratinib": 14, "hydroxychloroquine": 15, "metformin": 16, "baricitinib": 17}
    filas = [r for r in csv.DictReader((T2 / "candidatos.csv").open())
             if r["modo"] != "avoid" and r["estado_regulatorio"].startswith("approved")]
    filas.sort(key=lambda r: -float(r["score"]))
    v3 = {r["farmaco"]: i for i, r in enumerate(filas, 1)}
    motivo = {"quercetin": "not an approved medicine", "fisetin": "not an approved medicine",
              "nicotinamide riboside": "not an approved medicine", "ataluren": "read-through activity = assay artefact",
              "gentamicin": "nephrotoxicity vs nephrocalcinosis (safety)"}
    destaca = {"nicotinamide": (CRITICO, "RMS depends on NAMPT/NAD$^{+}$"), "sirolimus": (NARANJA, ""),
               "amlexanox": (GRIS, "UGA-A context + NMD")}
    fig, ax = plt.subplots(figsize=(ANCHO, 3.9))
    base = len(v3) + 1.2                      # los "retirados" empiezan bajo el último puesto real
    k = 0
    for d, r2 in v2.items():
        c, nota = destaca.get(d, (GRIS, ""))
        if d in v3:
            r3 = v3[d]
            ax.plot([0, 1], [r2, r3], color=c, linewidth=1.6 if d in destaca else 0.8, alpha=1 if d in destaca else 0.7,
                    solid_capstyle="round")
            ax.scatter([0, 1], [r2, r3], s=14, color=c, zorder=3, edgecolor="white", linewidth=0.5)
            ax.text(1.03, r3, f"{r3}  {d}" + (f"  · {nota}" if nota else ""), va="center", fontsize=7,
                    color=TEXTO if d in destaca else TEXTO2, fontweight="semibold" if d in destaca else "normal")
        else:
            yr = base + k
            k += 1
            ax.plot([0, 1], [r2, yr], color=GRIS_CLARO, linewidth=0.8, linestyle=(0, (3, 2)))
            ax.scatter([1], [yr], marker="x", s=16, color=TENUE, zorder=3, linewidth=0.8)
            ax.text(1.03, yr, f"removed · {d}: {motivo[d]}", va="center", fontsize=6.7, color=TENUE, style="italic")
        ax.text(-0.03, r2, f"{d}  {r2}", va="center", ha="right", fontsize=7, color=TEXTO if d in destaca else TEXTO2)
    # entradas NUEVAS: líneas que la propia revisión añadió (H7, considerada y declinada)
    nuevos = [(d, r) for d, r in v3.items() if d not in v2]
    for j, (d, r3) in enumerate(nuevos):
        ax.plot([0, 1], [base + k + 0.4 + j, r3], color=GRIS_CLARO, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.scatter([0], [base + k + 0.4 + j], marker="+", s=20, color=TENUE, zorder=3, linewidth=0.8)
        ax.text(-0.03, base + k + 0.4 + j, f"added · {d}", va="center", ha="right", fontsize=6.7, color=TENUE,
                style="italic")
        ax.scatter([1], [r3], s=14, color=GRIS, zorder=3, edgecolor="white", linewidth=0.5)
        ax.text(1.03, r3, f"{r3}  {d}  · declined (H7)", va="center", fontsize=7, color=TEXTO2)
    ax.axhline(len(v3) + 0.55, xmin=0.53, xmax=0.98, color=GRIS_CLARO, linewidth=0.6)
    ax.set_ylim(base + k + len(nuevos) + 0.9, 0.2)
    ax.set_xlim(-0.55, 1.75)
    ax.axis("off")
    ax.text(0, -0.4, "before the evidence reviews", ha="center", fontsize=7.8, fontweight="semibold", color=TEXTO)
    ax.text(1, -0.4, "after the evidence reviews", ha="center", fontsize=7.8, fontweight="semibold", color=TEXTO)
    guarda(fig, "t2_cambio")


def ejes():
    """Por línea (su mejor fármaco aprobado), coloreado por eje; y qué dos hechos subirían la línea correctora."""
    import yaml
    import candidatos
    cfg = yaml.safe_load((ROOT / "scripts" / "t2" / "hipotesis.yaml").read_text())
    filas = [r for r in csv.DictReader((T2 / "candidatos.csv").open())
             if r["modo"] != "avoid" and r["estado_regulatorio"].startswith("approved")]
    for r in filas:
        r["score"] = float(r["score"])
    filas.sort(key=lambda r: -r["score"])
    COL = {"correct": VERDE, "clear": AZUL, "tolerate": NARANJA}
    fig, (a, b) = plt.subplots(1, 2, figsize=(ANCHO, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.55))
    # a · mejor fármaco de cada línea
    lineas = []
    for h in cfg["hipotesis"]:
        fs = [r for r in filas if r["hipotesis"] == h["id"]]
        if fs:
            lineas.append((h["id"], h.get("eje", ""), fs[0]["farmaco"], fs[0]["score"], filas.index(fs[0]) + 1))
    lineas.sort(key=lambda x: -x[3])
    ys = []
    for i, (hid, eje, farm, sc, rk) in enumerate(lineas):
        y = len(lineas) - 1 - i
        ys.append((y, f"{hid} · {farm}"))
        a.barh(y, sc, color=COL[eje], height=0.62, edgecolor="white")
        a.text(sc + 0.004, y, f"{sc:.2f} · #{rk}", va="center", fontsize=6.8, color=TEXTO2)
    a.set_yticks([y for y, _ in ys], [s for _, s in ys], fontsize=7)
    a.tick_params(axis="y", length=0)
    a.spines["left"].set_visible(False)
    a.set_xlim(0.45, 0.79)
    a.set_xlabel("score of the line's best approved drug")
    a.set_title("a  Each line by its best approved drug", fontsize=8.4)
    a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, label=e) for e, c in COL.items()], loc="lower right",
             fontsize=6.5, handlelength=0.9, borderaxespad=0.1)
    # b · escenarios: puesto del mejor fármaco corrector
    cand = [dict(r) for r in filas]
    esc = candidatos.contrafactual(cfg, cand)
    nombres = ["today", "class A\nconfirmed", "in\nremission", "both"]
    for j, (nombre, fs) in enumerate(esc):
        top = [f for _, f in fs]
        rk = next(i + 1 for i, f in enumerate(top) if f["hipotesis"] == "H1")
        corr = top[rk - 1]["farmaco"]
        alto = len(top) + 1 - rk
        b.bar(j, alto, color=VERDE if rk <= 2 else lighter_hex(VERDE), width=0.62, edgecolor="white")
        b.text(j, alto + 0.4, f"#{rk}", ha="center", fontsize=7.4, color=TEXTO, va="bottom", fontweight="semibold")
        b.text(j, alto - 0.6, corr, ha="center", va="top", fontsize=6.4, color="white", rotation=90)
    b.set_xticks(range(len(nombres)), nombres, fontsize=6.6)
    b.set_yticks([])
    b.spines["left"].set_visible(False)
    b.set_ylim(0, len(cand) + 4.5)
    b.set_title("b  The NAD$^{+}$ line (H1): best drug's rank", fontsize=8.4)
    guarda(fig, "t2_ejes")


def lighter_hex(c, k=0.45):
    c = c.lstrip("#")
    r, g, bb = (int(c[i:i + 2], 16) for i in (0, 2, 4))
    return "#{:02x}{:02x}{:02x}".format(*(int(v + (255 - v) * k) for v in (r, g, bb)))


def nof1():
    """Potencia del n-de-1 sobre %MN-RET (σ estimada de los datos) y cuándo sirve contar más células."""
    import json
    r = json.loads((T2 / "nof1.json").read_text())
    fig, (a, b) = plt.subplots(1, 2, figsize=(ANCHO, 2.95), gridspec_kw=dict(width_ratios=[1.1, 1], wspace=0.38))
    estilos = [("reduction >= 25 %, 95 % certainty", ROJO, "≥25 % reduction, 95 % sure"),
               ("any reduction, 95 % certainty", AZUL, "any reduction, 95 % sure"),
               ("any reduction, 90 % certainty", AQUA, "any reduction, 90 % sure")]
    fp = r["fp_regla"]
    for clave, color, et in estilos:
        d = r["potencia_efecto"][clave]
        ef = [100 * float(k) for k in d]
        po = [100 * v[0] for v in d.values()]
        fpt = "< 1 %" if fp[clave][0] < 0.01 else f"{100 * fp[clave][0]:.0f} %"
        a.plot(ef, po, color=color, marker="o", markersize=3.2, markeredgecolor="white", markeredgewidth=0.5,
               label=f"{et} · {fpt} false pos.")
    a.axhline(80, color=GRIS, linewidth=0.6, linestyle=(0, (3, 3)))
    a.text(83, 75, "80 % power", fontsize=6.6, color=TEXTO2, ha="right")
    a.set_xlabel("true reduction in %MN-RET")
    a.set_ylabel("probability of declaring a responder")
    a.set_xlim(-3, 84)
    a.set_ylim(-4, 138)  # espacio por encima del 100 % para que la leyenda no tape las curvas
    a.set_yticks(range(0, 101, 20))
    a.set_title("a  Power, 4 visits per condition", fontsize=8.4)
    a.legend(loc="upper left", fontsize=5.8, handlelength=1.2, borderaxespad=0.2, bbox_to_anchor=(0.0, 1.0),
             title="decision rule · false positives (no effect, 8 visits)", title_fontsize=5.8, alignment="left")
    # b · rejilla basal × células
    basales = sorted({float(k.split("|")[0]) for k in r["rejilla"]})
    celdas = sorted({int(k.split("|")[1]) for k in r["rejilla"]})
    colores = RAMPA[::-1][:len(basales)]
    for c, bas in zip(colores, basales):
        po = [100 * r["rejilla"][f"{bas:.4f}|{n}"][0] for n in celdas]
        b.plot(celdas, po, color=c, marker="o", markersize=3, markeredgecolor="white", markeredgewidth=0.5,
               label=f"{100 * bas:g} %")
    b.set_xscale("log")
    b.set_xticks(celdas, [f"{n // 1000}k" for n in celdas])
    b.minorticks_off()
    b.set_xlabel("reticulocytes scored per sample")
    b.set_ylabel("power at a true 35 % reduction")
    b.set_xlim(celdas[0] * 0.8, celdas[-1] * 1.4)
    b.set_ylim(0, 60)
    b.legend(title="baseline %MN-RET", loc="lower right", fontsize=6.2, title_fontsize=6.4, handlelength=1.2,
             borderaxespad=0.2)
    b.set_title("b  More cells help only at low baselines", fontsize=8.4)
    guarda(fig, "t2_nof1")


if __name__ == "__main__":
    for f in (firma, ranking, cambio, ejes, nof1):
        f()
