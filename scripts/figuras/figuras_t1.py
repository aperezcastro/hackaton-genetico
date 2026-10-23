"""figuras_t1.py — figuras del informe del Track 1 (PDF vectorial en reports/latex/fig/).

Entradas: agregados de results/t1/ (embudo, top de Exomiser, SpliceAI del panel, equilibrio alélico) y datos
públicos (UniProt O60566, modelo AlphaFold). Ninguna figura muestra genotipos más allá de las variantes enviadas.

    .venv/bin/python scripts/figuras/figuras_t1.py
"""
import csv
import json
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from estilo import *  # noqa: E402,F403

R1 = ROOT / "results" / "t1"
aplica()


def embudo():
    e = json.loads((R1 / "embudo.json").read_text())
    pasos = [("Normalised, non-reference", e["0_total_normalizadas"]),
             ("FILTER = PASS", e["1_pass"]),
             ("gnomAD v4.1 AF < 1 %", e["2_pass_af<1%"]),
             ("HIGH / MODERATE / splice", e["3_raras_funcionales(HIGH/MOD/splice)"]),
             ("Recessive gene model (99 genes)", e["4_variantes_en_genes_recesivos_candidatos"]),
             ("MVA / SAC panel genes", 2)]
    fig, ax = plt.subplots(figsize=(ANCHO, 2.3))
    y = np.arange(len(pasos))[::-1]
    for i, ((et, n), yy) in enumerate(zip(pasos, y)):
        ax.barh(yy, n, height=0.62, color=RAMPA[min(i + 2, len(RAMPA) - 1)], edgecolor="white", linewidth=1.5)
        txt = f"{n:,}" if i < 5 else "2 variants · 1 gene: BUB1B"
        ax.text(n * 1.25, yy, txt, va="center", ha="left", fontsize=8, color=TEXTO,
                fontweight="semibold" if i == 5 else "normal")
    ax.set_yticks(y, [p[0] for p in pasos])
    ax.set_xscale("log")
    ax.set_xlim(1, 4e8)
    ax.set_xlabel("variants (log scale)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    guarda(fig, "t1_embudo")


def proteina():
    import Bio.PDB as bp
    ch = bp.PDBParser(QUIET=True).get_structure("b", ROOT / "refs/estructura/bubr1_af.pdb")[0]["A"]
    res = sorted((r.id[1], r["CA"].get_bfactor()) for r in ch if "CA" in r)
    x, pl = zip(*res)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(ANCHO, 2.75), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.0, 1.35], hspace=0.08))
    # pista 1: confianza AlphaFold
    a1.fill_between(x, pl, color=AZUL, alpha=0.18, linewidth=0)
    a1.plot(x, pl, color=AZUL, linewidth=0.8)
    for lim, et in ((90, "very high"), (70, "confident")):
        a1.axhline(lim, color=GRIS, linewidth=0.5, linestyle=(0, (3, 3)))
        a1.text(1058, lim, et, fontsize=6.8, color=TEXTO2, va="center")
    a1.set_ylim(20, 100)
    a1.set_ylabel("AlphaFold\npLDDT")
    a1.annotate("Asn1002: pLDDT 91\n(buried, RSA 0.16)", xy=(1002, 91), xytext=(860, 40), fontsize=7, color=TEXTO,
                arrowprops=dict(arrowstyle="-", color=TEXTO2, linewidth=0.6))
    a1.spines["bottom"].set_visible(False)
    a1.tick_params(axis="x", length=0)
    # pista 2: dominios (UniProt) + variantes
    a2.set_ylim(-1.7, 1.9)
    a2.add_patch(mpatches.FancyBboxPatch((1, -0.18), 1049, 0.36, boxstyle="round,pad=0,rounding_size=0.08",
                                         fc=GRIS_CLARO, ec="none"))
    for (s, e_, et, c) in ((62, 226, "BUB1-N", VIOLETA), (766, 1050, "kinase-like domain", AZUL)):
        a2.add_patch(mpatches.FancyBboxPatch((s, -0.3), e_ - s, 0.6, boxstyle="round,pad=0,rounding_size=0.1",
                                             fc=c, ec="white", lw=1, alpha=0.9))
        a2.text((s + e_) / 2, 0, et, ha="center", va="center", fontsize=7.4, color="white", fontweight="semibold")
    a2.add_patch(mpatches.Rectangle((224, -0.3), 8, 0.6, fc=TEXTO2, ec="none"))
    a2.text(228, -0.45, "D-box", ha="center", va="top", fontsize=6.6, color=TEXTO2)
    for pos, et, c, alto in ((737, "p.Leu737*\nnull (NMD)", NARANJA, 1.35), (1002, "p.Asn1002Lys\nhypomorph", AZUL, 1.35)):
        a2.plot([pos, pos], [0.3, alto - 0.25], color=c, linewidth=1.2)
        a2.scatter([pos], [alto - 0.25], s=34, color=c, zorder=3, edgecolor="white", linewidth=0.8)
        a2.text(pos, alto - 0.05, et, ha="center", va="bottom", fontsize=7.4, color=TEXTO, fontweight="semibold",
                linespacing=1.05)
    # Q921H no se dibuja: es un pasajero funcionalmente normal en cis con L844F (Suijkerbuijk 2010), no un alelo MVA
    conocidas = [(550, "R550Q"), (814, "R814H"), (844, "L844F"), (909, "I909T"), (1012, "L1012P")]
    for i, (pos, et) in enumerate(conocidas):
        y0 = -0.95 - 0.38 * (i % 2)
        a2.plot([pos, pos], [-0.3, y0 + 0.12], color=GRIS, linewidth=0.8)
        a2.scatter([pos], [y0 + 0.12], s=16, color=GRIS, zorder=3, edgecolor="white", linewidth=0.6)
        a2.text(pos, y0, et, ha="center", va="top", fontsize=6.6, color=TEXTO2)
    a2.text(1, -1.55, "published MVA1 missense alleles (UniProt; Q921H omitted, a passenger in cis with L844F)",
            fontsize=6.8, color=TEXTO2, style="italic")
    a2.set_yticks([])
    a2.spines["left"].set_visible(False)
    a2.set_xlim(0, 1110)
    a2.set_xlabel("BUBR1 residue (1,050 aa; UniProt O60566)")
    guarda(fig, "t1_proteina")


def exomiser():
    """Solo se nombran los 3 primeros (discutidos en el informe); los puestos 4-10 se agregan SIN nombres de gen ni
    modo de herencia, para no insinuar variantes del paciente en otros genes (posibles hallazgos incidentales)."""
    todos = json.loads((R1 / "exomiser_top.json").read_text())["A"]
    top = todos[:3]
    resto = todos[3:10]
    fig, ax = plt.subplots(figsize=(ANCHO_MEDIO + 0.2, 2.2))
    filas = [(f"{t['rank']}  {t['gene']} ({t['moi']})", t["combined"], t["gene"] == "BUB1B") for t in top]
    filas.append((f"ranks 4–10 (7 genes, not named)", max(t["combined"] for t in resto), False))
    y = np.arange(len(filas))[::-1]
    for (et, v, es), yy in zip(filas, y):
        ax.barh(yy, v, height=0.62, color=AZUL if es else GRIS, edgecolor="white", linewidth=1.2)
        txt = f"{v:.3f}" if not et.startswith("ranks") else f"≤ {v:.3f}"
        ax.text(v + 0.012, yy, txt, va="center", fontsize=7.2, color=TEXTO if es else TEXTO2,
                fontweight="semibold" if es else "normal")
    ax.set_yticks(y, [f[0] for f in filas], fontsize=7.4)
    for lab, f in zip(ax.get_yticklabels(), filas):
        if f[2]:
            lab.set_fontweight("semibold"); lab.set_color(TEXTO)
    ax.set_xlim(0, 0.72)
    ax.set_xlabel("Exomiser combined score")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title("Genome-wide ranking (HPO + variants)", fontsize=8.6)
    guarda(fig, "t1_exomiser")


def acmg():
    var = [("c.2210T>G\np.(Leu737*)", [("PVS1", 8, 4), ("PM2_Sup", 1, 1)],
            "Likely pathogenic · 9\n(PM3 not applied: partner is a VUS)"),
           ("c.3006T>G\np.(Asn1002Lys)", [("PM2_Sup", 1, 1), ("PP3", 1, 1), ("PM1_Sup", 1, 1), ("PM3_Sup", 1, 1), ("PP4", 1, 1)],
            "VUS · 4–5 (PP4 contestable; LP needs\nphase confirmed and PP4)"),
           ("c.2679-1026A>G\n(deep intronic)", [("PM2_Sup", 1, 1), ("BP4", -1, 0)], "VUS · 0")]
    fig, ax = plt.subplots(figsize=(5.4, 2.5))
    for x0, x1, et, c in ((-2, 0, "LB", "#f6f6f4"), (0, 6, "VUS", "#f6f6f4"), (6, 10, "Likely\npathogenic", "#eef4fc"),
                          (10, 12.5, "Pathogenic", "#e3eefb")):
        ax.axvspan(x0, x1, color=c, zorder=0, linewidth=0)
        if et:
            ax.text((x0 + x1) / 2, 2.62, et, ha="center", va="bottom", fontsize=6.8, color=TEXTO2)
    tono = {8: RAMPA[6], 1: RAMPA[1]}
    for i, (nombre, crit, clase) in enumerate(var):
        yy = 2 - i
        izq = 0
        for c, p, _ in crit:
            if p > 0:
                ax.barh(yy, p, left=izq, height=0.52, color=tono[8 if p >= 8 else 1], edgecolor="white", linewidth=1.5)
                ax.text(izq + p / 2, yy, c, ha="center", va="center", fontsize=6.6,
                        color="white" if p >= 8 else TEXTO)
                izq += p
            else:
                ax.barh(yy, p, left=0, height=0.52, color=NARANJA, edgecolor="white", linewidth=1.5)
                ax.text(p / 2, yy, c, ha="center", va="center", fontsize=6.6, color="white")
        ax.text(12.4, yy - 0.42, clase, ha="right", va="top", fontsize=6.6, color=TEXTO2, style="italic")
    ax.axvline(0, color=TEXTO2, linewidth=0.6)
    ax.set_yticks([2, 1, 0], [v[0] for v in var], fontsize=7.2)
    ax.set_xlim(-2, 12.5)
    ax.set_ylim(-0.95, 2.95)
    ax.set_xlabel("ACMG/AMP Bayesian points (Tavtigian et al. 2020)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    guarda(fig, "t1_acmg")


def fase():
    v1, v3, v2 = 40209701, 40216470, 40220612
    fig, ax = plt.subplots(figsize=(ANCHO, 1.45))
    ax.plot([v1 - 400, v2 + 400], [0, 0], color=GRIS, linewidth=3, solid_capstyle="round", zorder=1)
    for p, et, c in ((v1, "p.Leu737*", NARANJA), (v3, "c.2679-1026A>G", GRIS), (v2, "p.Asn1002Lys", AZUL)):
        ax.scatter([p], [0], s=46, color=c, zorder=3, edgecolor="white", linewidth=1)
        ax.text(p, 0.32, et, ha="center", fontsize=7.4, color=TEXTO, fontweight="semibold")
    for a, b in ((v1, v3), (v3, v2)):
        ax.annotate("", xy=(a + 80, -0.32), xytext=(b - 80, -0.32),
                    arrowprops=dict(arrowstyle="<->", color=TEXTO2, linewidth=0.7))
        ax.text((a + b) / 2, -0.45, f"{b - a:,} bp without any heterozygous site", ha="center", va="top",
                fontsize=7, color=TEXTO2)
    # par de lecturas a escala
    x0 = v1 + 600
    ax.add_patch(mpatches.Rectangle((x0, 0.62), 150, 0.07, color=VIOLETA))
    ax.add_patch(mpatches.Rectangle((x0 + 350, 0.62), 150, 0.07, color=VIOLETA))
    ax.plot([x0 + 150, x0 + 350], [0.655, 0.655], color=VIOLETA, linewidth=0.6, linestyle=(0, (2, 2)))
    ax.text(x0 + 560, 0.655, "one Illumina read pair (~500 bp insert), to scale", va="center", fontsize=6.8,
            color=TEXTO2)
    ax.set_ylim(-0.95, 0.9)
    ax.set_xlim(v1 - 500, v2 + 500)
    ax.axis("off")
    ax.text(v1 - 500, -0.9, "chr15 (GRCh38) · 10,911 bp between the two submitted alleles · 3 heterozygous PASS sites",
            fontsize=6.8, color=TENUE)
    guarda(fig, "t1_fase")


def spliceai():
    filas = list(csv.DictReader((R1 / "spliceai.tsv").open(), delimiter="\t"))
    mx = {}
    for r in filas:
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        mx[k] = max(mx.get(k, (0, ""))[0], float(r["DS_max"])), r["gen"]
    rng = np.random.default_rng(3)
    fig, ax = plt.subplots(figsize=(ANCHO_MEDIO + 0.2, 1.6))
    nombres = {("15", "40209701"): "p.Leu737*", ("15", "40220612"): "p.Asn1002Lys", ("15", "40216470"): "deep intronic"}
    for k, (d, g) in mx.items():
        yy = rng.uniform(-0.25, 0.25)
        es = (k[0], k[1]) in nombres
        ax.scatter([d], [yy], s=22 if es else 14, color=AZUL if es else GRIS, edgecolor="white", linewidth=0.6,
                   zorder=3 if es else 2)
    for x_, et in ((0.2, "possible"), (0.5, "likely"), (0.8, "high confidence")):
        ax.axvline(x_, color=GRIS, linewidth=0.5, linestyle=(0, (3, 3)))
        ax.text(x_ + 0.01, 0.42, et, fontsize=6.6, color=TEXTO2)
    ax.text(0.98, -0.62, f"{len(mx)} scorable rare PASS panel variants · max Δ ≤ 0.03\n"
            "BUB1B variants in blue · windows 50 and 2,000 nt", fontsize=6.6, color=TEXTO2, va="center", ha="right",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(-0.95, 0.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("SpliceAI maximum delta score")
    guarda(fig, "t1_spliceai")


def balance():
    d = json.loads((R1 / "balance_alelico.json").read_text())
    dg = json.loads((R1 / "balance_alelico_diagnostico.json").read_text())
    fig = plt.figure(figsize=(ANCHO, 3.6))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.05], hspace=0.75, wspace=0.35)
    a0 = fig.add_subplot(gs[0, :])
    cr = list(d["autosomas"])
    so = [d["autosomas"][c]["sobredispersion"] for c in cr]
    a0.bar(range(len(cr)), so, width=0.7, color=[NARANJA if s > 1.9 else GRIS for s in so], edgecolor="white")
    a0.axhline(1, color=TEXTO2, linewidth=0.6, linestyle=(0, (3, 3)))
    a0.set_xticks(range(len(cr)), cr)
    a0.set_ylim(0, 3.9)
    a0.set_ylabel("overdispersion")
    a0.set_title("Allelic balance of heterozygous SNVs per autosome (2.24 M sites)", fontsize=8.6)
    a0.text(19, 3.6, "excess on 9, 17, 20–22 …", fontsize=7, color=TEXTO, ha="center")
    sel = ["1", "9", "17", "20", "21", "22"]
    for i, c in enumerate(sel):
        ax = fig.add_subplot(gs[1, i])
        bins = dg["bins"][c]
        xs = [int(k.split("-")[0]) + 2.5 for k in bins]
        ys = [v["sobredisp"] for v in bins.values()]
        ax.axhspan(0, 2, color="#f3f6fb", linewidth=0)
        ax.plot(xs, ys, color=AZUL, linewidth=1.0)
        ax.scatter([x for x, y in zip(xs, ys) if y > 2], [y for y in ys if y > 2], s=12, color=NARANJA, zorder=3,
                   edgecolor="white", linewidth=0.5)
        ax.set_title(f"chr{c}", fontsize=7.6, loc="center")
        ax.set_ylim(0, 14)
        ax.tick_params(labelsize=6.6)
        if i:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel("overdispersion\n(5-Mb bins)", fontsize=7.2)
        ax.set_xlabel("Mb", fontsize=6.8, labelpad=1)
    fig.text(0.5, -0.04, "… is focal: confined to pericentromeric / segmental-duplication bins (orange), where depth is "
             "inflated and allele fraction depressed (collapsed paralogues).\nEverywhere else every bin is balanced: "
             "no detectable whole-chromosome or arm-level imbalance.", ha="center", fontsize=6.9, color=TEXTO2)
    guarda(fig, "t1_balance")


if __name__ == "__main__":
    for f in (embudo, proteina, exomiser, acmg, fase, spliceai, balance):
        f()
