#!/usr/bin/env python3
"""filtra.py — corre DENTRO de Modal sobre la salida de VEP (anotar.py). Produce tablas cortas y el embudo.

    python filtra.py <dir_anot> <dir_salida> <croms separados por comas>

Salidas (en <dir_salida>, que anotar.py baja a results/t1/):
  embudo.json            conteos de cada paso del filtrado (QC; no contiene variantes)
  raras_funcionales.tsv  PASS, AF_max(gnomAD exomas/genomas) < 1 %, impacto HIGH/MODERATE o splicing, o ClinVar P/LP
  genes_recesivos.tsv    genes con ≥2 variantes raras_funcionales (candidatos a heterocigoto compuesto) o ≥1 homocigota
  panel_mva.tsv          TODA variante (cualquier consecuencia: intrónica, UTR, promotor…) con AF_max < 5 % en los
                         genes del panel MVA/SAC — para no perder el alelo hipomórfico discreto

Decisiones (justificarlas en el informe): AF_max = máximo de gnomADe_AF y gnomADg_AF globales (v4.1); el
transcrito que se informa es el MANE Select si existe, si no el de mayor impacto; «PASS» según FILTER del VCF
(Sentieon + GATK VariantFiltration).
"""
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

PANEL = {
    # MVA descritos (OMIM/literatura)
    "BUB1B", "CEP57", "TRIP13", "CENATAC", "CCDC84", "BUB1", "BUB3",
    # candidatos de la llamada ERN-ITHACA y del mismo circuito (SAC, cinetocoro, RZZ, MCC, APC/C)
    "MAD1L1", "MAD2L1", "MAD2L1BP", "SLF2", "SMC5", "SMC6", "NSMCE2", "CEP192", "TTK", "KNL1", "ZWILCH", "ZW10",
    "KNTC1", "CDC20", "CENPE", "AURKB", "INCENP", "CDCA8", "BIRC5", "PLK1", "SGO1", "ESCO2", "NDC80", "SPC24",
    "SPC25", "NUF2", "SKA1", "SKA3", "PPP2R5C", "CENPF", "CENPJ", "PLK4", "STIL", "SASS6",
    # espliceosoma menor (vía de CENATAC)
    "RNU4ATAC", "RNPC3", "ZRSR2", "SNRNP48", "PDCD7",
}
IMP = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def flt(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def main():
    src, dst, croms = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    os.makedirs(dst, exist_ok=True)
    emb = Counter()
    raras, panel = [], []
    cols = None
    for c in croms:
        p = f"{src}/{c}.vcf.gz"
        if not os.path.exists(p):
            emb["croms_sin_fichero"] += 1
            continue
        with gzip.open(p, "rt") as f:
            for ln in f:
                if ln.startswith("##INFO=<ID=CSQ"):
                    cols = ln.split("Format: ")[1].rstrip('">\n').split("|")
                    continue
                if ln.startswith("#"):
                    continue
                F = ln.rstrip("\n").split("\t")
                emb["0_total_normalizadas"] += 1
                chrom, pos, _, ref, alt, qual, filt, info = F[:8]
                fmt = dict(zip(F[8].split(":"), F[9].split(":"))) if len(F) > 9 else {}
                gt = fmt.get("GT", "").replace("|", "/")
                if gt in ("0/0", "./.", "0", "."):
                    emb["x_genotipo_ref_o_nulo"] += 1
                    continue
                pas = filt == "PASS"
                emb["1_pass"] += pas
                csq = next((kv[4:] for kv in info.split(";") if kv.startswith("CSQ=")), "")
                anns = [dict(zip(cols, a.split("|"))) for a in csq.split(",")] if csq and cols else []
                af = max([flt(a.get("gnomADe_AF")) for a in anns] + [flt(a.get("gnomADg_AF")) for a in anns] + [0.0])
                mane = [a for a in anns if a.get("MANE_SELECT")]
                best = max(anns, key=lambda a: IMP.get(a.get("IMPACT"), 0)) if anns else {}
                rep = max(mane, key=lambda a: IMP.get(a.get("IMPACT"), 0)) if mane else best
                imp = best.get("IMPACT", "")
                cons_all = ",".join(sorted({x for a in anns for x in a.get("Consequence", "").split("&")}))
                splice = "splice" in cons_all
                clin = ",".join(sorted({a.get("CLIN_SIG", "") for a in anns if a.get("CLIN_SIG")} |
                                       {a.get("ClinVarNCBI_CLNSIG", "") for a in anns if a.get("ClinVarNCBI_CLNSIG")}))
                plp = "pathogenic" in clin.lower() and "conflicting" not in clin.lower()
                genes = sorted({a.get("SYMBOL") for a in anns if a.get("SYMBOL")})
                ad = fmt.get("AD", "")
                try:
                    a0, a1 = [int(x) for x in ad.split(",")[:2]]
                    vaf = round(a1 / (a0 + a1), 3) if a0 + a1 else None
                except ValueError:
                    vaf = None
                row = dict(chrom=chrom, pos=pos, ref=ref, alt=alt, filter=filt, qual=qual, gt=fmt.get("GT", ""),
                           ad=ad, vaf=vaf, dp=fmt.get("DP", ""), gq=fmt.get("GQ", ""), pgt=fmt.get("PGT", ""),
                           pid=fmt.get("PID", ""), gene=rep.get("SYMBOL", ""), genes=";".join(genes),
                           consequence=rep.get("Consequence", ""), impact_max=imp, mane=rep.get("MANE_SELECT", ""),
                           hgvsc=rep.get("HGVSc", ""), hgvsp=rep.get("HGVSp", ""), exon=rep.get("EXON", ""),
                           intron=rep.get("INTRON", ""), af_max=af, gnomADe_AF=rep.get("gnomADe_AF", ""),
                           gnomADg_AF=rep.get("gnomADg_AF", ""), am_class=rep.get("am_class", ""),
                           am_pathogenicity=rep.get("am_pathogenicity", ""), sift=rep.get("SIFT", ""),
                           polyphen=rep.get("PolyPhen", ""), clinvar=clin,
                           clinvar_stars=next((a.get("ClinVarNCBI_CLNREVSTAT", "") for a in anns
                                               if a.get("ClinVarNCBI_CLNREVSTAT")), ""),
                           existing=rep.get("Existing_variation", ""))
                if set(genes) & PANEL and af < 0.05:
                    panel.append(row)
                if not pas:
                    continue
                if af < 0.01:
                    emb["2_pass_af<1%"] += 1
                    if af < 0.001:
                        emb["2b_pass_af<0.1%"] += 1
                    func = IMP.get(imp, 0) >= 2 or splice
                    emb["3_raras_funcionales(HIGH/MOD/splice)"] += func
                    if func or plp:
                        raras.append(row)
                elif plp and af < 0.05:
                    raras.append(row)
    por_gen = defaultdict(list)
    for r in raras:
        por_gen[r["gene"] or r["genes"]].append(r)
    genes_rec = []
    for g, rs in por_gen.items():
        hom = [r for r in rs if r["gt"].replace("|", "/") == "1/1"]
        if len(rs) >= 2 or hom:
            pids = {r["pid"] for r in rs if r["pid"]}
            genes_rec.append(dict(gene=g, n_variantes=len(rs), n_hom=len(hom), en_panel=g in PANEL,
                                  mismo_bloque_fase=len(pids) == 1 and all(r["pid"] for r in rs),
                                  variantes=";".join(f'{r["chrom"]}:{r["pos"]}{r["ref"]}>{r["alt"]}({r["gt"]})' for r in rs)))
    genes_rec.sort(key=lambda d: (-d["en_panel"], -d["n_hom"], -d["n_variantes"]))
    emb["4_variantes_en_genes_recesivos_candidatos"] = sum(d["n_variantes"] for d in genes_rec)
    emb["5_genes_recesivos_candidatos"] = len(genes_rec)
    emb["6_genes_recesivos_en_panel"] = sum(d["en_panel"] for d in genes_rec)
    emb["7_variantes_panel_af<5%_cualquier_consecuencia"] = len(panel)

    def tsv(nombre, filas):
        if not filas:
            open(f"{dst}/{nombre}", "w").write("")
            return
        ks = list(filas[0].keys())
        with open(f"{dst}/{nombre}", "w") as f:
            f.write("\t".join(ks) + "\n")
            for r in filas:
                f.write("\t".join(str(r[k]) for k in ks) + "\n")
    tsv("raras_funcionales.tsv", raras)
    tsv("genes_recesivos.tsv", genes_rec)
    tsv("panel_mva.tsv", panel)
    json.dump(dict(sorted(emb.items())), open(f"{dst}/embudo.json", "w"), indent=1)
    for k, v in sorted(emb.items()):
        print(f"{k:55s} {v}")


if __name__ == "__main__":
    main()
