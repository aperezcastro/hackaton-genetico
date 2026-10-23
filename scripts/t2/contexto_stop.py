"""contexto_stop.py — ¿es el alelo nonsense del caso un buen candidato a «lectura a través» (read-through)?

La eficiencia del read-through depende sobre todo de (i) la identidad del codón de parada (UGA > UAG > UAA) y
(ii) el nucleótido +4, el que sigue al codón (preferencia C > U > A/G), además de la longitud de la 3'UTR
(Mangkalaphiban et al., Nat Commun 2024, PMID 38509072).

Este script descarga la CDS pública de referencia de BUB1B (Ensembl, dato público a nivel gen), aplica el cambio
c.2210T>G y calcula el contexto. NO lee ningún fichero del paciente.

    .venv/bin/python scripts/t2/contexto_stop.py [--cds c.2210T>G]

Salida: identidad del stop, base +4, tetranucleótido y veredicto cualitativo.
"""
import pathlib
import re
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2]
CDS = ROOT / "refs" / "estructura" / "bub1b_cds.txt"
URL = "https://rest.ensembl.org/sequence/id/ENST00000287598?type=cds;content-type=text/plain"
CAMBIO = sys.argv[sys.argv.index("--cds") + 1] if "--cds" in sys.argv else "c.2210T>G"
# Preferencia del nucleótido +4 (Mangkalaphiban 2024): C mejor, A/G peores
PREF_4 = {"C": "favourable", "T": "intermediate", "A": "unfavourable", "G": "unfavourable"}
PREF_STOP = {"TGA": "most permissive", "TAG": "intermediate", "TAA": "least permissive"}


def cds():
    if not CDS.exists():
        CDS.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(URL, CDS)
    return CDS.read_text().strip().replace("\n", "")


def main():
    s = cds()
    pos, ref, alt = re.match(r"c\.(\d+)([ACGT])>([ACGT])", CAMBIO).groups()
    pos = int(pos)
    if s[pos - 1] != ref:
        sys.exit(f"✗ la CDS tiene {s[pos - 1]} en c.{pos}, no {ref}: ¿transcrito equivocado?")
    mut = s[:pos - 1] + alt + s[pos:]
    cod_i = (pos - 1) // 3                      # índice 0-based del codón afectado
    stop, mas4 = mut[cod_i * 3:cod_i * 3 + 3], mut[cod_i * 3 + 3]
    if stop not in PREF_STOP:
        sys.exit(f"✗ {CAMBIO} no genera un codón de parada (codón resultante {stop})")
    total = len(s) // 3
    print(f"BUB1B ENST00000287598 · CDS {len(s)} nt · {total} codones")
    print(f"{CAMBIO}: codón {cod_i + 1} {s[cod_i * 3:cod_i * 3 + 3]} → {stop} (p.{cod_i + 1}Ter)")
    print(f"  contexto: {' '.join(mut[(i - 1) * 3:i * 3] for i in range(cod_i - 1, cod_i + 4))}")
    print(f"  identidad del stop : {stop} ({PREF_STOP[stop]})")
    print(f"  nucleótido +4      : {mas4} ({PREF_4[mas4]})")
    print(f"  tetranucleótido    : {stop}{mas4}")
    # NMD depende de la distancia al ÚLTIMO EMPALME exón-exón (regla de 50-55 nt), no al final de la CDS; aquí no
    # se cargan los exones, así que solo se informa de los codones restantes. El stop está en el exón 17 de 23
    # (VEP, Ensembl 116), muy por delante del último empalme → NMD esperado (ver el informe del T1).
    print(f"  codones hasta el final de la CDS: {total - cod_i - 1}")
    bueno = stop == "TGA" and PREF_4[mas4] == "favourable"
    print(f"\n  veredicto: {'contexto favorable' if bueno else 'contexto MIXTO/POBRE'} para read-through")


if __name__ == "__main__":
    main()
