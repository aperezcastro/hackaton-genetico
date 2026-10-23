"""estructura.py — contexto estructural de los alelos BUBR1 en el modelo AlphaFold público (AF-O60566-F1, v6).

Descarga la estructura (no envía ninguna variante a ningún servicio) y calcula, para el residuo del caso y los alelos
MVA1 conocidos (UniProt): pLDDT, accesibilidad relativa al disolvente (Shrake–Rupley; máx. de Tien et al. 2013),
contactos < 4,5 Å y distancias Cα.

    .venv/bin/python scripts/t2/estructura.py
"""
import pathlib
import statistics
import urllib.request

from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

ROOT = pathlib.Path(__file__).resolve().parents[2]
PDB = ROOT / "refs" / "estructura" / "bubr1_af.pdb"
URL = "https://alphafold.ebi.ac.uk/files/AF-O60566-F1-model_v6.pdb"
CASO, MVA1 = [737, 1002], [814, 844, 909, 921, 1012]      # residuos del caso; missense MVA1 de UniProt
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
          "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
          "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

if not PDB.exists():
    PDB.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(URL, PDB)
modelo = PDBParser(QUIET=True).get_structure("bubr1", PDB)[0]
ch = modelo["A"]
ShrakeRupley().compute(modelo, level="R")
ns = NeighborSearch(list(modelo.get_atoms()))

print("residuo  pLDDT  RSA   contactos<4.5Å (no contiguos)")
for i in CASO + MVA1:
    r = ch[i]
    cont = sorted({b.get_parent().id[1] for a in r for b in ns.search(a.coord, 4.5) if abs(b.get_parent().id[1] - i) > 2})
    print(f"{r.get_resname()}{i:<5} {r['CA'].get_bfactor():5.1f}  {r.sasa / MAXASA[r.get_resname()]:.2f}  {cont}")
print("\nCα N1002 → alelos MVA1:", {j: round(ch[1002]["CA"] - ch[j]["CA"], 1) for j in MVA1})
print("pLDDT medio dominio quinasa (766–1050):", round(statistics.mean(ch[i]["CA"].get_bfactor() for i in range(766, 1051)), 1))
