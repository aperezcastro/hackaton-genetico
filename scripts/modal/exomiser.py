"""exomiser.py — prioriza el genoma del caso con Exomiser 15.1 (datos 2512, hg38) en Modal, combinando variantes y
fenotipo (términos HPO del caso, results/t1/hpo.txt). Dos análisis:

  A · «estándar»: preset de genoma con priorityScoreFilter (hiPHIVE ≥ 0,501): el ranking clásico.
  B · «abierto»: igual pero SIN filtro de fenotipo, para no perder genes aún no asociados a la enfermedad
      (el sesgo de «buscar bajo la farola»).

Herencias: AR hom-alt (0,1 %), AR comp-het (2 %), AD (0,1 %), X (como la plantilla), MT 0,2 %.
Patogenicidad: REVEL, MVP, AlphaMissense, SpliceAI (incluidos en los datos 2512). Solo PASS.

    .venv/bin/python scripts/modal/exomiser.py [--max-usd 1] [--timeout 60]

Baja a data/derivados/ las salidas completas de Exomiser (exomiser_{A,B}.genes.tsv, .variants.tsv…): son tablas
a escala genómica, así que van a la carpeta que el agente no puede leer. `scripts/exomiser_top.py` extrae de ellas
los 12 primeros puestos (results/t1/exomiser_top.json) para la figura del informe.
"""
import pathlib
import sys
import time

import modal

if modal.is_local():
    import comun
    comun.credenciales()
arg = lambda n, d: type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
TIMEOUT_MIN = arg("--timeout", 60)
MAX_USD = arg("--max-usd", 1.0)
CPU, MEM = 4.0, 24
ROOT = pathlib.Path(__file__).resolve().parents[2] if modal.is_local() else pathlib.Path("/")
VCF = "WGS_EX2312012_HGWCNDSX7.vcf.gz"

exo = modal.Volume.from_name("mva-exomiser")
datos = modal.Volume.from_name("mva-data")
# Exomiser 15 está compilado para Java 21 (class file 65): Debian slim solo trae 17 → imagen Temurin 21
image = (modal.Image.from_registry("eclipse-temurin:21-jre-jammy", add_python="3.12")
         .apt_install("tabix").pip_install("pyyaml"))
app = modal.App("mva-exomiser", image=image)

PROPS = """
exomiser.data-directory=/ex/data
exomiser.hg38.data-version=2512
exomiser.phenotype.data-version=2512
exomiser.hg38.clin-var-data-version=2512
logging.level.org.monarchinitiative.exomiser=WARN
"""


def analisis(abierto):
    """Análisis en formato Exomiser ≥ 14: solo el bloque de análisis; muestra (fenopacket) y VCF van aparte."""
    pasos = [
        {"hiPhivePrioritiser": {}},
        {"geneBlacklistFilter": {}},
        {"failedVariantFilter": {}},
        {"regulatoryFeatureFilter": {}},
        {"frequencyFilter": {"maxFrequency": 2.0}},
        {"pathogenicityFilter": {"keepNonPathogenic": True}},
        {"inheritanceFilter": {}},
        {"omimPrioritiser": {}},
    ]
    if not abierto:
        pasos.insert(1, {"priorityScoreFilter": {"priorityType": "HIPHIVE_PRIORITY", "minPriorityScore": 0.501}})
    return {
        "inheritanceModes": {"AUTOSOMAL_DOMINANT": 0.1, "AUTOSOMAL_RECESSIVE_HOM_ALT": 0.1,
                             "AUTOSOMAL_RECESSIVE_COMP_HET": 2.0, "X_DOMINANT": 0.1, "X_RECESSIVE_HOM_ALT": 0.1,
                             "X_RECESSIVE_COMP_HET": 2.0, "MITOCHONDRIAL": 0.2},
        "analysisMode": "PASS_ONLY",
        "frequencySources": ["UK10K", "GNOMAD_E_AFR", "GNOMAD_E_AMR", "GNOMAD_E_ASJ", "GNOMAD_E_EAS", "GNOMAD_E_FIN",
                             "GNOMAD_E_NFE", "GNOMAD_E_SAS", "GNOMAD_G_AFR", "GNOMAD_G_AMR", "GNOMAD_G_ASJ",
                             "GNOMAD_G_EAS", "GNOMAD_G_FIN", "GNOMAD_G_NFE", "GNOMAD_G_SAS"],
        "pathogenicitySources": ["REVEL", "MVP", "ALPHA_MISSENSE", "SPLICE_AI"],
        "steps": pasos,
    }


def muestra(hpo):
    return {"id": "WGS_EX2312012", "subject": {"id": "WGS_EX2312012"},  # nombre de la muestra en el VCF
            "phenotypicFeatures": [{"type": {"id": h}} for h in hpo]}


@app.function(volumes={"/ex": exo, "/data": datos}, cpu=CPU, memory=MEM * 1024, timeout=TIMEOUT_MIN * 60)
def corre(hpo: list):
    import glob
    import os
    import subprocess
    import yaml
    t0 = time.time()
    # contigs sin «chr» en el VCF; Exomiser los acepta. Copia local (lectura rápida)
    subprocess.run(["cp", f"/data/{VCF}", "/tmp/caso.vcf.gz"], check=True)
    subprocess.run(["cp", f"/data/{VCF}.tbi", "/tmp/caso.vcf.gz.tbi"], check=True)
    open("/tmp/application.properties", "w").write(PROPS)
    os.makedirs("/tmp/out", exist_ok=True)
    logs, files = "", {}
    yaml.safe_dump(muestra(hpo), open("/tmp/muestra.yml", "w"), sort_keys=False)
    for nombre, abierto in (("exomiser_A", False), ("exomiser_B", True)):
        yaml.safe_dump(analisis(abierto), open(f"/tmp/{nombre}.yml", "w"), sort_keys=False)
        r = subprocess.run(["java", f"-Xmx{MEM - 4}g", "-Dspring.config.location=/tmp/application.properties", "-jar", "/ex/exomiser-cli-15.1.0/exomiser-cli-15.1.0.jar",
                            "analyse", "--sample", "/tmp/muestra.yml", "--analysis", f"/tmp/{nombre}.yml",
                            "--vcf", "/tmp/caso.vcf.gz", "--assembly", "GRCh38",
                            "--output-directory", "/tmp/out", "--output-filename", nombre,
                            "--output-format", "TSV_GENE,TSV_VARIANT"],
                           capture_output=True, text=True, cwd="/ex/exomiser-cli-15.1.0")
        o = r.stdout + r.stderr
        logs += f"== {nombre}: rc={r.returncode}\n" + (o if len(o) < 4000 else o[:2000] + "\n[…]\n" + o[-2000:])
    for f in glob.glob("/tmp/out/*"):
        files[os.path.basename(f)] = open(f, "rb").read()
    return {"log": logs, "files": files, "secs": time.time() - t0}


if __name__ == "__main__":
    hpo = [l.split("\t")[0] for l in (ROOT / "results" / "t1" / "hpo.txt").read_text().splitlines() if l.startswith("HP:")]
    comun.tope(CPU, MEM, TIMEOUT_MIN * 60, MAX_USD, label="exomiser A+B")
    comun.vigilante(TIMEOUT_MIN * 60 + 300, "exomiser")
    t0 = time.time()
    with modal.enable_output(), app.run():
        r = corre.remote(hpo)
    print(r["log"])
    dst = ROOT / "data" / "derivados"
    dst.mkdir(parents=True, exist_ok=True)
    for f, b in r["files"].items():
        (dst / f).write_bytes(b)
    usd = comun.coste(CPU, MEM, r["secs"])
    print(f"✓ {len(r['files'])} ficheros → data/derivados/ · {r['secs'] / 60:.1f} min · ~{usd:.3f} $")
    comun.apunta("exomiser A+B", CPU, MEM, r["secs"], time.time() - t0, usd, bool(r["files"]))
