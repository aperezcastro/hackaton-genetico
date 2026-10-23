"""compila_informes.py — regenera las figuras, el anexo del formulario de métodos y los dos informes en LaTeX.

    .venv/bin/python scripts/compila_informes.py [--sin-figuras]

Salida: reports/latex/track1.pdf y reports/latex/track2.pdf (LuaLaTeX, tipografía IBM Plex incluida en
reports/latex/fonts/, OFL). Las figuras T1 necesitan los agregados de results/t1/ (solo existen donde están los datos
del caso); las del T2 solo usan datos públicos cacheados.
"""
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
TEX = ROOT / "reports" / "latex"
PY = str(ROOT / ".venv" / "bin" / "python")


def paso(cmd, cwd=ROOT, etiqueta=""):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"✗ {etiqueta}: {(r.stdout + r.stderr)[-1500:]}")
    print(f"✓ {etiqueta}")


def main():
    if "--sin-figuras" not in sys.argv:
        paso([PY, "scripts/figuras/figuras_t1.py"], etiqueta="figuras T1")
        paso([PY, "scripts/figuras/figuras_t2.py"], etiqueta="figuras T2")
    paso([PY, "scripts/formulario_metodos.py"], etiqueta="formulario de métodos (anexo A)")
    for t in ("track1", "track2"):
        paso(["latexmk", "-lualatex", "-interaction=nonstopmode", "-halt-on-error", f"{t}.tex"], cwd=TEX,
             etiqueta=f"LaTeX {t}")
        log = (TEX / f"{t}.log").read_text(errors="ignore")
        faltan = sorted({l for l in log.splitlines() if "Missing character" in l})
        if faltan:
            print("  ⚠ glifos ausentes:", *faltan[:5], sep="\n    ")
        indef = [l for l in log.splitlines() if "undefined" in l.lower() and "Reference" in l]
        if indef:
            print("  ⚠ referencias indefinidas:", *indef[:5], sep="\n    ")


if __name__ == "__main__":
    main()
