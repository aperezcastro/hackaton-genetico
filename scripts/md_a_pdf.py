"""md_a_pdf.py — convierte un informe Markdown en PDF (Markdown → HTML con estilo → Chrome sin interfaz).

    .venv/bin/python scripts/md_a_pdf.py entrada.md salida.pdf [--compacto]

  --compacto  márgenes y letra más pequeños (fichas de una página, como la tarjeta de seguridad)
"""
import pathlib
import subprocess
import sys
import tempfile

import markdown

CSS = """
@page { size: A4; margin: 16mm 14mm; }
body { font-family: 'DejaVu Sans', Arial, sans-serif; font-size: 9.6pt; line-height: 1.42; color: #111; }
h1 { font-size: 17pt; border-bottom: 2px solid #1f4e79; padding-bottom: 4px; color: #1f4e79; }
h2 { font-size: 13pt; color: #1f4e79; margin-top: 18px; border-bottom: 1px solid #ccd; }
h3 { font-size: 11pt; color: #333; }
table { border-collapse: collapse; width: 100%; margin: 8px 0; font-size: 8.4pt; page-break-inside: auto; }
th, td { border: 1px solid #bbb; padding: 3px 5px; vertical-align: top; text-align: left; }
th { background: #e8eef5; }
tr { page-break-inside: avoid; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 8.4pt; background: #f3f3f3; padding: 0 2px; }
pre { background: #f6f6f6; border: 1px solid #ddd; padding: 6px; font-size: 7.8pt; white-space: pre-wrap; }
blockquote { border-left: 3px solid #1f4e79; margin: 6px 0; padding: 2px 10px; color: #333; background: #f5f8fb; }
hr { border: 0; border-top: 1px solid #ccc; }
"""

COMPACTO = """
@page { size: A4; margin: 10mm 10mm; }
body { font-size: 8.4pt; line-height: 1.3; }
h1 { font-size: 14pt; margin: 0 0 6px; } table { font-size: 7.6pt; margin: 5px 0; } th, td { padding: 2px 4px; }
p { margin: 4px 0; }
"""
src, dst = map(pathlib.Path, sys.argv[1:3])
if "--compacto" in sys.argv:
    CSS += COMPACTO
texto = src.read_text()
html = markdown.markdown(texto, extensions=["tables", "fenced_code", "sane_lists"])
# título del PDF = primer encabezado (si no, Chrome pone el nombre del HTML temporal)
titulo = next((l.lstrip("# ").strip().replace("*", "") for l in texto.splitlines() if l.startswith("#")), src.stem)
with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False) as f:
    f.write(f"<!doctype html><html><head><meta charset='utf-8'><title>{titulo}</title><style>{CSS}</style></head>"
            f"<body>{html}</body></html>")
subprocess.run(["google-chrome", "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                f"--print-to-pdf={dst.resolve()}", f"file://{f.name}"], check=True, capture_output=True)
print(f"✓ {dst}")
