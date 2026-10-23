"""estilo.py — estilo común de las figuras de los informes (matplotlib → PDF vectorial para LaTeX).

Paleta categórica validada (skill dataviz, `validate_palette.js --mode light`: todas las comprobaciones PASS; aqua
y magenta < 3:1 de contraste → siempre con etiqueta directa). Tipografía IBM Plex, la misma del informe y del vídeo.
Reglas: marcas finas, ejes y rejilla discretos, texto en tinta de texto (nunca en el color de la serie), etiquetas
directas en lugar de leyendas cuando hay ≤ 4 series, un solo eje Y por gráfico.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
FIG = ROOT / "reports" / "latex" / "fig"
FUENTES = ROOT / "reports" / "latex" / "fonts"

# paleta categórica (orden fijo, nunca ciclado) y tintas
AZUL, NARANJA, AQUA, AMARILLO, MAGENTA, VERDE, VIOLETA, ROJO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
TEXTO, TEXTO2, TENUE = "#0b0b0b", "#52514e", "#8a8984"
GRIS, GRIS_CLARO, SUPERFICIE = "#b9b8b2", "#e4e3df", "#fcfcfb"
CRITICO, BUENO, AVISO = "#d03b3b", "#0ca30c", "#fab219"          # estado: solo con icono/etiqueta
RAMPA = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95"]  # ordinal (≥ 250)

ANCHO = 6.7          # pulgadas: ancho de texto del informe A4
ANCHO_MEDIO = 3.25


def aplica():
    for f in FUENTES.glob("IBMPlexSans-*.ttf"):
        font_manager.fontManager.addfont(str(f))
    for f in FUENTES.glob("IBMPlexMono-*.ttf"):
        font_manager.fontManager.addfont(str(f))
    plt.rcParams.update({
        "font.family": "IBM Plex Sans", "font.size": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "semibold",
        "axes.titlelocation": "left", "axes.labelsize": 8.5, "axes.labelcolor": TEXTO2,
        "axes.edgecolor": GRIS, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": TEXTO2, "ytick.color": TEXTO2, "xtick.labelsize": 7.8, "ytick.labelsize": 7.8,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.grid": False, "grid.color": GRIS_CLARO, "grid.linewidth": 0.5,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "legend.frameon": False, "legend.fontsize": 7.8, "pdf.fonttype": 42, "svg.fonttype": "none",
        "text.color": TEXTO, "lines.linewidth": 1.4,
    })


def guarda(fig, nombre):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{nombre}.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(FIG / f"{nombre}.png", bbox_inches="tight", pad_inches=0.03, dpi=200)   # vista previa
    plt.close(fig)
    print(f"✓ fig/{nombre}.pdf")


def barras_redondeadas(ax, horizontal=True):
    """Extremos de datos ligeramente redondeados (spec de marcas): matplotlib no lo hace nativo; usamos joinstyle."""
    for p in ax.patches:
        p.set_joinstyle("round")
