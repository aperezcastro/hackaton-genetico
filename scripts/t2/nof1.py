"""nof1.py — ¿cómo sabríamos si un fármaco funciona EN ESTE NIÑO? Biomarcador y potencia de un ensayo n-de-1.

El problema que resuelve: todos los informes de reposicionamiento acaban en «habría que probarlo». Con <50 casos en
el mundo no hay ensayo aleatorizado posible, así que la pregunta operativa es: ¿existe una medida que cambie lo
bastante rápido para decidir en semanas, y cuántas medidas hacen falta?

Biomarcador propuesto: **% de reticulocitos micronucleados (MN-RET)** por citometría de flujo.
  · Un micronúcleo es un trozo de cromosoma que queda fuera del núcleo tras un reparto defectuoso: es el producto
    directo del defecto de esta enfermedad, no un marcador indirecto.
  · Los reticulocitos son eritrocitos recién salidos de la médula y circulan 1-2 días, así que la medida es una
    TASA DE ERRORES RECIENTE, no el acumulado de años que da el cariotipo de linfocitos.
  · Es sangre periférica (una gota), validado en humanos (Dertinger 2003, PMID 14644356; plataforma de citometría
    PMID 21164195) y usado en niños (EII pediátrica, PMID 32761646: límite superior del intervalo de tolerancia
    del 90 % en controles = 0,42 %).
  · Nunca se ha medido en MVA. Es, a la vez, medida de BENEFICIO (debería bajar si el fármaco reduce los errores)
    y de DAÑO (subiría con un fármaco que los añade, p. ej. la rapamicina, PMID 9914383).

Modelo de medida (todo público, ningún dato del paciente):
  visita j, condición x_j ∈ {0=off, 1=on}:  tasa_j = λ0 · exp(θ·x_j + b_j),  b_j ~ N(0, σ²)
  recuento y_j ~ Binomial(N_ret, tasa_j)
  σ = variabilidad biológica entre visitas (log); θ = log del cociente de tasas (efecto del fármaco).
Análisis: modelo normal en log(tasa) con σ ESTIMADA de los propios datos (posterior t de θ con priores planos;
la cuadratura de Gauss-Hermite que se usó antes se descartó porque inflaba la potencia con N grande).
Criterio de respuesta PREDECLARADO: P(reducción ≥ `umbral`) ≥ `certeza`.

    .venv/bin/python scripts/t2/nof1.py [--rapido]

Salida: reports/t2/nof1.md (+ datos para la figura t2_nof1).
"""
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "t2"
RAPIDO = "--rapido" in sys.argv

# ── Supuestos, todos declarados y variados en la sensibilidad ────────────────────────────────────────────────
NORMAL_ALTO = 0.0042      # 0,42 %: techo del 90 % en controles pediátricos (PMID 32761646)
BASAL = 0.02              # 2 %: basal SUPUESTO en MVA (~5x el techo normal). Nunca medido: se varía abajo.
SIGMA = 0.35              # variabilidad biológica entre visitas (CV log ~35 %); se varía abajo
N_RET = 20_000            # reticulocitos CD71+ por muestra, como el estudio pediátrico de referencia (PMID 32761646)
VISITAS = 4               # visitas por condición (diseño on/off/on/off equilibrado)
UMBRAL = 0.25             # regla estricta: «responde» = reducción relativa >= 25 % ...
CERTEZA = 0.95            # ... con probabilidad posterior >= 95 %
N_SIM = 1000 if RAPIDO else 4000
SEMILLA = 2026            # números aleatorios COMUNES para todas las celdas: las diferencias entre configuraciones
                          # no las crea el azar de cada celda (antes cada celda tenía su semilla)

# ── Modelo de análisis ──────────────────────────────────────────────────────────────────────────────────────
#   z_j = log(y_j / N) = α + θ·x_j + ε_j,   Var(ε_j) = σ² + 1/y_j ≈ τ² (común)
# σ NO se supone conocida (corregido el 4 oct tras la revisión independiente): se estima con los mismos datos.
# Con priores no informativos en (α, θ, log τ), el posterior de θ es una t de Student con n − 2 grados de
# libertad centrada en la estimación por mínimos cuadrados (resultado estándar). Tratar 1/y_j como parte de una
# varianza común es exacto cuando el recuento es grande y conservador cuando no lo es.


def posterior_theta(y, n, x):
    """Media, error estándar y grados de libertad del posterior t de θ. y: (sim, J); x: (J,)."""
    y = np.maximum(np.asarray(y), 1)
    z = np.log(y / n)
    on, off = x == 1, x == 0
    theta = z[..., on].mean(-1) - z[..., off].mean(-1)
    resid = np.concatenate([z[..., on] - z[..., on].mean(-1, keepdims=True),
                            z[..., off] - z[..., off].mean(-1, keepdims=True)], axis=-1)
    gl = x.size - 2
    tau2 = (resid ** 2).sum(-1) / gl
    se = np.sqrt(tau2 * (1 / on.sum() + 1 / off.sum()))
    return theta, se, gl


def p_reduccion(theta, se, gl, umbral=None):
    """P(reducción >= umbral) bajo el posterior t. umbral None = P(cualquier reducción)."""
    from scipy.stats import t as student
    corte = np.log(1 - umbral) if umbral else 0.0
    return student.cdf((corte - theta) / se, gl)


def _base(n_sim, J):
    """Ruido biológico y uniforme de recuento COMUNES a todas las celdas con el mismo número de visitas."""
    rng = np.random.default_rng(SEMILLA + J)
    return rng.standard_normal((n_sim, 2 * J)), rng.random((n_sim, 2 * J))


def simula(efecto, sigma=SIGMA, n_ret=N_RET, visitas=VISITAS, basal=BASAL, n_sim=None,
           umbral=UMBRAL, certeza=CERTEZA, devolver_se=False):
    """Fracción de simulaciones que declaran respuesta (efecto 0 → tasa de falsos positivos). Con devolver_se,
    también el error de Monte Carlo de esa fracción."""
    from scipy.stats import binom
    n_sim = n_sim or N_SIM
    x = np.array([0, 1] * visitas, float)
    zb, u = _base(n_sim, visitas)
    tasa = np.clip(basal * np.exp(np.log(1 - efecto) * x + sigma * zb), 1e-9, 0.5)
    y = binom.ppf(u, n_ret, tasa)                     # recuentos binomiales con la misma uniforme en cada celda
    theta, se, gl = posterior_theta(y, n_ret, x)
    p = float((p_reduccion(theta, se, gl, umbral) >= certeza).mean())
    return (p, (p * (1 - p) / n_sim) ** 0.5) if devolver_se else p


def efecto_minimo_detectable(sigma=SIGMA, n_ret=N_RET, visitas=VISITAS, basal=BASAL, umbral=UMBRAL,
                             certeza=CERTEZA, potencia=0.8):
    """Reducción verdadera más pequeña que se detecta con la potencia pedida (búsqueda binaria)."""
    lo, hi = 0.0, 0.95
    for _ in range(14):
        med = (lo + hi) / 2
        if simula(med, sigma, n_ret, visitas, basal, N_SIM, umbral, certeza) < potencia:
            lo = med
        else:
            hi = med
    return (lo + hi) / 2


def visitas_necesarias(efecto, sigma, umbral=UMBRAL, certeza=CERTEZA, potencia=0.8, n_ret=N_RET, maximo=30):
    """Visitas por condición (on/off) necesarias para detectar `efecto` con la potencia pedida. None si > maximo."""
    for v in range(2, maximo + 1):
        if simula(efecto, sigma=sigma, n_ret=n_ret, visitas=v, umbral=umbral, certeza=certeza) >= potencia:
            return v
    return None


REGLAS = [("reduction >= 25 %, 95 % certainty", UMBRAL, CERTEZA),
          ("any reduction, 95 % certainty", None, 0.95),
          ("any reduction, 90 % certainty", None, 0.90)]


def main():
    res = {"supuestos": {"normal_alto": NORMAL_ALTO, "basal_mva_supuesta": BASAL, "sigma": SIGMA, "n_ret": N_RET,
                         "visitas_por_condicion": VISITAS, "n_sim": N_SIM, "sigma_conocida": False}}
    print("· curva de potencia por regla (σ estimada de los datos)")
    efectos = [0.0, 0.15, 0.25, 0.35, 0.50, 0.65, 0.80]
    res["potencia_efecto"] = {regla: {f"{e:.2f}": simula(e, umbral=u, certeza=c, devolver_se=True) for e in efectos}
                              for regla, u, c in REGLAS}
    print("· mandos del diseño (35 % de reducción real, regla convencional, números aleatorios comunes)")
    conv = dict(umbral=None, certeza=0.95, devolver_se=True)
    res["potencia_visitas"] = {f"{v}": simula(0.35, visitas=v, **conv) for v in (2, 4, 6, 8, 12)}
    res["potencia_sigma"] = {f"{s:.2f}": simula(0.35, sigma=s, **conv) for s in (0.20, 0.35, 0.50)}
    res["potencia_regla"] = {regla: simula(0.35, umbral=u, certeza=c, devolver_se=True) for regla, u, c in REGLAS}
    print("· rejilla basal × células puntuadas (el «contar más células no ayuda» depende del basal)")
    res["rejilla"] = {f"{b:.4f}|{n}": simula(0.35, basal=b, n_ret=n, **conv)
                      for b in (0.001, 0.0042, 0.01, 0.02, 0.05) for n in (10_000, 20_000, 100_000, 500_000)}
    print("· efecto mínimo detectable (80 % de potencia)")
    res["mde"] = {f"{regla} | {v} visits": efecto_minimo_detectable(visitas=v, umbral=u, certeza=c)
                  for regla, u, c in REGLAS[:2] for v in (4, 8)}
    print("· tabla de diseño: visitas por condición según la σ que se mida en el peldaño 0")
    res["fp_regla"] = {regla: simula(0.0, visitas=8, umbral=u, certeza=c, devolver_se=True) for regla, u, c in REGLAS}
    res["diseno"] = {f"{regla} | sigma {s:.2f} | effect {e:.2f}": visitas_necesarias(e, s, u, c)
                     for regla, u, c in REGLAS for s in (0.20, 0.35, 0.50) for e in (0.35, 0.50)}
    for k in ("potencia_regla", "mde", "fp_regla"):
        print("  ", k, {kk: (round(vv[0], 3) if isinstance(vv, tuple) else round(vv, 3)) for kk, vv in res[k].items()})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "nof1.json").write_text(json.dumps(res, indent=1))
    escribe_md(res)
    print("✓ reports/t2/nof1.md")


def _pc(v):
    p, se = v
    return f"{100 * p:.1f} % (± {100 * se:.1f})"


def escribe_md(r):
    s = r["supuestos"]
    reglas = [x[0] for x in REGLAS]
    fp = r["fp_regla"]
    fp_txt = {k: ("< 1 %" if v[0] < 0.01 else f"{100 * v[0]:.0f} %") for k, v in fp.items()}
    l = ["# How we would know whether it works in this child", "",
         "Generated by `scripts/t2/nof1.py`. Public parameters only; no patient data. The endpoint is the "
         "**frequency of micronucleated reticulocytes (%MN-RET)** measured by flow cytometry on peripheral blood.", "",
         "## Why this endpoint", "",
         "| Property | Value for this disease |", "|---|---|",
         "| What it counts | micronuclei — chromosome fragments left outside the nucleus by a missegregation "
         "event, i.e. the direct product of the defect rather than a proxy for it |",
         "| Time window | reticulocytes circulate 1–2 days, so this is a **recent error rate**, unlike a "
         "lymphocyte karyotype, which accumulates over years |",
         "| Sample | 240 µL of whole blood; flow cytometry scores ≥ 20,000 CD71+ reticulocytes (PMID 32761646) |",
         "| Human validation | Dertinger 2003 (PMID 14644356); flow platform PMID 21164195; rose within 3–5 days "
         "of starting chemo- or radiotherapy in 9 of 10 adults (PMID 17059891) |",
         "| Paediatric reference | upper bound of the 90 % tolerance interval in paediatric controls = "
         f"**{100 * s['normal_alto']:.2f} %** (PMID 32761646) |",
         "| Measured in MVA | **never** — rung 0, worth doing whether or not any drug is ever given |",
         "| Direction of information | falls if a drug reduces missegregation; **rises if a drug adds it**, which is "
         "why it doubles as a safety monitor (rapamycin is aneugenic in human lymphocytes in vitro, PMID 9914383) |", "",
         "## Model and decision rules", "",
         f"Alternating on/off blocks; {s['n_ret']:,} reticulocytes per sample; assumed MVA baseline "
         f"{100 * s['basal_mva_supuesta']:.1f} % (about {s['basal_mva_supuesta'] / s['normal_alto']:.0f}× the paediatric "
         f"ceiling — **an assumption, never measured**), so the expected count is about "
         f"{s['basal_mva_supuesta'] * s['n_ret']:.0f} micronucleated cells per sample and log(y/N) is close to normal. "
         "**The between-visit variability σ is not assumed known**: it is estimated from the same data, so the "
         "posterior for the log rate ratio is a Student t with n − 2 degrees of freedom (flat priors on the mean, the "
         f"effect and log σ). True σ in the simulations = {s['sigma']:.2f}. Every cell uses the same random numbers "
         f"(common random numbers, {s['n_sim']:,} simulated trials), and each figure carries its Monte Carlo standard "
         "error.", "",
         "| Responder rule | False-positive rate (8 visits per condition) |", "|---|---|"]
    l += [f"| {k} | {_pc(fp[k])} |" for k in reglas]
    l += ["", "The conventional rule (any reduction with 95 % posterior certainty) is calibrated: its false-positive "
          "rate is the advertised 5 %. The strict rule buys an almost zero false-positive rate at a large cost in "
          "power, because its threshold sits close to the effects one can hope for.", "",
          "## Power", "",
          "| True reduction | " + " | ".join(reglas) + " |", "|---|" + "---|" * len(reglas)]
    for e in r["potencia_efecto"][reglas[0]]:
        l.append(f"| {100 * float(e):.0f} %{' (false positives)' if float(e) == 0 else ''} | " +
                 " | ".join(_pc(r["potencia_efecto"][k][e]) for k in reglas) + " |")
    l += ["", "**Minimum detectable effect (80 % power):**", ""]
    l += [f"- {k}: **{100 * v:.0f} % reduction**" for k, v in r["mde"].items()]
    l += ["", "## What the design is limited by (true 35 % reduction, conventional rule)", "",
          "| Knob | Setting | Power |", "|---|---|---|"]
    l += [f"| visits per condition | {k} | {_pc(v)} |" for k, v in r["potencia_visitas"].items()]
    l += [f"| between-visit variability σ | {k} | {_pc(v)} |" for k, v in r["potencia_sigma"].items()]
    l += [f"| decision rule | {k} | {_pc(v)} |" for k, v in r["potencia_regla"].items()]
    # rejilla basal × células
    basales = sorted({float(k.split("|")[0]) for k in r["rejilla"]})
    celdas = sorted({int(k.split("|")[1]) for k in r["rejilla"]})
    l += ["", "## Counting more cells only helps when the baseline is low", "",
          "Power at a true 35 % reduction (conventional rule, 4 visits per condition), by baseline %MN-RET and "
          "reticulocytes scored per sample. The counting term 1/y is negligible beside σ² = 0.12 when the expected "
          "count is in the hundreds, and matters only near or below the paediatric reference.", "",
          "| Baseline %MN-RET | " + " | ".join(f"{n:,} cells" for n in celdas) + " |", "|---|" + "---|" * len(celdas)]
    for b in basales:
        l.append(f"| {100 * b:.2f} % | " + " | ".join(_pc(r["rejilla"][f"{b:.4f}|{n}"]) for n in celdas) + " |")
    l += ["", "## Design table: what to do with the rung-0 result", "",
          "Visits per condition (on/off) needed for 80 % power, 20,000 CD71+ reticulocytes per sample. **σ is the "
          "between-visit variability measured at rung 0**; plan with the upper bound of its interval, not the point "
          "estimate. «—» means more than 30 visits per condition, i.e. not feasible.", "",
          "| Decision rule (false-positive rate) | σ measured | 35 % true reduction | 50 % true reduction |",
          "|---|---|---|---|"]
    for k in reglas:
        for sg in (0.20, 0.35, 0.50):
            f = lambda e: r["diseno"][f"{k} | sigma {sg:.2f} | effect {e:.2f}"]
            g = lambda v: "—" if v is None else str(v)
            l.append(f"| {k} ({fp_txt[k]}) | {sg:.2f} | {g(f(0.35))} | {g(f(0.50))} |")
    l += ["", "## Per-line biomarkers, so every hypothesis has a stopping rule", "",
          "| Line | Biomarker that decides it | Window |", "|---|---|---|",
          "| H1 NAD+/SIRT2 | BUBR1 immunoblot ± cycloheximide chase in patient cells — also decides allele "
          "scenario A vs B | days, before any dosing |",
          "| H3 mTORC1 | phospho-S6 in PBMC or muscle, lean mass, grip strength, **plus %MN-RET as a safety readout** "
          "| weeks |",
          "| H4 senolysis | p16INK4A in T cells, circulating SASP factors | weeks |",
          "| H6 interferon/JAK | interferon-stimulated gene score in blood — measure it first; no MVA value is "
          "published | weeks |",
          "| H7 proteostasis | folded/total BUBR1 **and** aneuploid fraction, which must not rise | cells first |",
          "| Every line | %MN-RET: benefit if it falls, harm if it rises | 2–4 weeks |"]
    (OUT / "nof1.md").write_text("\n".join(l) + "\n")


if __name__ == "__main__":
    main()
