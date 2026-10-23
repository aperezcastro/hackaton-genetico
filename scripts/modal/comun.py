"""comun.py — utilidades compartidas por los trabajos de Modal del pipeline: credenciales desde .env, tope de coste antes de lanzar, vigilante de reloj local y registro de
gasto en scripts/modal/gastos.csv.

Tarifas usadas (Modal, 2026): 0,0000131 $ por núcleo·s, 0,00000222 $ por GiB·s (CPU). GPU H100 ≈ 0,001097 $/s.
"""
import csv
import datetime
import os
import pathlib
import sys
import threading

ROOT = pathlib.Path(__file__).resolve().parents[2]
GASTOS = ROOT / "scripts" / "modal" / "gastos.csv"
CPU_USD_S, MEM_USD_S = 0.0000131, 0.00000222
GPU_USD_S = {"H100": 0.001097, "A100-80GB": 0.000694, "A100": 0.000583, "L40S": 0.000542, "A10G": 0.000306,
             "T4": 0.000164}


def credenciales():
    env = ROOT / ".env"
    for line in env.read_text().splitlines() if env.exists() else []:
        if "=" in line:
            k, v = [x.strip() for x in line.split("=", 1)]
            if k == "modal-token-id":
                os.environ["MODAL_TOKEN_ID"] = v
            if k == "modal-token-secret":
                os.environ["MODAL_TOKEN_SECRET"] = v


def coste(cpu, mem_gib, segundos, gpu=None, n=1):
    s = (cpu * CPU_USD_S + mem_gib * MEM_USD_S + (GPU_USD_S.get(gpu, 0) if gpu else 0)) * segundos * n
    return s


def tope(cpu, mem_gib, timeout_s, max_usd, gpu=None, n=1, label=""):
    peor = coste(cpu, mem_gib, timeout_s, gpu, n)
    print(f"Modal · {label} · {n} contenedor(es) · {cpu} núcleos · {mem_gib} GiB{' · ' + gpu if gpu else ''} · "
          f"tope {timeout_s / 60:.0f} min → máximo {peor:.2f} $")
    if peor > max_usd:
        sys.exit(f"✗ el máximo posible ({peor:.2f} $) pasa de --max-usd {max_usd}")
    return peor


def vigilante(max_s, label):
    def fire():
        print(f"\n✗ VIGILANTE: {label} supera {max_s / 60:.0f} min de reloj — se detiene la app", flush=True)
        os._exit(3)
    th = threading.Timer(max_s, fire)
    th.daemon = True
    th.start()


def apunta(tarea, cpu, mem, segundos_ocupados, segundos_reloj, usd, ok, proyecto="hackaton-genetico"):
    nuevo = not GASTOS.exists()
    with GASTOS.open("a", newline="") as f:
        w = csv.writer(f)
        if nuevo:
            w.writerow(["fecha", "tarea", "proyecto", "cpu", "mem_gib", "render_s", "total_s", "coste_usd_est", "ok"])
        w.writerow([f"{datetime.datetime.now():%Y-%m-%d %H:%M}", tarea, proyecto, cpu, mem, round(segundos_ocupados),
                    round(segundos_reloj), round(usd, 3), ok])


def arg(n, d):
    return type(d)(sys.argv[sys.argv.index(n) + 1]) if n in sys.argv else d
