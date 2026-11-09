#!/usr/bin/env python3
# ===========================================================================
# ANALYSIS: first Mars FDTD — plane P wave through TAYAK crust
#
# Compares arrivals in the synthetic traces against 1-D ray predictions
# for vertical incidence: direct P, Ps conversions at 77.4 / 10 / 1 km,
# PPmP Moho multiple. Primary metric: relative delays (phase - direct P),
# measured on envelope maxima; the wavelet-center bias cancels because all
# phases share the source wavelet. Stdlib only. Writes picks.csv +
# section.png (hand-rolled PNG writer).
# ===========================================================================
import csv, math, struct, zlib, os

# ---- TAYAK layers, grid-quantised (identical to config/env_mars.jl) -------
LAYERS = [  # (d_top_km, d_bot_km, vp, vs)
    (0.0,   1.0,   3677.71, 1739.80),
    (1.0,   10.0,  4952.25, 2780.97),
    (10.0,  77.4,  5846.66, 3281.16),
    (77.4,  80.0,  7400.90, 4246.33),
    (80.0,  100.0, 7413.64, 4247.35),
    (100.0, 110.0, 7436.12, 4180.19),
    (110.0, 120.0, 7446.65, 4181.17),
]
Z_REC = 0.1      # receiver depth in km (iz = 2)
Z_SRC = 100.0    # source depth in km

def travel(z0, z1, mode):
    """1-D travel time between depths z0 (deeper) and z1 (shallower)."""
    t = 0.0
    for top, bot, vp, vs in LAYERS:
        a = max(min(z1, bot), top)      # shallow end, clipped to layer
        b = min(max(z0, top), bot)      # deep end, clipped to layer
        if b > a:
            t += (b - a) * 1000.0 / (vp if mode == "P" else vs)
    return t

# ---- predictions -----------------------------------------------------------
t_P = travel(Z_SRC, Z_REC, "P")
def ps_delay(d):    # Ps at interface depth d: S above, P below
    return travel(d, Z_REC, "S") - travel(d, Z_REC, "P")
def ppm_delay(d):   # free-surface reflection + reflection off interface d
    return 2.0 * travel(d, Z_REC, "P")

PHASES = [
    ("P",      t_P,                    "direct P from 100 km"),
    ("Ps10",   t_P + ps_delay(10.0),   "Ps conversion at 10 km"),
    ("Ps77",   t_P + ps_delay(77.4),   "Ps conversion at Moho (77.4 km)"),
    ("PPmP10", t_P + ppm_delay(10.0),  "multiple off 10 km interface"),
    ("PPmP77", t_P + ppm_delay(77.4),  "PPmP Moho multiple"),
]

# ---- read traces -----------------------------------------------------------
base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mars_run")
with open(os.path.join(base, "traces.csv")) as f:
    rdr = csv.reader(f)
    header = next(rdr)
    rows = [r for r in rdr]
meta = {}
for line in open(os.path.join(base, "meta.jl")):
    if "=" in line:
        k, v = line.split("=", 1)
        meta[k.strip().replace("const ", "").strip()] = float(v.strip())
dt, f0 = meta["dt"], meta["f0"]
t = [float(r[0]) for r in rows]
traces = [[float(v) for v in c] for c in zip(*rows)][1:]   # nrec x nt
nrec = len(traces)

# ---- envelope picking ------------------------------------------------------
PULSE_SPAN = 4.0 + 1.5     # source wavelet length + ringdown (s)

def envelope(tr):
    """Moving-average of |trace| over ~0.5 s."""
    w = max(1, int(round(0.5 / dt)))
    ab = [abs(v) for v in tr]
    out, s = [], 0.0
    for i in range(len(ab)):
        s += ab[i]
        if i >= w: s -= ab[i - w]
        out.append(s / min(i + 1, w))
    return out

def peak_pick(env, t_pred, half_win=4.0):
    """Envelope-maximum time within +/-half_win of the predicted onset."""
    lo = next((i for i, tv in enumerate(t) if tv >= t_pred - half_win), None)
    hi = next((i for i, tv in enumerate(t) if tv >= t_pred + half_win), len(t))
    if lo is None or hi <= lo: return None, 0.0
    im = max(range(lo, hi), key=lambda i: env[i])
    return t[im], env[im]

def onset_pick(env, t_peak, frac=0.10, half_win=4.0):
    """First time the envelope exceeds frac of the phase's window peak."""
    lo = next((i for i, tv in enumerate(t) if tv >= t_peak - half_win), None)
    hi = next((i for i, tv in enumerate(t) if tv >= t_peak + half_win), len(t))
    if lo is None or hi is None or hi <= lo: return None
    pk = max(env[lo:hi])
    for i in range(lo, hi):
        if env[i] >= frac * pk:
            return t[i]
    return None

print(f"receivers: {nrec}, nt: {len(t)}, dt: {dt:.5f} s, f0: {f0} Hz\n")
envs = [envelope(tr) for tr in traces]

# ---- cross-correlation picking ---------------------------------------------
# All phases share the source wavelet, so aligning each phase window
# against the direct-P envelope estimates the onset delay without the
# envelope-peak / wavelet-shape bias.
IDX = lambda tv: min(max(int(round(tv / dt)), 0), len(t) - 1)

def xcorr_pick(env, t_ref0, t_search0, t_search1):
    """Lag maximizing normalised correlation of env with the P reference."""
    k0, k1 = IDX(t_ref0), IDX(t_ref0 + 5.5)     # P reference wavelet
    r = [env[k] for k in range(k0, k1)]
    rr = sum(v * v for v in r) ** 0.5
    if rr == 0.0: return None, 0.0
    s0, s1 = IDX(t_search0), IDX(t_search1) - (k1 - k0)
    best, best_sc = None, -1.0
    for L in range(s0, max(s0 + 1, s1 + 1)):
        seg = [env[k] for k in range(L, L + (k1 - k0))]
        ss = sum(v * v for v in seg) ** 0.5
        if ss == 0.0: continue
        sc = sum(r[i] * seg[i] for i in range(len(r))) / (rr * ss)
        if sc > best_sc: best_sc, best = sc, L
    if best is None: return None, 0.0
    pk = max(env[best:best + (k1 - k0)])
    # onset = search start + lag, corrected for the reference offset
    return t[best] - (t_ref0 - t_P), pk

tP_pred = t_P
print(f"{'phase':8s} {'pred s':>7s} {'onset s':>8s} {'err s':>7s} {'delay pred':>10s} "
      f"{'delay meas':>10s} {'peakenv':>8s}  note")
picks_rows = [["phase", "pred_onset_s", "meas_onset_s", "onset_err_s",
               "pred_delay_s", "meas_delay_s", "peak_env", "note"]]
for name, tp, note in PHASES:
    ons, pks = [], []
    for env in envs:
        to, pk = xcorr_pick(env, tP_pred - 1.0, tp - 2.5, tp + 8.0)
        if to is not None:
            ons.append(to); pks.append(pk)
    on = sum(ons) / len(ons) if ons else None
    pkv = max(pks) if pks else 0.0
    pred_d = tp - tP_pred
    meas_d = (on - tP_pred) if on is not None else None
    if name != "P" and tp < tP_pred + PULSE_SPAN:
        flag = "window contaminated by direct-P coda (unresolved)"
    elif meas_d is not None:
        flag = note + (f", err {meas_d - pred_d:+.2f} s" if name != "P" else "")
    else:
        flag = note
    fmt = lambda v: f"{v:.2f}" if v is not None else "  -  "
    print(f"{name:8s} {tp:7.2f} {fmt(on):>8s} {fmt((on - tp) if on is not None else None):>7s} "
          f"{pred_d:10.2f} {fmt(meas_d):>10s} {pkv:8.1e}  {flag}")
    picks_rows.append([name, f"{tp:.3f}",
                       f"{on:.3f}" if on is not None else "",
                       f"{on - tp:.3f}" if on is not None else "",
                       f"{pred_d:.3f}", f"{meas_d:.3f}" if meas_d is not None else "",
                       f"{pkv:.3e}", flag])

with open(os.path.join(base, "picks.csv"), "w", newline="") as f:
    csv.writer(f).writerows(picks_rows)

# ---- record-section PNG (hand-rolled, zlib only) ----------------------------
SPP = 16           # px per second
RH  = 34           # px per trace row
L, R = 6, 0        # left margin for labels
W = L + int(t[-1] * SPP) + R
H = nrec * RH + 40
RED, GRN, BLU, BG = (60, 90, 200), (0, 170, 60), (200, 90, 60), (24, 24, 30)

img = [[BG] * W for _ in range(H)]
gmax = max(max(abs(v) for v in tr) for tr in traces)
clip = gmax / 3.0

def put(x, y, c):
    if 0 <= x < W and 0 <= y < H: img[y][x] = c

# stacked, individually-normalised wiggles + variable-density fill
for j, tr in enumerate(traces):
    y0 = j * RH + RH // 2
    amp = max(abs(v) for v in tr)
    prev_y = y0
    for i, tv in enumerate(t):
        x = L + int(tv * SPP)
        if abs(tr[i]) > amp * 0.15:
            c = RED if tr[i] > 0 else BLU
            sgn = 1 if tr[i] > 0 else -1
            h = min(RH // 2 - 1, int(abs(tr[i]) / amp * (RH // 2 - 1)))
            for yy in range(min(y0, y0 + sgn * h), max(y0, y0 + sgn * h) + 1):
                put(x, yy, c)
        y = int(round(y0 - tr[i] / clip * (RH // 2 - 1)))
        y = max(y0 - (RH // 2 - 1), min(y0 + (RH // 2 - 1), y))
        lo2, hi2 = min(prev_y, y), max(prev_y, y)
        for yy in range(lo2, hi2 + 1):
            put(x, yy, (235, 235, 235) if abs(tr[i]) > amp * 0.02 else BG)
        prev_y = y

# predicted arrival curves (green, dashed)
for name, tp, _ in PHASES:
    x = L + int(tp * SPP)
    for y in range(nrec * RH):
        if (y // 3) % 2 == 0: put(x, y, GRN)

# receiver labels + time axis (tiny hand-drawn bitmap font)
DIG = {
    "0":[(1,1),(2,1),(3,1),(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(2,5),(1,5),(3,5),(1,6),(3,6)],
    "1":[(3,1),(2,1),(3,2),(3,3),(3,4),(3,5),(3,6),(2,4),(1,3)],
    "2":[(1,1),(2,1),(3,2),(3,3),(2,4),(1,5),(1,6),(2,6),(3,6)],
    "3":[(1,1),(2,1),(3,1),(3,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
    "4":[(1,1),(1,2),(1,3),(2,3),(3,1),(3,2),(3,3),(3,4),(3,5),(3,6)],
    "5":[(1,1),(2,1),(3,1),(1,2),(1,3),(2,3),(3,4),(3,5),(1,5),(2,6),(3,6)],
    "6":[(2,1),(3,2),(1,2),(1,3),(1,4),(2,4),(3,4),(3,5),(3,6),(1,5),(1,6),(2,6)],
    "7":[(1,1),(2,1),(3,1),(3,2),(3,3),(2,4),(2,5),(1,6)],
    "8":[(2,1),(1,1),(3,1),(1,2),(3,2),(2,3),(1,4),(3,4),(1,5),(3,5),(2,6),(1,6),(3,6)],
    "9":[(2,1),(1,1),(3,1),(1,2),(3,2),(1,3),(2,3),(3,4),(3,5),(1,4),(3,6),(2,6)],
    ".":[(2,5),(2,6)],
    "k":[(1,1),(2,1),(1,2),(2,2),(1,3),(2,3),(2,4),(3,5),(3,6)],
    "m":[(1,1),(2,1),(1,2),(1,3),(2,3),(1,4),(1,5),(3,2),(3,5),(1,6),(3,6)],
    "s":[(2,1),(3,2),(1,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
}
def draw_text(s, x0, y0, col):
    for ci, ch in enumerate(s):
        for dx, dy in DIG.get(ch, []):
            put(x0 + ci * 4 + dx - 1, y0 + dy - 1, col)

for j in range(nrec):
    draw_text(f"{14 + j}km", 0, j * RH + RH // 2 - 4, (150, 150, 160))
for sec in range(0, int(t[-1]) + 1, 10):
    draw_text(f"{sec}s", L + sec * SPP - 4, H - 14, (150, 150, 160))
    for y in range(nrec * RH): put(L + sec * SPP, y, (60, 60, 70))

raw = b"".join(b"\x00" + bytes(c for px in row for c in px) for row in img)
def chunk(tag, data):
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff)
png = (b"\x89PNG\r\n\x1a\n"
       + chunk(b"IHDR", struct.pack(">IIBBBBB", W, H, 8, 2, 0, 0, 0))
       + chunk(b"IDAT", zlib.compress(raw, 6))
       + chunk(b"IEND", b""))
with open(os.path.join(base, "section.png"), "wb") as f:
    f.write(png)
print(f"\nsection.png: {W}x{H}px, picks.csv written")
