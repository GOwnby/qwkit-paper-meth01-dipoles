"""METH-01: numbers quoted in the note, from scan.json -> numbers.log"""
import json, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(here, "scan.json")))
L = []
P = lambda s: (L.append(s), print(s))
M = ("raw", "raw+100", "norm", "GS", "L-bound", "1st", "1st-all", "L-b2", "L-b10", "L-E0.5", "L-E1.0", "L-E2.0", "L-all")
for r in R:
    P(f"== {r['name']}{' matched K' if r['matched'] else ''} padding {r['pad']:.0f} A: {r['nloc']} bound (edge {r['edge']:.3f} eV), "
      f"E {np.round(r['E'], 4).tolist()}, max one-band overlap {r['S_max']:.4f}, two-band overlap {r['kane_overlap']:.1e}, "
      f"light-particle weight {np.round(r['light_weight'], 4).tolist()}, states in 0.5/1/2 eV buffers {list(r['nE'].values())}")
    P("   two-band dipoles (A): " + ", ".join(f"z{k} {v:.3f}" for k, v in r["kane"].items()))
    for m in M:
        P(f"   {m:8s} " + " ".join(f"z{k} {100 * v:+6.2f}%" for k, v in r["err"][m].items())
          + (f" | peak chi2 / two-band {r['chi'][m]:.4f}" if m in r["chi"] else ""))
for matched in (False, True):
    P(f"== summary{' (matched K)' if matched else ''}: largest |error| (%) over allowed transitions and paddings")
    for name in ("ADQW", "TW", "SQW"):
        rr = [r for r in R if r["name"] == name and r["matched"] == matched]
        if not rr:
            continue
        P(f"{name} (paddings {sorted(int(r['pad']) for r in rr)}): " + ", ".join(
            f"{m} {max(100 * max(abs(v) for v in r['err'][m].values()) for r in rr):.2f}" for m in M))
        if rr[0]["chi"]:
            P(f"   peak chi2 error (%): " + ", ".join(f"{m} {min(100 * (r['chi'][m] - 1) for r in rr):+.2f}..{max(100 * (r['chi'][m] - 1) for r in rr):+.2f}" for m in M))
P("== |1st - L-bound| (percentage points), largest; 1st-all errors")
P(str(max(abs(100 * (r["err"]["1st"][k] - r["err"]["L-bound"][k])) for r in R for k in r["err"]["1st"])))
for name in ("ADQW", "TW", "SQW"):
    r = [r for r in R if r["name"] == name and r["pad"] == 200.0 and not r["matched"]][0]
    P(f"{name} 200 A: 1st-all " + " ".join(f"z{k} {100 * v:+.2f}%" for k, v in r["err"]["1st-all"].items()))
open(os.path.join(here, "numbers.log"), "w").write("\n".join(L) + "\n")
