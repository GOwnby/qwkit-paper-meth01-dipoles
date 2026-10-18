"""METH-01: dipoles of the bound states of an energy-dependent-mass (one-band Kane) Hamiltonian in a finite box.
The one-band states of different energies are not orthogonal; their two-band Kane completions are (core.kane_valence,
discretized exactly as the solver, so the two-component overlaps vanish to the eigenvalue accuracy). Fixed (not
self-consistent) potentials, so that the box size changes nothing but the box: In0.53Ga0.47As/In0.52Al0.48As
64/10/29 A double well (ADQW) at 0 kV/cm, 58/14/34/14/58 A triple well (TW) at 55 kV/cm, 60 A single well (SQW) at
0 kV/cm; barrier padding 100-800 A (half on each side); 120 states solved. Bound states: at least 80% of the
probability within 25 A of the wells AND energy below the lower barrier edge of the active region (resonances above
the barrier are excluded however localized they look). Procedures compared with the two-band dipoles z^K:
  raw      normalized one-band states (origin at the active-region centre; raw+100: origin moved by +100 A)
  norm     two-band norm N_n but no light-particle dipole term
  GS       Gram-Schmidt of the bound states from the lowest up
  L-bound  Loewdin over the bound states only
  1st      first-order Loewdin over the bound states only: z - (S z + z S)/2 with the bound-state overlaps
  1st-all  the same with all 120 states as partners (then truncated)
  L-b{n}   Loewdin over the bound states plus the next n box states (n = 2, 10)
  L-E{x}   Loewdin over all states up to x eV above the topmost bound state (x = 0.5, 1, 2)
  L-all    Loewdin over all 120 solved states
A variant with the barrier gap chosen so that E_g/m* is the same in well and barrier (K = hbar^2 E_g / 2 m* uniform)
tests the origin of the residual of L-all. Reported: relative errors of the allowed dipoles among the lowest four
bound states and of the peak second-order susceptibility (1-2-3 double resonance, 10 meV) against z^K. Writes
scan.json."""
import json, os
import numpy as np
from multiprocessing import Pool
from qwkit import core as q

here = os.path.dirname(os.path.abspath(__file__))
CASES = {"ADQW": ([(64, "well"), (10, "barrier"), (29, "well")], 0.0),
         "TW": ([(58, "well"), (14, "barrier"), (34, "well"), (14, "barrier"), (58, "well")], 55.0),
         "SQW": ([(60, "well")], 0.0)}
PADS = (100.0, 150.0, 200.0, 300.0, 400.0, 800.0)
NS = 120


def chi_peak(E, Z):
    hw = np.arange(0.5, 1.5, 0.002) * (E[1] - E[0])
    return max(abs(q.chi2_shg(E, Z, h, 0.010, 1e18)) for h in hw)


def job(args):
    name, pad, matched = args
    lay, F = CASES[name]
    q.add_material("barrier", 0.52, 0.075, 0.75 * 0.075 / 0.043 if matched else 1.45)   # matched: E_g / m* as in the well
    z, V0, m, Eg, act = q.build(lay, pad=pad / 2, dz=0.5)
    V = q.potential(z, V0, act, F)
    dz = z[1] - z[0]; zc = 0.5 * (act[0] + act[1])
    E, psi = q.solve_states(z, V, m, Eg, nstates=NS, tol=1e-11, iters=60)
    psi = psi / np.sqrt(np.sum(psi ** 2, axis=0) * dz)
    edge = min(np.interp(act[0] - 1.0, z, V), np.interp(act[1] + 1.0, z, V))     # lower barrier edge just outside the wells
    loc = np.where(q.localized(z, psi, act) & (E < edge))[0]
    nl = len(loc); n = min(nl, 4); top = int(loc[-1])
    P = psi[:, loc]
    S = P.T @ P * dz
    Zraw = (P.T * (z - zc)) @ P * dz
    ZK, ovK = q.kane_dipoles(z, P, E[loc], V, m, Eg, zc)
    zh, phi = q.kane_valence(z, P, E[loc], V, m, Eg)
    Nn = np.sqrt(np.diag(S) + np.sum(phi ** 2, axis=0) * dz)
    meth = {"raw": Zraw, "raw+100": (P.T * (z - zc + 100.0)) @ P * dz, "norm": Zraw / np.outer(Nn, Nn)}
    Q, R = np.linalg.qr(P * np.sqrt(dz)); Q = Q * np.sign(np.diag(R))[None, :] / np.sqrt(dz)
    meth["GS"] = (Q.T * (z - zc)) @ Q * dz
    meth["L-bound"] = q.dipoles(z - zc, q.orthonormalize(z, P))
    D = S - np.eye(nl)
    meth["1st"] = Zraw - 0.5 * (D @ Zraw + Zraw @ D)
    Sall = psi.T @ psi * dz - np.eye(NS); Zall = (psi.T * (z - zc)) @ psi * dz
    meth["1st-all"] = (Zall - 0.5 * (Sall @ Zall + Zall @ Sall))[np.ix_(loc, loc)]
    nE = {}
    for nb in (2, 10):
        meth[f"L-b{nb}"] = q.dipoles(z - zc, q.orthonormalize(z, psi[:, :top + 1 + nb])[:, loc], orthonormal=False)
    for de in (0.5, 1.0, 2.0):
        k = int(np.searchsorted(E, E[top] + de)); assert k < NS - 1, "solve more states"
        meth[f"L-E{de}"] = q.dipoles(z - zc, q.orthonormalize(z, psi[:, :k])[:, loc], orthonormal=False); nE[str(de)] = k - top - 1
    meth["L-all"] = q.dipoles(z - zc, q.orthonormalize(z, psi)[:, loc], orthonormal=False)
    big = 0.02 * np.max(np.abs(ZK[:n, :n] - np.diag(np.diag(ZK[:n, :n]))))
    pr = [(i, j) for i in range(n) for j in range(i + 1, n) if abs(ZK[i, j]) > big]
    res = dict(name=name, pad=pad, matched=matched, nloc=nl, E=E[loc][:n].tolist(), edge=float(edge), S_max=float(np.max(np.abs(D))),
               kane_overlap=ovK, light_weight=(np.sum(phi ** 2, axis=0) * dz).tolist(), nE=nE,
               kane={f"{i + 1}{j + 1}": float(ZK[i, j]) for i, j in pr}, err={}, chi={})
    ck = chi_peak(E[loc][:3], ZK[:3, :3]) if n >= 3 else None
    for k_, Z in meth.items():
        res["err"][k_] = {f"{i + 1}{j + 1}": float(abs(Z[i, j]) / abs(ZK[i, j]) - 1) for i, j in pr}
        if ck and name != "SQW":
            res["chi"][k_] = float(chi_peak(E[loc][:3], Z[:3, :3]) / ck)
    line = " | ".join(f"{k_} " + " ".join(f"{100 * v:+.2f}" for v in res["err"][k_].values())
                      for k_ in ("raw", "norm", "L-bound", "L-b10", "L-E0.5", "L-all"))
    print(f"{name}{' (matched K)' if matched else ''} pad {pad:.0f}: {nl} bound (edge {edge:.3f} eV), max overlap "
          f"{res['S_max']:.3f}, two-band overlap {ovK:.1e} | {list(res['kane'])}: {line}", flush=True)
    return res


if __name__ == "__main__":
    jobs = [(c, p, False) for c in CASES for p in PADS] + [(c, 200.0, True) for c in CASES] + [(c, 400.0, True) for c in CASES]
    with Pool(int(os.environ.get("NPROC", "2"))) as pool:
        out = pool.map(job, jobs, chunksize=1)
    json.dump(out, open(os.path.join(here, "scan.json"), "w"), indent=1)
