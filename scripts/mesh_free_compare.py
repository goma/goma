"""Mesh-independent comparison of two Goma Exodus outputs (2D, TRI3).

Intended for adaptive (Omega_h) runs where the two meshes differ, e.g. a
serial run against a parallel run, so node-by-node exodiff is not meaningful.

For each time step reports:
  * level set region F<0: area, centroid and F=0 interface length (exact for P1)
  * relative L2 differences of VX, VY, P and max abs difference of F, sampled on
    a common nx x ny grid (P excludes points within 2h of either interface).

Requires numpy, scipy and the SEACAS exodus.py module (found via SEACAS_DIR or
ACCESS, i.e. $SEACAS_DIR/lib, or already on PYTHONPATH).

usage: mesh_free_compare.py A.exoII B.exoII [nx ny]

Each Exodus file holds a single mesh, so compare matching adapt segments
(out.exoII vs out.exoII, out.exoII-s.0 vs out.exoII-s.0, ...).
"""
import os
import sys

for _env in ("SEACAS_DIR", "ACCESS"):
    if os.environ.get(_env):
        sys.path.insert(0, os.path.join(os.environ[_env], "lib"))
import numpy as np
from scipy.spatial import cKDTree
import exodus  # noqa: E402


def load(path):
    e = exodus.exodus(path, "r", array_type="numpy")
    x, y, _ = e.get_coords()
    conn = []
    for bid in e.get_ids("EX_ELEM_BLOCK"):
        c, ne, npe = e.get_elem_connectivity(bid)
        assert npe == 3, f"expected TRI3, got {npe} nodes/elem"
        conn.append(np.asarray(c).reshape(ne, npe) - 1)
    conn = np.vstack(conn)
    names = e.get_node_variable_names()
    times = np.asarray(e.get_times())
    fields = {
        n: [np.asarray(e.get_node_variable_values(n, k + 1)) for k in range(len(times))]
        for n in ("VX", "VY", "P", "F")
        if n in names
    }
    gnames = e.get_global_variable_names()
    mv = [e.get_global_variable_value("MESH_VOLUME", k + 1) for k in range(len(times))] if "MESH_VOLUME" in gnames else None
    e.close()
    return dict(xy=np.column_stack([x, y]), conn=conn, times=times, fields=fields, mesh_vol=mv,
                nnodes=len(x), nelem=len(conn))


def neg_region(xy, conn, f):
    """Area, centroid of {f<0} and length of {f=0} on a P1 mesh (exact)."""
    area = 0.0
    mom = np.zeros(2)
    length = 0.0
    for tri in conn:
        p = xy[tri]
        v = f[tri]
        neg = v < 0
        k = neg.sum()
        if k == 0:
            continue
        if k == 3:
            d1, d2 = p[1] - p[0], p[2] - p[0]
            a = 0.5 * abs(d1[0] * d2[1] - d1[1] * d2[0])
            area += a
            mom += a * p.mean(axis=0)
            continue
        # clip triangle: build polygon of negative part
        poly, cut = [], []
        for i in range(3):
            j = (i + 1) % 3
            if v[i] < 0:
                poly.append(p[i])
            if (v[i] < 0) != (v[j] < 0):
                t = v[i] / (v[i] - v[j])
                q = p[i] + t * (p[j] - p[i])
                poly.append(q)
                cut.append(q)
        poly = np.array(poly)
        # shoelace
        xs, ys = poly[:, 0], poly[:, 1]
        cr = xs * np.roll(ys, -1) - np.roll(xs, -1) * ys
        a = 0.5 * cr.sum()
        if abs(a) > 0:
            cx = ((xs + np.roll(xs, -1)) * cr).sum() / (6 * a)
            cy = ((ys + np.roll(ys, -1)) * cr).sum() / (6 * a)
            area += abs(a)
            mom += abs(a) * np.array([cx, cy])
        if len(cut) == 2:
            length += np.linalg.norm(cut[0] - cut[1])
    return area, (mom / area if area > 0 else np.full(2, np.nan)), length


class Locator:
    def __init__(self, xy, conn):
        self.xy, self.conn = xy, conn
        self.cent = xy[conn].mean(axis=1)
        self.tree = cKDTree(self.cent)

    def interp_weights(self, pts, k=12):
        _, cand = self.tree.query(pts, k=k)
        elem = np.full(len(pts), -1)
        lam = np.zeros((len(pts), 3))
        for i, pt in enumerate(pts):
            for e in cand[i]:
                p = self.xy[self.conn[e]]
                T = np.column_stack([p[1] - p[0], p[2] - p[0]])
                l12 = np.linalg.solve(T, pt - p[0])
                l = np.array([1 - l12.sum(), *l12])
                if l.min() >= -1e-10:
                    elem[i], lam[i] = e, l
                    break
        return elem, lam

    @staticmethod
    def apply(conn, elem, lam, f):
        return (f[conn[elem]] * lam).sum(axis=1)


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    nx = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    ny = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    A, B = load(a_path), load(b_path)
    print(f"A: {a_path}\n   nodes={A['nnodes']} elems={A['nelem']} steps={len(A['times'])}")
    print(f"B: {b_path}\n   nodes={B['nnodes']} elems={B['nelem']} steps={len(B['times'])}")
    if not np.allclose(A["times"], B["times"]):
        print("WARNING: time arrays differ", A["times"], B["times"])
    lo = np.maximum(A["xy"].min(0), B["xy"].min(0))
    hi = np.minimum(A["xy"].max(0), B["xy"].max(0))
    pad = 1e-6 * (hi - lo)
    gx, gy = np.meshgrid(np.linspace(lo[0] + pad[0], hi[0] - pad[0], nx),
                         np.linspace(lo[1] + pad[1], hi[1] - pad[1], ny))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    # meshes are static between adapt steps but may change; Exodus stores one mesh,
    # so a single locator per file suffices.
    LA, LB = Locator(A["xy"], A["conn"]), Locator(B["xy"], B["conn"])
    eA, wA = LA.interp_weights(pts)
    eB, wB = LB.interp_weights(pts)
    ok = (eA >= 0) & (eB >= 0)
    pts, eA, wA, eB, wB = pts[ok], eA[ok], wA[ok], eB[ok], wB[ok]
    h = np.sqrt((hi - lo).prod() / A["nelem"])  # rough element size
    print(f"   sample points used: {ok.sum()} / {len(ok)}   (~h = {h:.3g})")

    hdr = (f"{'t':>6} | {'areaA':>9} {'areaB':>9} {'dA%':>7} | {'dCx':>8} {'dCy':>8} | "
           f"{'lenA':>8} {'lenB':>8} | {'VX L2':>8} {'VY L2':>8} {'P L2*':>8} {'F Linf':>8}")
    print(hdr)
    print("-" * len(hdr))
    worst = dict(area=0, cent=0, vx=0, vy=0, p=0, f=0)
    for k, t in enumerate(A["times"]):
        aA, cA, lA = neg_region(A["xy"], A["conn"], A["fields"]["F"][k])
        aB, cB, lB = neg_region(B["xy"], B["conn"], B["fields"]["F"][k])
        s = {}
        for n in ("VX", "VY", "P", "F"):
            s[n] = (Locator.apply(A["conn"], eA, wA, A["fields"][n][k]),
                    Locator.apply(B["conn"], eB, wB, B["fields"][n][k]))
        band = (np.abs(s["F"][0]) > 2 * h) & (np.abs(s["F"][1]) > 2 * h)

        def rel(n, m=None):
            a, b = s[n]
            if m is not None:
                a, b = a[m], b[m]
            return np.linalg.norm(a - b) / max(np.linalg.norm(a), 1e-30)

        dA = 100 * (aB - aA) / aA if aA else np.nan
        dc = cB - cA
        r = dict(vx=rel("VX"), vy=rel("VY"), p=rel("P", band),
                 f=np.abs(s["F"][0] - s["F"][1]).max())
        print(f"{t:6.2f} | {aA:9.4f} {aB:9.4f} {dA:7.3f} | {dc[0]:8.4f} {dc[1]:8.4f} | "
              f"{lA:8.4f} {lB:8.4f} | {r['vx']:8.2e} {r['vy']:8.2e} {r['p']:8.2e} {r['f']:8.2e}")
        worst["area"] = max(worst["area"], abs(dA))
        worst["cent"] = max(worst["cent"], np.linalg.norm(dc))
        for key in ("vx", "vy", "p", "f"):
            worst[key] = max(worst[key], r[key])
    if A["mesh_vol"] and B["mesh_vol"]:
        print(f"MESH_VOLUME A={A['mesh_vol'][-1]:.6g} B={B['mesh_vol'][-1]:.6g}")
    print("WORST: " + " ".join(f"{k}={v:.3g}" for k, v in worst.items()))
    print("(P L2* excludes points within 2h of either interface; F Linf is absolute)")


if __name__ == "__main__":
    main()
