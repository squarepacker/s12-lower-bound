"""Weighted point certificates for s(12): LP re-weighting with cutting planes from Daniel's
verifier (TIGHT_DUMP) and column generation of new points priced by LP duals.

Floats propose, the verifier decides: every certificate written here is only a candidate until a
complete sweep of Daniel's `verify` (and an independent checker) accepts it.

Inspired by jlevy/squares devtools.s12_reweight (MIT) and Daniel's search/tighten.py; written
from scratch.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

import numpy as np
import highspy
from shapely.geometry import Polygon, box as sbox

SCR = Path("/tmp/claude-0/-home-claude/827af079-8944-5ba5-b8f5-20e5625f6a0c/scratchpad")
VERIFY = SCR / "evsp/s12/verify/target/release/verify"
DANIEL = SCR / "evsp/s12/certificates/s12_lower_3.9686.txt"
WSCALE = 10_000_000


# ----------------------------------------------------------------------------- certificates
def read_cert(path):
    toks = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].strip()
        if line:
            toks.append(line.split())
    s_num, s_den = map(int, toks[0])
    D = int(toks[1][0]); W = int(toks[2][0]); m = int(toks[3][0])
    pts = [tuple(int(v) for v in t) for t in toks[4:4 + m]]
    return s_num, s_den, D, W, pts


def write_cert(path, U, D, pts_w, W=WSCALE):
    lines = [f"{U} {D}", f"{D}", f"{W}", f"{len(pts_w)}"]
    lines += [f"{x} {y} {w}" for x, y, w in pts_w]
    Path(path).write_text("\n".join(lines) + "\n")


def images(X, Y, U):
    return sorted({(X, Y), (U - X, Y), (X, U - Y), (U - X, U - Y),
                   (Y, X), (U - Y, X), (Y, U - X), (U - Y, U - X)})


def canon(X, Y, U):
    """Representative of the D4 orbit in the triangle 0 <= Y <= X <= U/2 (doubled test)."""
    best = None
    for (a, b) in images(X, Y, U):
        if 2 * a <= U and 2 * b <= U and b <= a:
            if best is None or (a, b) < best:
                best = (a, b)
    return best


@dataclass
class Geometry:
    U: int            # container side in grid units
    D: int            # grid denominator; side s = U / D
    reps: list = field(default_factory=list)      # orbit representatives (X, Y)
    members: list = field(default_factory=list)   # per orbit: list of (X, Y) images

    @property
    def s(self):
        return self.U / self.D

    def add_orbit(self, X, Y):
        r = canon(X, Y, self.U)
        if r is None:
            raise ValueError("no canonical rep")
        if r in self._repset:
            return None
        self._repset.add(r)
        self.reps.append(r)
        self.members.append(images(r[0], r[1], self.U))
        return len(self.reps) - 1

    def __post_init__(self):
        self._repset = set()

    def points_array(self):
        """(n_points, 2) float array of all points and orbit index per point."""
        P, O = [], []
        for o, mem in enumerate(self.members):
            for (x, y) in mem:
                P.append((x, y)); O.append(o)
        return np.array(P, dtype=np.int64), np.array(O, dtype=np.int64)


def daniel_geometry(U, D):
    s_num, s_den, D0, W0, pts = read_cert(DANIEL)
    U0 = s_num * D0 // s_den
    assert Fraction(s_num, s_den) * D0 == U0
    g = Geometry(U, D)
    w0 = {}
    for (x, y, w) in pts:
        r = canon(x, y, U0)
        if r not in w0:
            w0[r] = w
    for r in sorted(w0):
        # scale the representative, then build exact images
        X = round(Fraction(r[0] * U, U0)); Y = round(Fraction(r[1] * U, U0))
        g.add_orbit(X, Y)
    return g


# ----------------------------------------------------------------------------- verifier runs
def bin_count(N):
    K = 0
    while (K + N) ** 2 < 2 * N * N:
        K += 1
    return K


def run_verify(cert, N, threads=1, bins=None, dump=None, thresh=None, topk=0, wit=None, overflow=False, cap_per_bin=3000):
    env = dict(os.environ)
    if bins is not None:
        env["VERIFY_BINS"] = f"{bins[0]}:{bins[1]}"
    if dump is not None:
        env["TIGHT_DUMP"] = str(dump); env["TIGHT_THRESH"] = str(thresh); env["TIGHT_MAX"] = str(cap_per_bin * (bin_count(N) + 1))
    exe = str(VERIFY) if not overflow else str(VERIFY) + "_ovf"
    args = [exe, str(cert), "12", str(N), str(threads), str(topk)]
    if wit is not None:
        args.append(str(wit))
    p = subprocess.run(args, env=env, capture_output=True, text=True)
    out = p.stdout + p.stderr
    least = None
    for line in out.splitlines():
        if line.startswith("min covered weight"):
            least = int(line.split("=")[1].split("/")[0])
    return out, least


# ----------------------------------------------------------------------------- dump parsing
@dataclass
class Frame:
    k: int; den: int; cn: int; sn: int; g: int; hh: int; sg_n: int; sg_d: int; wm_n: int; wm_d: int


def parse_dump(text):
    frames, cells = {}, []
    for line in text.splitlines():
        if line.startswith("bin "):
            f = line.split()
            k = int(f[1])
            frames[k] = Frame(k, int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]),
                              int(f[11]), int(f[12]))
        elif line.startswith("c "):
            f = line.split()
            area = float(f[9])
            cells.append((int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), area))
    return frames, cells


# ----------------------------------------------------------------------------- rows (poses)
@dataclass
class Pose:
    k: int
    cx: float; cy: float      # centre, container units
    c: float; s_: float       # cos, sin of theta_k
    hx: float                 # half-widths of the CORE: points inside the tested square for every
    hy: float                 # centre of the (clipped) cell, in the bin's rotated frame


def pose_members_float(pose, P_float):
    dx = P_float[:, 0] - pose.cx; dy = P_float[:, 1] - pose.cy
    u0 = pose.c * dx + pose.s_ * dy
    u1 = -pose.s_ * dx + pose.c * dy
    return (np.abs(u0) <= pose.hx) & (np.abs(u1) <= pose.hy)


class BinCtx:
    """Points of the current geometry in one bin's rotated integer frame (doubled), for the exact
    midpoint test, plus the rotated admissible box polygon for clipping."""

    def __init__(self, geo, P_int, fr):
        k1 = 2 * fr.sg_d * fr.wm_d
        self.fr = fr
        X = [int(v) for v in P_int[:, 0]]; Y = [int(v) for v in P_int[:, 1]]
        self.q0 = [2 * (fr.cn * x + fr.sn * y) * k1 for x, y in zip(X, Y)]
        self.q1 = [2 * (-fr.sn * x + fr.cn * y) * k1 for x, y in zip(X, Y)]
        self.qf = np.array([self.q0, self.q1], dtype=np.float64).T
        self.reach = 2 * fr.hh
        # admissible centre box in container coords: [L, s-L]^2, L = wm/2
        s = geo.U / geo.D
        L = fr.wm_n / fr.wm_d / 2.0
        c = fr.cn / fr.g; sn = fr.sn / fr.g
        verts = [(L, L), (s - L, L), (s - L, s - L), (L, s - L)]
        # u = R(-theta) v
        self.poly = Polygon([(c * x + sn * y, -sn * x + c * y) for x, y in verts])
        self.c = c; self.sn = sn
        self.h = fr.sg_n / fr.sg_d / 2.0

    def members(self, a, b, c0, c1):
        m0, m1 = a + b, c0 + c1
        loose = self.reach * (1 + 1e-9) + 4.0
        near = np.nonzero((np.abs(self.qf[:, 0] - m0) <= loose) & (np.abs(self.qf[:, 1] - m1) <= loose))[0]
        out = []
        for i in near:
            if abs(self.q0[i] - m0) <= self.reach and abs(self.q1[i] - m1) <= self.reach:
                out.append(int(i))
        return out

    def pose(self, a, b, c0, c1):
        den = self.fr.den
        if not hasattr(self, "pb"):
            self.pb = self.poly.bounds
        x0 = max(a / den, self.pb[0]); x1 = min(b / den, self.pb[2])
        y0 = max(c0 / den, self.pb[1]); y1 = min(c1 / den, self.pb[3])
        if x1 <= x0 or y1 <= y0:
            return None
        u0, u1 = (x0 + x1) / 2, (y0 + y1) / 2
        cx = self.c * u0 - self.sn * u1
        cy = self.sn * u0 + self.c * u1
        return Pose(self.fr.k, cx, cy, self.c, self.sn, self.h - (x1 - x0) / 2, self.h - (y1 - y0) / 2)


# ----------------------------------------------------------------------------- the LP
class Master:
    def __init__(self, geo, margin=2e-6, ipm=False, defer_cols=False):
        self.geo = geo
        self.margin = margin
        self.rows = []          # list of (Pose, dict orbit->count)
        self.keys = set()
        self.h = highspy.Highs()
        self.h.setOptionValue("output_flag", False)
        if ipm:
            self.h.setOptionValue("solver", "ipm")
            self.h.setOptionValue("run_crossover", "off")
        else:
            self.h.setOptionValue("solver", "simplex")
        self.ncols = 0
        self.inf = highspy.kHighsInf
        if not defer_cols:
            self.sync_cols()

    def _arrays(self):
        n = len(self.rows)
        if getattr(self, "_arr_n", -1) != n:
            P = [r[0] for r in self.rows]
            self._A = np.array([[p.cx, p.cy, p.c, p.s_, p.hx, p.hy] for p in P], dtype=np.float64).reshape(-1, 6)
            self._arr_n = n
        return self._A

    def col_counts(self, mem):
        A = self._arrays()
        cnt = np.zeros(len(A), dtype=np.int64)
        for (px, py) in mem:
            dx = px - A[:, 0]; dy = py - A[:, 1]
            u0 = A[:, 2] * dx + A[:, 3] * dy
            u1 = -A[:, 3] * dx + A[:, 2] * dy
            cnt += ((np.abs(u0) <= A[:, 4]) & (np.abs(u1) <= A[:, 5]))
        return cnt

    def sync_cols(self):
        n = len(self.geo.reps)
        while self.ncols < n:
            o = self.ncols
            size = len(self.geo.members[o])
            idx = np.zeros(0, dtype=np.int32); val = np.zeros(0)
            if self.rows:
                mem = np.array(self.geo.members[o], dtype=np.float64) / self.geo.D
                cnt = self.col_counts(mem)
                nz = np.nonzero(cnt)[0]
                idx = nz.astype(np.int32); val = cnt[nz].astype(np.float64)
            self.h.addCol(float(size), 0.0, self.inf, len(idx), idx, val)
            self.ncols += 1

    def add_rows(self, items):
        """items: (pose, idx int array of orbit ids sorted, cnt int array); returns number added."""
        starts, idxs, vals, n = [], [], [], 0
        tot = 0
        for pose, idx, cnt in items:
            key = hash((idx.tobytes(), cnt.tobytes()))
            if key in self.keys:
                continue
            self.keys.add(key)
            self.rows.append((pose, key))
            starts.append(tot); idxs.append(idx.astype(np.int32)); vals.append(cnt.astype(np.float64))
            tot += len(idx); n += 1
        if n:
            self.h.addRows(n, np.full(n, 1.0 + self.margin), np.full(n, self.inf), tot, np.array(starts, dtype=np.int32),
                           np.concatenate(idxs), np.concatenate(vals))
        return n

    def purge(self, max_slack):
        """Delete rows whose slack at the current solution exceeds max_slack (they are re-added by
        later dumps if they become near-tight again)."""
        sol = self.h.getSolution()
        act = np.array(sol.row_value)
        slack = act - (1.0 + self.margin)
        drop = np.nonzero(slack > max_slack)[0]
        if len(drop) == 0:
            return 0
        self.h.deleteRows(len(drop), drop.astype(np.int32))
        keep = np.ones(len(self.rows), dtype=bool); keep[drop] = False
        for r in drop:
            k = self.rows[r][1]
            if k is not None:
                self.keys.discard(k)
        self.rows = [r for r, kp in zip(self.rows, keep) if kp]
        self._arr_n = -1
        return len(drop)

    def key_of(self, idx, cnt):
        return hash((idx.tobytes(), cnt.tobytes()))

    def state(self):
        lp = self.h.getLp()
        a = lp.a_matrix_
        return dict(poses=[r[0] for r in self.rows], keys=[r[1] for r in self.rows], start=np.array(a.start_), index=np.array(a.index_),
                    value=np.array(a.value_), ncols=self.ncols, nrows=len(self.rows))

    def load_state(self, st):
        """Rebuild the model from a saved CSC matrix (columns must equal the geometry's first ncols)."""
        assert self.ncols == st["ncols"] or self.ncols == 0
        nr = st["nrows"]
        # rows first (empty), then columns with their CSC entries
        h = self.h
        if self.ncols:
            raise RuntimeError("load_state needs an empty master")
        h.addRows(nr, np.full(nr, 1.0 + self.margin), np.full(nr, self.inf), 0, np.zeros(nr, dtype=np.int32),
                  np.zeros(0, dtype=np.int32), np.zeros(0))
        self.rows = [(p, None) for p in st["poses"]]
        if "keys" in st:
            self.rows = [(p, k) for p, k in zip(st["poses"], st["keys"])]
            self.keys = set(k for k in st["keys"] if k is not None)
        start = st["start"]; index = st["index"]; value = st["value"]
        for j in range(st["ncols"]):
            a, b = start[j], start[j + 1]
            size = len(self.geo.members[j])
            h.addCol(float(size), 0.0, self.inf, int(b - a), index[a:b].astype(np.int32), value[a:b].astype(np.float64))
        self.ncols = st["ncols"]

    def solve(self):
        self.h.run()
        st = self.h.getModelStatus()
        if st not in (highspy.HighsModelStatus.kOptimal,):
            # IPM without crossover may stop "imprecise"; accept if primal is (nearly) feasible
            info = self.h.getInfo()
            if not (str(st).endswith("kUnknown") or "Imprecise" in str(st)) or info.max_primal_infeasibility > 1e-7:
                raise RuntimeError(f"LP status {st}")
        sol = self.h.getSolution()
        x = np.array(sol.col_value); y = np.array(sol.row_dual)
        obj = self.h.getInfo().objective_function_value
        return x, y, obj


def weights_to_int(x, eps=1e-12):
    return [max(0, math.ceil(float(v) * WSCALE - 0.05)) for v in x]


def cert_points(geo, wint):
    out = []
    for o, mem in enumerate(geo.members):
        for (X, Y) in mem:
            out.append((X, Y, wint[o]))
    return out


# ----------------------------------------------------------------------------- pricing
def price_grid(geo, master, y, delta, top=60, min_sep=0.004, tol=1e-6):
    """Reduced-cost map on the quadrant grid; returns candidate representatives (X, Y)."""
    s = geo.s
    M = int(math.floor((s / 2) / delta))
    xs = np.arange(M + 1) * delta
    B = np.zeros((M + 1, M + 1), dtype=np.float64)   # B[j, i] at (x_i, y_j)
    pos = np.nonzero(y > 1e-12)[0]
    for r in pos:
        pose, _ = master.rows[r]
        w = y[r]
        for fx in (False, True):
            for fy in (False, True):
                cx = s - pose.cx if fx else pose.cx
                cy = s - pose.cy if fy else pose.cy
                # orientation of the reflected square: reflecting flips the angle sign; a square is
                # symmetric under that so axes stay {(c, s'), (-s', c)} with s' = +-s
                c = pose.c; sn = -pose.s_ if (fx ^ fy) else pose.s_
                hx, hy = pose.hx, pose.hy
                if hx <= 0 or hy <= 0:
                    continue
                e1 = np.array([c, sn]); e2 = np.array([-sn, c])
                V = np.array([[cx, cy]]) + np.array([[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]]) @ np.vstack([e1, e2])
                ylo, yhi = V[:, 1].min(), V[:, 1].max()
                if yhi < 0 or ylo > s / 2 or V[:, 0].max() < 0 or V[:, 0].min() > s / 2:
                    continue
                j0 = max(0, int(math.ceil(ylo / delta))); j1 = min(M, int(math.floor(yhi / delta)))
                if j1 < j0:
                    continue
                yy = xs[j0:j1 + 1]
                xmin = np.full(len(yy), np.inf); xmax = np.full(len(yy), -np.inf)
                for e in range(4):
                    p0 = V[e]; p1 = V[(e + 1) % 4]
                    y0, y1 = p0[1], p1[1]
                    if y0 == y1:
                        msk = np.abs(yy - y0) < 1e-15
                        xmin[msk] = np.minimum(xmin[msk], min(p0[0], p1[0])); xmax[msk] = np.maximum(xmax[msk], max(p0[0], p1[0]))
                        continue
                    t = (yy - y0) / (y1 - y0)
                    msk = (t >= 0) & (t <= 1)
                    xv = p0[0] + t * (p1[0] - p0[0])
                    xmin = np.where(msk, np.minimum(xmin, xv), xmin)
                    xmax = np.where(msk, np.maximum(xmax, xv), xmax)
                ok = xmax >= xmin
                i0 = np.clip(np.ceil(xmin / delta - 1e-9), 0, M + 1).astype(np.int64)
                i1 = np.clip(np.floor(xmax / delta + 1e-9), -1, M).astype(np.int64)
                rows_j = np.arange(j0, j1 + 1)
                ok &= i1 >= i0
                for jj, a, b in zip(rows_j[ok], i0[ok], i1[ok]):
                    B[jj, a] += w
                    if b + 1 <= M:
                        B[jj, b + 1] -= w
    B = np.cumsum(B, axis=1)
    F = B + B.T          # f_8 at (x_i, y_j)
    # fundamental triangle y <= x
    iu = np.tril_indices(M + 1)          # j >= i  -> (j, i) with y_j >= x_i ; we want y <= x
    F_tri = np.where(np.tri(M + 1, dtype=bool).T, F, -np.inf)   # keep j <= i
    order = np.argsort(F_tri, axis=None)[::-1]
    picks = []
    sep = int(round(min_sep / delta))
    for flat in order[: 200000]:
        j, i = divmod(int(flat), M + 1)
        val = F_tri[j, i]
        if val <= 8.0 + tol:
            break
        if any(abs(i - pi) <= sep and abs(j - pj) <= sep for (pi, pj, _) in picks):
            continue
        picks.append((i, j, val))
        if len(picks) >= top:
            break
    out = []
    for (i, j, val) in picks:
        X = int(round(xs[i] * geo.D)); Y = int(round(xs[j] * geo.D))
        out.append((X, Y, float(val)))
    return out, float(np.nanmax(F_tri))
