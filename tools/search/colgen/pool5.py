"""Row-pool master, CSR storage (float32 data, int32 indices), batched column additions, and
compact save/load so that a restart does not lose the cutting planes."""
from __future__ import annotations

import numpy as np
import highspy
import scipy.sparse as sp


class PoolMaster:
    def __init__(self, geo, margin=2e-6, ipm=False):
        self.geo = geo
        self.margin = margin
        self.lo = 1.0 + margin
        self.inf = highspy.kHighsInf
        self.h = highspy.Highs()
        self.h.setOptionValue("output_flag", False)
        if ipm:
            self.h.setOptionValue("solver", "ipm")
            self.h.setOptionValue("run_crossover", "off")
        else:
            self.h.setOptionValue("solver", "simplex")
        self.ipm = ipm
        self.P = np.zeros((0, 7))            # pool poses: k, cx, cy, c, s, hx, hy
        self.keys = set()
        self.A = sp.csr_matrix((0, 0), dtype=np.float32)
        self.active = np.zeros(0, dtype=np.int64)   # model row -> pool row
        self.is_active = np.zeros(0, dtype=bool)
        self.ncols = 0
        self.sync_cols()

    # ------------------------------------------------------------------ pool
    def n_pool(self):
        return len(self.P)

    def pose_row(self, r):
        return self.P[r]

    def add_pool_rows(self, items):
        rows, cols, vals, poses = [], [], [], []
        n = 0
        for pose, idx, cnt in items:
            key = hash((idx.tobytes(), cnt.tobytes()))
            if key in self.keys:
                continue
            self.keys.add(key)
            poses.append((pose.k, pose.cx, pose.cy, pose.c, pose.s_, pose.hx, pose.hy))
            rows.append(np.full(len(idx), n, dtype=np.int32)); cols.append(idx.astype(np.int32)); vals.append(cnt.astype(np.float32))
            n += 1
        if n:
            B = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                              shape=(n, self.ncols), dtype=np.float32)
            self.A = sp.vstack([self.A, B], format="csr", dtype=np.float32) if self.A.shape[0] else B
            self.P = np.vstack([self.P, np.array(poses)])
            self.is_active = np.concatenate([self.is_active, np.zeros(n, dtype=bool)])
        return n

    def _counts(self, mem):
        A = self.P
        cnt = np.zeros(len(A), dtype=np.int32)
        for (px, py) in mem:
            dx = px - A[:, 1]; dy = py - A[:, 2]
            u0 = A[:, 3] * dx + A[:, 4] * dy
            u1 = -A[:, 4] * dx + A[:, 3] * dy
            cnt += ((np.abs(u0) <= A[:, 5]) & (np.abs(u1) <= A[:, 6]))
        return cnt

    def sync_cols(self):
        n = len(self.geo.reps)
        if self.ncols >= n:
            return
        new = list(range(self.ncols, n))
        npool = len(self.P)
        rr, cc, vv = [], [], []
        act = self.active
        for j, o in enumerate(new):
            size = len(self.geo.members[o])
            if npool:
                mem = np.array(self.geo.members[o], dtype=np.float64) / self.geo.D
                cnt = self._counts(mem)
                nz = np.nonzero(cnt)[0]
                rr.append(nz.astype(np.int32)); cc.append(np.full(len(nz), j, dtype=np.int32)); vv.append(cnt[nz].astype(np.float32))
                if len(act):
                    va = cnt[act]
                    m = np.nonzero(va)[0]
                    self.h.addCol(float(size), 0.0, self.inf, len(m), m.astype(np.int32), va[m].astype(np.float64))
                else:
                    self.h.addCol(float(size), 0.0, self.inf, 0, np.zeros(0, dtype=np.int32), np.zeros(0))
            else:
                self.h.addCol(float(size), 0.0, self.inf, 0, np.zeros(0, dtype=np.int32), np.zeros(0))
        if npool:
            C = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(npool, len(new)), dtype=np.float32)
            self.A = sp.hstack([self.A, C], format="csr", dtype=np.float32)
        else:
            self.A = sp.csr_matrix((0, n), dtype=np.float32)
        self.ncols = n

    def _activate(self, prow):
        B = self.A[prow]
        n = len(prow)
        self.h.addRows(n, np.full(n, self.lo), np.full(n, self.inf), B.nnz, B.indptr[:-1].astype(np.int32),
                       B.indices.astype(np.int32), B.data.astype(np.float64))
        self.active = np.concatenate([self.active, np.asarray(prow, dtype=np.int64)])
        self.is_active[prow] = True

    def warm_activate(self, x0, slack=0.002, cap=20000):
        """Activate the pool rows that are tight (within slack) or violated at x0."""
        if not self.A.shape[0]:
            return 0
        act = np.asarray(self.A @ np.asarray(x0, dtype=np.float32), dtype=np.float64)
        cand = np.nonzero((act < self.lo + slack) & ~self.is_active)[0]
        cand = cand[np.argsort(act[cand])][:cap]
        if len(cand):
            self._activate(cand)
        return len(cand)

    # ------------------------------------------------------------------ solving
    def x(self):
        return np.array(self.h.getSolution().col_value)

    def solve(self, max_add=4000, tol=1e-9, max_iter=60):
        it = 0
        while True:
            if len(self.active):
                self.h.run()
                st = self.h.getModelStatus()
                if st != highspy.HighsModelStatus.kOptimal:
                    info = self.h.getInfo()
                    if not self.ipm or info.max_primal_infeasibility > 1e-7:
                        raise RuntimeError(f"LP status {st} pinf {info.max_primal_infeasibility}")
                x = np.maximum(self.x(), 0.0)
            else:
                x = np.zeros(self.ncols)
            act = self.A @ x.astype(np.float32) if self.A.shape[0] else np.zeros(0)
            act = np.asarray(act, dtype=np.float64)
            viol = np.nonzero((act < self.lo - 1e-6 * 0 - tol) & ~self.is_active)[0]
            if len(viol) == 0 or it >= max_iter:
                break
            order = viol[np.argsort(act[viol])][:max_add]
            self._activate(order)
            it += 1
        y = np.array(self.h.getSolution().row_dual) if len(self.active) else np.zeros(0)
        obj = float(np.dot([len(m) for m in self.geo.members[: self.ncols]], x))
        return x, y, obj

    def purge(self, max_slack):
        if not len(self.active):
            return 0
        act = np.array(self.h.getSolution().row_value)
        drop = np.nonzero(act - self.lo > max_slack)[0]
        if len(drop) == 0:
            return 0
        self.h.deleteRows(len(drop), drop.astype(np.int32))
        keep = np.ones(len(self.active), dtype=bool); keep[drop] = False
        self.is_active[self.active[drop]] = False
        self.active = self.active[keep]
        return len(drop)

    def prune_pool(self, x, max_slack):
        act = np.asarray(self.A @ x.astype(np.float32), dtype=np.float64)
        keep = (act - self.lo <= max_slack) | self.is_active
        if keep.all():
            return 0
        newidx = np.cumsum(keep) - 1
        self.A = self.A[keep]
        self.P = self.P[keep]
        self.active = newidx[self.active]
        self.is_active = self.is_active[keep]
        self.keys = set()
        return int((~keep).sum())

    def drop_cols(self, thr):
        """Delete zero-weight columns whose reduced cost exceeds thr; returns the keep mask."""
        sol = self.h.getSolution()
        x = np.array(sol.col_value); rc = np.array(sol.col_dual)
        drop = np.nonzero((x <= 1e-12) & (rc > thr))[0]
        keep = np.ones(self.ncols, dtype=bool)
        if len(drop) == 0:
            return keep
        self.h.deleteCols(len(drop), drop.astype(np.int32))
        keep[drop] = False
        self.A = self.A[:, keep].tocsr()
        self.ncols = int(keep.sum())
        self.keys = set()
        return keep

    # ------------------------------------------------------------------ persistence
    def save(self, path):
        A = self.A
        np.savez(path, data=A.data, indices=A.indices, indptr=A.indptr, shape=np.array(A.shape), P=self.P)

    def load(self, path):
        z = np.load(path)
        shape = tuple(int(v) for v in z["shape"])
        if shape[1] != self.ncols:
            raise ValueError(f"pool has {shape[1]} columns, master {self.ncols}")
        self.A = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)
        self.P = z["P"]
        self.is_active = np.zeros(shape[0], dtype=bool)
        self.active = np.zeros(0, dtype=np.int64)
