"""Driver: cutting-plane re-weighting (+ optional column generation) at side U/D, net N."""
from __future__ import annotations

import argparse
import os
import json
import math
import pickle
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pool5 import PoolMaster
from push2 import Pose
from push2 import (Geometry, Master, BinCtx, Pose, bin_count, canon, cert_points, daniel_geometry,
                  parse_dump, price_grid, read_cert, run_verify, weights_to_int, write_cert, WSCALE,
                  DANIEL, images)


def bin_rows(geo, P_int, orbit_of, n_orb, text, known=None, per_bin=60, PU=None, OU=None):
    """Unique (pose, counts) rows from one bin's dump text."""
    frames, cells = parse_dump(text)
    out = []
    by_bin = {}
    for c in cells:
        by_bin.setdefault(c[0], []).append(c)
    for k, cl in by_bin.items():
        fr = frames[k]
        ctx = BinCtx(geo, P_int, fr)
        reach = float(ctx.reach)
        tol = 1e-9 * reach
        seen = {}
        cnt_all = {}
        bs = max(50, int(3_000_000 // max(1, len(ctx.qf))))
        for start in range(0, len(cl), bs):
            blk = cl[start:start + bs]
            m0 = np.array([float(c[1]) + float(c[2]) for c in blk])
            m1 = np.array([float(c[3]) + float(c[4]) for c in blk])
            d0 = np.abs(ctx.qf[None, :, 0] - m0[:, None])
            d1 = np.abs(ctx.qf[None, :, 1] - m1[:, None])
            inside = (d0 <= reach - tol) & (d1 <= reach - tol)
            amb = ((np.abs(d0 - reach) <= tol) & (d1 <= reach + tol)) | ((np.abs(d1 - reach) <= tol) & (d0 <= reach + tol))
            if amb.any():
                for ci, pi in zip(*np.nonzero(amb)):
                    aa, bb, cc0, cc1 = blk[ci][1:5]
                    M0, M1 = aa + bb, cc0 + cc1
                    inside[ci, pi] = abs(ctx.q0[pi] - M0) <= ctx.reach and abs(ctx.q1[pi] - M1) <= ctx.reach
            cnt = np.zeros((len(blk), n_orb), dtype=np.int32)
            ci, pi = np.nonzero(inside)
            np.add.at(cnt, (ci, orbit_of[pi]), 1)
            for r in range(len(blk)):
                key = cnt[r].tobytes()
                if key in seen:
                    continue
                seen[key] = start + r
                cnt_all[start + r] = cnt[r].copy()
        cand = sorted(seen.values(), key=lambda r: cl[r][5])   # lowest captured weight first
        taken = 0
        for r in cand:
            if taken >= per_bin:
                break
            cr = cnt_all[r]
            nzk = np.nonzero(cr)[0]
            kk = tuple((int(o), int(cr[o])) for o in nzk)
            area = cl[r][6]
            if area <= 0:
                continue
            pose = ctx.pose(*cl[r][1:5])
            if pose is None:
                continue
            full = cr.astype(np.int16).copy()
            if PU is not None and len(PU) and pose.hx > 0 and pose.hy > 0:
                dx = PU[:, 0] - pose.cx; dy = PU[:, 1] - pose.cy
                u0 = pose.c * dx + pose.s_ * dy; u1 = -pose.s_ * dx + pose.c * dy
                hit = np.nonzero((np.abs(u0) <= pose.hx) & (np.abs(u1) <= pose.hy))[0]
                if len(hit):
                    full += np.bincount(OU[hit], minlength=n_orb).astype(np.int16)
            nzf = np.nonzero(full)[0].astype(np.int32)
            vals = full[nzf]
            if known is not None and hash((nzf.tobytes(), vals.tobytes())) in known:
                continue
            out.append((pose, nzf, vals))
            taken += 1
    return out


_G = {}


def _rows_worker(text):
    g = _G
    return bin_rows(g["geo"], g["P"], g["O"], g["n"], text, known=g["known"], per_bin=g["per_bin"], PU=g["PU"], OU=g["OU"])


def rows_parallel(geo, P_int, orbit_of, texts, known, per_bin, workers, PU=None, OU=None):
    import multiprocessing as mp
    _G.update(geo=geo, P=P_int, O=orbit_of, n=len(geo.reps), known=known, per_bin=per_bin, PU=PU, OU=OU)
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        return pool.map(_rows_worker, texts, chunksize=4)


def dump_bins(cert, N, bins, thresh, workers, tmp):
    def one(k):
        f = Path(tmp) / f"d{k}.txt"
        out, least = run_verify(cert, N, threads=1, bins=(k, k), dump=f, thresh=thresh)
        text = f.read_text() if f.exists() else ""
        if f.exists():
            f.unlink()
        return k, text, least
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(one, bins))


SCAN = Path(__file__).parent / "tools" / "indep_scan"


def full_scan(cert, N, K, thresh, workers=2):
    """Bins k < K whose exact minimum (indep_scan, exact sigma) is below thresh."""
    import subprocess
    cuts = [round(K * i / workers) for i in range(workers + 1)]
    procs = []
    for i in range(workers):
        procs.append(subprocess.Popen([str(SCAN), str(cert), str(N), "1000000000000000", "12", str(thresh),
                                       str(cuts[i + 1]), str(cuts[i])], stdout=subprocess.PIPE, text=True))
    bad = []
    least = None
    for p in procs:
        out, _ = p.communicate()
        for line in out.splitlines():
            if line.startswith("BIN "):
                _, k, v = line.split()
                bad.append((int(v), int(k)))
            if line.startswith("min captured weight"):
                v = int(line.split("=")[1].split("/")[0])
                least = v if least is None else min(least, v)
    return sorted(bad), least



DUMP = Path(__file__).parent / "tools" / "indep_dump"
QGRID = 10 ** 15


def indep_bins(cert, N, bins, thresh, cap, workers):
    """Run indep_dump on each bin (sigma rounded down to 1e-6 like Daniel's verifier)."""
    import subprocess
    def one(k):
        p = subprocess.run([str(DUMP), str(cert), str(N), str(QGRID), "1", str(thresh), str(k), str(k + 1), str(cap)],
                           capture_output=True, text=True)
        best = None
        for line in p.stdout.splitlines():
            if line.startswith("B "):
                best = int(line.split()[6])
                break
        return k, p.stdout, best
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(one, bins))


def indep_scan_all(cert, N, K, workers=2):
    import subprocess
    cuts = [round(K * i / workers) for i in range(workers + 1)]
    procs = [subprocess.Popen([str(DUMP), str(cert), str(N), str(QGRID), "1", "-1", str(cuts[i]), str(cuts[i + 1]), "0"],
                              stdout=subprocess.PIPE, text=True) for i in range(workers)]
    bad, least = [], None
    for p in procs:
        out, _ = p.communicate()
        for line in out.splitlines():
            if line.startswith("B "):
                f = line.split(); k = int(f[1]); v = int(f[6])
                least = v if least is None else min(least, v)
                if v < WSCALE:
                    bad.append((v, k))
    return sorted(bad), least


def bin_rows_indep(geo, P_int, orbit_of, n_orb, text, known=None, per_bin=60, PU=None, OU=None):
    out = []
    frame = None
    seen = set()
    taken = 0
    for line in text.splitlines():
        if line.startswith("B "):
            f = line.split()
            k = int(f[1]); cn = int(f[2]); sn = int(f[3]); g = int(f[4]); hq = int(f[5])
            frame = (k, cn / g, sn / g, hq / QGRID)
            continue
        if not line.startswith("C ") or taken >= per_bin:
            continue
        f = line.split()
        a, b, ylo, yhi = int(f[2]), int(f[3]), int(f[4]), int(f[5])
        n = int(f[6]); idx = np.array(f[7:7 + n], dtype=np.int64)
        cnt = np.bincount(orbit_of[idx], minlength=n_orb).astype(np.int16) if n else np.zeros(n_orb, dtype=np.int16)
        kb = cnt.tobytes()
        if kb in seen:
            continue
        seen.add(kb)
        k, c, s, h = frame
        u0 = (a + b) / 2 / QGRID; u1 = (ylo + yhi) / 2 / QGRID
        pose = Pose(k, c * u0 - s * u1, s * u0 + c * u1, c, s, h, h)   # pose row: the tested square at the cell centre
        full = cnt.copy()
        if PU is not None and len(PU) and pose.hx > 0 and pose.hy > 0:
            dx = PU[:, 0] - pose.cx; dy = PU[:, 1] - pose.cy
            uu0 = c * dx + s * dy; uu1 = -s * dx + c * dy
            hit = np.nonzero((np.abs(uu0) <= pose.hx) & (np.abs(uu1) <= pose.hy))[0]
            if len(hit):
                full += np.bincount(OU[hit], minlength=n_orb).astype(np.int16)
        nzf = np.nonzero(full)[0].astype(np.int32)
        vals = full[nzf]
        if known is not None and hash((nzf.tobytes(), vals.tobytes())) in known:
            continue
        out.append((pose, nzf, vals))
        taken += 1
    return out


def _rows_worker4(text):
    g = _G
    return bin_rows_indep(g["geo"], g["P"], g["O"], g["n"], text, known=g["known"], per_bin=g["per_bin"], PU=g["PU"], OU=g["OU"])


def rows_parallel4(geo, P_int, orbit_of, texts, known, per_bin, workers, PU=None, OU=None):
    import multiprocessing as mp
    _G.update(geo=geo, P=P_int, O=orbit_of, n=len(geo.reps), known=known, per_bin=per_bin, PU=PU, OU=OU)
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        return pool.map(_rows_worker4, texts, chunksize=8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--U", type=int, required=True)
    ap.add_argument("--D", type=int, required=True)
    ap.add_argument("--N", type=int, default=48000)
    ap.add_argument("--rounds", type=int, default=20)
    ap.add_argument("--stride", type=int, default=48)
    ap.add_argument("--slack", type=float, default=0.004)
    ap.add_argument("--margin", type=float, default=2e-6)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--colgen", type=int, default=0, help="new columns per round (0 = off)")
    ap.add_argument("--delta", type=float, default=0.0005)
    ap.add_argument("--start", type=str, default="", help="state pickle to resume / warm start")
    ap.add_argument("--state", type=str, required=True)
    ap.add_argument("--out", type=str, required=True)
    ap.add_argument("--log", type=str, required=True)
    ap.add_argument("--focus", type=str, default="")
    ap.add_argument("--stop_total", type=float, default=12.0)
    ap.add_argument("--scan_cap", type=int, default=300)
    ap.add_argument("--scan_margin", type=float, default=0.0)
    ap.add_argument("--ipm", type=int, default=1)
    ap.add_argument("--per_bin", type=int, default=60)
    ap.add_argument("--pool", type=str, default="", help="comma-separated certificate files whose points join the column pool")
    ap.add_argument("--pool_top", type=int, default=150)
    ap.add_argument("--purge", type=float, default=0.003)
    ap.add_argument("--price_rounds", type=int, default=1)
    ap.add_argument("--prune", type=float, default=0.01)
    ap.add_argument("--stride0", type=int, default=4)
    ap.add_argument("--tr0", type=float, default=0.0, help="initial trust-region half-width on orbit weights (0 = off)")
    ap.add_argument("--tr_growth", type=float, default=1.6)
    ap.add_argument("--drop_rc", type=float, default=0.0, help="drop zero-weight columns with reduced cost above this (0 = off)")
    ap.add_argument("--always_scan", type=int, default=0, help="scan every bin each round and dump up to this many violated bins")
    ap.add_argument("--slack0", type=float, default=0.05)
    ap.add_argument("--per_bin0", type=int, default=40)
    args = ap.parse_args()

    mstate = None
    if args.start and Path(args.start).exists():
        st = pickle.load(open(args.start, "rb"))
        geo = st["geo"]
        if geo.U != args.U or geo.D != args.D:
            g2 = Geometry(args.U, args.D)
            for (X, Y) in geo.reps:
                g2.add_orbit(round(X * args.U / geo.U), round(Y * args.U / geo.U))
            geo = g2
            w0 = st["w"][:len(geo.reps)]
            rows0 = []
        else:
            w0 = st["w"]
            mstate = st.get("master")
            rows0 = []
            for pose, cnt in st.get("rows", []):
                ks = np.array(sorted(cnt), dtype=np.int32)
                rows0.append((pose, ks, np.array([cnt[k] for k in ks], dtype=np.int16)))
    else:
        geo = daniel_geometry(args.U, args.D)
        s_num, s_den, D0, W0, pts = read_cert(DANIEL)
        U0 = s_num * D0 // s_den
        wd = {}
        for (x, y, w) in pts:
            wd.setdefault(canon(x, y, U0), w)
        w0 = [wd[r] for r in sorted(wd)]
        rows0 = []

    pool = []
    if args.pool:
        from fractions import Fraction as _F
        seenp = set()
        for fname in args.pool.split(","):
            lines = [l.split("#")[0].split() for l in open(fname)]
            lines = [l for l in lines if l]
            sn0, sd0 = map(int, lines[0]); D0 = int(lines[1][0]); m0 = int(lines[3][0])
            U0 = _F(sn0, sd0) * D0
            nb = len(pool)
            for tline in lines[4:4 + m0]:
                xx, yy, ww = map(int, tline)
                if ww <= 0:
                    continue
                r = canon(round(_F(xx) * geo.U / U0), round(_F(yy) * geo.U / U0), geo.U)
                if r is not None and r not in seenp:
                    seenp.add(r); pool.append(r)
            print(f"pool {fname}: +{len(pool) - nb} candidate orbits", flush=True)
    pool_alive = np.array([r not in geo._repset for r in pool], dtype=bool)
    pool_index = {r: i for i, r in enumerate(pool)}
    print(f"pool: {len(pool)} candidates, {int(pool_alive.sum())} inactive", flush=True)
    PP, PO, PS = [], [], []
    for i, r in enumerate(pool):
        im = images(r[0], r[1], geo.U)
        PS.append(len(im))
        for (a, b) in im:
            PP.append((a / geo.D, b / geo.D)); PO.append(i)
    PP = np.array(PP, dtype=np.float64).reshape(-1, 2); PO = np.array(PO, dtype=np.int64); PS = np.array(PS, dtype=np.float64)

    def price_pool(master, y, top):
        """Reduced costs of the inactive pool orbits by the core test at positive-dual rows."""
        if not len(pool):
            return []
        f = np.zeros(len(pool))
        for r in np.nonzero(y > 1e-12)[0]:
            _k, pcx, pcy, pc, ps_, phx, phy = master.P[master.active[r]]
            dx = PP[:, 0] - pcx; dy = PP[:, 1] - pcy
            u0 = pc * dx + ps_ * dy; u1 = -ps_ * dx + pc * dy
            hit = (np.abs(u0) <= phx) & (np.abs(u1) <= phy)
            if hit.any():
                f += y[r] * np.bincount(PO[hit], minlength=len(pool))
        rc = PS - f
        rc[~pool_alive] = np.inf
        order = np.argsort(rc)
        return [int(i) for i in order[:top] if rc[i] < -1e-7]

    print(f"[{time.strftime('%H:%M:%S')}] building master", flush=True)
    master = PoolMaster(geo, margin=args.margin, ipm=bool(args.ipm))
    if args.start and Path(args.start + ".pool.npz").exists() and mstate is None and not rows0:
        try:
            master.load(args.start + ".pool.npz")
            print(f"loaded pool of {master.n_pool()} rows", flush=True)
            x0 = np.array(list(w0) + [0] * (master.ncols - len(w0)), dtype=np.float64)[: master.ncols] / WSCALE
            na = master.warm_activate(x0, slack=0.002)
            print(f"warm-activated {na} rows", flush=True)
        except ValueError as e:
            print("pool not loaded:", e, flush=True)
    elif rows0:
        master.add_pool_rows(rows0)
    wint = list(w0) + [0] * (len(geo.reps) - len(w0))

    K = bin_count(args.N)
    focus = sorted({int(v) for v in args.focus.split(",") if v.strip()})
    tmp = tempfile.mkdtemp(dir=str(Path(args.state).parent))
    cert_path = Path(tmp) / "cand.txt"
    logf = open(args.log, "a")
    best_total = None
    pending = []
    for rnd in range(args.rounds):
        t0 = time.time()
        P_all, O_all = geo.points_array()
        wpos = np.array([wint[o] > 0 for o in range(len(geo.reps))])
        mask = wpos[O_all]
        P_int, orbit_of = P_all[mask], O_all[mask]
        PU = P_all[~mask].astype(np.float64) / geo.D; OU = O_all[~mask]
        write_cert(cert_path, geo.U, geo.D, [(int(X), int(Y), wint[o]) for (X, Y), o in zip(P_int, orbit_of)])
        off = (rnd * 7) % args.stride
        bins = sorted(set(range(off, K, args.stride if rnd > 0 else args.stride0)) | {0} | {k for k in focus if k < K} | set(pending))
        if args.always_scan > 0:
            write_cert(cert_path, geo.U, geo.D, [(int(X), int(Y), wint[o]) for (X, Y), o in zip(P_int, orbit_of)])
            _ts = time.time()
            badall, sl_all = indep_scan_all(cert_path, args.N, K, workers=args.workers)
            print(f"[{time.strftime('%H:%M:%S')}] scan: {len(badall)} bad bins, least {sl_all}, {time.time() - _ts:.0f}s", flush=True)
            logf.write(json.dumps(dict(round=rnd, scan_bad=len(badall), scan_least=sl_all, t_scan=round(time.time() - _ts, 1))) + "\n"); logf.flush()
            if not badall and sl_all is not None and sl_all >= WSCALE:
                print("CLEAN: every bin passes", flush=True)
                write_cert(args.out, geo.U, geo.D, [p for p in cert_points(geo, wint) if p[2] > 0])
                break
            ks = sorted(k for _, k in badall)
            step = max(1, len(ks) // args.always_scan)
            bins = sorted(set(bins) | set(ks[::step][: args.always_scan]) | {k for _, k in badall[:50]})
        pending = []
        thresh = int(WSCALE * (1 + (args.slack if rnd > 0 else args.slack0)))
        pb = args.per_bin if rnd > 0 else args.per_bin0
        print(f"[{time.strftime('%H:%M:%S')}] round {rnd}: dumping {len(bins)} bins, written pts {len(P_int)}", flush=True)
        res = indep_bins(cert_path, args.N, bins, thresh, 3 * pb, args.workers)
        print(f"[{time.strftime('%H:%M:%S')}] dumps done", flush=True)
        least = min((r[2] for r in res if r[2] is not None), default=None)
        worst_bins = sorted([(r[2], r[0]) for r in res if r[2] is not None])[:5]
        t1 = time.time()
        added = 0
        ncells = 0
        allrows = rows_parallel4(geo, P_int, orbit_of, [r[1] for r in res], master.keys, pb, args.workers, PU=PU, OU=OU)
        del res
        for rows in allrows:
            ncells += len(rows)
            added += master.add_pool_rows(rows)
        t2 = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] rows done: pool {master.n_pool()}, active {len(master.active)}; solving", flush=True)
        if args.tr0 > 0:
            delta = args.tr0 * (args.tr_growth ** rnd)
            nc = master.ncols
            x0 = np.array(wint[:nc] + [0] * (nc - len(wint[:nc])), dtype=np.float64) / WSCALE
            if delta < 1.0:
                lo = np.maximum(0.0, x0 - delta); hi = np.full(nc, master.inf)
            else:
                lo = np.zeros(nc); hi = np.full(nc, master.inf)
            master.h.changeColsBounds(nc, np.arange(nc, dtype=np.int32), lo, hi)
        x, y, obj = master.solve()
        t3 = time.time()
        print(f"[{time.strftime('%H:%M:%S')}] LP {obj:.6f} active {len(master.active)}", flush=True)
        newcols = 0; fmax = None
        if pool:
            for _rep in range(args.price_rounds):
                pick = price_pool(master, y, args.pool_top)
                print(f"[{time.strftime('%H:%M:%S')}] priced: {len(pick)} new", flush=True)
                if not pick:
                    break
                for i in pick:
                    geo.add_orbit(*pool[i]); pool_alive[i] = False
                newcols += len(pick)
                master.sync_cols()
                x, y, obj = master.solve()
        if args.colgen > 0:
            cands, fmax = price_grid(geo, master, y, args.delta, top=args.colgen)
            for (X, Y, val) in cands:
                if geo.add_orbit(X, Y) is not None:
                    newcols += 1
            if newcols:
                master.sync_cols()
                x, y, obj = master.solve()
        t4 = time.time()
        npurged = master.purge(args.purge) if args.purge > 0 else 0
        ndrop = 0
        if args.drop_rc > 0:
            keep = master.drop_cols(args.drop_rc)
            if not keep.all():
                ndrop = int((~keep).sum())
                g2 = Geometry(geo.U, geo.D)
                for o, kp in enumerate(keep):
                    if kp:
                        g2.add_orbit(*geo.reps[o])
                    elif geo.reps[o] in pool_index:
                        pool_alive[pool_index[geo.reps[o]]] = True
                geo = g2
                master.geo = geo
                x = x[keep]
        wint = weights_to_int(x)
        total = sum(wint[o] * len(geo.members[o]) for o in range(len(geo.reps))) / WSCALE
        rec = dict(round=rnd, off=off, bins=len(bins), least_before=least, worst=worst_bins,
                   uniq=ncells, added=added, rows=len(master.active), pool=master.n_pool(), cols=len(geo.reps), newcols=newcols,
                   fmax=fmax, lp=obj, total=total, t_dump=round(t1 - t0, 1), t_rows=round(t2 - t1, 1),
                   t_lp=round(t3 - t2, 1), t_price=round(t4 - t3, 1), purged=npurged, dropped=ndrop)
        print(json.dumps(rec), flush=True)
        logf.write(json.dumps(rec) + "\n"); logf.flush()
        write_cert(args.out, geo.U, geo.D, [p for p in cert_points(geo, wint) if p[2] > 0])
        nprune = master.prune_pool(x, args.prune) if args.prune > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] pruned {nprune} pool rows, pool {master.n_pool()}", flush=True)
        with open(args.state + ".tmp", "wb") as fh:
            pickle.dump({"geo": geo, "w": wint, "pool_alive": pool_alive}, fh)
        os.replace(args.state + ".tmp", args.state)
        master.save(args.state + ".pool.tmp.npz")
        os.replace(args.state + ".pool.tmp.npz", args.state + ".pool.npz")
        if least is not None and (least >= WSCALE or added < 120) and newcols == 0:
            # strided screen clean: full exact scan of every bin with the new weights
            write_cert(cert_path, geo.U, geo.D, [p for p in cert_points(geo, wint) if p[2] > 0])
            ts = time.time()
            bad, sl = indep_scan_all(cert_path, args.N, K, workers=args.workers)
            rec2 = dict(round=rnd, scan_bad=len(bad), scan_least=sl, scan_worst=bad[:5], t_scan=round(time.time() - ts, 1))
            print(json.dumps(rec2), flush=True); logf.write(json.dumps(rec2) + "\n"); logf.flush()
            if not bad and sl is not None and sl >= WSCALE:
                break
            # spread the pending bins over the whole bad set (worst first, then evenly)
            ks = sorted(k for _, k in bad)
            step = max(1, len(ks) // args.scan_cap)
            pending = sorted(set(ks[::step][: args.scan_cap]) | {k for _, k in bad[:20]})
        if total >= args.stop_total and rnd >= 3:
            pass
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
