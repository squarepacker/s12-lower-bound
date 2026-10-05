"""Driver: cutting-plane re-weighting (+ optional column generation) at side U/D, net N."""
from __future__ import annotations

import argparse
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
from push import (Geometry, Master, BinCtx, Pose, bin_count, canon, cert_points, daniel_geometry,
                  parse_dump, price_grid, read_cert, run_verify, weights_to_int, write_cert, WSCALE,
                  DANIEL, images)


def bin_rows(geo, P_int, orbit_of, n_orb, text, known=None, per_bin=60):
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
            if known is not None and kk in known:
                continue
            area = cl[r][6]
            if area <= 0:
                continue
            pose = ctx.pose(*cl[r][1:5])
            if pose is None:
                continue
            out.append((pose, {o: c for o, c in kk}))
            taken += 1
    return out


_G = {}


def _rows_worker(text):
    g = _G
    return bin_rows(g["geo"], g["P"], g["O"], g["n"], text, known=g["known"], per_bin=g["per_bin"])


def rows_parallel(geo, P_int, orbit_of, texts, known, per_bin, workers):
    import multiprocessing as mp
    _G.update(geo=geo, P=P_int, O=orbit_of, n=len(geo.reps), known=known, per_bin=per_bin)
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
    args = ap.parse_args()

    if args.start and Path(args.start).exists():
        st = pickle.load(open(args.start, "rb"))
        geo = st["geo"]
        if geo.U != args.U or geo.D != args.D:
            # rescale the geometry (warm start at another side): rescale representatives
            g2 = Geometry(args.U, args.D)
            for (X, Y) in geo.reps:
                g2.add_orbit(round(X * args.U / geo.U), round(Y * args.U / geo.U))
            geo = g2
            w0 = st["w"][:len(geo.reps)]
            rows0 = []
        else:
            w0 = st["w"]
            rows0 = st.get("rows", [])
    else:
        geo = daniel_geometry(args.U, args.D)
        # Daniel's weights per orbit as the first probe
        s_num, s_den, D0, W0, pts = read_cert(DANIEL)
        U0 = s_num * D0 // s_den
        wd = {}
        for (x, y, w) in pts:
            wd.setdefault(canon(x, y, U0), w)
        w0 = [wd[r] for r in sorted(wd)]
        rows0 = []

    master = Master(geo, margin=args.margin, ipm=bool(args.ipm))
    master.add_rows(rows0)
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
        P_int, orbit_of = geo.points_array()
        write_cert(cert_path, geo.U, geo.D, cert_points(geo, wint))
        off = (rnd * 7) % args.stride
        bins = sorted(set(range(off, K, args.stride)) | {0} | {k for k in focus if k < K} | set(pending))
        pending = []
        thresh = int(WSCALE * (1 + args.slack))
        res = dump_bins(cert_path, args.N, bins, thresh, args.workers, tmp)
        least = min((r[2] for r in res if r[2] is not None), default=None)
        worst_bins = sorted([(r[2], r[0]) for r in res if r[2] is not None])[:5]
        t1 = time.time()
        added = 0
        ncells = 0
        allrows = rows_parallel(geo, P_int, orbit_of, [r[1] for r in res], master.keys, args.per_bin, args.workers)
        del res
        for rows in allrows:
            ncells += len(rows)
            added += master.add_rows(rows)
        t2 = time.time()
        x, y, obj = master.solve()
        t3 = time.time()
        newcols = 0; fmax = None
        if args.colgen > 0:
            cands, fmax = price_grid(geo, master, y, args.delta, top=args.colgen)
            for (X, Y, val) in cands:
                if geo.add_orbit(X, Y) is not None:
                    newcols += 1
            if newcols:
                master.sync_cols()
                x, y, obj = master.solve()
        t4 = time.time()
        wint = weights_to_int(x)
        total = sum(wint[o] * len(geo.members[o]) for o in range(len(geo.reps))) / WSCALE
        rec = dict(round=rnd, off=off, bins=len(bins), least_before=least, worst=worst_bins,
                   uniq=ncells, added=added, rows=len(master.rows), cols=len(geo.reps), newcols=newcols,
                   fmax=fmax, lp=obj, total=total, t_dump=round(t1 - t0, 1), t_rows=round(t2 - t1, 1),
                   t_lp=round(t3 - t2, 1), t_price=round(t4 - t3, 1))
        print(json.dumps(rec), flush=True)
        logf.write(json.dumps(rec) + "\n"); logf.flush()
        write_cert(args.out, geo.U, geo.D, cert_points(geo, wint))
        pickle.dump({"geo": geo, "w": wint, "rows": master.rows}, open(args.state, "wb"))
        if least is not None and (least >= WSCALE or added < 120) and newcols == 0:
            # strided screen clean: full exact scan of every bin with the new weights
            write_cert(cert_path, geo.U, geo.D, cert_points(geo, wint))
            ts = time.time()
            bad, sl = full_scan(cert_path, args.N, K, int(WSCALE * (1 + args.scan_margin)), workers=args.workers)
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
