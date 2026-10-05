# s(12) ≥ 7943/2000 = 3.9715

DOI: https://doi.org/10.5281/zenodo.23106581 (all versions; v1.0: https://doi.org/10.5281/zenodo.23106582)

**Claim.** Twelve unit squares cannot be packed — interiors pairwise disjoint, each square
rotated freely — into any square of side smaller than

    7943/2000 = 3.9715

that is, **s(12) ≥ 7943/2000**.  The previous best lower bound known to us is
s(12) ≥ 15680000/3949423 = 3.97020020… (Joshua Levy's Squares Project, October 2026,
[jlevy/squares](https://github.com/jlevy/squares), T-079); the improvement is exactly
10266889/7898846000 ≈ 0.0013.  The conjectured value is s(12) = 4.

This is a computer-assisted result.  It has **not** been peer reviewed.  A write-up with the
verification lemmas and their proofs, the computations and a discussion of the limits is in
[`paper/s12_lower_3.9715.pdf`](paper/s12_lower_3.9715.pdf) (LaTeX source next to it).

## What is new, and what is not

* **The points are Evan Daniel's.**  The certificate keeps the 1736 points of Daniel's
  `s12/certificates/s12_lower_3.9686.txt` (evand/square-packing, commit `7d6f46d`), moved to the
  larger container by the dilation with factor (7943/2000)/(15680/3951) = 31382793/31360000 ≈
  1.000727 and rounded to the grid 1/4000000.  One representative of each of the 223 orbits of
  the 8 symmetries of the square was rounded and the orbit rebuilt exactly, so the set is still
  symmetric; every coordinate is within 11/39200000 of the exact image.
* **The weights are new.**  They were recomputed by linear programming.  The idea of keeping
  Daniel's points and solving for new weights on a fine angle net is Route B of Joshua Levy's
  Squares Project
  ([research note](https://github.com/jlevy/squares/blob/32cd041c54c0e6ecebb7b3cf06e1fa0d3492c74f/docs/project/research/research-2026-10-02-s12-beyond-rescaling.md),
  [review](https://github.com/jlevy/squares/blob/4b3c0573ace10d7fb8c6175ff1b677dd714ebac3/docs/project/reviews/review-2026-10-02-s12-reweighted-certificate.md)),
  which gave 15680000/3949423.  Route B takes its rows from every 96th angle bin, with the
  offset rotating between rounds, and checked its weights by complete sweeps (the first refused
  six bins, which three focused rounds repaired).  Here the complete exact scan is part of every
  cycle of the search: it runs whenever the sampled bins are clean or a round adds few rows, and
  the violated bins it finds feed the next round until a scan is clean.  In the run that
  produced the certificate the scans found 695, 675, 67 and finally 0 violated bins, although
  the sampled bins had been clean before the first and the third of them
  (`search/lp_3.9715.jsonl`; Section 4.1 of the paper).
* **How the weights were found plays no role in the proof**, which consists of the verification
  below.

## The certificate

`s12_lower_3.9715.txt`, in Daniel's format (`s12/certificates/FORMAT.md` in his repository):

```
15886000 4000000   # container [0, 15886000/4000000]^2 = [0, 7943/2000]^2
4000000            # coordinate denominator D
10000000           # weight denominator W
1736               # number of points
X Y w              # 1736 lines, integers: point (X/D, Y/D) with weight w/W
```

Total weight `119974808/10^7 = 11.9974808 < 12`; all 1736 weights are positive; the weighted set
is invariant under the 8 symmetries of the square; all points lie in the container.  SHA-256 in
`SHA256SUMS`.

## Why the certificate proves the bound

The file asserts: *every closed unit square Q contained in [0, s]² (any centre, any angle)
contains certificate points of total weight ≥ 1.*  Given that, suppose twelve unit squares with
disjoint interiors fit in a square of side s′ < s.  Scale the picture by s/s′ > 1: the
container becomes [0, s]² and each square becomes a square of side s/s′ > 1, still with
disjoint interiors.  The concentric closed unit square of each lies in its open interior, so
these twelve closed unit squares are pairwise disjoint and lie in [0, s]².  Each captures
weight ≥ 1 and no point is counted twice, so 12 ≤ total weight = 11.9974808 — a contradiction.
(This is Daniel's reduction; it is formalised in Lean in his repository,
`s12/lean/Sqpack/Basic.lean`.)

The covering property is checked by computer over the whole continuum of positions and angles:

* **Angle bins.** θ_k = 2·arctan(k/N), so cos θ_k and sin θ_k are rational.  A unit square at
  any angle in [θ_k, θ_{k+1}] contains the concentric square of side σ_k = 1/(cos δ + sin δ) at
  angle θ_k, δ = θ_{k+1} − θ_k (in the frame of the unit square the smaller one is tilted by
  φ ∈ [0, δ] and has half-extent (σ/2)(cos φ + sin φ) ≤ 1/2, as cos φ + sin φ increases on
  [0°, 45°]).
* **Admissible centres.** A unit square at angle θ ∈ [0°, 90°] lies in [0, s]² iff its centre
  is in [w/2, s − w/2]², w = cos θ + sin θ; w is concave, so over a bin the union of these
  boxes is the box for min(w(θ_k), w(θ_{k+1})).
* **Fixed angle.** At angle θ_k a point is captured iff the centre lies in a closed square
  around it, so the captured weight is a step function of the centre and its exact minimum over
  the admissible box is found by sweeping the arrangement.

## Verification

| checker | angle range | N | minimum captured weight | verdict | log (`logs/3.9715/`) |
|---|---|---|---|---|---|
| Daniel's `verify` (commit `7d6f46d`), overflow checks on | [0°,45°] + exact D4 check | 24000 | 6737611/10⁷ | rejected | `daniel_verify_ovf_N24000.log` |
| Daniel's `verify`, overflow checks on | [0°,45°] + exact D4 check | **96000** | **10000050/10⁷ = 1.0000050** | **VERIFIED** | `daniel_verify_ovf_N96000.log` |
| `tools/indep_check.cpp` | **[0°,90°], no symmetry assumed** | 24000 | 6737611/10⁷ | rejected | `indep_check_N24000.log` |
| `tools/indep_check.cpp` | [0°,90°] | 48000 | 9288526/10⁷ | rejected | `indep_check_N48000.log` |
| `tools/indep_check.cpp` | [0°,90°] | **96000** | **10000050/10⁷** | **VERIFIED** | `indep_check_N96000.log` |
| `tools/indep_check.cpp` | [0°,90°] | 192000 | 10000050/10⁷ | VERIFIED | `indep_check_N192000.log` |

The bound rests on the two verified runs at N = 96000; either one alone is a complete check.
Every minimum in the table is attained at θ = 0, and on the nets where both programs were run
they report the same minimum.  On N = 96000 the reported minimum is 10000050/10⁷ in every one
of the 39765 bins covering [0°,45°] (`logs/3.9715/tightscan_N96000.txt.gz`); in the bins we
inspected (k = 0, 12000, 30000, 39000; the cell dumps are not included) it is attained by
squares in a corner of the container, which capture exactly the 58 points of the certificate
in [0,1]² (total weight 10000050/10⁷).  Daniel's verifier took about 49 minutes on two threads,
the independent checker about 8 minutes on one core (`.time` files).  Daniel's verifier prints
its verdict (`VERIFIED` / `NOT VERIFIED`); read that line, since it exits with status 0 in both
cases (as the Squares Project's reviews noted).

The coarse nets fail because of rows of points that the shrunk test squares can slip past; see
"How far these points go" below.  A pass on one net is a complete check; a failure on a coarser
net says nothing against it.

**Controls** (inputs and outputs in `controls/3.9715/`; both must be rejected, and are):

| control | change | `indep_check`, N = 96000 | Daniel's `verify`, single bins |
|---|---|---|---|
| `control1_lowered_orbit.txt` | heaviest orbit (8 points, two of them in [0,1]²) lowered by 100/10⁷ | rejected, 9999850/10⁷ at θ = 0 | 9999850/10⁷ at k = 0 and at k = 30000 (θ ≈ 34.7°) |
| `control2_scaled_further.txt` | the same integers over D = 3999600 (side 3.97189…) | rejected, 6737611/10⁷ at θ = 0 | 6737611/10⁷ at k = 0; 10000050/10⁷ at k = 30000 |

A single-bin run of Daniel's verifier says "PARTIAL RUN … cannot verify anything"; it is used here
only to show that it rejects.

### About `tools/indep_check.cpp`

Written from the statement above, not from Daniel's code; unchanged since v1.0.  Differences
from his `verify`: the full range [0°, 90°] is swept without using the symmetry of the point set;
σ_k is used exactly rather than rounded down to 10⁻⁶; the sweep is a segment tree over open
elementary intervals rather than a sliding window.  All positions are mapped to an integer grid
of pitch 1/Q with **conservative** rounding (points floored, half-sides floored and reduced by one
grid unit, centre box rounded outward), so rounding can only lower the reported minimum.  All
integer arithmetic is `__int128` with overflow checks.  It shares with Daniel's verifier only the
elementary framework above (angle net, shrink lemma, admissible box).  It is the producer's own
checker, written with the help of an AI system and not reviewed by anyone else, so it is a second
implementation of the same method rather than an independent method.  In v1.1 the search for the
weights stopped when `tools/indep_scan.cpp` (the same checker restricted to a range of bins)
found no violation on [0°,45°], so for those bins its acceptance was built into the search;
Daniel's verifier is the check the search did not anticipate.

## Reproduce

```sh
# Daniel's verifier
git clone https://github.com/evand/square-packing
cd square-packing && git checkout 7d6f46d && cd s12/verify
CARGO_PROFILE_RELEASE_OVERFLOW_CHECKS=true cargo build --release
./target/release/verify /path/to/s12_lower_3.9715.txt 12 96000 2 0     # ~50 min on 2 threads

# independent checker
g++ -O2 -o indep_check tools/indep_check.cpp
./indep_check s12_lower_3.9715.txt 96000                               # ~8 min, one core
```

`logs/3.9715/` contains the output of every run listed above.

## How far these points go

Daniel's points include rows on the lines x = 1 − 3/3951 and x = 2(1 − 3/3951) and their images
(in this certificate at x = 0.999967 and x = 1.9999338).  The check shrinks the tested square of
the bin θ = 0 by 1 − σ₀ ≈ 2/N, and once the container is large enough the shrunk square fits
between the rows, or between a row and the side of the container, and captures neither row.  On
N = 96000 this starts at about s = 3.97155 (on N = 24000 and 48000 earlier, which is why the
certificate is rejected there).  Consistently, the linear program for these points exceeds 12
already at s = 3.97155 (12.0052, a lower bound from 2244 constraints, for weights that respect the
symmetries) and reaches 12.0758 at s = 3.9716.  Finer nets move the threshold towards
560/141 = 3.97163…, where actual unit squares start to fit between the rows; we did not try.
Attempts to add new points (column generation) were inconclusive on our hardware.  Details:
Section 6 of the paper, logs in `search/`.

## Earlier results

| lower bound for s(12) | where |
|---|---|
| 15680/3951 = 3.9686155 | Evan Daniel, [evand/square-packing](https://github.com/evand/square-packing), August 2026 |
| 31360/7901 = 3.9691178 | this repository, **v1.0**: Daniel's certificate rescaled by 7902/7901 |
| 1568000/395039 = 3.9692284 | Joshua Levy's Squares Project, Route A: Daniel's certificate rescaled further, N = 96000 |
| 15680000/3949423 = 3.9702002 | Joshua Levy's Squares Project, Route B: Daniel's points with new weights |
| **7943/2000 = 3.9715** | this repository, **v1.1** |

### v1.0: s(12) ≥ 31360/7901 ≈ 3.969118

Released 2 October 2026 (Zenodo [10.5281/zenodo.23106582](https://doi.org/10.5281/zenodo.23106582));
its certificate, controls, logs and tools are unchanged in this repository (README.md and
SHA256SUMS were rewritten and LICENSE was extended for v1.1).  The Squares Project replayed it
and records it as [T-078](https://jlevy.github.io/squares/all-results.html#t-078)
(jlevy/squares#309).

The certificate `s12_lower_3.969118.txt` is **Evan Daniel's 1736-point weighted certificate**
`s12/certificates/s12_lower_3.9686.txt` (commit `7d6f46d`) with every coordinate and the
container multiplied by **7902/7901**; no point is added, moved independently or reweighted.
Line by line it is the original with `X' = 2X`, `Y' = 2Y`, `D' = 7901`, `s' = 31360/7901`,
since `(X/3951)·(7902/7901) = 2X/7901` and `(15680/3951)·(7902/7901) = 31360/7901` exactly.
Total weight `119738036/10^7 = 11.9738036 < 12`.  The improvement over 15680/3951 is exactly
15680/31216851 ≈ 0.000502.

The new observation was about verification.  Daniel's verifier checks the covering condition on
a rational angle net of N bins and, inside each bin, shrinks the square by the factor
`1/(cos δ + sin δ)` (δ = bin width); a coarse net therefore hides some slack.  The scaled
certificate is **rejected on the nets N = 6000 and N = 12000** (the nets the original was checked
on) but **verifies on N = 24000** (and, with the independent checker, also on N = 48000 and
N = 96000).

| checker | angle range | N | minimum captured weight | verdict |
|---|---|---|---|---|
| Daniel's `verify` (commit `7d6f46d`) | [0°,45°] + exact D4 check | 6000 | 9849809/10⁷ = 0.9849809 | rejected |
| Daniel's `verify` | [0°,45°] + exact D4 check | 12000 | 9867834/10⁷ = 0.9867834 | rejected |
| Daniel's `verify` | [0°,45°] + exact D4 check | **24000** | **10000056/10⁷ = 1.0000056** | **VERIFIED** |
| Daniel's `verify`, rebuilt with `-C overflow-checks=on` | [0°,45°] + exact D4 check | 24000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` | **[0°,90°], no symmetry assumed** | 6000 | 9849809/10⁷ | rejected |
| `tools/indep_check.cpp` | [0°,90°] | 12000 | 9867834/10⁷ | rejected |
| `tools/indep_check.cpp` | [0°,90°] | **24000** | **10000056/10⁷** | **VERIFIED** |
| `tools/indep_check.cpp`, grid Q = 2⁴⁰ | [0°,90°] | 24000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` | [0°,90°] | 48000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` | [0°,90°] | 96000 | 10000056/10⁷ | VERIFIED |

The two programs report the same minimum on every net, including the two nets on which the
certificate fails (where the minimum is attained at about 18.5° and 3.5° respectively).  As the
Squares Project's review noted, the bin numbers reported there differ by one between the two
programs; the values agree.

**Controls** (`tools/indep_check.cpp`; inputs in `controls/`, outputs in `logs/control_*.log`):

* Daniel's original certificate (side 15680/3951) at N = 6000: minimum 10000056/10⁷, the value
  he reports.
* Daniel's 56-point uniform set (`s12_56points_3.8.txt`): accepted when scaled to its critical
  container 760/199 and rejected at 1520/397 with an axis-parallel square capturing 4/5, as
  documented in his `FORMAT.md`.
* The certificate scaled one step further (container 31360/7900): rejected at N = 24000
  (minimum 9849809/10⁷) and at N = 96000 (minimum 9867834/10⁷, near 3.05°;
  `logs/control_scaled_further_31360_7900_N96000.log`, added in v1.1).
* 30 000 random admissible poses evaluated exactly (`tools/spot_check.py`): none below 1.

*Correction to the v1.0 text.*  v1.0 called the certificate "essentially critical" because the
next scaling step is rejected.  A rejection only says that the shrunk test squares capture less
than 1; it does not show that an actual unit square does.  For the rejecting bin at N = 96000,
actual unit squares at both end angles of the bin, on a 401 × 401 grid of centres around the
reported pose, capture at least 10048468/10⁷ = 1.0048468, confirmed in exact arithmetic at the
minimising centre (`tools/exact_pose.py`, `logs/exact_pose_31360_7900_N96000_k2557.log` and
`logs/exact_pose_31360_7900_N96000_k2558.log`).  How far that rescaled set is a cover was not
determined.  v1.1 goes further with new weights instead.

## Files

| path | what | since |
|---|---|---|
| `s12_lower_3.9715.txt` | certificate for 7943/2000 | v1.1 |
| `s12_lower_3.969118.txt` | certificate for 31360/7901 | v1.0 |
| `paper/` | write-up of the v1.1 result (PDF and LaTeX source) | v1.1 |
| `logs/3.9715/` | outputs of the v1.1 verification runs, with start/end times, and the per-bin scan | v1.1 |
| `logs/*.log` | outputs of the v1.0 runs (three added in v1.1: the 31360/7900 control at N = 96000 and the two exact-pose checks) | v1.0 |
| `controls/3.9715/` | altered certificates that must be rejected, with outputs | v1.1 |
| `controls/*.txt` | the v1.0 controls | v1.0 |
| `search/` | logs of the linear-programming runs (JSON lines; `lp` = optimum, `total` = rounded total weight, `rows`, `scan_bad`/`scan_least` = full exact scan, least weight times 10⁷); `colgen_*` were run on N = 48000, the others on N = 96000.  The 3.9715 run was restarted once, keeping its rows, so its log has two segments whose round numbers both start at 0 | v1.1 |
| `tools/indep_check.cpp` | the independent checker | v1.0 |
| `tools/spot_check.py`, `tools/exact_pose.py` | exact evaluation of actual unit squares | v1.0, v1.1 |
| `tools/indep_scan.cpp`, `tools/indep_dump.cpp` | the checker restricted to a bin range; a lister of near-tight cells | v1.1 |
| `tools/search/` | the re-weighting scripts (`run.py`, `push.py`) and the column-generation scripts (`colgen/`, final versions; the logged runs used earlier versions).  They contain hard-coded paths of the machine they ran on and are included for documentation | v1.1 |
| `SHA256SUMS` | checksums of the certificates, controls, logs, search logs and tools (all files except `README.md`, `LICENSE`, `.zenodo.json`, `paper/` and itself) | v1.0, extended in v1.1 |
| `.zenodo.json` | archive metadata for Zenodo | v1.1 |

## Credits

* **Evan Daniel** — the 1736 points used by every certificate in this repository, the verifier,
  the reduction and its Lean formalisation
  ([evand/square-packing](https://github.com/evand/square-packing), MIT licence; see `LICENSE`).
* **Joshua Levy's Squares Project** ([jlevy/squares](https://github.com/jlevy/squares)) — the
  idea of keeping Daniel's points and re-weighting them by linear programming on a fine angle net
  (Route B, s(12) ≥ 15680000/3949423).  Its tool `packing/devtools/s12_reweight.py` (on the branch
  of the research note above, MIT licence) and Daniel's `s12/search/tighten.py` guided
  `tools/search/`, which was written from scratch.  The project also replayed and reviewed v1.0.
* The weighted-certificate method: S. Burns and G. Massaccesi (n = 17, 2026); the unavoidable-set
  method: F. Göbel (1979), W. Stromquist (2003, s(11) ≥ 3.7889, the bound for s(12) before 2026);
  see Daniel's `s12/CREDITS.md` for the full lineage.

**Author:** Ryu Sungjoon (@squarepacker). The rescaling, the re-weighting, the verification runs and the tools in `tools/` were prepared with the help of Claude (Anthropic).

## Limitations

* Not peer reviewed; the proof rests on computation.  Two separately written checkers agree,
  but they share the elementary framework above, and the v1.1 search used a version of one of
  them (see above).
* The improvement is small and comes from re-weighting Daniel's points, not from a new
  mathematical idea.  Weighted point covers cannot reach s = 4 at all: Daniel exhibits, in exact
  arithmetic, a fractional packing of mass 12.0282 at s = 3.99, so no cover of weight < 12 exists
  at any side ≥ 3.99.  Proving s(12) = 4 needs a different argument.
* No Lean proof of these certificates.  Daniel's Lean development, which kernel-checks his
  15680/3951 certificate (`s12_ge_15680_3951` in `s12/lean/Sqpack/S12HLower.lean`), could in
  principle be applied to them.
