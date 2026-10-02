# s(12) ≥ 31360/7901 ≈ 3.969118

**Claim.** Twelve unit squares cannot be packed — interiors pairwise disjoint, each square
rotated freely — into any square of side smaller than

    31360/7901 = 3.96911783318567...

that is, **s(12) ≥ 31360/7901**.  The previous best lower bound known to us is
s(12) ≥ 15680/3951 = 3.96861554… (Evan Daniel, August 2026,
[evand/square-packing](https://github.com/evand/square-packing)); the improvement is exactly
15680/31216851 ≈ 0.000502.  The conjectured value is s(12) = 4.

This is a computer-assisted result.  It has **not** been peer reviewed.

## What is new, and what is not

The certificate in this repository is **Evan Daniel's 1736-point weighted certificate**
`s12/certificates/s12_lower_3.9686.txt` (evand/square-packing, commit `7d6f46d`), with every
coordinate and the container multiplied by **7902/7901**.  No point is added, moved
independently or reweighted.

The only new observation is about verification.  Daniel's verifier checks the covering
condition on a rational angle net of N bins and, inside each bin, shrinks the square by the
factor `1/(cos δ + sin δ)` (δ = bin width); a coarse net therefore hides some slack.  The scaled
certificate is **rejected on the nets N = 6000 and N = 12000** (the nets the original was checked
on) but **verifies on N = 24000** (and, with the independent checker below, also on N = 48000 and
N = 96000).  Scaling one step further (container 31360/7900 = 3.969620) is rejected at
N = 24000, so this certificate is essentially critical; going further needs new weights or a
new idea.

## The certificate

`s12_lower_3.969118.txt`, in Daniel's format (`certificates/FORMAT.md` in his repository):

```
31360 7901      # container [0, 31360/7901]^2
7901            # coordinate denominator D
10000000        # weight denominator W
1736            # number of points
X Y w           # 1736 lines, integers
```

Line by line it is the original with `X' = 2X`, `Y' = 2Y`, `D' = 7901`, `s' = 31360/7901`,
since `(X/3951)·(7902/7901) = 2X/7901` and `(15680/3951)·(7902/7901) = 31360/7901` exactly.
Total weight `119738036/10^7 = 11.9738036 < 12`; all weights are non-negative; all points lie in
the container.  SHA-256 in `SHA256SUMS`.

## Why the certificate proves the bound

The file asserts: *every closed unit square Q contained in [0, s]² (any centre, any angle)
contains certificate points of total weight ≥ 1.*  Given that, suppose twelve unit squares with
disjoint interiors fit in a square of side s′ < s.  Scale the picture by s/s′ > 1: the
container becomes [0, s]² and each square becomes a square of side s/s′ > 1, still with
disjoint interiors.  The concentric closed unit square of each lies in its open interior, so
these twelve closed unit squares are pairwise disjoint and lie in [0, s]².  Each captures
weight ≥ 1 and no point is counted twice, so 12 ≤ total weight = 11.9738036 — a contradiction.
(This is Daniel's reduction; it is formalised in Lean in his repository, `lean/Sqpack/Basic.lean`.)

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

| checker | angle range | N | minimum captured weight | verdict |
|---|---|---|---|---|
| Daniel's `verify` (commit `7d6f46d`) | [0°,45°] + exact D4 check | 6000 | 9849809/10⁷ = 0.9849809 | rejected |
| Daniel's `verify` | [0°,45°] + exact D4 check | 12000 | 9867834/10⁷ = 0.9867834 | rejected |
| Daniel's `verify` | [0°,45°] + exact D4 check | **24000** | **10000056/10⁷ = 1.0000056** | **VERIFIED** |
| Daniel's `verify`, rebuilt with `-C overflow-checks=on` | [0°,45°] + exact D4 check | 24000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` (this repository) | **[0°,90°], no symmetry assumed** | 6000 | 9849809/10⁷ | rejected |
| `tools/indep_check.cpp` | [0°,90°] | 12000 | 9867834/10⁷ | rejected |
| `tools/indep_check.cpp` | [0°,90°] | **24000** | **10000056/10⁷** | **VERIFIED** |
| `tools/indep_check.cpp`, grid Q = 2⁴⁰ | [0°,90°] | 24000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` | [0°,90°] | 48000 | 10000056/10⁷ | VERIFIED |
| `tools/indep_check.cpp` | [0°,90°] | 96000 | 10000056/10⁷ | VERIFIED |

The two programs agree to the last integer on every net, including the two nets on which the
certificate fails (where the minimum is attained at about 18.5° and 3.5° respectively).
A run that verifies on any single net is already a complete check; the finer nets are there to
show that the result does not hinge on one discretisation.

**Controls** (`tools/indep_check.cpp`; inputs in `controls/`, outputs in `logs/control_*.log`):

* Daniel's original certificate (side 15680/3951) at N = 6000: minimum 10000056/10⁷, the value
  he reports.
* Daniel's 56-point uniform set (`s12_56points_3.8.txt`): accepted when scaled to its critical
  container 760/199 and rejected at 1520/397 with an axis-parallel square capturing 4/5, as
  documented in his `FORMAT.md`.
* The present certificate scaled one step further (container 31360/7900): rejected at
  N = 24000 (minimum 9849809/10⁷).
* 30 000 random admissible poses evaluated exactly (`tools/spot_check.py`): none below 1.

### About `tools/indep_check.cpp`

Written from the statement above, not from Daniel's code.  Differences from his `verify`:
the full range [0°, 90°] is swept without using the symmetry of the point set; σ_k is used
exactly rather than rounded down to 10⁻⁶; the sweep is a segment tree over open elementary
intervals rather than a sliding window.  All positions are mapped to an integer grid of pitch
1/Q with **conservative** rounding (points floored, half-sides floored and reduced by one grid
unit, centre box rounded outward), so rounding can only lower the reported minimum.  All
integer arithmetic is `__int128` with overflow checks.  It shares with Daniel's verifier only
the elementary framework above (angle net, shrink lemma, admissible box).

## Reproduce

```sh
# Daniel's verifier
git clone https://github.com/evand/square-packing && cd square-packing/s12/verify
cargo build --release
./target/release/verify /path/to/s12_lower_3.969118.txt 12 24000 2 0      # ~8 min on 2 cores

# independent checker
g++ -O2 -o indep_check tools/indep_check.cpp
./indep_check s12_lower_3.969118.txt 24000                                  # ~1-2 min, one core
```

`logs/` contains the output of every run listed above.

## Credits

* **Evan Daniel** — the certificate this file is derived from, the verifier, the reduction and
  its Lean formalisation ([evand/square-packing](https://github.com/evand/square-packing),
  MIT licence; see `LICENSE`).
* The weighted-certificate method: S. Burns and G. Massaccesi (n = 17, 2026); the unavoidable-set
  method: F. Göbel (1979), W. Stromquist (2003, s(11) ≥ 3.7889, the bound for s(12) before 2026);
  see Daniel's `CREDITS.md` for the full lineage.

**Author:** Ryu Sungjoon (@squarepacker). The rescaling, the verification runs and `tools/indep_check.cpp` were prepared with the help of Claude (Anthropic).

## Limitations

* Not peer reviewed; the proof rests on computation (two independent checkers agree).
* The improvement is small and comes from recovering discretisation slack, not from a new
  mathematical idea.  Weighted point covers cannot reach s = 4 at all: Daniel shows by LP duality
  that none exists at any side ≥ 3.99.
* No Lean proof of this specific certificate yet.  Daniel's Lean pose-box-tree checker, which
  kernel-checks his 15680/3951 certificate, could in principle be applied to this one.
