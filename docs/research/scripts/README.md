# Verification scripts

These scripts check numerical claims made in the research notes. They were written during the
research session of 27–28 October 2026 and copied here from the session's scratch folder. The
only changes are the import paths, which are now relative, and a header comment in each file.
`identities_check.py` was called `check.py`.

On 28 October the triangle functions (`analyse`, `softFit`, the shape sphere) moved from
`js/procrustes.js` to `js/triangle.js` (commit 8ff5a17). So `rank_rule.mjs` imports `analyse`
from `'../../../js/triangle.js'`. The two jump tests import both modules and merge them into one
namespace, `Q`, which leaves their code unchanged.

The scratch helpers that only fetched or printed API records were not copied, and neither was a
debugging variant of the rank test.

Run every script from the repository root. The JavaScript files need Node with ES modules
(tested with Node 24.22.1). The Python files need numpy (tested with Python 3.12 and numpy 2.4.2).
Nothing was installed for these runs.

The outputs below are from runs on 28 October 2026 against `js/procrustes.js` and
`js/triangle.js` as of commit 8ff5a17. In that code the scale field is called `beta` (formerly
`rho`). None of the scripts reads that field.
The `rho` in `identities_check.py` is Dryden's Procrustes distance ρ, not the scale.

| script | checks | notes file |
|---|---|---|
| [rank_rule.mjs](rank_rule.mjs) | the side-length rule picks the best of the 12 fits | [novelty_jump_locus.md](../research_notes/Procrustes%20and%20Barbour%20best%20matching/novelty_jump_locus.md), claim B |
| [jump_test.mjs](jump_test.mjs) | where the hard fit jumps, and the soft fit at a mirror wall | [correspondence_soft_alignment.md](../research_notes/Procrustes%20and%20Barbour%20best%20matching/correspondence_soft_alignment.md), novelty claim A |
| [jump_test2.mjs](jump_test2.mjs) | the jump positions along a line, mirrors off and on | novelty claim A |
| [holonomy_check.py](holonomy_check.py) | chained Procrustes fits round a loop rotate by Ω/2 | [shape_sphere_bridge.md](../research_notes/Procrustes%20and%20Barbour%20best%20matching/shape_sphere_bridge.md), §4 |
| [polygon_check.py](polygon_check.py) | Ω/2 of the sampled loop, computed independently | shape_sphere_bridge.md, §4 |
| [identities_check.py](identities_check.py) | the Procrustes distance and singular-value identities | [procrustes_history_mathematics.md](../research_notes/Procrustes%20and%20Barbour%20best%20matching/procrustes_history_mathematics.md), §4 Inferences |

## rank_rule.mjs

    node docs/research/scripts/rank_rule.mjs

The claim: with mirrors allowed, the best of the 12 candidate fits pairs the vertices of ABC and
DEF by the rank of the opposite side length, and mirrors exactly when ABC and the relabelled DEF
have opposite signed area. The script draws 20,000 random pairs of triangles (vertices uniform in
[−5, 5]²) and compares the rule with exhaustive search, without and with scaling. It also asks how
often ranking alone picks the best fit when mirrors are off.

    scaling=false: mirrors on, rule = best in 20000/20000 (same rms in 20000/20000, worst rms gap 0.00e+0); mirrors off, ranking alone = best in 12169/20000
    scaling=true: mirrors on, rule = best in 20000/20000 (same rms in 20000/20000, worst rms gap 0.00e+0); mirrors off, ranking alone = best in 11921/20000

The rule agreed in 40,000 of 40,000 trials. With mirrors off, ranking alone was best in about 60%.
There is no random seed, so the mirrors-off counts change a little from run to run. A run before
the module split gave 11913 and 12015. The jump tests below gave the same output before and after
the split.

## jump_test.mjs

    node docs/research/scripts/jump_test.mjs

Each line moves one triangle through a symmetric shape, from −10⁻⁴ to +10⁻⁴, and reports the
change of best candidate and the largest jump of a fitted vertex. The last three lines put ABC
exactly on a collinear shape and print the two largest soft-fit weights.

    ABC crosses isosceles: best id 0 -> 10; max vertex-partner jump = 0.1169
    ABC crosses collinear: best id 6 -> 0; max vertex-partner jump = 2.6742
    DEF crosses isosceles: best id 9 -> 2; max vertex-partner jump = 0.0000
    DEF crosses collinear: best id 3 -> 9; max vertex-partner jump = 0.0000
    collinear ABC, sigma=0.3: top weights 0:0.306 6:0.306; area(soft fit)=-0.0000 vs area(DEF)=3.0800
    collinear ABC, sigma=1: top weights 0:0.113 6:0.113; area(soft fit)=-0.0000 vs area(DEF)=3.0800
    collinear ABC, sigma=3: top weights 0:0.086 6:0.086; area(soft fit)=0.0000 vs area(DEF)=3.0800

When ABC crosses an isosceles or collinear shape, the fit jumps. When DEF crosses one, the best
candidate changes but the fitted triangle does not move: only the labels change. With ABC exactly
collinear, a fit and its mirror image get equal weight for every σ, and the soft fit is flat
(area 0, against 3.08 for DEF).

## jump_test2.mjs

    node docs/research/scripts/jump_test2.mjs

A = (0, 0) and B = (3, 0) are fixed. The apex C of ABC runs along y = 1.7 from x = −1 to x = 4 in
4,000 steps. DEF is fixed and scalene. Every change of best candidate is listed with the size of
the jump.

    mirror=false: apex x from -1 to 4 at height 1.7 (isosceles at x=1.5):
      x≈0.134 id 2->0 jump 0.433
      x≈2.172 id 0->1 jump 0.408
    mirror=true: apex x from -1 to 4 at height 1.7 (isosceles at x=1.5):
      x≈0.529 id 9->0 jump 0.329
      x≈1.500 id 0->10 jump 0.118
      x≈2.473 id 10->1 jump 0.329

With mirrors on, all three switches are where ABC is isosceles: |BC| = |AB| at
x = 3 − √6.11 ≈ 0.528, |CA| = |CB| at x = 1.5, and |CA| = |AB| at x = √6.11 ≈ 2.472. The
printed x is the first grid point past each crossing. With mirrors off, the two switches are at
other places, which depend on DEF.

## holonomy_check.py

    python3 docs/research/scripts/holonomy_check.py

A triangle is carried round a circular loop on the shape sphere. Each shape is fitted to the
previous fitted one by the best mass-weighted rotation. The net rotation after the loop is
compared with half the solid angle enclosed, Ω/2 = π(1 − cos α) for a cap of angular radius α.
Loops are centred on the pole (the equilateral shape when the masses are equal) and on a point of
the equator (a collinear shape). The masses are (1, 1, 1) and (1, 2, 5), and the loop is cut into
N = 12, 100 and 2000 steps. "one-shot" fits the last shape directly onto the first.

    m=(1, 1, 1) loop=pole    alpha=0.3 Omega/2=+0.14031 chained(N=12,100,2000)=+0.13439, +0.14023, +0.14031  one-shot=-1.4e-17
    m=(1, 1, 1) loop=pole    alpha=1.0 Omega/2=+1.44418 chained(N=12,100,2000)=+1.41630, +1.44379, +1.44418  one-shot=+8.7e-18
    m=(1, 1, 1) loop=pole    alpha=2.0 Omega/2=-1.83423 chained(N=12,100,2000)=-1.80902, -1.83387, -1.83423  one-shot=-5.1e-18
    m=(1, 1, 1) loop=equator alpha=0.3 Omega/2=+0.14031 chained(N=12,100,2000)=-0.13439, -0.14023, -0.14031  one-shot=-1.5e-33
    m=(1, 1, 1) loop=equator alpha=1.0 Omega/2=+1.44418 chained(N=12,100,2000)=-1.41630, -1.44379, -1.44418  one-shot=+6.2e-33
    m=(1, 1, 1) loop=equator alpha=2.0 Omega/2=-1.83423 chained(N=12,100,2000)=+1.80902, +1.83387, +1.83423  one-shot=-3.1e-33
    m=(1, 2, 5) loop=pole    alpha=0.3 Omega/2=+0.14031 chained(N=12,100,2000)=+0.13439, +0.14023, +0.14031  one-shot=-8.7e-18
    m=(1, 2, 5) loop=pole    alpha=1.0 Omega/2=+1.44418 chained(N=12,100,2000)=+1.41630, +1.44379, +1.44418  one-shot=+2.9e-18
    m=(1, 2, 5) loop=pole    alpha=2.0 Omega/2=-1.83423 chained(N=12,100,2000)=-1.80902, -1.83387, -1.83423  one-shot=-8.2e-18
    m=(1, 2, 5) loop=equator alpha=0.3 Omega/2=+0.14031 chained(N=12,100,2000)=-0.13439, -0.14023, -0.14031  one-shot=+3.1e-33
    m=(1, 2, 5) loop=equator alpha=1.0 Omega/2=+1.44418 chained(N=12,100,2000)=-1.41630, -1.44379, -1.44418  one-shot=-3.1e-33
    m=(1, 2, 5) loop=equator alpha=2.0 Omega/2=-1.83423 chained(N=12,100,2000)=+1.80902, +1.83387, +1.83423  one-shot=+3.1e-33

The chained rotation converges to Ω/2, and with N = 2000 it agrees to five decimals. The
equator loops give −Ω/2 because the script runs them in the opposite sense. The result is the
same for both sets of masses. The one-shot fit of a closed loop is the identity, so the whole
rotation is holonomy. These are the numbers quoted in shape_sphere_bridge.md.

## polygon_check.py

    python3 docs/research/scripts/polygon_check.py

This computes Ω/2 for the geodesic polygon through the same N sampled points of the pole loops,
by summing the solid angles of triangles (the Van Oosterom–Strackee formula). It does not fit
anything.

    alpha=0.3 N=12: half solid angle of geodesic polygon = +0.13439
    alpha=0.3 N=100: half solid angle of geodesic polygon = +0.14023
    alpha=1.0 N=12: half solid angle of geodesic polygon = +1.41630
    alpha=1.0 N=100: half solid angle of geodesic polygon = +1.44379
    alpha=2.0 N=12: half solid angle of geodesic polygon = -1.80902
    alpha=2.0 N=100: half solid angle of geodesic polygon = -1.83387

These equal the chained values of holonomy_check.py for N = 12 and N = 100 to all printed
digits. So a chain of N exact Procrustes fits rotates by half the solid angle of the geodesic
polygon through the sampled shapes. This holds for each N, not only in the limit.

## identities_check.py

    python3 docs/research/scripts/identities_check.py

This checks, on 2,000 random configurations of 3 to 7 planar points for each identity:

- the complex full fit β̂e^{iθ̂} = w\*z/(w\*w) is no worse than a search over 4,001 turns, and the
  residual divided by ‖z‖² equals d_F²;
- for A = ZᵀW, |w\*z| = σ₁ + sign(det A)·σ₂, |wᵀz| = σ₁ − sign(det A)·σ₂ and
  det A = (|w\*z|² − |wᵀz|²)/4, so the Kabsch reflection correction fires exactly when the
  mirror fit is better;
- d_P = 2 sin(ρ/2) and d_F = sin ρ;
- on random triangles, cos²ρ = (1 + n·n′)/2 for the shape-sphere points n, n′ used by
  `shapePoint` in `js/triangle.js` (the script recomputes them itself), and reflecting a triangle
  negates the height n₃.

Output (about 14 seconds):

    all identities hold: True
    equilateral -> [-0.  0.  1.]  collinear -> [-0.684211  0.729285  0.      ]

The equilateral triangle maps to the pole and a collinear one to the equator.
