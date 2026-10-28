# Procrustes analysis

An interactive page fitting one triangle onto another with complex numbers. It is built on
the GeoGebra construction
[Procrustes analysis, reflections optimized](https://www.geogebra.org/m/M6drf6Ff).

The construction centres both triangles on their centroids K and J and writes the vertices
as complex numbers. For each pairing of vertices it forms s = Σ v_i·u_π(i) with v the
conjugates of ABC, and turns DEF by −θ. The page does the same and shows, side by side,
where least squares differs:

- **atan(y/x) loses the half-plane.** The correct turn is −arg s = −atan2(y, x). With atan,
  the fit is half a turn off whenever Re s < 0.
- **The smallest |θ| is not the best fit.** The best pairing is the one with the largest
  |s|, that is the smallest residual.
- **The sixth sum pairs u₁ twice** (v₁u₁ + v₂u₃ + v₃u₁). The intended sum is v₁u₁ + v₂u₃ + v₃u₂.
- **Reversing the vertex order is not a reflection.** It changes the order of travel, not the shape.
  A mirror image is fitted by conjugating u, which gives 12 candidates in all.

The page also draws Kendall's shape sphere, enlarged from its usual radius ½ to 1. There
the best pairing is the image of DEF nearest to ABC. With mirror images allowed, the best
pairing changes exactly when either triangle passes through a flat or an isosceles shape,
but the fitted triangle jumps only when ABC does. When DEF does, the tied pairings place the
same triangle and only the names of its vertices change.

[docs/reflections.md](docs/reflections.md) explains why fitted triangles jump, especially
around reflections, and how the page animates them without jumping.

## Library

`js/procrustes.js` has no dependencies and no DOM. Points are `[x, y]` pairs.

```js
import { fit, distance } from 'procrustes-analysis'

const f = fit(target, source, {
  pairings: 'identity', // or 'cyclic' (closed polygons), 'all' (n ≤ 8), or an array of pairings
  mirror: false,        // also try the mirror image of source
  scaling: false,       // also fit the size
})
f.best.rms     // root mean square residual of the best candidate
f.best.perm    // which source point went with each target point
f.points       // the fitted source points, in target order

distance(a, b, { mirror: true }) // full Procrustes distance, 0 … 1, free of position, turn and size
```

The triangle functions (`analyse`, `shapePoint`, `softFit`, …) are the special case n = 3
that the page uses.

## Command line

```bash
npx procrustes fit target.json source.json --pairings cyclic --mirror --scaling
npx procrustes distance a.json b.json --mirror
echo '[[0,0],[2,0],[0,1]]' | npx procrustes fit - '[[5,5],[5,7],[4,5]]'
```

Point sets are JSON, `[[x, y], …]` or `{"points": [[x, y], …]}`, given as a file, inline or
`-` for standard input. The output is JSON: `fit` prints the rms, the pairing, the turn, the
scale, the transform x ↦ a·x + b (a·x̄ + b when mirrored) and the fitted points. The exit
status is 0 on success, 1 for bad input and 2 for bad usage.

## Running

```bash
npm run serve        # python3 -m http.server 8138
open http://localhost:8138
```

## Tests

```bash
npm test               # unit, randomised property, command-line and package tests
npm run test:coverage  # the same, failing below 98 % of lines or 90 % of branches
npm run test:types     # the published types accept real use and reject misuse
npm run test:e2e       # the page in Chromium (npx playwright install chromium first)
```

The property tests draw 300 cases each from a fixed seed and check the fit against an
independent real-matrix (Kabsch) solution. GitHub Actions runs the first three when the
library, the command or their tests change, and the browser tests only when the page or
the library changes.

## Files

| | |
|---|---|
| `js/procrustes.js` | complex arithmetic, fitting n points, the 12 triangle candidates, the GeoGebra rule, paths for animation; no DOM |
| `js/procrustes.d.ts` | types |
| `bin/procrustes.js` | the command line |
| `js/app.js` | the page |
| `index.html` | layout and styles |

## License

MIT
