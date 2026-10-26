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

[docs/reflections.md](docs/reflections.md) explains why fitted triangles jump, especially
around reflections, and how the page animates them without jumping.

## Running

```bash
npm run serve        # python3 -m http.server 8138
open http://localhost:8138
```

## Tests

```bash
npm test
```

## Files

| | |
|---|---|
| `js/procrustes.js` | complex arithmetic, the 12 candidates, the GeoGebra rule, paths for animation; no DOM |
| `js/app.js` | the page |
| `index.html` | layout and styles |
