import { test } from 'node:test'
import assert from 'node:assert/strict'
import {
  analyse, placed, ggbPlaced, mapAt, relativeMap, applyMap, add, sub, mul, conj, expi, abs,
  candidate, centre, GGB_PERMS, withHysteresis,
} from '../js/procrustes.js'

// The GeoGebra file's starting triangles.
const ABC = [[-14.9845, 12.5254], [-11.0133, 3.7193], [-8.158, 4.4858]]
const DEF = [[17.6394, 3.7687], [14.2589, 10.6926], [22.1603, 5.4386]]

const close = (a, b, eps = 1e-9) => abs(sub(a, b)) < eps
const moved = (pts, a, t) => pts.map((p) => add(mul(p, a), t))

test('a turned and moved copy is fitted back exactly', () => {
  const copy = moved(ABC, expi(2.1), [30, -4])
  const an = analyse(ABC, copy, {})
  assert.ok(an.best.rms < 1e-9)
  const fit = placed(an, an.best.map)
  an.best.perm.forEach((j, i) => assert.ok(close(fit[j], ABC[i])))
})

test('relabelling the copy is found among the six pairings', () => {
  const copy = moved(ABC, expi(-0.7), [5, 5])
  const shuffled = [copy[2], copy[0], copy[1]]
  const an = analyse(ABC, shuffled, {})
  assert.ok(an.best.rms < 1e-9)
  assert.deepEqual(an.best.perm, [1, 2, 0])
})

test('a mirror image is fitted only when mirror images are allowed', () => {
  const mirrored = moved(ABC.map(conj), expi(0.4), [10, 0])
  assert.ok(analyse(ABC, mirrored, { mirror: false }).best.rms > 0.5)
  const an = analyse(ABC, mirrored, { mirror: true })
  assert.ok(an.best.rms < 1e-9)
  assert.equal(an.best.mirror, true)
})

test('reversing the vertex order is not a reflection: it cannot fit a mirror image', () => {
  const mirrored = moved(ABC.map(conj), expi(0.4), [10, 0])
  const reversed = analyse(ABC, mirrored, { mirror: false }).all.filter((c) => [3, 4, 5].includes(c.id))
  for (const c of reversed) assert.ok(c.rms > 0.5)
})

test('atan(y/x) turns a half-turned copy the wrong way, atan2 does not', () => {
  const half = moved(ABC, [-1, 0], [20, 0])
  const an = analyse(ABC, half, {})
  assert.ok(an.best.rms < 1e-9)
  const ggb = ggbPlaced(an)
  const err = Math.max(...ABC.map((p, i) => abs(sub(p, ggb[i]))))
  assert.ok(err > 3, `GeoGebra result off by ${err}`)
  assert.ok(Math.abs(an.ggbPick.thetaAtan) < 1e-9)
})

test('the smallest |θ| is not the smallest error', () => {
  // ABC turned by 50° and relabelled: the true pairing needs 50°, another needs less.
  const turned = moved(ABC, expi((50 * Math.PI) / 180), [25, 0])
  const an = analyse(ABC, [turned[1], turned[2], turned[0]], {})
  assert.ok(an.best.rms < 1e-9)
  assert.ok(an.ggbPick.rmsAtan > 1, `GeoGebra pick error ${an.ggbPick.rmsAtan}`)
})

test('the sixth GeoGebra sum pairs D twice', () => {
  assert.deepEqual(GGB_PERMS[5], [0, 2, 0])
  const z = centre(ABC)
  const u = centre(DEF)
  const typo = candidate(z, u, GGB_PERMS[5])
  const meant = candidate(z, u, [0, 2, 1])
  assert.ok(abs(sub(typo.s, meant.s)) > 1)
})

test('scaling fits a shrunken copy', () => {
  const small = moved(ABC, mul(expi(1), [0.3, 0]), [0, 0])
  assert.ok(analyse(ABC, small, { scaling: false }).best.rms > 1)
  assert.ok(analyse(ABC, small, { scaling: true }).best.rms < 1e-9)
})

test('the flip path ends at the mirror map and moves in small steps', () => {
  const map = { a: mul(expi(1.3), [1.2, 0]), mirror: true }
  const w = [2, -1.5]
  assert.ok(close(mapAt(map, 1, w), applyMap(map, w)))
  assert.ok(close(mapAt(map, 0, w), w))
  let last = w
  for (let i = 1; i <= 200; i++) {
    const p = mapAt(map, i / 200, w)
    assert.ok(abs(sub(p, last)) < 0.1)
    last = p
  }
})

test('the turn path takes the shorter arc', () => {
  const map = { a: expi(3.0), mirror: false }
  const mid = mapAt(map, 0.5, [1, 0])
  assert.ok(Math.abs(Math.atan2(mid[1], mid[0]) - 1.5) < 1e-12)
  const [x, y] = mapAt({ a: expi(3.5), mirror: false }, 0.5, [1, 0])
  assert.ok(Math.abs(Math.atan2(y, x) - (3.5 - 2 * Math.PI) / 2) < 1e-12)
})

test('the relative map takes one fit onto the other, mirror or not', () => {
  const from = { a: mul(expi(0.5), [1.1, 0]), mirror: false }
  for (const to of [{ a: expi(2), mirror: false }, { a: mul(expi(-1), [0.9, 0]), mirror: true }]) {
    const rel = relativeMap(from, to)
    const w = [1.7, 0.4]
    assert.ok(close(applyMap(rel, applyMap(from, w)), applyMap(to, w)))
  }
})

test('hysteresis keeps the current pairing near a tie', () => {
  const an = analyse(ABC, DEF, { mirror: true })
  const second = [...an.all].sort((a, b) => a.E - b.E)[1]
  const gap = (second.E - an.best.E) / second.E
  assert.equal(withHysteresis(an, second.id, gap + 0.01).id, second.id)
  assert.equal(withHysteresis(an, second.id, gap / 2).id, an.best.id)
})

test('on the shape sphere |s|² is (1 + n·m)/2 of the norms, for every pairing', async () => {
  const { shapePoint, shapeImages, dot3 } = await import('../js/procrustes.js')
  const z = centre(ABC)
  const u = centre(DEF)
  const norm2 = (p) => p.reduce((t, q) => t + abs(q) ** 2, 0)
  const an = analyse(ABC, DEF, { mirror: true })
  shapeImages(DEF, true).forEach((m, k) => {
    const lhs = abs(an.all[k].s) ** 2 / (norm2(z) * norm2(u))
    assert.ok(Math.abs(lhs - (1 + dot3(shapePoint(ABC), m)) / 2) < 1e-12)
  })
})

test('the best pairing is the image of DEF nearest to ABC on the sphere', async () => {
  const { shapePoint, shapeImages, dot3 } = await import('../js/procrustes.js')
  const n = shapePoint(ABC)
  for (const opts of [{ mirror: false }, { mirror: true }, { mirror: true, scaling: true }]) {
    const an = analyse(ABC, DEF, opts)
    const dots = shapeImages(DEF, opts.mirror).map((m) => dot3(n, m))
    assert.equal(dots.indexOf(Math.max(...dots)), an.best.id)
  }
})

test('relabelling turns the sphere, mirroring reflects it through the equator', async () => {
  const { shapePoint } = await import('../js/procrustes.js')
  const [x, y, h] = shapePoint(DEF)
  const cyc = shapePoint([DEF[1], DEF[2], DEF[0]])
  const turn = Math.atan2(cyc[1], cyc[0]) - Math.atan2(y, x)
  assert.ok(Math.abs(Math.cos(turn) - Math.cos((2 * Math.PI) / 3)) < 1e-12)
  assert.ok(Math.abs(cyc[2] - h) < 1e-12)
  const mir = shapePoint(DEF.map(conj))
  assert.ok(Math.abs(mir[0] - x) + Math.abs(mir[1] - y) + Math.abs(mir[2] + h) < 1e-12)
  const eq = shapePoint([[0, 0], [1, 0], [0.5, Math.sqrt(3) / 2]])
  assert.ok(Math.abs(eq[2] - 1) < 1e-12)
})

const place = (pts, t, a) => pts.map((p) => add(mul(p, expi(a)), t))
const gap = (A, D, opts) => {
  const e = analyse(A, D, opts).all.map((c) => c.E).sort((x, y) => x - y)
  return (e[1] - e[0]) / e[1]
}

test('with mirror images allowed, the best fit ties exactly when ABC is isosceles or flat', () => {
  const isosceles = place([[-2, 0], [2, 0], [0, 3]], [4, -1], 0.8)
  const flat = place([[0, 0], [3, 0], [-1.2, 0]], [-2, 5], 2.1)
  for (const D of [DEF, [[1, 2], [5, -1], [0, 7]], [[0, 0], [1, 3], [4, 1]]]) {
    assert.ok(gap(isosceles, D, { mirror: true }) < 1e-12)
    assert.ok(gap(flat, D, { mirror: true }) < 1e-12)
    assert.ok(gap(ABC, D, { mirror: true }) > 1e-3)
  }
})

test('without mirror images an isosceles ABC is no tie: the edges then move with DEF', () => {
  const isosceles = place([[-2, 0], [2, 0], [0, 3]], [4, -1], 0.8)
  assert.ok(gap(isosceles, DEF, { mirror: false }) > 1e-3)
})

test('the soft fit becomes the hard fit as σ goes to zero', async () => {
  const { softFit } = await import('../js/procrustes.js')
  const an = analyse(ABC, DEF, { mirror: true })
  const hard = placed(an, an.best.map)
  softFit(an, 1e-3).pts.forEach((p, i) => assert.ok(close(p, hard[an.best.perm[i]], 1e-9)))
})

test('across a tie the soft fit moves in small steps while the hard fit leaps', async () => {
  const { softFit } = await import('../js/procrustes.js')
  // move C across the isosceles shape (apex C on the perpendicular bisector of AB)
  const frames = Array.from({ length: 201 }, (_, i) => [[-2, 0], [2, 0], [-0.5 + i / 200, 3]])
  const D = [[10, 0], [11, 5], [18, 1]]
  const step = (fit) => Math.max(...frames.slice(1).map((f, i) => Math.max(...[0, 1, 2].map((v) => abs(sub(fit(f)[v], fit(frames[i])[v]))))))
  const hardFit = (f) => {
    const an = analyse(f, D, { mirror: true })
    const p = placed(an, an.best.map)
    return [0, 1, 2].map((i) => p[an.best.perm[i]])
  }
  const soft = (f) => softFit(analyse(f, D, { mirror: true }), 0.6).pts
  assert.ok(step(hardFit) > 1, `hard ${step(hardFit)}`)
  assert.ok(step(soft) < 0.15, `soft ${step(soft)}`)
})

test('the fitted triangle jumps when ABC crosses an isosceles or flat shape, not when DEF does', () => {
  const A0 = [[0, 0], [4, 1], [1, 3]]
  const D0 = [[10, 0], [11, 5], [18, 1]]
  const fitted = (A, D) => {
    const an = analyse(A, D, { mirror: true })
    return placed(an, an.best.map)
  }
  // the largest move of the fitted triangle as a set of points, between neighbouring frames
  const setStep = (frames) => {
    let worst = 0
    for (let i = 1; i < frames.length; i++) {
      const [P, R] = [fitted(...frames[i - 1]), fitted(...frames[i])]
      const best = Math.min(...[[0, 1, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0], [1, 0, 2], [0, 2, 1]]
        .map((pm) => Math.max(...pm.map((j, k) => abs(sub(P[k], R[j]))))))
      worst = Math.max(worst, best)
    }
    return worst
  }
  const sweep = (make) => Array.from({ length: 401 }, (_, i) => make(i / 400))
  const isosceles = (t) => [[-2, 0], [2, 0], [-0.5 + t, 3]]
  const flat = (t) => [[0, 0], [4, 0], [1.5, -0.5 + t]]
  for (const shape of [isosceles, flat]) {
    assert.ok(setStep(sweep((t) => [shape(t), D0])) > 1)
    assert.ok(setStep(sweep((t) => [A0, shape(t)])) < 0.01)
  }
})
