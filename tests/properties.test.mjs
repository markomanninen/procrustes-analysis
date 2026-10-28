// Randomised checks of the fit against what least squares requires, with a fixed seed so a
// failure is repeatable. Each test states the seed and case it failed on.
import { test } from 'node:test'
import assert from 'node:assert/strict'
import { fit, distance, mapAt, applyMap, relativeMap, centre, add, sub, mul, conj, expi, abs, scale } from '../js/procrustes.js'
import { analyse, shapePoint, shapeImages, dot3, softFit } from '../js/triangle.js'

const SEED = 20261028
const CASES = 300

// mulberry32
function rng(seed) {
  let a = seed >>> 0
  return () => {
    a = (a + 0x6d2b79f5) >>> 0
    let t = a
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}
const cases = (salt, body) => {
  const r = rng(SEED + salt)
  const between = (lo, hi) => lo + (hi - lo) * r()
  const int = (lo, hi) => Math.floor(between(lo, hi + 1))
  const points = (n) => Array.from({ length: n }, () => [between(-10, 10), between(-10, 10)])
  for (let k = 0; k < CASES; k++) body({ k, between, int, points }, `seed ${SEED + salt}, case ${k}`)
}

const sumsq = (pts) => pts.reduce((t, p) => t + abs(p) ** 2, 0)
const err = (a, b) => a.reduce((t, p, i) => t + abs(sub(p, b[i])) ** 2, 0)
const similarity = (pts, a, t, mirror) => pts.map((p) => add(mul(a, mirror ? conj(p) : p), t))

// The same problem with 2×2 real matrices, as in Kabsch: with H = Σ u_i z_iᵀ the best
// rotation maximises tr(R H) and the best reflection tr(F H). Derived without complex numbers.
function realMatrixBest(z, u, { mirror, scaling }) {
  let [h11, h12, h21, h22] = [0, 0, 0, 0]
  z.forEach((p, i) => {
    h11 += u[i][0] * p[0]; h12 += u[i][0] * p[1]
    h21 += u[i][1] * p[0]; h22 += u[i][1] * p[1]
  })
  const fits = []
  const t = Math.atan2(h12 - h21, h11 + h22) // R = [[c, −s], [s, c]]
  fits.push(u.map(([x, y]) => [Math.cos(t) * x - Math.sin(t) * y, Math.sin(t) * x + Math.cos(t) * y]))
  if (mirror) {
    const f = Math.atan2(h12 + h21, h11 - h22) // F = [[c, s], [s, −c]]
    fits.push(u.map(([x, y]) => [Math.cos(f) * x + Math.sin(f) * y, Math.sin(f) * x - Math.cos(f) * y]))
  }
  return Math.min(...fits.map((v) => {
    const k = scaling ? v.reduce((a, p, i) => a + p[0] * z[i][0] + p[1] * z[i][1], 0) / sumsq(v) : 1
    return err(z, v.map((p) => scale(p, k)))
  }))
}

test('the fit agrees with the real-matrix (Kabsch) solution', () => {
  cases(1, ({ int, points, r }, at) => {
    const n = int(2, 12)
    const [A, B] = [points(n), points(n)]
    for (const mirror of [false, true]) for (const scaling of [false, true]) {
      const f = fit(A, B, { mirror, scaling })
      const want = realMatrixBest(f.z, f.u, { mirror, scaling })
      assert.ok(Math.abs(f.best.E - want) <= 1e-9 * (1 + want), `${at}, mirror ${mirror}, scaling ${scaling}: ${f.best.E} vs ${want}`)
    }
  })
})

test('no nearby turn, size or shift does better than the fit', () => {
  cases(2, ({ int, points, between }, at) => {
    const n = int(2, 10)
    const [A, B] = [points(n), points(n)]
    const f = fit(A, B, { scaling: true })
    for (let j = 0; j < 8; j++) {
      const a = mul(f.best.map.a, mul(expi(between(-0.05, 0.05)), [1 + between(-0.05, 0.05), 0]))
      const shift = [between(-0.05, 0.05), between(-0.05, 0.05)]
      const other = f.u.map((w) => add(add(f.K, shift), mul(a, w)))
      assert.ok(err(other, A) >= f.best.E - 1e-9, at)
    }
  })
})

test('a similar copy is recovered whatever its turn, size, place, labels and handedness', () => {
  cases(3, ({ int, points, between }, at) => {
    const n = int(3, 9)
    const A = points(n)
    const mirror = between(0, 1) < 0.5
    const roll = int(0, n - 1)
    const reverse = between(0, 1) < 0.5
    let B = similarity(A, mul(expi(between(-Math.PI, Math.PI)), [between(0.2, 5), 0]), [between(-50, 50), between(-50, 50)], mirror)
    B = [...B.slice(roll), ...B.slice(0, roll)]
    if (reverse) B.reverse()
    const f = fit(A, B, { pairings: 'cyclic', mirror: true, scaling: true })
    assert.ok(f.best.rms < 1e-8, `${at}: rms ${f.best.rms}`)
    assert.equal(f.best.mirror, mirror, at)
    f.points.forEach((p, i) => assert.ok(abs(sub(p, A[i])) < 1e-8, at))
  })
})

test('more pairings, mirror images or scaling never fit worse', () => {
  cases(4, ({ int, points }, at) => {
    const n = int(3, 6)
    const [A, B] = [points(n), points(n)]
    const E = (o) => fit(A, B, o).best.E
    const id = E({}); const cyc = E({ pairings: 'cyclic' }); const all = E({ pairings: 'all' })
    assert.ok(all <= cyc + 1e-9 && cyc <= id + 1e-9, at)
    assert.ok(E({ mirror: true }) <= id + 1e-9, at)
    assert.ok(E({ scaling: true }) <= id + 1e-9, at)
  })
})

test('the fitted points are the best candidate applied to the source', () => {
  cases(5, ({ int, points }, at) => {
    const n = int(3, 5)
    const [A, B] = [points(n), points(n)]
    const f = fit(A, B, { pairings: 'all', mirror: true })
    assert.deepEqual([...f.best.perm].sort((a, b) => a - b), [...Array(n).keys()], at)
    assert.ok(Math.abs(err(f.points, A) - f.best.E) < 1e-9, at)
    assert.equal(f.best.E, Math.min(...f.all.map((c) => c.E)), at)
  })
})

test('distance: invariant, symmetric, between 0 and 1, and its arcsine is a metric', () => {
  cases(6, ({ int, points, between }, at) => {
    const n = int(3, 10)
    const [A, B, C] = [points(n), points(n), points(n)]
    const d = distance(A, B)
    assert.ok(d >= 0 && d <= 1, at)
    assert.ok(Math.abs(d - distance(B, A)) < 1e-12, at)
    const A2 = similarity(A, mul(expi(between(-3, 3)), [between(0.1, 10), 0]), [between(-9, 9), 0], false)
    assert.ok(Math.abs(d - distance(A2, B)) < 1e-9, at)
    assert.ok(distance(A, B, { mirror: true }) <= d + 1e-12, at)
    const rho = (x, y) => Math.asin(Math.min(1, distance(x, y)))
    assert.ok(rho(A, C) <= rho(A, B) + rho(B, C) + 1e-9, at)
    // d² is the residual of the scaled fit relative to the target's spread
    const f = fit(A, B, { scaling: true })
    assert.ok(Math.abs(d * d - f.best.E / sumsq(f.z)) < 1e-9, at)
  })
})

test('triangles: the best pairing is the nearest image on the shape sphere', () => {
  cases(7, ({ points }, at) => {
    const [A, B] = [points(3), points(3)]
    for (const mirror of [false, true]) {
      const an = analyse(A, B, { mirror })
      const dots = shapeImages(B, mirror).map((m) => dot3(shapePoint(A), m))
      const top = Math.max(...dots)
      assert.ok(top - dots[an.best.id] < 1e-12, at)
    }
  })
})

test('soft weights are a distribution and the soft fit lies within the candidates', () => {
  cases(8, ({ points, between }, at) => {
    const an = analyse(points(3), points(3), { mirror: true })
    const { weights, pts } = softFit(an, between(0.01, 5))
    assert.ok(Math.abs(weights.reduce((a, b) => a + b, 0) - 1) < 1e-12, at)
    assert.ok(weights.every((w) => w >= 0), at)
    assert.ok(weights[an.best.id] === Math.max(...weights), at)
    assert.ok(pts.every((p) => abs(sub(p, an.K)) <= Math.sqrt(sumsq(an.u)) + 1e-9), at)
  })
})

test('animation paths start at the identity, end at the map and never leap', () => {
  cases(9, ({ between }, at) => {
    const map = { a: mul(expi(between(-Math.PI, Math.PI)), [between(0.3, 3), 0]), mirror: between(0, 1) < 0.5 }
    const w = [between(-5, 5), between(-5, 5)]
    assert.ok(abs(sub(mapAt(map, 0, w), w)) < 1e-9, at)
    assert.ok(abs(sub(mapAt(map, 1, w), applyMap(map, w))) < 1e-9, at)
    let last = w
    for (let i = 1; i <= 100; i++) {
      const p = mapAt(map, i / 100, w)
      assert.ok(abs(sub(p, last)) < 0.6, at)
      last = p
    }
    const to = { a: mul(expi(between(-3, 3)), [between(0.3, 3), 0]), mirror: between(0, 1) < 0.5 }
    assert.ok(abs(sub(applyMap(relativeMap(map, to), applyMap(map, w)), applyMap(to, w))) < 1e-9, at)
  })
})

test('degenerate input: two points, collinear points, a collapsed source', () => {
  const two = fit([[0, 0], [2, 0]], [[5, 5], [5, 8]], { scaling: true })
  assert.ok(two.best.rms < 1e-12)
  const line = [[0, 0], [1, 0], [3, 0], [4, 0]]
  assert.ok(fit(line, line.map(([x]) => [0, 2 * x]), { scaling: true }).best.rms < 1e-12)
  const collapsed = fit(line, line.map(() => [7, 7]), { scaling: true })
  assert.ok(Number.isFinite(collapsed.best.E))
  assert.ok(Math.abs(collapsed.best.E - sumsq(centre(line))) < 1e-12)
  collapsed.points.forEach((p) => assert.ok(abs(sub(p, [2, 0])) < 1e-12))
})

test("'best' finds the same best pairing as trying every ordering, for any shapes", () => {
  cases(10, ({ int, points, between }, at) => {
    const n = int(2, 7)
    const A = points(n)
    const kind = int(0, 2)
    // a noisy relabelled copy, a symmetric regular polygon (many ties), or unrelated points
    const B = kind === 0
      ? [...similarity(A, expi(between(-3, 3)), [4, -2], false).map((p) => add(p, [between(-0.3, 0.3), between(-0.3, 0.3)]))].sort(() => between(-1, 1))
      : kind === 1 ? Array.from({ length: n }, (_, i) => expi((2 * Math.PI * i) / n)) : points(n)
    for (const mirror of [false, true]) for (const scaling of [false, true]) {
      const all = fit(A, B, { pairings: 'all', mirror, scaling }).best.E
      const f = fit(A, B, { pairings: 'best', mirror, scaling })
      assert.ok(Math.abs(f.best.E - all) <= 1e-9 * (1 + all), `${at}, kind ${kind}, mirror ${mirror}, scaling ${scaling}`)
      assert.ok(f.all.length <= (mirror ? 2 : 1) * Math.max(1, n * (n - 1)), `${at}: ${f.all.length} vertices`)
    }
  })
})

test("'best' recovers a shuffled, turned, noisy copy of many points", () => {
  const r = rng(SEED + 11)
  const n = 40
  const A = Array.from({ length: n }, () => [20 * r() - 10, 20 * r() - 10])
  const order = [...Array(n).keys()].sort(() => r() - 0.5)
  const B = order.map((i) => add(mul(expi(2.2), A[i]), [30 + 0.05 * r(), -5 + 0.05 * r()]))
  const f = fit(A, B, { pairings: 'best' })
  // target point i is paired with the source point that came from it
  f.best.perm.forEach((j, i) => assert.equal(order[j], i))
  assert.ok(f.best.rms < 0.05)
})

test('triangles: the side-length rule picks the best pairing without trying any', async () => {
  const { trianglePairing } = await import('../js/triangle.js')
  cases(12, ({ points }, at) => {
    const [A, D] = [points(3), points(3)]
    for (const scaling of [false, true]) {
      const an = analyse(A, D, { mirror: true, scaling })
      const { perm, mirror } = trianglePairing(A, D)
      const c = an.all.find((x) => x.mirror === mirror && String(x.perm) === String(perm))
      assert.ok(c.E - an.best.E <= 1e-9 * (1 + an.best.E), `${at}, scaling ${scaling}`)
    }
  })
})

test('holonomy: carried round a loop of shapes, a triangle turns by minus half the solid angle', async () => {
  const { transport } = await import('../js/procrustes.js')
  const { solidAngle } = await import('../js/triangle.js')
  cases(13, ({ k, between }, at) => {
    if (k >= 20) return // each case carries 1500 frames
    const [cx, cy, R, B] = [between(-0.5, 0.5), between(0.5, 1.5), between(0.1, 1.4), between(0.5, 2)]
    const N = 1500
    const dir = between(0, 1) < 0.5 ? 1 : -1
    const frames = Array.from({ length: N + 1 }, (_, i) => {
      const t = (dir * 2 * Math.PI * i) / N
      return [[-B, 0], [B, 0], [cx + R * Math.cos(t), cy + R * Math.sin(t)]]
    })
    const carried = transport(frames)
    const [first, last] = [carried[0], carried[N]]
    const turn = Math.atan2(...[1, 0].map((k) => first.reduce((s, p, i) => s + mul(conj(p), last[i])[k], 0)))
    const omega = solidAngle(frames.slice(0, N).map(shapePoint))
    const diff = Math.atan2(Math.sin(turn + omega / 2), Math.cos(turn + omega / 2))
    assert.ok(Math.abs(diff) < 1e-4, `${at}: turn ${turn}, −Ω/2 ${-omega / 2}`)
    // and the carried frames keep their shapes (distance takes a square root of rounding)
    assert.ok(distance(last, frames[N]) < 1e-6, at)
  })
})
