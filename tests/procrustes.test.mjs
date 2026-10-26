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
