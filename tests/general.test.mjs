import { test } from 'node:test'
import assert from 'node:assert/strict'
import {
  fit, distance, permutations, cyclicPairings, analyse, shapePoint, dot3,
  add, sub, mul, conj, expi, abs,
} from '../js/procrustes.js'

const close = (a, b, eps = 1e-9) => abs(sub(a, b)) < eps
const moved = (pts, a, t) => pts.map((p) => add(mul(p, a), t))
const HEX = [[3, 0], [1.5, 2.2], [-1.4, 2.7], [-3.1, 0.2], [-1.2, -2.4], [1.9, -2.5]]

test('n points: a turned, scaled and moved copy is fitted back exactly', () => {
  const copy = moved(HEX, mul(expi(-2.4), [0.6, 0]), [7, 3])
  assert.ok(fit(HEX, copy).best.rms > 1)
  const f = fit(HEX, copy, { scaling: true })
  assert.ok(f.best.rms < 1e-9)
  f.points.forEach((p, i) => assert.ok(close(p, HEX[i])))
})

test('n points: a mirror image needs mirror, a relabelled polygon needs cyclic pairings', () => {
  const mirrored = moved(HEX.map(conj), expi(1), [0, 0])
  assert.ok(fit(HEX, mirrored).best.rms > 0.5)
  assert.ok(fit(HEX, mirrored, { mirror: true }).best.rms < 1e-9)
  const rolled = [...HEX.slice(4), ...HEX.slice(0, 4)].reverse()
  assert.ok(fit(HEX, rolled).best.rms > 0.5)
  const f = fit(HEX, rolled, { pairings: 'cyclic' })
  assert.ok(f.best.rms < 1e-9)
  f.points.forEach((p, i) => assert.ok(close(p, HEX[i])))
})

test('the best rotation beats every other rotation', () => {
  const D = [[0, 1], [4, 2], [3, -1], [-2, -3], [-1, 4], [2, 2]]
  const f = fit(HEX, D)
  const u = f.u
  for (let k = 0; k < 360; k++) {
    const pts = u.map((w) => add(f.K, mul(expi((k * Math.PI) / 180), w)))
    const E = pts.reduce((t, p, i) => t + abs(sub(p, HEX[i])) ** 2, 0)
    assert.ok(E >= f.best.E - 1e-9)
  }
})

test('pairings: n! orderings, 2n cyclic ones, and the triangle matches analyse()', () => {
  assert.equal(permutations(5).length, 120)
  assert.equal(new Set(permutations(4).map(String)).size, 24)
  assert.throws(() => permutations(9), RangeError)
  const cyc = cyclicPairings(5)
  assert.equal(cyc.length, 10)
  assert.equal(new Set(cyc.map(String)).size, 10)
  const T = [[0, 0], [4, 1], [1, 3]]
  const D = [[5, 5], [6, 9], [9, 6]]
  const a = analyse(T, D, { mirror: true }).all.map((c) => c.E)
  const b = fit(T, D, { pairings: 'cyclic', mirror: true }).all.map((c) => c.E)
  a.forEach((e, i) => assert.ok(Math.abs(e - b[i]) < 1e-12))
})

test('distance: zero for similar shapes, symmetric, at most 1, and the shape sphere chord for triangles', () => {
  const copy = moved(HEX, mul(expi(0.9), [2.5, 0]), [-4, 1])
  assert.ok(distance(HEX, copy) < 1e-7)
  assert.ok(distance(HEX, HEX.map(conj)) > 0.1)
  assert.ok(distance(HEX, HEX.map(conj), { mirror: true }) < 1e-7)
  const D = [[0, 1], [4, 2], [3, -1], [-2, -3], [-1, 4], [2, 2]]
  assert.ok(Math.abs(distance(HEX, D) - distance(D, HEX)) < 1e-12)
  assert.ok(distance(HEX, D) <= 1)
  const T = [[0, 0], [4, 1], [1, 3]]
  const U = [[5, 5], [6, 9], [9, 6]]
  const chord = Math.sqrt((1 - dot3(shapePoint(T), shapePoint(U))) / 2)
  assert.ok(Math.abs(distance(T, U) - chord) < 1e-12)
})

test('fit rejects mismatched or degenerate input', () => {
  assert.throws(() => fit(HEX, HEX.slice(1)), RangeError)
  assert.throws(() => fit(HEX, HEX, { pairings: 'nope' }), RangeError)
  assert.throws(() => distance([[1, 1], [1, 1]], [[0, 0], [1, 0]]), RangeError)
})
