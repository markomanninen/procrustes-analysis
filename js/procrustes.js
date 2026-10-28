// Procrustes fitting with complex numbers: n points onto n points, by moving, turning and,
// if allowed, scaling and mirroring. Points and complex numbers are both [re, im]. The
// targets are z, the fitted points u. Triangles in detail are in triangle.js.

export const add = (a, b) => [a[0] + b[0], a[1] + b[1]]
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1]]
export const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
export const conj = (a) => [a[0], -a[1]]
export const scale = (a, k) => [a[0] * k, a[1] * k]
export const abs = (a) => Math.hypot(a[0], a[1])
export const arg = (a) => Math.atan2(a[1], a[0])
export const expi = (t) => [Math.cos(t), Math.sin(t)]
export const div = (a, b) => scale(mul(a, conj(b)), 1 / (b[0] * b[0] + b[1] * b[1]))

export function centroid(pts) {
  return scale(pts.reduce(add, [0, 0]), 1 / pts.length)
}

export const centre = (pts) => {
  const c = centroid(pts)
  return pts.map((p) => sub(p, c))
}

export const wrap = (t) => {
  const r = (((t + Math.PI) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI) - Math.PI
  return r === -Math.PI ? Math.PI : r
}

// The similarity map w ↦ a·w (proper) or w ↦ a·w̄ (mirror), on points centred at J.
export const applyMap = (map, w) => mul(map.a, map.mirror ? conj(w) : w)

// Σ|z_i − map(u_perm(i))|², summed point by point: Σ|z|² + β²Σ|u|² − 2β|s| says the same
// but cancels badly when the fit is close.
export const residual = (z, u, perm, map) => perm.reduce((t, j, i) => t + abs(sub(z[i], applyMap(map, u[j]))) ** 2, 0)

// Fits u (centred) onto z (centred) with the pairing perm: z_i goes with u_perm(i).
// Minimises Σ|z_i − β e^{iφ} u_perm(i)|². With s = Σ z̄_i u_perm(i) the best angle is −arg s.
// β is the scale, named as in Dryden and Mardia; ρ there is the Riemannian shape distance.
export function candidate(z, u, perm, { mirror = false, scaling = false } = {}) {
  const uu = mirror ? u.map(conj) : u
  let s = [0, 0]
  perm.forEach((j, i) => (s = add(s, mul(conj(z[i]), uu[j]))))
  const Su = uu.reduce((t, p) => t + p[0] ** 2 + p[1] ** 2, 0)
  const phi = -arg(s)
  // A source collapsed to one point fits equally badly at every size; keep it at 1
  const beta = scaling && Su > 0 ? abs(s) / Su : 1
  const map = { a: scale(expi(phi), beta), mirror }
  const E = residual(z, u, perm, map)
  return { perm, mirror, s, phi, beta, E, rms: Math.sqrt(E / z.length), map }
}

// The map that takes points fitted by `from` to points fitted by `to`, y ↦ c·y or y ↦ c·ȳ.
export function relativeMap(from, to) {
  if (from.mirror === to.mirror) return { a: div(to.a, from.a), mirror: false }
  return { a: div(to.a, conj(from.a)), mirror: true }
}

// A continuous path for a map, t from 0 to 1. A proper map turns along the shorter arc and
// scales geometrically. A mirror map y ↦ c·ȳ is one reflection, about the line at angle
// arg(c)/2, and is drawn as a half turn in space about that line: the part of each point
// across the line shrinks by cos(πt), passes through zero and comes back reversed.
export function mapAt(map, t, w) {
  const beta = abs(map.a)
  const phi = arg(map.a)
  const k = beta ** t
  if (!map.mirror) return scale(mul(expi(wrap(phi) * t), w), k)
  const alpha = phi / 2
  const q = mul(w, expi(-alpha))
  return scale(mul(expi(alpha), [q[0], q[1] * Math.cos(Math.PI * t)]), k)
}

// ---- pairings ------------------------------------------------------------------------

// Every ordering of 0 … n−1, in lexicographic order. n! of them, so n is capped.
export function permutations(n) {
  if (n > 8) throw new RangeError(`permutations(${n}): ${n}! pairings is too many; pass pairings yourself`)
  const out = []
  const walk = (prefix, rest) => {
    if (!rest.length) return out.push(prefix)
    rest.forEach((x, k) => walk([...prefix, x], [...rest.slice(0, k), ...rest.slice(k + 1)]))
  }
  walk([], [...Array(n).keys()])
  return out
}

// The pairings of two closed polygons that keep neighbours together: the n starting
// vertices, and with `reverse` the same n with the direction of travel reversed. For a
// triangle these are all six.
export function cyclicPairings(n, { reverse = true } = {}) {
  const forward = [...Array(n).keys()].map((k) => [...Array(n).keys()].map((i) => (i + k) % n))
  if (!reverse) return forward
  const backward = forward.map((_, k) => [...Array(n).keys()].map((i) => (((n - 1 - k - i) % n) + 2 * n) % n))
  return [...forward, ...backward]
}

const PAIRINGS = {
  identity: (n) => [[...Array(n).keys()]],
  cyclic: (n) => cyclicPairings(n),
  all: (n) => permutations(n),
}

// Fits `source` onto `target`, n points each, over the given pairings and, with `mirror`,
// their mirror images. `pairings` is 'identity' (the default: point i goes with point i),
// 'cyclic', 'all' or an array of pairings. Returns the candidates, the best one and the
// fitted source points in target order.
export function fit(target, source, { pairings = 'identity', mirror = false, scaling = false } = {}) {
  if (target.length !== source.length) throw new RangeError('target and source need the same number of points')
  if (target.length < 2) throw new RangeError('at least two points are needed')
  const list = typeof pairings === 'string' ? PAIRINGS[pairings]?.(target.length) : pairings
  if (!list) throw new RangeError(`unknown pairings '${pairings}'`)
  const K = centroid(target)
  const J = centroid(source)
  const z = centre(target)
  const u = centre(source)
  const mirrors = mirror ? [false, true] : [false]
  const all = mirrors.flatMap((m) => list.map((p) => candidate(z, u, p, { mirror: m, scaling })))
  all.forEach((c, i) => (c.id = i))
  const best = all.reduce((a, b) => (b.E < a.E ? b : a))
  const points = best.perm.map((j) => add(K, applyMap(best.map, u[j])))
  return { K, J, z, u, all, best, points }
}

// The full Procrustes distance between two shapes, √(1 − |⟨z, u⟩|² / (‖z‖²‖u‖²)): the sine
// of the angle between them on the preshape sphere, free of position, turn and size. It is
// 0 for similar shapes and 1 at most. With `mirror`, a mirror image counts as the same shape.
// For a triangle it is half the chord on the shape sphere of triangle.js: √((1 − n·m)/2).
export function distance(a, b, { mirror = false, pairing } = {}) {
  const z = centre(a)
  const u = centre(b)
  const perm = pairing ?? [...z.keys()]
  const nz = z.reduce((t, p) => t + abs(p) ** 2, 0)
  const nu = u.reduce((t, p) => t + abs(p) ** 2, 0)
  if (nz === 0 || nu === 0) throw new RangeError('a shape whose points all coincide has no direction')
  const s2 = (m) => abs(candidate(z, u, perm, { mirror: m }).s) ** 2
  const best = Math.max(s2(false), mirror ? s2(true) : 0)
  return Math.sqrt(Math.max(0, 1 - best / (nz * nu)))
}
