// Procrustes fitting with complex numbers: n points onto n points, and one triangle onto
// another in detail. Points and complex numbers are both [re, im]. In the triangle parts ABC
// are the targets z and DEF are fitted, u.

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

// Which vertex of DEF each of A, B, C is paired with: perm[i] indexes D, E, F.
export const PERMS = [
  [0, 1, 2], [1, 2, 0], [2, 0, 1], // the three cyclic orders keep the direction of travel
  [2, 1, 0], [1, 0, 2], [0, 2, 1], // the three others reverse it
]

// The six sums of the GeoGebra construction, s_1 … s_6, as written there. The sixth repeats
// u₁ (v₁u₁ + v₂u₃ + v₃u₁); the intended pairing is [0, 2, 1].
export const GGB_PERMS = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0], [1, 0, 2], [0, 2, 0]]

export const wrap = (t) => {
  const r = (((t + Math.PI) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI) - Math.PI
  return r === -Math.PI ? Math.PI : r
}

// The similarity map w ↦ a·w (proper) or w ↦ a·w̄ (mirror), on points centred at J.
export const applyMap = (map, w) => mul(map.a, map.mirror ? conj(w) : w)

// Fits DEF (u, centred) onto ABC (z, centred) with the pairing perm.
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
  // GeoGebra: θ = atan(y/x), which is −arg s only when Re s > 0, and off by half a turn otherwise
  const thetaAtan = s[0] === 0 ? Math.sign(s[1]) * (Math.PI / 2) : Math.atan(s[1] / s[0])
  // Residuals summed point by point: Σ|z|² + ρ²Σ|u|² − 2ρ|s| says the same but cancels badly
  const residual = (a) => perm.reduce((t, j, i) => t + abs(sub(z[i], mul(a, uu[j]))) ** 2, 0)
  const E = residual(scale(expi(phi), beta))
  const Eatan = residual(scale(expi(-thetaAtan), beta))
  return {
    perm, mirror, s, phi, beta, E, rms: Math.sqrt(E / z.length),
    thetaAtan, rmsAtan: Math.sqrt(Eatan / z.length),
    map: { a: scale(expi(phi), beta), mirror },
  }
}

export function analyse(abc, def, opts = {}) {
  const { K, J, z, u, all, best } = fit(abc, def, { pairings: 'all', mirror: opts.mirror, scaling: opts.scaling })
  const ggb = GGB_PERMS.map((p) => candidate(z, u, p, { scaling: false }))
  const ggbPick = ggb.reduce((a, b) => (Math.abs(b.thetaAtan) < Math.abs(a.thetaAtan) ? b : a))
  return { K, J, z, u, all, best, ggb, ggbPick }
}

// Points of DEF after a map, placed at ABC's centroid.
export function placed(an, map) {
  return an.u.map((w) => add(an.K, applyMap(map, w)))
}

// GeoGebra's result: DEF rotated as a whole by −θ about its own centroid, then moved to K.
export function ggbPlaced(an) {
  return placed(an, { a: expi(-an.ggbPick.thetaAtan), mirror: false })
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

// Keeps the current pairing until another is better by more than `margin` (a fraction of
// the current error). Ties between pairings are where fits jump; this stops them flickering.
export function withHysteresis(an, currentId, margin) {
  const current = an.all.find((c) => c.id === currentId)
  if (!current) return an.best
  return an.best.E < current.E * (1 - margin) - 1e-12 ? an.best : current
}

export const LETTERS = ['D', 'E', 'F']
export const pairing = (perm) => perm.map((j, i) => `${'ABC'[i]}${LETTERS[j]}`).join(' ')

// ---- Kendall's shape sphere ----------------------------------------------------------
//
// A triangle's shape, free of position, turn and size, as a unit vector. The poles are the
// two equilateral triangles, the equator the collinear ones, and the height is
// (4/√3)·area / Σ|p − K|². This is Kendall's shape sphere S²(½) enlarged to radius 1, so the
// angle between two points here is twice the Riemannian shape distance ρ. For any pairing, |s|²/(‖z‖²‖u‖²) = (1 + n_ABC · n_DEF′)/2, so the
// best pairing is the image of DEF nearest to ABC, and two pairings tie on the great circle
// halfway between their images. With mirror images allowed those circles are fixed: the
// equator and three circles of isosceles shapes. The best pairing changes when either
// triangle crosses one, but the fitted figure jumps only when ABC does; a flat or isosceles
// DEF is its own relabelled mirror image, so the two tied pairings place the same figure.
export function shapePoint(tri) {
  const w1 = scale(sub(tri[1], tri[0]), 1 / Math.SQRT2)
  const w2 = scale(sub(scale(tri[2], 2), add(tri[0], tri[1])), 1 / Math.sqrt(6))
  const n = abs(w1) ** 2 + abs(w2) ** 2
  if (n === 0) return [0, 0, 1]
  const c = mul(conj(w1), w2)
  return [(abs(w1) ** 2 - abs(w2) ** 2) / n, (2 * c[0]) / n, (2 * c[1]) / n]
}

// DEF as each candidate sees it: relabelled so its i-th vertex is paired with ABC's i-th,
// and mirrored for a mirror candidate. Ordered as analyse() numbers its candidates.
export function shapeImages(def, mirror) {
  const mirrors = mirror ? [false, true] : [false]
  return mirrors.flatMap((m) => PERMS.map((p) => shapePoint(p.map((j) => (m ? conj(def[j]) : def[j])))))
}

export const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

// ---- soft fit ------------------------------------------------------------------------
//
// Instead of the single best pairing, every pairing is weighted by e^{−(rms² − rms_best²)/σ²}
// and each vertex of ABC gets the weighted mean of where the pairings put its partner.
// The weights change continuously with the vertices, so the fit does too: at a tie the two
// pairings share the weight and the fit lies between them. σ → 0 is the hard fit. This is
// the softassign idea of Rangarajan, Chui and Bookstein; lowering σ step by step is
// deterministic annealing.
export function softFit(an, sigma) {
  const best = an.best.rms ** 2
  const raw = an.all.map((c) => (sigma > 0 ? Math.exp(-(c.rms ** 2 - best) / (sigma * sigma)) : c.id === an.best.id ? 1 : 0))
  const total = raw.reduce((a, b) => a + b, 0)
  const weights = raw.map((v) => v / total)
  const partner = (c, i) => add(an.K, applyMap(c.map, an.u[c.perm[i]]))
  const pts = [0, 1, 2].map((i) => an.all.reduce((acc, c, k) => add(acc, scale(partner(c, i), weights[k])), [0, 0]))
  return { pts, weights }
}

// ---- any number of points ------------------------------------------------------------
//
// Everything above that does not name D, E, F works for n points: candidate() takes any
// pairing of the same length as the points. The functions below choose the pairings.

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
// triangle these are all six, in the order of PERMS.
export function cyclicPairings(n, { reverse = true } = {}) {
  const forward = [...Array(n).keys()].map((k) => [...Array(n).keys()].map((i) => (i + k) % n))
  if (!reverse) return forward
  const backward = forward.map((_, k) => [...Array(n).keys()].map((i) => (((n - 1 - k - i) % n) + 2 * n) % n))
  return [...forward, ...backward]
}

const PAIRINGS = {
  identity: (n) => [[...Array(n).keys()]],
  cyclic: (n) => cyclicPairings(n),
  all: (n) => (n === 3 ? PERMS : permutations(n)),
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
// For a triangle it is the chord length on the shape sphere over two: √((1 − n·m)/2).
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
