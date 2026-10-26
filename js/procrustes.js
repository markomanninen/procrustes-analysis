// Procrustes fitting of one triangle onto another with complex numbers.
// Points and complex numbers are both [re, im]. ABC are the targets z, DEF are fitted, u.

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
// Minimises Σ|z_i − ρ e^{iφ} u_perm(i)|². With s = Σ z̄_i u_perm(i) the best angle is −arg s.
export function candidate(z, u, perm, { mirror = false, scaling = false } = {}) {
  const uu = mirror ? u.map(conj) : u
  let s = [0, 0]
  perm.forEach((j, i) => (s = add(s, mul(conj(z[i]), uu[j]))))
  const Su = uu.reduce((t, p) => t + p[0] ** 2 + p[1] ** 2, 0)
  const phi = -arg(s)
  const rho = scaling ? abs(s) / Su : 1
  // GeoGebra: θ = atan(y/x), which is −arg s only when Re s > 0, and off by half a turn otherwise
  const thetaAtan = s[0] === 0 ? Math.sign(s[1]) * (Math.PI / 2) : Math.atan(s[1] / s[0])
  // Residuals summed point by point: Σ|z|² + ρ²Σ|u|² − 2ρ|s| says the same but cancels badly
  const residual = (a) => perm.reduce((t, j, i) => t + abs(sub(z[i], mul(a, uu[j]))) ** 2, 0)
  const E = residual(scale(expi(phi), rho))
  const Eatan = residual(scale(expi(-thetaAtan), rho))
  return {
    perm, mirror, s, phi, rho, E, rms: Math.sqrt(E / z.length),
    thetaAtan, rmsAtan: Math.sqrt(Eatan / z.length),
    map: { a: scale(expi(phi), rho), mirror },
  }
}

export function analyse(abc, def, opts = {}) {
  const K = centroid(abc)
  const J = centroid(def)
  const z = centre(abc)
  const u = centre(def)
  const mirrors = opts.mirror ? [false, true] : [false]
  const all = mirrors.flatMap((m) => PERMS.map((p) => candidate(z, u, p, { mirror: m, scaling: opts.scaling })))
  all.forEach((c, i) => (c.id = i))
  const best = all.reduce((a, b) => (b.E < a.E ? b : a))
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
  const rho = abs(map.a)
  const phi = arg(map.a)
  const k = rho ** t
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
// (4/√3)·area / Σ|p − K|². For any pairing, |s|²/(‖z‖²‖u‖²) = (1 + n_ABC · n_DEF′)/2, so the
// best pairing is the image of DEF nearest to ABC, and two pairings tie on the great circle
// halfway between their images.
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
