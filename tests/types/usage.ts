// Compiled, never run: the published types must accept real use and reject misuse.
import { fit, distance, cyclicPairings, type Point, type Fit } from 'procrustes-analysis'
import { analyse, softFit, type TriangleAnalysis } from 'procrustes-analysis/triangle'

const a: Point[] = [[0, 0], [1, 0], [0, 1]]
const f: Fit = fit(a, a, { pairings: 'cyclic', mirror: true, scaling: true })
const rms: number = f.best.rms
const turned: Point = f.points[0]
const d: number = distance(a, a, { mirror: true })
const custom = fit(a, a, { pairings: cyclicPairings(3, { reverse: false }) })
const an: TriangleAnalysis = analyse(a, a, { mirror: true })
const sphere: number = an.best.thetaAtan
const weights: number[] = softFit(an, 0.5).weights

// @ts-expect-error unknown pairing rule
fit(a, a, { pairings: 'some' })
// @ts-expect-error points are pairs
fit([[0, 0, 0]], a)
// @ts-expect-error distance has no scaling
distance(a, a, { scaling: true })
// @ts-expect-error GeoGebra's angle belongs to the triangle analysis only
f.best.thetaAtan

export { rms, turned, d, custom, sphere, weights }
