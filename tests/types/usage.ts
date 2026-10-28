// Compiled, never run: the published types must accept real use and reject misuse.
import { fit, distance, analyse, cyclicPairings, type Point, type Fit } from 'procrustes-analysis'

const a: Point[] = [[0, 0], [1, 0], [0, 1]]
const f: Fit = fit(a, a, { pairings: 'cyclic', mirror: true, scaling: true })
const rms: number = f.best.rms
const turned: Point = f.points[0]
const d: number = distance(a, a, { mirror: true })
const custom = fit(a, a, { pairings: cyclicPairings(3, { reverse: false }) })
const sphere: number = analyse(a, a).best.phi

// @ts-expect-error unknown pairing rule
fit(a, a, { pairings: 'some' })
// @ts-expect-error points are pairs
fit([[0, 0, 0]], a)
// @ts-expect-error distance has no scaling
distance(a, a, { scaling: true })

export { rms, turned, d, custom, sphere }
