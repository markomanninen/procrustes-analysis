import type { Candidate, Pairing, Point, SimilarityMap } from './procrustes.js'

/** A candidate as analyse() returns it, with GeoGebra's angle beside the least-squares one. */
export interface TriangleCandidate extends Candidate {
  /** GeoGebra's angle atan(Im s / Re s). */
  thetaAtan: number
  /** The rms error after turning by −thetaAtan. */
  rmsAtan: number
}

export interface TriangleAnalysis {
  K: Point
  J: Point
  z: Point[]
  u: Point[]
  all: TriangleCandidate[]
  best: TriangleCandidate
  /** The six sums of the GeoGebra construction, as written there. */
  ggb: TriangleCandidate[]
  /** The sum GeoGebra picks: the smallest |θ|. */
  ggbPick: TriangleCandidate
}

export type Vec3 = [number, number, number]

export const PERMS: Pairing[]
export const GGB_PERMS: Pairing[]
export const LETTERS: string[]
export function pairing(perm: Pairing): string

export function analyse(abc: Point[], def: Point[], opts?: { mirror?: boolean; scaling?: boolean }): TriangleAnalysis
export function placed(an: { K: Point; u: Point[] }, map: SimilarityMap): Point[]
export function ggbPlaced(an: TriangleAnalysis): Point[]
export function withHysteresis<C extends Candidate>(an: { all: C[]; best: C }, currentId: number, margin: number): C

export function shapePoint(tri: Point[]): Vec3
export function shapeImages(def: Point[], mirror: boolean): Vec3[]
export function dot3(a: Vec3, b: Vec3): number
/** The best pairing with mirror images allowed, by the order of the side lengths: no candidates tried. */
export function trianglePairing(abc: Point[], def: Point[]): { perm: Pairing; mirror: boolean }
/** Signed solid angle enclosed by a closed loop on the unit sphere, counterclockwise from outside positive. */
export function solidAngle(loop: Vec3[]): number
export function softFit(an: TriangleAnalysis, sigma: number): { pts: Point[]; weights: number[] }
