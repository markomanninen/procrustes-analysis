/** A point, which is also a complex number: [re, im]. */
export type Complex = [number, number]
export type Point = Complex
/** pairing[i] is the index of the source point paired with target point i. */
export type Pairing = number[]

/** The similarity w ↦ a·w, or w ↦ a·w̄ when mirror is set, on centred points. */
export interface SimilarityMap { a: Complex; mirror: boolean }

export interface Candidate {
  id: number
  perm: Pairing
  mirror: boolean
  /** Σ z̄_i u_perm(i); the best turn is −arg s. */
  s: Complex
  phi: number
  /** The scale; 1 without scaling. */
  beta: number
  /** Sum of squared residuals. */
  E: number
  rms: number
  /** GeoGebra's angle atan(Im s / Re s). */
  thetaAtan: number
  rmsAtan: number
  map: SimilarityMap
}

export interface FitOptions {
  pairings?: 'identity' | 'cyclic' | 'all' | Pairing[]
  mirror?: boolean
  scaling?: boolean
}

export interface Fit {
  /** Centroid of the target. */
  K: Point
  /** Centroid of the source. */
  J: Point
  z: Point[]
  u: Point[]
  all: Candidate[]
  best: Candidate
  /** The fitted source points, in target order. */
  points: Point[]
}

export interface TriangleAnalysis {
  K: Point
  J: Point
  z: Point[]
  u: Point[]
  all: Candidate[]
  best: Candidate
  ggb: Candidate[]
  ggbPick: Candidate
}

export function add(a: Complex, b: Complex): Complex
export function sub(a: Complex, b: Complex): Complex
export function mul(a: Complex, b: Complex): Complex
export function div(a: Complex, b: Complex): Complex
export function conj(a: Complex): Complex
export function scale(a: Complex, k: number): Complex
export function abs(a: Complex): number
export function arg(a: Complex): number
export function expi(t: number): Complex
export function wrap(t: number): number
export function centroid(pts: Point[]): Point
export function centre(pts: Point[]): Point[]

export const PERMS: Pairing[]
export const GGB_PERMS: Pairing[]
export const LETTERS: string[]
export function pairing(perm: Pairing): string

export function applyMap(map: SimilarityMap, w: Complex): Complex
export function candidate(z: Point[], u: Point[], perm: Pairing, opts?: { mirror?: boolean; scaling?: boolean }): Candidate

export function permutations(n: number): Pairing[]
export function cyclicPairings(n: number, opts?: { reverse?: boolean }): Pairing[]
export function fit(target: Point[], source: Point[], opts?: FitOptions): Fit
export function distance(a: Point[], b: Point[], opts?: { mirror?: boolean; pairing?: Pairing }): number

export function analyse(abc: Point[], def: Point[], opts?: { mirror?: boolean; scaling?: boolean }): TriangleAnalysis
export function placed(an: { K: Point; u: Point[] }, map: SimilarityMap): Point[]
export function ggbPlaced(an: TriangleAnalysis): Point[]
export function relativeMap(from: SimilarityMap, to: SimilarityMap): SimilarityMap
export function mapAt(map: SimilarityMap, t: number, w: Complex): Complex
export function withHysteresis(an: { all: Candidate[]; best: Candidate }, currentId: number, margin: number): Candidate

export type Vec3 = [number, number, number]
export function shapePoint(tri: Point[]): Vec3
export function shapeImages(def: Point[], mirror: boolean): Vec3[]
export function dot3(a: Vec3, b: Vec3): number
export function softFit(an: TriangleAnalysis, sigma: number): { pts: Point[]; weights: number[] }
