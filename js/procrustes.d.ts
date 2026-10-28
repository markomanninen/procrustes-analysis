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
  map: SimilarityMap
}

export interface FitOptions {
  /** 'all' lists every ordering (n ≤ 8); 'best' finds the best for any n. */
  pairings?: 'identity' | 'cyclic' | 'all' | 'best' | Pairing[]
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

export function applyMap(map: SimilarityMap, w: Complex): Complex
export function residual(z: Point[], u: Point[], perm: Pairing, map: SimilarityMap): number
export function candidate(z: Point[], u: Point[], perm: Pairing, opts?: { mirror?: boolean; scaling?: boolean }): Candidate

export function permutations(n: number): Pairing[]
export function cyclicPairings(n: number, opts?: { reverse?: boolean }): Pairing[]
/** The pairings at the vertices of the permutation polygon of centred z and u: the only ones that can be best. */
export function polygonPairings(z: Point[], u: Point[]): Pairing[]
export function fit(target: Point[], source: Point[], opts?: FitOptions): Fit
export function distance(a: Point[], b: Point[], opts?: { mirror?: boolean; pairing?: Pairing }): number

export function relativeMap(from: SimilarityMap, to: SimilarityMap): SimilarityMap
export function mapAt(map: SimilarityMap, t: number, w: Complex): Complex
