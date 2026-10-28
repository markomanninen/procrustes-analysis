#!/usr/bin/env node
// Command line front to js/procrustes.js. Reads point sets as JSON, prints JSON.
import { readFileSync, realpathSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { parseArgs } from 'node:util'
import { fit, distance, add, sub, mul, conj, arg, abs } from '../js/procrustes.js'

const USAGE = `Usage:
  procrustes fit <target> <source> [--pairings identity|cyclic|all|best] [--mirror] [--scaling]
  procrustes distance <a> <b> [--mirror]

<target>, <source>, <a>, <b> are JSON files holding [[x, y], …] or {"points": [[x, y], …]},
inline JSON such as '[[0,0],[1,0],[0,1]]', or - for standard input.

fit prints the best fit: rms and squared error, the pairing (pairing[i] is the source point
that goes with target point i), whether it mirrors, the turn, the scale, the transform
x ↦ a·x + b (or a·x̄ + b when mirrored) on points as complex numbers, and the fitted points
in target order. distance prints the full Procrustes distance, 0 … 1.

Exit status: 0 done, 1 bad input, 2 bad usage.`

class UsageError extends Error {}

function readPoints(spec, name, stdin) {
  let text
  if (spec === '-') text = stdin()
  else if (/^\s*[[{]/.test(spec)) text = spec
  else {
    try { text = readFileSync(spec, 'utf8') } catch (e) { throw new Error(`${name}: cannot read ${spec}: ${e.code ?? e.message}`) }
  }
  let data
  try { data = JSON.parse(text) } catch (e) { throw new Error(`${name}: not JSON: ${e.message}`) }
  const pts = Array.isArray(data) ? data : data?.points
  const ok = Array.isArray(pts) && pts.every((p) => Array.isArray(p) && p.length === 2 && p.every(Number.isFinite))
  if (!ok) throw new Error(`${name}: expected [[x, y], …] of finite numbers`)
  return pts
}

const round = (x) => (Math.abs(x) < 1e-12 ? 0 : Number(x.toPrecision(12)))
const roundAll = (v) => {
  if (Array.isArray(v)) return v.map(roundAll)
  if (typeof v === 'number') return round(v)
  if (v && typeof v === 'object') return Object.fromEntries(Object.entries(v).map(([k, x]) => [k, roundAll(x)]))
  return v
}

export function run(argv, { stdin = () => readFileSync(0, 'utf8') } = {}) {
  let parsed
  try {
    parsed = parseArgs({
      args: argv,
      allowPositionals: true,
      options: {
        pairings: { type: 'string', default: 'identity' },
        mirror: { type: 'boolean', default: false },
        scaling: { type: 'boolean', default: false },
        help: { type: 'boolean', short: 'h', default: false },
      },
    })
  } catch (e) {
    return { code: 2, err: `${e.message}\n\n${USAGE}` }
  }
  const { values: o, positionals: [cmd, ...files] } = parsed
  if (o.help || !cmd) return { code: o.help ? 0 : 2, out: o.help ? USAGE : undefined, err: o.help ? undefined : USAGE }
  try {
    if (!['fit', 'distance'].includes(cmd)) throw new UsageError(`unknown command '${cmd}'`)
    if (files.length !== 2) throw new UsageError(`${cmd} takes two point sets, got ${files.length}`)
    if (!['identity', 'cyclic', 'all', 'best'].includes(o.pairings)) throw new UsageError(`--pairings must be identity, cyclic, all or best, not '${o.pairings}'`)
    if (cmd === 'distance' && (o.scaling || o.pairings !== 'identity')) throw new UsageError('distance takes only --mirror')
    if (files.filter((f) => f === '-').length > 1) throw new UsageError('only one point set can come from standard input')
    const [a, b] = files.map((f, i) => readPoints(f, ['first', 'second'][i], stdin))
    if (cmd === 'distance') return { code: 0, out: JSON.stringify({ distance: round(distance(a, b, { mirror: o.mirror })) }) }
    const f = fit(a, b, { pairings: o.pairings, mirror: o.mirror, scaling: o.scaling })
    const { map, perm, mirror, rms, E, beta } = f.best
    // x ↦ K + a·(x − J) is x ↦ a·x + (K − a·J), and the same with conjugates when mirrored
    const shift = sub(f.K, mul(map.a, mirror ? conj(f.J) : f.J))
    return {
      code: 0,
      out: JSON.stringify(roundAll({
        rms, E, pairing: perm, mirror,
        turn_degrees: (arg(map.a) * 180) / Math.PI,
        scale: beta,
        transform: { a: map.a, b: shift, mirror },
        points: f.points,
      }), null, 2),
    }
  } catch (e) {
    return { code: e instanceof UsageError ? 2 : 1, err: e instanceof UsageError ? `${e.message}\n\n${USAGE}` : e.message }
  }
}

const isMain = () => {
  try { return realpathSync(process.argv[1]) === fileURLToPath(import.meta.url) } catch { return false }
}

if (isMain()) {
  const { code, out, err } = run(process.argv.slice(2))
  if (out) process.stdout.write(out + '\n')
  if (err) process.stderr.write(err + '\n')
  process.exitCode = code
}
