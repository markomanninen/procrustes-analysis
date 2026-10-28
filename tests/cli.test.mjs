import { test } from 'node:test'
import assert from 'node:assert/strict'
import { spawnSync } from 'node:child_process'
import { writeFileSync, mkdtempSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { fileURLToPath } from 'node:url'
import { run } from '../bin/procrustes.js'

const BIN = fileURLToPath(new URL('../bin/procrustes.js', import.meta.url))
const json = (r) => JSON.parse(r.out)

test('fit prints the transform, and the transform maps the source onto the fitted points', () => {
  const target = [[0, 0], [2, 0], [0, 1], [1, 3]]
  const source = [[5, 5], [5, 7], [4, 5], [2, 6]].map(([x, y]) => [x * 1.5, y * 1.5])
  const r = run(['fit', JSON.stringify(target), JSON.stringify(source), '--scaling'])
  assert.equal(r.code, 0)
  const f = json(r)
  assert.equal(f.rms, 0)
  assert.equal(f.scale, Number((1 / 1.5).toPrecision(12)))
  const { a, b } = f.transform
  source.forEach(([x, y], i) => {
    const p = [a[0] * x - a[1] * y + b[0], a[0] * y + a[1] * x + b[1]]
    assert.ok(Math.hypot(p[0] - target[i][0], p[1] - target[i][1]) < 1e-9)
  })
})

test('fit with --mirror and --pairings cyclic finds a reversed mirror image', () => {
  const target = [[0, 0], [4, 0], [5, 2], [1, 3]]
  const source = target.map(([x, y]) => [x + 10, -y]).reverse()
  const f = json(run(['fit', JSON.stringify(target), JSON.stringify(source), '--mirror', '--pairings', 'cyclic']))
  assert.equal(f.rms, 0)
  assert.equal(f.mirror, true)
  assert.deepEqual(f.pairing, [3, 2, 1, 0])
  f.points.forEach((p, i) => assert.deepEqual(p, target[i]))
})

test('points come from files, {"points": …} objects and standard input', () => {
  const dir = mkdtempSync(join(tmpdir(), 'procrustes-'))
  const file = join(dir, 'a.json')
  writeFileSync(file, JSON.stringify({ points: [[0, 0], [1, 0], [0, 1]] }))
  const r = run(['distance', file, '-', '--mirror'], { stdin: () => '[[0,0],[1,0],[0,-1]]' })
  assert.deepEqual(json(r), { distance: 0 })
})

test('fit --pairings best matches points given in any order', () => {
  const target = [[0, 0], [4, 0], [5, 2], [1, 3], [2, 5]]
  const source = [target[3], target[0], target[4], target[1], target[2]].map(([x, y]) => [-y + 9, x])
  const f = json(run(['fit', JSON.stringify(target), JSON.stringify(source), '--pairings', 'best']))
  assert.equal(f.rms, 0)
  assert.deepEqual(f.pairing, [1, 3, 4, 0, 2])
})

test('bad input exits 1, bad usage exits 2, and both say why', () => {
  const tri = '[[0,0],[1,0],[0,1]]'
  const cases = [
    [['fit', tri, '[[0,0],[1,0]]'], 1, /same number of points/],
    [['fit', tri, '[[0,0],[1,"x"],[0,1]]'], 1, /finite numbers/],
    [['fit', tri, '/no/such/file.json'], 1, /cannot read/],
    [['fit', tri, '{oops'], 1, /not JSON/],
    [['fit', tri], 2, /two point sets/],
    [['fit', tri, tri, '--pairings', 'some'], 2, /--pairings/],
    [['distance', tri, tri, '--scaling'], 2, /only --mirror/],
    [['fit', '-', '-'], 2, /standard input/],
    [['spin', tri, tri], 2, /unknown command/],
    [['fit', tri, tri, '--frobnicate'], 2, /Unknown option/],
    [[], 2, /Usage/],
  ]
  for (const [args, code, msg] of cases) {
    const r = run(args)
    assert.equal(r.code, code, args.join(' '))
    assert.match(r.err, msg, args.join(' '))
  }
  assert.equal(run(['--help']).code, 0)
})

test('the installed command runs as a process and sets the exit status', () => {
  const ok = spawnSync(process.execPath, [BIN, 'distance', '-', '[[0,0],[1,0],[0,1]]'], { input: '[[0,0],[2,0],[0,2]]', encoding: 'utf8' })
  assert.equal(ok.status, 0)
  assert.deepEqual(JSON.parse(ok.stdout), { distance: 0 })
  const bad = spawnSync(process.execPath, [BIN, 'fit'], { encoding: 'utf8' })
  assert.equal(bad.status, 2)
  assert.match(bad.stderr, /Usage/)
})
