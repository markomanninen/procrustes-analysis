import { test } from 'node:test'
import assert from 'node:assert/strict'
import { execFileSync } from 'node:child_process'
import { readFileSync } from 'node:fs'

const pkg = JSON.parse(readFileSync(new URL('../package.json', import.meta.url)))

test('the package ships the library, its types, the command and nothing else', () => {
  const npm = process.platform === 'win32' ? 'npm.cmd' : 'npm'
  const [packed] = JSON.parse(execFileSync(npm, ['pack', '--dry-run', '--json', '--ignore-scripts'], { encoding: 'utf8' }))
  const files = packed.files.map((f) => f.path).sort()
  assert.deepEqual(files, ['LICENSE', 'README.md', 'bin/procrustes.js', 'js/procrustes.d.ts', 'js/procrustes.js', 'package.json'])
})

test('every path named in package.json is one that ships', async () => {
  const named = [pkg.exports['.'].default, pkg.exports['.'].types, pkg.types, ...Object.values(pkg.bin)]
  for (const p of named) assert.ok(pkg.files.includes(p.replace(/^\.\//, '')), p)
  const lib = await import(new URL(`../${pkg.exports['.'].default}`, import.meta.url))
  const declared = readFileSync(new URL(`../${pkg.types}`, import.meta.url), 'utf8')
  for (const name of Object.keys(lib)) assert.match(declared, new RegExp(`export (function|const) ${name}\\b`), `${name} has no type`)
})
