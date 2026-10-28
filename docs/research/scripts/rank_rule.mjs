// Side-length rule test (novelty_jump_locus.md, claim B). 20,000 random triangle pairs, with
// and without scaling.
// Run from the repo root: node docs/research/scripts/rank_rule.mjs
// Prediction: with mirrors on, the best of the 12 candidates pairs vertices by the rank of the
// opposite side length, and mirrors exactly when the two triangles have opposite orientation.
import { analyse } from '../../../js/triangle.js'

const rnd = () => [Math.random() * 10 - 5, Math.random() * 10 - 5]
const tri = () => [rnd(), rnd(), rnd()]
const opp = (t) => [0, 1, 2].map((i) => Math.hypot(t[(i + 1) % 3][0] - t[(i + 2) % 3][0], t[(i + 1) % 3][1] - t[(i + 2) % 3][1]))
const area = (t) => (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[2][0] - t[0][0]) * (t[1][1] - t[0][1])
const ranks = (v) => { const o = [0, 1, 2].sort((a, b) => v[a] - v[b]); const r = []; o.forEach((i, k) => (r[i] = k)); return r }

function rule(abc, def) {
  const ra = ranks(opp(abc)), rd = ranks(opp(def))
  const byRank = []; rd.forEach((k, j) => (byRank[k] = j))
  const perm = ra.map((k) => byRank[k]); return { perm, mirror: Math.sign(area(abc)) !== Math.sign(area(perm.map((j) => def[j]))) }
}

for (const scaling of [false, true]) {
  let agree = 0, agreeRms = 0, n = 20000, worst = 0
  let noMirrorAgree = 0
  for (let t = 0; t < n; t++) {
    const abc = tri(), def = tri()
    const an = analyse(abc, def, { mirror: true, scaling })
    const r = rule(abc, def)
    const same = an.best.mirror === r.mirror && an.best.perm.join() === r.perm.join()
    if (same) agree++
    const ruleCand = an.all.find((c) => c.mirror === r.mirror && c.perm.join() === r.perm.join())
    const gap = ruleCand.rms - an.best.rms
    if (gap < 1e-9) agreeRms++
    worst = Math.max(worst, gap)
    // Mirrors off: does ranking alone still pick the best proper candidate?
    const an0 = analyse(abc, def, { mirror: false, scaling })
    const r0 = rule(abc, def)
    if (an0.best.perm.join() === r0.perm.join()) noMirrorAgree++
  }
  console.log(`scaling=${scaling}: mirrors on, rule = best in ${agree}/${n} (same rms in ${agreeRms}/${n}, worst rms gap ${worst.toExponential(2)}); mirrors off, ranking alone = best in ${noMirrorAgree}/${n}`)
}
