// Scan for jumps of the hard fit: ABC's apex runs along a horizontal line through the
// isosceles position, with mirrors off and on, and every change of best candidate is listed
// with the size of the jump.
// Run from the repo root: node docs/research/scripts/jump_test2.mjs
import * as P from '../../../js/procrustes.js'
import * as T from '../../../js/triangle.js'
const Q = { ...P, ...T } // complex helpers from procrustes.js; analyse, softFit from triangle.js
const partners = (an, c) => [0,1,2].map(i => Q.add(an.K, Q.applyMap(c.map, an.u[c.perm[i]])))
const maxd = (P, R) => Math.max(...P.map((p,i)=>Q.abs(Q.sub(p,R[i]))))
const generic = [[0,0],[3,0.4],[1.1,2.2]]
// Proper-only fits (mirror off): scan ABC's apex along a horizontal line through the isosceles position, report every switch
for (const mirror of [false, true]) {
  let prev = null, out = []
  for (let k = 0; k <= 4000; k++) {
    const x = -1 + 5 * k / 4000, abc = [[0,0],[3,0],[x,1.7]]
    const an = Q.analyse(abc, generic, {mirror}); const pts = partners(an, an.best)
    if (prev && an.best.id !== prev.id) out.push(`x≈${x.toFixed(3)} id ${prev.id}->${an.best.id} jump ${maxd(prev.pts, pts).toFixed(3)}`)
    prev = {id: an.best.id, pts}
  }
  console.log(`mirror=${mirror}: apex x from -1 to 4 at height 1.7 (isosceles at x=1.5):\n  ` + out.join('\n  '))
}
