// Where the hard fit jumps, and what the soft fit does there.
// Moves ABC or DEF through an isosceles and a collinear triangle and reports the change of
// best candidate and the largest jump of a fitted vertex. Then puts ABC exactly on the
// collinear line and prints the two largest soft-fit weights for three values of sigma.
// Run from the repo root: node docs/research/scripts/jump_test.mjs
import * as P from '../../../js/procrustes.js'
import * as T from '../../../js/triangle.js'
const Q = { ...P, ...T } // complex helpers from procrustes.js; analyse, softFit from triangle.js
const partners = (an, c) => [0,1,2].map(i => Q.add(an.K, Q.applyMap(c.map, an.u[c.perm[i]])))
const hard = (abc, def) => { const an = Q.analyse(abc, def, {mirror:true}); return {an, pts: partners(an, an.best), id: an.best.id} }
const maxd = (P, R) => Math.max(...P.map((p,i)=>Q.abs(Q.sub(p,R[i]))))
const area = (t) => ((t[1][0]-t[0][0])*(t[2][1]-t[0][1]) - (t[2][0]-t[0][0])*(t[1][1]-t[0][1]))/2
function sweep(name, mk, fixedIsABC) {
  const eps = 1e-4
  for (const side of [-eps, eps]) {}
  const [a1, d1] = mk(-eps), [a2, d2] = mk(eps)
  const h1 = hard(a1, d1), h2 = hard(a2, d2)
  console.log(`${name}: best id ${h1.id} -> ${h2.id}; max vertex-partner jump = ${maxd(h1.pts, h2.pts).toFixed(4)}`)
}
const generic = [[0,0],[3,0.4],[1.1,2.2]]          // scalene, not isosceles
// ABC passes through an isosceles triangle (apex x = 1.5 on the axis of A=(0,0), B=(3,0))
sweep('ABC crosses isosceles', t => [[[0,0],[3,0],[1.5+t,1.7]], generic])
// ABC passes through a collinear triangle
sweep('ABC crosses collinear', t => [[[0,0],[3,0],[1.1,t]], generic])
// DEF passes through isosceles / collinear, ABC generic
const abc = [[0,0],[2.6,-0.3],[0.7,1.9]]
sweep('DEF crosses isosceles', t => [abc, [[0,0],[3,0],[1.5+t,1.7]]])
sweep('DEF crosses collinear', t => [abc, [[0,0],[3,0],[1.1,t]]])
// Soft fit exactly at a collinear ABC: equal weights for a proper/mirror pair?
for (const sig of [0.3, 1, 3]) {
  const an = Q.analyse([[0,0],[3,0],[1.1,0]], generic, {mirror:true})
  const s = Q.softFit(an, sig)
  const top = s.weights.map((w,k)=>[w,k]).sort((a,b)=>b[0]-a[0]).slice(0,2)
  console.log(`collinear ABC, sigma=${sig}: top weights ${top.map(([w,k])=>k+':'+w.toFixed(3)).join(' ')}; area(soft fit)=${area(s.pts).toFixed(4)} vs area(DEF)=${area(generic).toFixed(4)}`)
}
