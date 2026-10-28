// Carried round a loop: the triangle on the left changes shape as C circles; the one on the
// right is the same triangle best matched step by step, so no step turns it, and yet after
// the loop it has turned by −Ω/2, half the solid angle the loop encloses on the shape sphere.
import { transport, add, mul, conj, arg, scale } from './procrustes.js'
import { shapePoint, solidAngle } from './triangle.js'

const $ = (id) => document.getElementById(id)
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim()
const deg = (r) => (r * 180) / Math.PI
const N = 720

const frame = (t, R) => [[-1, 0], [1, 0], [0.2 + R * Math.cos(t), 0.9 + R * Math.sin(t)]]
const turnBetween = (a, b) => arg(a.reduce((s, p, i) => add(s, mul(conj(p), b[i])), [0, 0]))

let loop = null
let at = N
let playing = null

function build() {
  const R = Number($('loop-size').value)
  const frames = Array.from({ length: N + 1 }, (_, i) => frame((2 * Math.PI * i) / N, R))
  const carried = transport(frames)
  const omega = solidAngle(frames.slice(0, N).map(shapePoint))
  loop = { R, frames, carried, omega, turn: turnBetween(carried[0], carried[N]) }
  $('out-loop-size').textContent = R.toFixed(2)
  $('loop-turn').textContent = `${deg(loop.turn).toFixed(2)}°`
  $('loop-omega').textContent = `${deg(-loop.omega / 2).toFixed(2)}°`
}

function draw() {
  const c = $('loop')
  const dpr = window.devicePixelRatio || 1
  const w = c.clientWidth
  const h = c.clientHeight
  c.width = Math.round(w * dpr)
  c.height = Math.round(h * dpr)
  const ctx = c.getContext('2d')
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
  ctx.clearRect(0, 0, w, h)
  const unit = Math.min(w / 4, h) / 3.4
  const poly = (pts, centre, colour, { dashed = false, fill = null, width = 2 } = {}) => {
    ctx.beginPath()
    pts.forEach(([x, y], i) => ctx[i ? 'lineTo' : 'moveTo'](centre[0] + x * unit, centre[1] - y * unit))
    ctx.closePath()
    if (fill) { ctx.fillStyle = fill; ctx.fill() }
    ctx.setLineDash(dashed ? [6, 5] : [])
    ctx.strokeStyle = colour
    ctx.lineWidth = width
    ctx.stroke()
    ctx.setLineDash([])
  }
  const left = [w / 4, h / 2 + unit * 0.4]
  const right = [(3 * w) / 4, h / 2]
  // the loop C runs round
  ctx.beginPath()
  ctx.arc(left[0] + 0.2 * unit, left[1] - 0.9 * unit, loop.R * unit, 0, 2 * Math.PI)
  ctx.strokeStyle = css('--muted')
  ctx.setLineDash([3, 5])
  ctx.lineWidth = 1
  ctx.stroke()
  ctx.setLineDash([])
  poly(loop.frames[at], left, css('--abc'), { fill: css('--abc-fill') })
  poly(loop.carried[0], right, css('--muted'), { dashed: true, width: 1.5 })
  poly(loop.carried[at], right, css('--def'), { fill: css('--def-fill') })
  // vertex A marks the turn
  const a0 = scale(loop.carried[0][0], unit)
  const a = scale(loop.carried[at][0], unit)
  for (const [p, colour] of [[a0, css('--muted')], [a, css('--def')]]) {
    ctx.beginPath()
    ctx.arc(right[0] + p[0], right[1] - p[1], 4, 0, 2 * Math.PI)
    ctx.fillStyle = colour
    ctx.fill()
  }
  ctx.fillStyle = css('--muted')
  ctx.font = `13px ${css('--body')}`
  ctx.fillText('shape: C goes round the dotted loop', 12, 20)
  ctx.fillText(`carried: turned ${deg(turnBetween(loop.carried[0], loop.carried[at])).toFixed(1)}° from the dashed start`, w / 2 + 12, 20)
}

function tick(now) {
  if (!playing) return
  at = Math.min(N, Math.round(((now - playing) / 4000) * N))
  draw()
  if (at < N) requestAnimationFrame(tick)
  else playing = null
}

$('loop-size').addEventListener('input', () => { build(); at = N; draw() })
$('loop-play').addEventListener('click', () => { at = 0; playing = performance.now(); requestAnimationFrame(tick) })
window.addEventListener('resize', draw)
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw)
build()
draw()
