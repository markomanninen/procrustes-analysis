import * as Q from './procrustes.js'

const $ = (id) => document.getElementById(id)
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim()
const deg = (r) => (r * 180) / Math.PI
const NAMES = ['A', 'B', 'C', 'D', 'E', 'F']

const BASE = [[-14.9845, 12.5254], [-11.0133, 3.7193], [-8.158, 4.4858]]
const moved = (pts, a, t) => pts.map((p) => Q.add(Q.mul(p, a), t))
const PRESETS = {
  geogebra: () => [...BASE, [17.6394, 3.7687], [14.2589, 10.6926], [22.1603, 5.4386]],
  mirror: () => [...BASE, ...moved(BASE.map(Q.conj), Q.expi(0.9), [26, 12]).map((p, i) => Q.add(p, [[0.4, -0.2], [-0.3, 0.3], [0.2, 0.1]][i]))],
  halfturn: () => [...BASE, ...moved(BASE, [-1, 0], [8, 15])],
  relabel: () => {
    const t = moved(BASE, Q.expi((50 * Math.PI) / 180), [24, 12])
    return [...BASE, t[1], t[2], t[0]]
  },
  tie: () => {
    const eq = [0, 1, 2].map((k) => Q.add([-9, 7], Q.scale(Q.expi(Math.PI / 2 + (k * 2 * Math.PI) / 3), 5)))
    eq[0] = Q.add(eq[0], [0.6, 0])
    return [...eq, ...moved(eq, Q.expi(1.05), [14, -2])]
  },
}

const state = {
  pts: PRESETS.geogebra(),
  view: null,
  chosen: null,
  shown: null,
  transition: null,
  t: 1,
  playing: null,
  history: [],
  trail: [],
  drag: -1,
}

function opts() {
  return {
    mirror: $('mirror').checked,
    scaling: $('scaling').checked,
    ggb: $('show-ggb').checked,
    smooth: $('smooth').checked,
    margin: +$('margin').value / 100,
  }
}

// ---- view --------------------------------------------------------------------------

function fitView() {
  const c = $('board')
  const xs = state.pts.map((p) => p[0])
  const ys = state.pts.map((p) => p[1])
  const pad = 4
  const [x0, x1, y0, y1] = [Math.min(...xs) - pad, Math.max(...xs) + pad, Math.min(...ys) - pad, Math.max(...ys) + pad]
  const k = Math.min(c.clientWidth / (x1 - x0), c.clientHeight / (y1 - y0))
  state.view = { k, cx: (x0 + x1) / 2, cy: (y0 + y1) / 2 }
}

const toPx = ([x, y]) => {
  const c = $('board')
  const { k, cx, cy } = state.view
  return [c.clientWidth / 2 + (x - cx) * k, c.clientHeight / 2 - (y - cy) * k]
}
const toWorld = ([px, py]) => {
  const c = $('board')
  const { k, cx, cy } = state.view
  return [cx + (px - c.clientWidth / 2) / k, cy - (py - c.clientHeight / 2) / k]
}

function prepare(canvas) {
  const dpr = window.devicePixelRatio || 1
  const w = canvas.clientWidth
  const h = canvas.clientHeight
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr)
    canvas.height = Math.round(h * dpr)
  }
  const ctx = canvas.getContext('2d')
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
  ctx.clearRect(0, 0, w, h)
  return { ctx, w, h }
}

// ---- what is shown -----------------------------------------------------------------

// The fitted triangle on screen. During a switch between pairings it travels along the
// relative map from the previous fit; while the fit is being played it travels from DEF.
function shownPoints(an, now) {
  const target = state.chosen.map
  if (state.t < 1) return an.u.map((w) => Q.add(Q.add(an.J, Q.scale(Q.sub(an.K, an.J), state.t)), Q.mapAt(target, state.t, w)))
  const tr = state.transition
  if (tr) {
    const s = Math.min(1, (now - tr.start) / tr.ms)
    const e = s < 0.5 ? 2 * s * s : 1 - (-2 * s + 2) ** 2 / 2
    if (s >= 1) state.transition = null
    else return tr.from.map((y) => Q.add(an.K, Q.mapAt(tr.rel, e, Q.sub(y, an.K))))
  }
  return Q.placed(an, target)
}

function choose(an, o) {
  const prev = state.chosen
  const next = prev ? Q.withHysteresis(an, prev.id, o.margin) : an.best
  const fresh = an.all.find((c) => c.id === next.id)
  if (prev && fresh.id !== prev.id) {
    state.switches = (state.switches || 0) + 1
    if (o.smooth && state.t >= 1) {
      state.transition = { from: Q.placed(an, prev.map), rel: Q.relativeMap(prev.map, fresh.map), start: performance.now(), ms: 650 }
    }
  }
  state.chosen = fresh
  return prev && fresh.id !== prev.id
}

// ---- drawing -----------------------------------------------------------------------

function tri(ctx, pts, { stroke, fill, dash, width = 2 }) {
  ctx.save()
  ctx.setLineDash(dash || [])
  ctx.beginPath()
  pts.map(toPx).forEach((p, i) => (i ? ctx.lineTo(...p) : ctx.moveTo(...p)))
  ctx.closePath()
  if (fill) {
    ctx.fillStyle = fill
    ctx.fill()
  }
  ctx.strokeStyle = stroke
  ctx.lineWidth = width
  ctx.stroke()
  ctx.restore()
}

function dot(ctx, p, color, label, r = 5, labelColor) {
  const [x, y] = toPx(p)
  ctx.fillStyle = color
  ctx.beginPath()
  ctx.arc(x, y, r, 0, 2 * Math.PI)
  ctx.fill()
  if (label) {
    ctx.fillStyle = labelColor || css('--ink')
    ctx.font = '600 13px "Atkinson Hyperlegible", sans-serif'
    ctx.fillText(label, x + 8, y - 8)
  }
}

function grid(ctx, w, h) {
  const [x0, y1] = toWorld([0, 0])
  const [x1, y0] = toWorld([w, h])
  const step = state.view.k > 25 ? 1 : state.view.k > 8 ? 2 : 5
  ctx.lineWidth = 1
  for (let x = Math.ceil(x0 / step) * step; x <= x1; x += step) {
    ctx.strokeStyle = x === 0 ? css('--rule') : css('--grid')
    const px = Math.round(toPx([x, 0])[0]) + 0.5
    ctx.beginPath()
    ctx.moveTo(px, 0)
    ctx.lineTo(px, h)
    ctx.stroke()
  }
  for (let y = Math.ceil(y0 / step) * step; y <= y1; y += step) {
    ctx.strokeStyle = y === 0 ? css('--rule') : css('--grid')
    const py = Math.round(toPx([0, y])[1]) + 0.5
    ctx.beginPath()
    ctx.moveTo(0, py)
    ctx.lineTo(w, py)
    ctx.stroke()
  }
}

function drawBoard(an, o, now) {
  const { ctx, w, h } = prepare($('board'))
  grid(ctx, w, h)
  const abc = state.pts.slice(0, 3)
  const def = state.pts.slice(3)
  tri(ctx, abc, { stroke: css('--abc'), fill: css('--abc-fill') })
  tri(ctx, def, { stroke: css('--def'), fill: css('--def-fill') })
  if (o.ggb) tri(ctx, Q.ggbPlaced(an), { stroke: css('--ggb'), dash: [3, 4], width: 1.6 })

  const fit = shownPoints(an, now)
  const c = state.chosen
  if (c.mirror && state.t >= 1 && !state.transition) {
    // the line the mirror fit reflects about, through K at angle arg(a)/2
    const alpha = Q.arg(c.map.a) / 2
    const d = Q.scale(Q.expi(alpha), 30)
    ctx.save()
    ctx.strokeStyle = css('--axis')
    ctx.setLineDash([8, 5])
    ctx.lineWidth = 1.3
    ctx.beginPath()
    ctx.moveTo(...toPx(Q.sub(an.K, d)))
    ctx.lineTo(...toPx(Q.add(an.K, d)))
    ctx.stroke()
    ctx.restore()
  }
  tri(ctx, fit, { stroke: css('--def'), dash: [7, 4], width: 2.4 })
  if (state.t >= 1 && !state.transition) {
    ctx.save()
    ctx.strokeStyle = css('--resid')
    ctx.lineWidth = 1.5
    c.perm.forEach((j, i) => {
      ctx.beginPath()
      ctx.moveTo(...toPx(abc[i]))
      ctx.lineTo(...toPx(fit[j]))
      ctx.stroke()
    })
    ctx.restore()
  }
  fit.forEach((p, j) => dot(ctx, p, css('--def'), `${Q.LETTERS[j]}′`, 3.5, css('--muted')))
  dot(ctx, an.K, css('--ink'), 'K', 3.5)
  dot(ctx, an.J, css('--ink'), 'J', 3.5)
  state.pts.forEach((p, i) => dot(ctx, p, i < 3 ? css('--abc') : css('--def'), NAMES[i], state.drag === i ? 8 : 6))
}

function drawHistory() {
  const { ctx, w, h } = prepare($('history'))
  const H = state.history
  const box = { x: 40, y: 10, w: w - 50, h: h - 30 }
  ctx.strokeStyle = css('--rule')
  ctx.strokeRect(box.x + 0.5, box.y + 0.5, box.w, box.h)
  ctx.fillStyle = css('--muted')
  ctx.font = '11px "JetBrains Mono", monospace'
  if (H.length < 2) {
    ctx.fillText('Drag a vertex to record how every pairing’s error changes.', box.x + 10, box.y + box.h / 2)
    return
  }
  const top = Math.max(...H.flatMap((r) => r.rms)) * 1.05 || 1
  const X = (i) => box.x + (i / (H.length - 1)) * box.w
  const Y = (v) => box.y + box.h - (v / top) * box.h
  ctx.textAlign = 'right'
  for (const v of [0, top / 2, top]) ctx.fillText(v.toFixed(1), box.x - 5, Y(v) + 4)
  ctx.textAlign = 'left'
  const n = H[0].rms.length
  for (let k = 0; k < n; k++) {
    ctx.strokeStyle = k < 6 ? css('--abc-soft') : css('--axis-soft')
    ctx.lineWidth = 1
    ctx.beginPath()
    H.forEach((r, i) => (i ? ctx.lineTo(X(i), Y(r.rms[k])) : ctx.moveTo(X(i), Y(r.rms[k]))))
    ctx.stroke()
  }
  ctx.strokeStyle = css('--def')
  ctx.lineWidth = 2.4
  ctx.beginPath()
  H.forEach((r, i) => (i ? ctx.lineTo(X(i), Y(r.chosen)) : ctx.moveTo(X(i), Y(r.chosen))))
  ctx.stroke()
  ctx.fillStyle = css('--resid')
  H.forEach((r, i) => r.switched && ctx.fillRect(X(i) - 1, box.y, 2, box.h))
}

// ---- the shape sphere ----------------------------------------------------------------

const rgbCache = {}
function rgb(color) {
  if (!rgbCache[color]) {
    const c = document.createElement('canvas')
    c.width = c.height = 1
    const x = c.getContext('2d', { willReadFrequently: true })
    x.fillStyle = color
    x.fillRect(0, 0, 1, 1)
    rgbCache[color] = Array.from(x.getImageData(0, 0, 1, 1).data)
  }
  return rgbCache[color]
}
const mix = (a, b, t) => [0, 1, 2].map((i) => Math.round(a[i] * (1 - t) + b[i] * t))
const imageLabel = (c) => c.perm.map((j) => Q.LETTERS[j]).join('') + (c.mirror ? '*' : '')

// Both halves of the sphere seen from outside: the upper half from above, the lower from
// below, so the equator is the rim of each disc. Every pixel takes the colour of the image
// of DEF nearest to it; where two colours meet, two pairings tie.
function drawSphere(an, o) {
  const { ctx, w, h } = prepare($('sphere'))
  const dpr = window.devicePixelRatio || 1
  const R = Math.max(40, Math.min(w / 4 - 22, h / 2 - 30))
  const nA = Q.shapePoint(state.pts.slice(0, 3))
  const imgs = Q.shapeImages(state.pts.slice(3), o.mirror)
  const panel = rgb(css('--panel'))
  const fills = imgs.map((_, k) => mix(panel, rgb(k < 6 ? css('--abc') : css('--axis')), k % 2 ? 0.1 : 0.2))
  fills[an.best.id] = mix(panel, rgb(css('--def')), 0.32)
  const edge = rgb(css('--ink'))
  const size = Math.ceil(2 * R * dpr)
  const sides = [1, -1]
  const project = (p, side) => ({ x: side * p[0], y: p[1], front: p[2] * side >= 0 })
  sides.forEach((side, d) => {
    const cx = (w / 4) * (1 + 2 * d)
    const cy = h / 2 - 8
    const off = document.createElement('canvas')
    off.width = off.height = size
    const octx = off.getContext('2d')
    const img = octx.createImageData(size, size)
    const cell = new Int8Array(size * size).fill(-1)
    for (let py = 0; py < size; py++) {
      for (let px = 0; px < size; px++) {
        const x = ((px + 0.5) / size) * 2 - 1
        const y = 1 - ((py + 0.5) / size) * 2
        const r2 = x * x + y * y
        if (r2 > 1) continue
        const P3 = [side * x, y, side * Math.sqrt(1 - r2)]
        let best = 0
        let bd = -2
        imgs.forEach((m, k) => {
          const v = Q.dot3(P3, m)
          if (v > bd) [best, bd] = [k, v]
        })
        cell[py * size + px] = best
      }
    }
    for (let i = 0; i < size * size; i++) {
      const k = cell[i]
      if (k < 0) continue
      const right = i % size < size - 1 ? cell[i + 1] : k
      const down = i + size < size * size ? cell[i + size] : k
      const border = (right >= 0 && right !== k) || (down >= 0 && down !== k)
      const c = border ? edge : fills[k]
      img.data.set([c[0], c[1], c[2], border ? 170 : 255], i * 4)
    }
    octx.putImageData(img, 0, 0)
    ctx.drawImage(off, cx - R, cy - R, 2 * R, 2 * R)
    ctx.strokeStyle = css('--rule')
    ctx.lineWidth = 1.5
    ctx.beginPath()
    ctx.arc(cx, cy, R, 0, 2 * Math.PI)
    ctx.stroke()
    const toPx = (q) => [cx + q.x * R, cy - q.y * R]

    // where ABC has been while you dragged, then where it is
    ctx.strokeStyle = css('--abc')
    ctx.lineWidth = 1.5
    ctx.beginPath()
    let pen = false
    for (const n of state.trail) {
      const q = project(n, side)
      if (!q.front) {
        pen = false
        continue
      }
      const [a, b] = toPx(q)
      if (pen) ctx.lineTo(a, b)
      else ctx.moveTo(a, b)
      pen = true
    }
    ctx.stroke()
    ctx.font = '10px "JetBrains Mono", monospace'
    imgs.forEach((m, k) => {
      const q = project(m, side)
      if (!q.front) return
      const [a, b] = toPx(q)
      const shown = k === state.chosen.id
      ctx.fillStyle = shown ? css('--def') : k < 6 ? css('--abc') : css('--axis')
      ctx.beginPath()
      ctx.arc(a, b, shown ? 5 : 3, 0, 2 * Math.PI)
      ctx.fill()
      ctx.fillStyle = css('--ink')
      ctx.fillText(imageLabel(an.all[k]), a + 6, b - 5)
    })
    const q = project(nA, side)
    if (q.front) {
      const [a, b] = toPx(q)
      ctx.strokeStyle = css('--abc')
      ctx.lineWidth = 2.5
      ctx.beginPath()
      ctx.arc(a, b, 7, 0, 2 * Math.PI)
      ctx.stroke()
      ctx.fillStyle = css('--ink')
      ctx.font = '600 12px "Atkinson Hyperlegible", sans-serif'
      ctx.fillText('ABC', a + 10, b + 14)
    }
    ctx.fillStyle = css('--muted')
    ctx.font = '12px "Atkinson Hyperlegible", sans-serif'
    ctx.textAlign = 'center'
    ctx.fillText(side > 0 ? 'upper half, seen from above' : 'lower half, seen from below', cx, cy + R + 20)
    ctx.textAlign = 'left'
  })
}

function tables(an) {
  const best = an.best
  const rows = [...an.all].sort((a, b) => a.E - b.E).map((c) => {
    const ggb = !c.mirror && c.perm.join() === an.ggbPick.perm.join()
    const marks = [c.id === best.id ? '<span class="chip best">least squares</span>' : '', ggb ? '<span class="chip ggb">GeoGebra pick</span>' : '', c.id === state.chosen.id ? '<span class="chip shown">shown</span>' : ''].join('')
    return `<tr class="${c.id === state.chosen.id ? 'on' : ''}">
      <td class="l">${Q.pairing(c.perm)}</td><td>${c.mirror ? 'mirror' : ''}</td>
      <td>${c.s[0].toFixed(2)} ${c.s[1] < 0 ? '−' : '+'} ${Math.abs(c.s[1]).toFixed(2)}i</td>
      <td>${Q.abs(c.s).toFixed(2)}</td><td>${deg(-c.thetaAtan).toFixed(2)}°</td><td>${deg(c.phi).toFixed(2)}°</td>
      <td>${c.rms.toFixed(3)}</td><td class="l">${marks}</td></tr>`
  })
  $('candidates').innerHTML = `<thead><tr><th class="l">Pairing</th><th></th><th>s = Σ z̄ᵢ u</th><th>|s|</th><th>−atan(y/x)</th><th>−arg s</th><th>rms</th><th></th></tr></thead><tbody>${rows.join('')}</tbody>`
  const g = an.ggb.map((c, k) => {
    const pick = c === an.ggbPick
    const typo = k === 5
    return `<tr class="${pick ? 'on' : ''}"><td class="l">s${'₁₂₃₄₅₆'[k]}</td><td class="l">${c.perm.map((j, i) => `v${'₁₂₃'[i]}u${'₁₂₃'[j]}`).join(' + ')}${typo ? ' <span class="chip warn">u₁ twice</span>' : ''}</td>
      <td>${deg(c.thetaAtan).toFixed(2)}°</td><td>${typo ? '–' : c.rmsAtan.toFixed(3)}</td><td>${typo ? '–' : c.rms.toFixed(3)}</td><td class="l">${pick ? '<span class="chip ggb">smallest |θ|</span>' : ''}</td></tr>`
  })
  $('ggb-table').innerHTML = `<thead><tr><th class="l">Sum</th><th class="l">As written</th><th>θ = atan(y/x)</th><th>rms, turned by −θ</th><th>rms, turned by −arg s</th><th></th></tr></thead><tbody>${g.join('')}</tbody>`
  const sorted = [...an.all].sort((a, b) => a.E - b.E)
  const gap = sorted.length > 1 ? (sorted[1].E - sorted[0].E) / Math.max(sorted[1].E, 1e-12) : 1
  $('gap').textContent = `${(100 * gap).toFixed(1)} %`
  $('gap').className = gap < 0.1 ? 'warn' : ''
  $('rms').textContent = state.chosen.rms.toFixed(3)
  $('turn').textContent = state.chosen.mirror ? `mirror about ${deg(Q.arg(state.chosen.map.a) / 2).toFixed(1)}°` : `${deg(state.chosen.phi).toFixed(2)}°`
  $('scale').textContent = state.chosen.rho.toFixed(3)
  $('switches').textContent = state.switches || 0
}

// ---- loop --------------------------------------------------------------------------

let frame = 0
function render(recordHistory = false) {
  const o = opts()
  const an = Q.analyse(state.pts.slice(0, 3), state.pts.slice(3), o)
  const switched = choose(an, o)
  if (recordHistory) {
    state.history.push({ rms: [...an.all].sort((a, b) => a.id - b.id).map((c) => c.rms), chosen: state.chosen.rms, switched })
    if (state.history.length > 400) state.history.shift()
    state.trail.push(Q.shapePoint(state.pts.slice(0, 3)))
    if (state.trail.length > 400) state.trail.shift()
  }
  const now = performance.now()
  if (state.playing) {
    state.t = Math.min(1, (now - state.playing) / 2200)
    $('t').value = state.t
    if (state.t >= 1) state.playing = null
  }
  drawBoard(an, o, now)
  drawHistory()
  drawSphere(an, o)
  tables(an)
  $('out-t').textContent = state.t.toFixed(2)
  $('out-margin').textContent = `${$('margin').value} %`
  if ((state.transition || state.playing) && !frame) frame = requestAnimationFrame(() => ((frame = 0), render()))
}

function resetAll(name) {
  state.pts = PRESETS[name]()
  state.chosen = null
  state.transition = null
  state.history = []
  state.trail = []
  state.switches = 0
  fitView()
  render()
}

// ---- input -------------------------------------------------------------------------

const board = $('board')
board.addEventListener('pointerdown', (e) => {
  const r = board.getBoundingClientRect()
  const m = [e.clientX - r.left, e.clientY - r.top]
  let best = -1
  let bd = 18
  state.pts.forEach((p, i) => {
    const d = Math.hypot(...Q.sub(toPx(p), m))
    if (d < bd) [best, bd] = [i, d]
  })
  if (best < 0) return
  state.drag = best
  board.setPointerCapture(e.pointerId)
  e.preventDefault()
})
board.addEventListener('pointermove', (e) => {
  if (state.drag < 0) return
  const r = board.getBoundingClientRect()
  state.pts[state.drag] = toWorld([e.clientX - r.left, e.clientY - r.top])
  render(true)
})
const end = () => {
  state.drag = -1
  render()
}
board.addEventListener('pointerup', end)
board.addEventListener('pointercancel', end)

// keyboard: pick a vertex with 1–6, move it with the arrow keys
let keyVertex = 0
board.addEventListener('keydown', (e) => {
  if (/^[1-6]$/.test(e.key)) keyVertex = +e.key - 1
  const step = e.shiftKey ? 1 : 0.2
  const d = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, step], ArrowDown: [0, -step] }[e.key]
  if (!d) return
  e.preventDefault()
  state.pts[keyVertex] = Q.add(state.pts[keyVertex], d)
  render(true)
})

$('preset').addEventListener('change', () => resetAll($('preset').value))
for (const id of ['mirror', 'scaling', 'show-ggb', 'smooth', 'margin']) $(id).addEventListener('input', () => render())
$('t').addEventListener('input', () => {
  state.playing = null
  state.t = +$('t').value
  render()
})
$('play').addEventListener('click', () => {
  state.transition = null
  state.t = 0
  state.playing = performance.now()
  render()
})
$('fit-view').addEventListener('click', () => {
  fitView()
  render()
})
window.addEventListener('resize', () => {
  fitView()
  render()
})
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => render())
new MutationObserver(() => render()).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] })

fitView()
render()
