import { test, expect } from '@playwright/test'

// Every test fails on an uncaught error or console error in the page.
let problems
test.beforeEach(async ({ page }) => {
  problems = []
  page.on('pageerror', (e) => problems.push(e.message))
  page.on('console', (m) => m.type() === 'error' && problems.push(m.text()))
  await page.goto('/')
  await expect(page.locator('#rms')).not.toBeEmpty()
})
test.afterEach(() => expect(problems).toEqual([]))

const rms = async (page) => Number(await page.locator('#rms').textContent())
const rows = (page) => page.locator('#candidates tbody tr')

test('the three canvases are drawn', async ({ page }) => {
  for (const id of ['board', 'history', 'sphere']) {
    const inked = await page.locator(`#${id}`).evaluate((c) => {
      const { data } = c.getContext('2d').getImageData(0, 0, c.width, c.height)
      let n = 0
      for (let i = 3; i < data.length; i += 4) if (data[i] > 0) n++
      return n
    })
    expect(inked, id).toBeGreaterThan(500)
  }
})

test('exact copies fit with zero error, and the table lists every pairing', async ({ page }) => {
  for (const preset of ['halfturn', 'relabel']) {
    await page.selectOption('#preset', preset)
    await expect.poll(() => rms(page)).toBe(0)
  }
  await expect(rows(page)).toHaveCount(12)
  await page.uncheck('#mirror')
  await expect(rows(page)).toHaveCount(6)
})

test('the page shows the error the library computes, with and without mirror images', async ({ page }) => {
  await page.selectOption('#preset', 'mirror')
  // the preset is a mirror image of ABC with a little noise added to each vertex
  const expected = await page.evaluate(async () => {
    const Q = await import('/js/procrustes.js')
    const BASE = [[-14.9845, 12.5254], [-11.0133, 3.7193], [-8.158, 4.4858]]
    const D = BASE.map(Q.conj).map((p) => Q.add(Q.mul(p, Q.expi(0.9)), [26, 12])).map((p, i) => Q.add(p, [[0.4, -0.2], [-0.3, 0.3], [0.2, 0.1]][i]))
    return [true, false].map((mirror) => Q.analyse(BASE, D, { mirror }).best.rms.toFixed(3))
  })
  await expect(page.locator('#rms')).toHaveText(expected[0])
  await expect(page.locator('#turn')).toContainText('mirror about')
  await page.uncheck('#mirror')
  await expect(page.locator('#rms')).toHaveText(expected[1])
  await expect(page.locator('#turn')).not.toContainText('mirror')
})

test('moving a vertex from the keyboard changes the fit', async ({ page }) => {
  const before = await page.locator('#rms').textContent()
  await page.locator('#board').focus()
  await page.keyboard.press('4')
  for (let i = 0; i < 15; i++) await page.keyboard.press('ArrowRight')
  await expect(page.locator('#rms')).not.toHaveText(before)
})

test('the soft fit shows its weights', async ({ page }) => {
  await expect(page.locator('#weights')).toHaveText('–')
  await page.locator('#sigma').fill('1')
  await expect(page.locator('#weights')).toContainText('%')
})
