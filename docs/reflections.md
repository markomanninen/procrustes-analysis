# Why fitted triangles jump, and how to animate reflections

Animating a Procrustes fit of triangles, the fitted triangle sometimes leaps across the
screen, most of all near reflections. There are three separate causes. Only one of them is
a property of the problem itself.

## 1. atan instead of atan2

For a pairing, let s = Σ z̄_i·u_π(i), with z the centred vertices of ABC and u those of DEF.
Turning DEF by φ leaves the squared error

    E(φ) = Σ|z|² + Σ|u|² − 2·Re(e^{iφ} s)

which is smallest at φ = −arg s. The GeoGebra construction uses θ = atan(y/x) of s. That
gives the same value for s and −s, so when Re s < 0 it turns DEF to the *largest* error
instead of the smallest. As s crosses the imaginary axis while a vertex is dragged, atan
jumps from +90° to −90°, and the triangle leaps by half a turn.

**Fix.** Use atan2(y, x). This cause disappears entirely.

## 2. The best pairing changes

There are 6 pairings of vertices, 12 with mirror images. The fit shown is the pairing with
the least error, which is a minimum over 12 smooth functions of the vertex positions. Such a
minimum is continuous, but the *argument* of the minimum is not. Where two pairings tie,
the best fitted triangle is one triangle on one side of the tie and a different one on the
other. Near-equilateral triangles have near ties everywhere, since three of the turns are
about 120° apart and fit almost equally well.

That jump is the correct answer changing. Nothing in the formula can make it continuous.
What can be controlled is how it is shown.

**Fix.**
- **Hysteresis.** Keep the current pairing until another is better by a margin, so that a
  vertex resting on a tie does not make the display flicker.
- **A path between the two fits.** When the pairing does change, move from the old fitted
  triangle to the new one along the map between them.

That map is again a similarity. With old fit y = a_o·w (or a_o·w̄) and new fit a_n·w (or
a_n·w̄), it is

    y ↦ (a_n / a_o)·y         when both are turns, or both are mirror fits
    y ↦ (a_n / ā_o)·ȳ         when one is a turn and the other a mirror fit

The first is turned along the shorter arc. The second is a reflection, handled as in §3.

## 3. A mirror image cannot be turned into

Turns in the plane have determinant +1 and reflections −1, and no continuous path of plane
maps joins them. So an animation that tries to interpolate a turn into a mirror fit must
pass through a degenerate map or jump. In the GeoGebra construction a reflection was
approached by reversing the order of vertices. That changes which vertex goes where. The
triangle's shape is unchanged, and so is its handedness. A genuine mirror image is fitted
only by conjugating u.

The way out is to leave the plane. Every mirror fit z ↦ a·z̄ is a single reflection about the
line through the centroid at angle arg(a)/2, with no extra turn needed. In space, that
reflection is the end of a half turn about the same line. Writing a point relative to the
line, q = w·e^{−iα} with α = arg(a)/2, the path

    w(t) = |a|^t · e^{iα} · (Re q + i·Im q·cos πt),   0 ≤ t ≤ 1

is continuous. It flattens the triangle onto the line at t = ½ and arrives at the mirror image
at t = 1. It is what a sheet of paper does when it is turned over. On screen it reads as a
flip, not a jump, and it never needs a turn mixed in.

## Summary

| cause | inherent? | remedy |
|---|---|---|
| atan(y/x) | no | atan2 |
| tie between pairings | yes | hysteresis and a path between fits |
| turn versus reflection | yes, within the plane | a half turn in space about the mirror line |
