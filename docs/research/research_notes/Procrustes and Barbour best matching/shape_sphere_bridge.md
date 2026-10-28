# The shape sphere as the bridge between Kendall/Procrustes shape theory and N-body relational mechanics

Scope: (a) the geometry of the triangle shape sphere; (b) how celestial mechanics uses it; (c) Edward Anderson's identification of Kendall's shape spaces with Barbour's reduced configuration spaces; (d) the gauge-theory ("falling cat") account of separating rotation from shape change, with the holonomy formula that separates chained infinitesimal best matching from one-shot Procrustes fitting. Literature is covered up to October 2026. Notation: q = (q1, q2, q3) in C^3 are planar vertex positions; m = (m1, m2, m3) are masses or weights (Kendall: all equal to 1); I = Σ m_i |q_i − q_cm|^2 is the moment of inertia, the squared size. Unless noted, a "Procrustes fit" means minimising a (mass-)weighted squared distance over translations and rotations, with scaling or reflection added only where stated.

---

## 1. Shape sphere geometry: Kendall's S²(1/2), the Hopf map, where special triangles lie, and the S3 × Z2 action

### Takeaway
After translations are removed, a planar three-point configuration is a vector Z = (Z1, Z2) in C² (normalised, mass-weighted Jacobi coordinates). Rotations act by Z ↦ e^{iθ}Z. The Hopf-type quadratic map Z ↦ w ∈ R³ realises the quotient by rotations. Size is |w| = I/2. Shape is the direction of w, a point on a round sphere which, with the quotient (Procrustes) metric, has radius 1/2. This holds for **any** positive masses. Collinear triangles fill the equator (w3 ∝ signed area). The poles are the triangles of maximal |area|/I, and these are equilateral **only when the masses are equal**, which is always the case in Kendall's unit-weight setting. When labels and reflections are ignored, the 12-element group S3 × Z2 tiles the sphere with 12 congruent spherical triangles (angles π/3, π/2, π/2). One tile is Kendall's "spherical blackboard".

### Cited Findings
- **Quotient construction (Montgomery).** Translations are removed by the complex-linear projection to normalised Jacobi coordinates Z1 = μ1(q1 − q2) and Z2 = μ2(q3 − (m1q1 + m2q2)/(m1 + m2)), defined as Z_i = ⟨q, E_i⟩/‖E_i‖ in the mass inner product. These coordinates diagonalise the mass metric, so that I = |Z1|² + |Z2|² (his eq. 24). A rotation acts as (Z1, Z2) ↦ (e^{iθ}Z1, e^{iθ}Z2) — [Montgomery, AMM preprint arXiv:1402.0841, §8.2–8.3](https://arxiv.org/abs/1402.0841). **Typo warning:** the preprint prints "μ1² = ‖E1‖² = 1/m1 + 1/m2". Its own definition Z_i = ⟨q, E_i⟩/‖E_i‖ instead gives the reduced-mass factors μ1 = √(m1m2/(m1 + m2)) and μ2 = √(m3(m1 + m2)/(m1 + m2 + m3)). My calculation; the numerical check in §4 uses the latter.
- **Hopf map.** The map π_rot(Z1, Z2) = (½(|Z1|² − |Z2|²), Re(Z1 Z̄2), Im(Z1 Z̄2)) = (w1, w2, w3) "realizes R³ as the quotient space of C² by the rotation group S¹". Its restriction to a sphere "is the famous Hopf map from the three-sphere to the two-sphere" — [Montgomery, arXiv:1402.0841, Prop. 3 and eq. (33)](https://arxiv.org/abs/1402.0841)
- **Properties of the map (Montgomery's Theorem 1):**
  - (A) Two triangles are oriented-congruent if and only if they have the same w.
  - (C) Triple collision maps to the origin.
  - (D) w3 equals the signed area times a mass-dependent constant.
  - (E) The collinear triangles form the plane w3 = 0.
  - (F) Reflection acts as (w1, w2, w3) ↦ (w1, w2, −w3).
  - (G) w1² + w2² + w3² = (½ I)².

  — [Montgomery, arXiv:1402.0841, Thm 1](https://arxiv.org/abs/1402.0841)
- **The metric is literally a Procrustes metric.** Montgomery defines the shape-space metric as the quotient metric d([q1], [q2]) = inf over g1, g2 ∈ G of ‖g1q1 − g2q2‖, where G is the group of rigid motions and the norm comes from the mass metric (his eq. 23). He proves that this distance agrees with the Riemannian distance of ds² = (dw1² + dw2² + dw3²)/(2√(w1² + w2² + w3²)). Its geodesics are the projections of horizontal lines in C³, and on every plane through the origin the metric is dr² + ¼r²dθ², a cone of angle factor c = 1/2 — [Montgomery, arXiv:1402.0841, eqs. (23), (51), Thm 3](https://arxiv.org/abs/1402.0841)
- **Radius 1/2.** The space of congruence classes is "a cone over the space S+ of similarity classes … isometric to the closed upper hemisphere of radius one-half. The equator represents collinear configurations", with dilation parameter λ = √I. For oriented triangles, S "is isometric to the two-sphere of radius one-half". The planar quotient is "isometric as Riemannian submersions to the standard Hopf fibration S¹ → S³(1) → S²(1/2)" — [Montgomery 1996, Nonlinearity 9:1341, arXiv:dg-ga/9510005, Thms 2–3 and §3.1](https://arxiv.org/abs/dg-ga/9510005)
- **Height formula.** On the shape sphere of radius 1/2 the normalised height above the equator is z1 = 4·√(m1m2m3/(m1 + m2 + m3))·Δ/I, where Δ is the oriented area — [Montgomery 1996, arXiv:dg-ga/9510005, §1.1 and Thm 3](https://arxiv.org/abs/dg-ga/9510005). Montgomery credits this formula, and the idea of passing to oriented triangles, to Hsiang's 1994 Berkeley preprint PAM-620 — [same, §1.3 and §3.2](https://arxiv.org/abs/dg-ga/9510005)
- **Kendall's shape spaces.** Kendall showed that the shape spaces of N points are spheres S^{N−2} in 1-d and complex projective spaces CP^{N−2} in 2-d — [Anderson, arXiv:1307.1923 (Stud. Hist. Phil. Mod. Phys. 51:1–8)](https://arxiv.org/abs/1307.1923). For triangles this gives CP¹ ≅ S². Anderson's triangleland papers use the same "sphere of radius 1/2", rescaled to a unit sphere when convenient — [Anderson, Triangleland I, arXiv:0809.1168](https://arxiv.org/abs/0809.1168). Kendall's founding paper is "Shape manifolds, Procrustean metrics, and complex projective spaces", Bull. London Math. Soc. 16 (1984) 81–121 — [Crossref record, DOI 10.1112/blms/16.2.81](https://doi.org/10.1112/blms/16.2.81). Its full text was not read this session.
- **Location of special triangles (equal masses).** The three isosceles loci are "three great circles on the shape sphere, each circle passing through the two Lagrange points which are represented as the North and South Poles". Each circle ISOSC_i meets the collinear plane in two points: an Euler point EUL_i, where mass i is at the midpoint of the other two, and a binary-collision point B_i — [Montgomery, arXiv:1402.0841, §14](https://arxiv.org/abs/1402.0841). The three central-configuration families of Euler (collinear) and the two of Lagrange (equilateral) are the five critical points marked on the sphere — [same, §5 and Fig. 3](https://arxiv.org/abs/1402.0841)
- **Unequal masses need care.** Montgomery 2002 deliberately uses "equal-mass Hopf coordinates" for arbitrary masses, so that the poles z = ±1 are the Lagrange (equilateral) points. In those coordinates the true mass metric is only conformal to the round one: ds²_m = λ(φ, θ)²(dφ² + cos²φ dθ²) with λ given by his eq. (5.7) — [Montgomery, "Infinitely many syzygies", ARMA 164 (2002) 311–340, §5, author PDF](https://people.ucsc.edu/~rmont/papers/InfinitelyManySyzygies.pdf). Moeckel and Montgomery note that a Möbius transformation of CP¹ placing the three binary collisions at the cube roots of unity "automatically" moves the equilateral points to the poles. These are their "equilateral coordinates" — [Moeckel & Montgomery, Pacific J. Math. 262 (2013) 129–189, arXiv:1202.0972, Example 4](https://arxiv.org/abs/1202.0972)
- **Reflection.** Reflection conjugates (Z1, Z2) and flips w3. Unoriented congruence classes are therefore the closed half-space w3 ≥ 0, with the collinear triangles as its boundary. Montgomery uses oriented congruence because "the collinear triangles are smooth points of shape space" only in that case — [Montgomery, arXiv:1402.0841, §7, §8.4](https://arxiv.org/abs/1402.0841)
- **Relabelling acts by isometries (equal masses).** Swapping masses 1 and 3 (σ2) "is a linear isometric involution … The effect of σ2 on the shape sphere is a half-twist about the line connecting EUL2 and B2". Composing σ2 with the reflection about the triangle's symmetry axis gives "reflection about the plane ISOSC2" — [Montgomery, arXiv:1402.0841, §14](https://arxiv.org/abs/1402.0841)
- **Anderson's four quotients of the triangle shape sphere:**
  - the full sphere;
  - mirror images identified: a hemisphere;
  - labels indistinguishable: "a lune of width π/3";
  - both identifications: "the Leibniz space demilune of width π/3: Kendall's spherical blackboard".

  Elsewhere he calls this smallest cell "the 1/6th hemisphere bounded by arcs of isoscelesness and collinearity" and "the Leibs(3,2) 1/12-sphere" — [Anderson, "The Smallest Shape Spaces III", arXiv:1711.10115](https://arxiv.org/abs/1711.10115); [Anderson, "Alice in Triangleland", arXiv:1711.11492](https://arxiv.org/abs/1711.11492). **Inconsistency in the source:** the introduction of arXiv:1711.10115 instead says "a lune covering 1/3 of the sphere for the indistinguishable shapes" and "an isosceles spherical triangle covering 1/6 of the shape sphere" for the blackboard. arXiv:1307.1923 Fig. 2 says "1/3 of a hemisphere" and draws half of it. See the Inferences below for the group-theoretic count.
- **Uniform measure on the sphere.** On the shape sphere, right triangles form "three kissing cap-circles", and "a shape-theoretic answer to … the probability that a triangle is obtuse very readily follows: 3/4" — [Anderson, arXiv:1711.10115 abstract](https://arxiv.org/abs/1711.10115)
- **Kendall's distributional work.** The 1985 paper on triangle shapes is Kendall, "Exact distributions for shapes of random triangles in convex sets", Adv. Appl. Probab. 17 (1985) 308–329 — [Crossref record, DOI 10.2307/1427143](https://doi.org/10.2307/1427143). It was not read this session. Anderson credits the "spherical blackboard" picture to Kendall's shape-statistics work — [arXiv:1711.10115](https://arxiv.org/abs/1711.10115)
- **Hsiang & Straume** published "Kinematic geometry of triangles and the study of the three-body problem" in Lobachevskii J. Math. 25 (2007) 9–130 (MR2312347, Zbl 1118.53011) — [Math-Net.ru record](https://www.mathnet.ru/eng/ljm2). The monograph does "Riemannian geometric reduction … to the congruence moduli space and the shape space". It derives "Euler equations for mass triangle motions … generalizing the classical Euler equations for rigid body motions" — [arXiv:math-ph/0608060](https://arxiv.org/abs/math-ph/0608060)
- **Iwai 1987** brought the fibre-bundle picture into molecular and three-body dynamics:
  - "A geometric setting for classical molecular dynamics", Ann. Inst. H. Poincaré Phys. Théor. 47 (1987) 199–219, with keywords including "Eckart frame" and "Coriolis force" — [EuDML](https://eudml.org/doc/76377), [Numdam](https://numdam.org/item/AIHPA_1987__47_2_199_0/)
  - "A gauge theory for the quantum planar three-body problem", J. Math. Phys. 28 (1987) 964–974 — [DOI 10.1063/1.527588](https://doi.org/10.1063/1.527588)
  - Anderson attributes the half-space (mirror-identified) monopole to Iwai — [arXiv:1711.10115](https://arxiv.org/abs/1711.10115)

### Inferences
- **The S3 × Z2 action (equal masses / Kendall weights).** This combines Montgomery's half-twist statement with the 120° symmetry of the equal-mass sphere.
  - Cyclic relabellings rotate the sphere by ±120° about the polar (equilateral) axis.
  - Transpositions are half-turns about equatorial axes through B_k and EUL_k.
  - So S3 acts as the rotation group D3 (order 6).
  - The reflection z ↦ −z adds the mirror in the equator.
  - The full group of order 12 is generated by mirrors in the equator and in the three isosceles great circles.
  - The four great circles cut the sphere into 12 congruent spherical triangles. Each has vertices at a pole (equilateral, angle π/3), an Euler point (angle π/2) and a binary-collision point (angle π/2). Its area is π/3 on the unit sphere, i.e. 1/12 of the sphere.

  So the label-only quotient has area 1/6 of the sphere (a π/3-wide lune). The label-and-mirror quotient (Kendall's spherical blackboard / Anderson's Leibniz space) has area 1/12. Anderson's "π/3 lune" and "1/12-sphere" statements fit this count; his "1/3" and "1/6 of the sphere" phrasings do not.
- **Kendall's Riemannian shape distance between triangles.** The cone factor is c = 1/2 (Montgomery Thm 3), so the geodesic (Procrustes) distance between two shapes on S²(1/2) is ρ = ½·∠(ŵ, ŵ′), the half-angle between their w-directions. Antipodal shapes, for example a triangle and its mirror image when these sit at opposite poles, are at the maximal distance ρ = π/2.
- **Where special triangles sit for unequal masses.** Apply Montgomery's height formula to an equilateral triangle of side s. Then Δ = (√3/4)s² and I = s²(m1m2 + m2m3 + m3m1)/M, so z1 = √(3·m1m2m3·M)/(m1m2 + m2m3 + m3m1). This is ≤ 1, with equality if and only if m1 = m2 = m3. So in **mass-metric** coordinates the Lagrange points sit off the poles for unequal masses. I checked the formula numerically this session: for masses (1, 2, 5), z1 = 4√(m1m2m3/M)·Δ/I reproduced w3/|w| from the Hopf map to about 10⁻¹⁶. Equator = collinear holds for all masses. This is why Montgomery 2002 and Moeckel–Montgomery 2013 switch to equal-mass or "equilateral" coordinates, which are only conformal to the mass metric. For Kendall (unit weights) the poles are always the equilateral triangles.
- **Discontinuities in the fit.** Collinear and isosceles triangles lie on the mirror circles of S3 × Z2. They are fixed points of a reflection or relabelling, so two of the 12 candidate relabelled or reflected copies coincide there. This is where minimising over the finite group can become non-unique and "jump". This matches the behaviour noted in this repository's recent commit messages ("the fit jumps where either triangle is flat or isosceles").

### Gaps
- I did not read Kendall 1984 or Kendall 1985 in full this session. The claim "Σ³₂ = S²(1/2), poles equilateral, equator collinear, spherical blackboard = 1/12 of the sphere" rests on Montgomery, Anderson and group theory, not on Kendall's own text.
- I did not read Lemaître's original papers. Moeckel–Montgomery credit him with the degree-4 "octahedral covering map of the shape sphere, branched over the binary collision points", and Murnaghan (1936) with a symmetric reduction — [arXiv:1202.0972](https://arxiv.org/abs/1202.0972). Those claims are secondary.

---

## 2. Celestial mechanics: how the three-body problem uses the shape sphere

### Takeaway
In celestial mechanics the shape sphere is the base of the rotation (and scaling) bundle. Zero-angular-momentum motions are exactly the horizontal motions, and the reduced dynamics is a particle moving on shape space, attracted by the three binary-collision rays in the collinear plane. Several results come straight out of this picture:
- infinitely many syzygies (equator crossings);
- the figure-eight, built from 12 arcs across the 12 fundamental domains;
- hyperbolic "pants" for the 1/r² problem (the sphere minus the three collision points);
- symbolic dynamics by syzygy sequences;
- Montgomery's reconstruction formula for the overall rotation.

### Cited Findings
- **Montgomery's AMM exposition.** "The shape sphere is the moduli space of oriented similarity classes of planar triangles and lies inside shape space, a Euclidean 3-space parametrizing oriented congruence classes of triangles" — [Montgomery, Amer. Math. Monthly 122 (2015) 299–321; preprint arXiv:1402.0841](https://arxiv.org/abs/1402.0841); [DOI 10.4169/amer.math.monthly.122.04.299](https://doi.org/10.4169/amer.math.monthly.122.04.299). The arXiv text is a "Preliminary Version" dated 15 Nov 2011, so the published wording may differ.
- **Zero momentum = orthogonal to the orbit.** A curve is orthogonal to the oriented-congruence class at a point if and only if its linear and angular momentum vanish there (Lemma 1). The extremals of the shape-space action principle "are precisely those solutions to Newton's equations whose linear and angular momentum are zero" (Prop. 4) — [Montgomery, arXiv:1402.0841, §9](https://arxiv.org/abs/1402.0841)
- **Saari decomposition.** Velocities split into vertical (rigid motion) and horizontal (shape) parts: K = ½|P|²/M + ½J²/I + ½|ẇ|²/I, so K_shape = ½|ẇ|²/(2|w|). Montgomery explicitly calls this the vertical/horizontal splitting of bundle theory — [Montgomery, arXiv:1402.0841, Prop. 5, eqs. (40)–(45)](https://arxiv.org/abs/1402.0841)
- **Infinitely many syzygies.** "Any bounded zero-angular-momentum solution of the Newtonian three-body problem suffers infinitely many syzygies (collinearities) provided that it does not suffer a triple collision … The proof relies on the conformal geometry of the shape sphere." Syzygy sequences in the letters 1, 2, 3 encode free homotopy types of curves that are periodic modulo rotation — [Montgomery, ARMA 164 (2002) 311–340, author PDF](https://people.ucsc.edu/~rmont/papers/InfinitelyManySyzygies.pdf). In shape-space terms, w is "always attracted toward the collinear plane" by the three binary-collision rays, so it must cross w3 = 0 infinitely often — [Montgomery, arXiv:1402.0841, §13, Thm 4](https://arxiv.org/abs/1402.0841). A follow-up is Montgomery, "The zero angular momentum, three-body problem: all but one solution has syzygies", ETDS 27 (2007) 1933–1946 — [DOI 10.1017/S0143385707000338](https://doi.org/10.1017/S0143385707000338). Only the record was seen.
- **The figure-eight (Chenciner & Montgomery 2000).**
  - The orbit has zero angular momentum, and "visits in turns every 'Euler configuration'" — [arXiv:math/0011268; Ann. of Math. 152 (2000) 881–901](https://arxiv.org/abs/math/0011268)
  - Construction on the shape sphere: minimise the shape-space action over paths from EUL1 to the isosceles circle ISOSC2. The minimiser meets both orthogonally and is extended by half-twists and reflections, which gives "twelve solution arcs each congruent to the original arc" and "a smooth solution … periodic in shape space with period 12T". "With some extra work … the solution is actually periodic of period 12T in the original inertial space" — [Montgomery, arXiv:1402.0841, §14, Thm 5](https://arxiv.org/abs/1402.0841)
  - Moore found the curve numerically in 1993 — [same, §14 "History"](https://arxiv.org/abs/1402.0841); [Moore, PRL 70 (1993) 3675, DOI 10.1103/PhysRevLett.70.3675](https://doi.org/10.1103/PhysRevLett.70.3675)
- **Hyperbolic pants.** "Modulo symmetries, the dynamics of the bounded zero-angular-momentum solutions [of the 1/r² three-body problem] is equivalent to a geodesic flow on the thrice-punctured sphere, or 'pair of pants'. The sphere is the shape sphere. The punctures are the binary collisions. The metric … is the Jacobi-Maupertuis metric." For equal masses the curvature is negative except at the two Lagrange points. Corollaries are the uniqueness of the 1/r² figure-eight and a complete symbolic dynamics — [Montgomery, ETDS 25 (2005) 921–947, arXiv:math/0405014](https://arxiv.org/abs/math/0405014). Krishnaswami & Senapati extend the curvature analysis to C², R³ and S³ in "Hopf" coordinates — [J. Math. Phys. 57 (2016) 102901, arXiv:1606.05091](https://arxiv.org/abs/1606.05091)
- **Free homotopy classes.** Reduced free homotopy classes "have a simple description when projected onto the shape sphere. They are coded by syzygy sequences", and every such class is realised by a reduced periodic solution of the Newtonian problem with nonzero angular momentum — [Moeckel & Montgomery, arXiv:1412.2263](https://arxiv.org/abs/1412.2263)
- **Moeckel.** Montgomery notes that the shape sphere "plays a central role in Moeckel's work on the three-body problem", citing Moeckel, "Some qualitative features of the three-body problem", Contemp. Math. 81 (1988) 1–22 — [Montgomery 1996, arXiv:dg-ga/9510005, §1.1 and refs](https://arxiv.org/abs/dg-ga/9510005). Moeckel's lecture notes "Central Configurations" appear in Central Configurations, Periodic Orbits, and Hamiltonian Systems (Birkhäuser, 2015), pp. 105–167 — [DOI 10.1007/978-3-0348-0933-7_2](https://doi.org/10.1007/978-3-0348-0933-7_2). Only the record was seen.
- **Regularisation.** In the planar problem "the map from the configuration space to the shape sphere is a Hopf fibration, a nontrivial circle bundle". Lemaître's regularising map "behaves like the complex squaring map at each of the binary collision points on the shape sphere" — [Moeckel & Montgomery 2013, arXiv:1202.0972](https://arxiv.org/abs/1202.0972)
- **Reconstruction formula (Montgomery 1996).** "Given that the initial and final triangles of a three-body motion are similar, what is the rotation R which relates the two triangles?" Answer: J0·Δθ = ∫ω_{J0} dt (dynamic phase) + ∫∫_D Ω_{J0} (geometric phase). In the planar case the geometric integrand is "simply the area form J0 dz1 ∧ dθ1 on the shape sphere, normalized so that the total area of this sphere is J0·2π" — [Montgomery, Nonlinearity 9 (1996) 1341–1360, arXiv:dg-ga/9510005, eq. (2) and §1.2](https://arxiv.org/abs/dg-ga/9510005). Kuang (2021) extends the formula to the rotation angle of a single particle without requiring the initial and final configurations to be similar — [arXiv:2106.15106](https://arxiv.org/abs/2106.15106)
- **Recent physics use of "Kendall's shape sphere".** Garaud & Niemi (July 2026) build a shape-space theory of near-equilateral atomic trimers using "the Higgs oscillator on Kendall's shape sphere" — [arXiv:2607.16846](https://arxiv.org/abs/2607.16846)

### Inferences
- The celestial-mechanics shape sphere and Kendall's sphere are the same Riemannian object. Montgomery's eq. (23) is a mass-weighted, scale-free Procrustes distance, and his Thm 3 gives S²(1/2), as Kendall found for unit weights. The two literatures differ in emphasis: statistics uses the quotient metric and measures on it, while mechanics also uses the bundle structure (horizontal lifts, holonomy) that statistics usually discards.
- The 12-arc construction of the figure-eight is a dynamical use of the same S3 × Z2 fundamental domain that Kendall used as his "spherical blackboard".

### Gaps
- I did not obtain the published AMM version, so the section and equation numbers above refer to the arXiv preliminary version.
- I did not read the 2007 ETDS syzygy paper or Moeckel's 1988 and 2015 texts. Their content here is from records, titles and Montgomery's description.

---

## 3. Edward Anderson: Kendall's shape theory as the geometry of Barbour's relational particle mechanics

### Takeaway
Anderson's central claims are three:
- Kendall's shape spaces are the reduced configuration spaces of Barbour's scale-free relational particle mechanics (RPM).
- The cones over them are those of Barbour–Bertotti's scaled RPM.
- Building mechanics directly on Kendall's spaces (Jacobi–Synge construction) gives the same theory as reducing Barbour's best-matched actions. He calls this the "Direct = Best-Matched Theorem".

He works throughout in mass-weighted relative Jacobi coordinates. Mirror-image identification (Z2) and label indistinguishability (S3) are optional extra quotients, and he catalogues the spaces and monopoles that result. He is the only author found who mentions best matching and Procrustes in one breath, and he says they extremise *different* objects.

### Cited Findings
- **The core claim.** "I previously showed that Kendall's work on shape geometry is in fact also the geometrical description of Barbour's relational mechanics' reduced configuration spaces (alias shape spaces)." Highlights: "Kendall's shape geometry is Barbour's mechanics' reduced configuration space"; "Kendall's shape statistics is a classical timeless records theory" — [Anderson, Stud. Hist. Phil. Mod. Phys. 51 (2015) 1–8, arXiv:1307.1923](https://arxiv.org/abs/1307.1923); [DOI 10.1016/j.shpsb.2015.05.002](https://doi.org/10.1016/j.shpsb.2015.05.002)
- **Which theory matches which space.**
  - Kendall's "underlying notion of shape coincides with … that used in Barbour's pure-shape RPM".
  - "the cone … over Kendall's pure-shape geometry turns out to be the reduced configuration space for the scaled RPM".
  - The relational space R(N, d) := Q(N, d)/Eucl(d) "coincides with scaled RPM's configuration space … the cone over the above shape space".

  — [Anderson, arXiv:1307.1923, §1 and §3](https://arxiv.org/abs/1307.1923). The two Barbour theories are Barbour–Bertotti 1982 (scaled) and Barbour 2003 (pure shape) — [Anderson, arXiv:1111.1472 abstract](https://arxiv.org/abs/1111.1472); [Barbour & Bertotti, Proc. R. Soc. Lond. A 382 (1982) 295–306](https://doi.org/10.1098/rspa.1982.0102); [Barbour, CQG 20 (2003) 1543–1570, arXiv:gr-qc/0211021](https://arxiv.org/abs/gr-qc/0211021)
- **The "Direct = Best-Matched" theorem.** "Apply the Jacobi–Synge approach … to Kendall's shape spaces and to the cone over these. I previously showed that these procedures coincide …; I term this the Direct = Best-Matched Theorem." "Kendall's work – and its straightforward extension by coning – amounts to having already derived the detailed topological and metric structure of the reduced configuration spaces for Barbour's RPM's." He adds that "this coincidence of procedures was by no means guaranteed. For instance, Bookstein gave a distinct metric geometry on the space of shapes" — [Anderson, arXiv:1307.1923, §4](https://arxiv.org/abs/1307.1923). The original proof is in "Foundations of relational particle dynamics", §4 "Equivalence of my direct construction and the Barbour-type indirect construction": "Kendall's work is (a good start on) the configuration space study for Barbour's mechanics" — [Anderson, CQG 25 (2008) 025003, arXiv:0706.3934](https://arxiv.org/abs/0706.3934)
- **Best matching vs Procrustes.** Anderson defines best matching as "bringing two configurations into minimum incongruence with each other by application of G's group action". In the caption to Fig. 3 he writes: "Matching is a procedure that bears much conceptual similarity to the Procrustes procedure of Shape Geometry [10]; however these each involve extremizing a different object" — [Anderson, arXiv:1307.1923, §2 and Fig. 3 caption](https://arxiv.org/abs/1307.1923)
- **Mass weighting.** In triangleland: "I use mass-weighted Jacobi variables as the most succinct 'original variables of the problem'". Reduction proceeds by "Step 1) passing to mass-weighted Jacobi coordinates. Step 2) eliminating the rotational auxiliary …" — [Anderson, Triangleland I, CQG 26 (2009) 135020, arXiv:0809.1168](https://arxiv.org/abs/0809.1168). In 1307.1923 the configuration-space radius is ρ = √I, and the pure-shape coordinates are a relative angle Φ and Θ = 2 arctan(ρ2/ρ1) of the mass-weighted Jacobi magnitudes — [arXiv:1307.1923, Fig. 3](https://arxiv.org/abs/1307.1923). "The Smallest Shape Spaces III" rests on "3-body problem techniques including mass-weighted relative Jacobi coordinates" and "Hopf mathematics" — [arXiv:1711.10115](https://arxiv.org/abs/1711.10115)
- **Mirror images and labelling.** "Shapes in the sense of constellations of points with various automorphisms quotiented out: continuous translations, rotations and dilations, and also discrete mirror image identification and labelling indistinguishability" — [Anderson, "The Smallest Shape Spaces I", arXiv:1711.10054](https://arxiv.org/abs/1711.10054). For triangles he gives the full sphere, the hemisphere, the π/3-lune and the Leibniz space (Kendall's spherical blackboard), and attaches a monopole to each:
  - "one due to Dirac" (relational space R³);
  - "one to Iwai" (half-space);
  - two new ones.

  — [arXiv:1711.10115](https://arxiv.org/abs/1711.10115). In "Monopoles of Twelve Types in 3-Body Problems" there are "8 cases … in 2-d, and 4 distinct ones in 3-d; these reflect counts of non-equivalent subgroup actions of S3 × C2 and S3 respectively. The S3 acts on particle labels; the extra C2 corresponds to the purely 2-d option of whether or not to identify mirror images" — [arXiv:1802.03465](https://arxiv.org/abs/1802.03465)
- **Topological shape theory.** "Kendall's Similarity Shape Theory … was subsequently shown to reside within (Classical and Quantum) Mechanics' Shape-and-Scale Theory, in which the points are interpreted as particles and the carrier space plays the role of absolute space." The generalisations cover affine, projective, conformal and supersymmetric versions and other carrier spaces. "Rubber shapes" have shape spaces that "are just graphs" — [Anderson, arXiv:1803.11126](https://arxiv.org/abs/1803.11126)
- **Six new mechanics.** The paper extends "-1) Kendall's (metric) shape space with his shape statistics and Barbour's mechanics thereupon" and "0) Leibnizian relational space alias metric scale-and-shape space to which corresponds Barbour-Bertotti mechanics". The new cases are conformal (angle) shape mechanics, affine (area-ratio) shape mechanics, affine scale-and-shape (area) mechanics, and supersymmetric counterparts — [Anderson, IJMPD 25 (2016) 1650044, arXiv:1505.00488](https://arxiv.org/abs/1505.00488)
- **Quadrilaterals.** Triangleland's shape space is "S² with isometry group SO(3)". For the quadrilateral "the shape space now being CP² with isometry group SU(3)/Z3", solving the model requires "Geometry, Shape Statistics and Particle Physics" — [Anderson, Relational Quadrilateralland I, IJMPD 23 (2014) 1450014, arXiv:1202.4186](https://arxiv.org/abs/1202.4186). Part II (quantum) is by Anderson & Kneller, IJMPD 23 (2014) 1450052 — [Crossref DOI 10.1142/S0218271814500527](https://doi.org/10.1142/S0218271814500527); [arXiv:1303.5645](https://arxiv.org/abs/1303.5645)
- **Anderson's book.** He also wrote The Problem of Time (Springer, Fundamental Theories of Physics, 2017) — [DOI 10.1007/978-3-319-58848-3](https://doi.org/10.1007/978-3-319-58848-3). His review arXiv:1111.1472 names Kendall ("relational models only: pure-shape configuration spaces") among the principal sources of his Problem-of-Time approach — [arXiv:1111.1472 abstract](https://arxiv.org/abs/1111.1472)

### Inferences
- What Anderson claims is a reduced-space and Lagrangian-level equivalence (same quotient manifold, same metric, same reduced mechanics). It is not a claim that the finite best-matching step and Procrustes fitting minimise the same functional. His "different object" remark plausibly refers to two differences:
  - Barbour's best matching extremises a Jacobi-type action ∫√((E − V)·Σ m dx·dx), weighted by the potential, along a path.
  - Procrustes minimises a squared distance between two given configurations.

  In the potential-free (geodesic) case the two coincide infinitesimally; see §4.
- Mass weighting is harmless to the identification. In mass-weighted Jacobi coordinates the mass metric becomes the Euclidean metric on C², so the quotient geometry is Kendall's for any masses (Montgomery's S²(1/2) result in §1). What changes with masses is where named triangles sit on the sphere (Lagrange points off the poles, collision points not 120° apart), and S3 stops acting by isometries unless masses are equal.
- On reflections, Kendall's shape spaces quotient by SO(d), as do Montgomery's shape sphere and Anderson's main RPMs. Identifying mirror images (O(2)) gives the hemisphere with a collinear boundary, so "allow reflections" in a Procrustes fit is the Z2 quotient.

### Gaps
- I could not read the full text of Anderson's 2017 book. I cannot confirm whether it restates the Kendall–Barbour equivalence or uses the word "Procrustes".
- I found no explicit statement by Anderson that "mass weighting can be absorbed into Kendall's unit-weight theory". The inference above is mine, from the Jacobi-coordinate argument.
- "Anderson & Kneller" co-authored only Quadrilateralland II in the material checked. Quadrilateralland I is Anderson alone.

---

## 4. Gauge theory of deformable bodies: the zero-angular-momentum connection, best matching, and the holonomy formula

### Takeaway
The mechanical (Guichardet–Iwai–Shapere–Wilczek) connection declares velocities with zero angular momentum to be horizontal, i.e. mass-metric-orthogonal to the rotation orbits. This is exactly Barbour's infinitesimal best-matching condition. Chaining infinitesimal best matches along a path of shapes is therefore parallel transport (a horizontal lift).

For the planar three-body problem in the Hopf trivialisation, the connection is A = ½ z dθ and its curvature is half the solid-angle form. Hence:
- **Holonomy rule:** a closed loop on the shape sphere enclosing solid angle Ω (measured on the unit sphere) produces a net rotation Δθ = Ω/2 (mod 2π, sign set by orientation).
- **Equivalently:** Δθ equals twice the area enclosed on Kendall's S²(1/2).
- **With J ≠ 0:** the dynamic term ∫J/I dt is added.

A one-shot Procrustes fit from start to end ignores this rotation, so the chained and one-shot fits differ by exactly this holonomy.

### Cited Findings
- **Who found the connection.** "The one-form α_{J0} is the component of the 'natural mechanical connection' A … The connection A first appears explicitly in Guichardet. It can be argued that it was discovered by Smale … It was later used by Iwai and rediscovered by Shapere and Wilczek. I used it … in studying the Falling Cat Problem" — [Montgomery 1996, arXiv:dg-ga/9510005, §1.3](https://arxiv.org/abs/dg-ga/9510005). Guichardet's paper is "On rotation and vibration motions of molecules", Ann. Inst. H. Poincaré Phys. Théor. 40 (1984) 329–342 — [Numdam PDF](https://www.numdam.org/item/AIHPA_1984__40_3_329_0.pdf). Gomes & Gryb add that Guichardet "dealt with the kinematics of deformable bodies such as molecules and falling cats", and that "the same discovery was made independently by Shapere and Wilczek (1989; 1989b)" — [Gomes & Gryb, arXiv:2011.01693](https://arxiv.org/abs/2011.01693)
- **Formula and meaning.** The connection is A(q) = 𝕀(q)⁻¹ ∘ J(q), with 𝕀 the inertia tensor and J = Σ m_a q_a × dq_a. "The horizontal space is the space of infinitesimal deformations with zero angular momentum" — [Montgomery 1996, arXiv:dg-ga/9510005, §2](https://arxiv.org/abs/dg-ga/9510005)
- **The planar case is the Hopf connection.** Restricted to planar configurations, "the restriction of A … is Γê3 where Γ is the canonical connection for the Hopf fibration. One can choose a local section s … such that s*Γ = ½ z1 dθ1 … It follows that s*A = (½ z1 dθ1) ê3". A is basic under dilations, σλ*A = A, so the holonomy depends only on the loop projected to the shape sphere — [Montgomery 1996, arXiv:dg-ga/9510005, §3.2 eqs. (6)–(7)](https://arxiv.org/abs/dg-ga/9510005)
- **Reconstruction formula, planar case.** Δθ = ∫ (J0/I(t)) dt + ∫∫_D dz1 ∧ dθ1, with the area form "normalized so that the total area of this sphere is … 2π". In the spatial case the geometric term lives on a two-sphere bundle S_{J0} over the shape sphere, with Ω_{J0} = J0{½ d(z1z2) ∧ dθ1 + dz2 ∧ dθ2} — [Montgomery 1996, eqs. (2)–(4), §1.2](https://arxiv.org/abs/dg-ga/9510005); [DOI 10.1088/0951-7715/9/5/014](https://doi.org/10.1088/0951-7715/9/5/014)
- **Monopole language.** Anderson describes the curvature on the cone R³ (relational space) as "the Dirac monopole's realization in the R³ relational space, and Iwai's monopole in the corresponding half-space", the latter for mirror-identified triangles — [Anderson, arXiv:1711.10115](https://arxiv.org/abs/1711.10115); [arXiv:1802.03465](https://arxiv.org/abs/1802.03465)
- **Littlejohn & Reinsch (Rev. Mod. Phys. 69 (1997) 213–276).** They review separating rotations from internal motions in n-body systems. The reduced configuration space carries non-Abelian gauge fields, those of "the 'falling cat' problem", which appear as gauge potentials (Coriolis effects) in the internal dynamics. Among the new results are Kaluza–Klein-type field equations and "geometrical interpretations of the Eckart frame" — [APS abstract page](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.69.213). Only the abstract summary was seen. Anderson cites them for internal/shape coordinates — [arXiv:1711.10115 refs.](https://arxiv.org/abs/1711.10115)
- **Eckart frame = one-shot Procrustes.** Kudin & Dymarsky show that "the rotation matrix minimizing the root-mean-square deviation between two molecular conformations" (Kabsch 1976) also satisfies the Eckart axis conditions — [J. Chem. Phys. 122 (2005) 224105, DOI 10.1063/1.1929739](https://doi.org/10.1063/1.1929739). Only the abstract highlight was seen.
- **Falling cat and isoholonomic problems.**
  - Montgomery, "Isoholonomic problems and some applications", Commun. Math. Phys. 128 (1990) 565–592 — [DOI 10.1007/BF02096874](https://doi.org/10.1007/BF02096874)
  - "How much does the rigid body rotate? A Berry's phase from the 18th century", Am. J. Phys. 59 (1991) 394–398 — [DOI 10.1119/1.16514](https://doi.org/10.1119/1.16514)
  - "Gauge theory of the falling cat", Fields Inst. Commun. 1 (1993) 193–218, in Enos (ed.), Dynamics and Control of Mechanical Systems: The Falling Cat and Related Problems (AMS). This is a secondary citation, recorded in [Wikipedia, Falling cat problem](https://en.wikipedia.org/wiki/Falling_cat_problem).
- **Shapere & Wilczek.**
  - "Self-propulsion at low Reynolds number", PRL 58 (1987) 2051–2054 — [DOI](https://doi.org/10.1103/PhysRevLett.58.2051)
  - "Geometry of self-propulsion at low Reynolds number", J. Fluid Mech. 198 (1989) 557–585 — [DOI](https://doi.org/10.1017/S002211208900025X)
  - "Gauge kinematics of deformable bodies", Am. J. Phys. 57 (1989) 514–518 — [DOI](https://doi.org/10.1119/1.15986). Records were seen via Crossref.
- **Best matching as a connection (Gryb 2010).** In best-matching theories, D_ω q^a = q̇^a + ω̇^α t_α^a_b q^b "is the covariant derivative of q with connection ω̇ along a trial curve γ … the pullback onto γ of a connection on the principal G-bundle". Solving the constraints gives corrected coordinates obtained "by subtracting the action of the open-path holonomy (P exp …) … This subtracts all vertical motion of the q's along the fiber bundle". For rotations the linear constraint "is the vanishing of total angular momentum" — [Gryb, CQG 27 (2010) 215018, arXiv:1003.1973, eqs. (12), (19), (43)](https://arxiv.org/abs/1003.1973)
- **Best matching as orthogonal projection (Gomes & Gryb 2021).**
  - Their bundle projection, with the connection defined by orthogonality to the gauge orbits in the kinetic metric, "is formally similar to the best-matching construction of Barbour (cf. Barbour & Bertotti (1982) and Mercati (2017, Ch. II.5))".
  - The work builds on Littlejohn & Reinsch.
  - "the rotational bundle is curved … curvature prohibits the existence of a horizontal section".
  - Barbour–Bertotti "requires the angular momentum of the universe to be zero".

  They then extend relationalism to nonzero angular momentum with an SO(3) charge — [Gomes & Gryb, Stud. Hist. Phil. Sci. 88 (2021) 138–155, arXiv:2011.01693](https://arxiv.org/abs/2011.01693)
- **Barbour's own infinitesimal statement.** He defines ds_bm := min over orbits of √((E − V)·Σ ½ m_a dx_a·dx_a). "The connecting 'strut' is orthogonal with respect to the supermetric on Q^N in the best-matched position. Best matching brings the centres of mass to coincidence and reduces the net rotation to zero" — [Barbour, "Shape Dynamics. An Introduction", in Quantum Field Theory and Gravity (Birkhäuser 2012) 257–297, arXiv:1105.0183, §3.2, eq. (9) and figure caption](https://arxiv.org/abs/1105.0183)
- **The finite best-matched distance.** d_BM(q^i, q^f) = inf over Ω ∈ SO(3) and θ of (Σ_a ‖r_a^i − Ω r_a^f − θ‖²)^{1/2}. Translations are fixed by passing to barycentric coordinates, and stationarity under Ω → (I + δω×)Ω gives the rotational condition — [Mercati, "A Shape Dynamics Tutorial", arXiv:1409.0105, eqs. (8)–(15)](https://arxiv.org/abs/1409.0105)
- **Philosophy literature on the holonomy.** Batterman discusses falling cats and holonomy-based explanation in "Falling cats, parallel parking, and polarized light", Stud. Hist. Phil. Mod. Phys. 34 (2003) 527–557 — [DOI 10.1016/S1355-2198(03)00062-5](https://doi.org/10.1016/S1355-2198(03)00062-5). Only the record was seen. Pooley & Brown (2002) treat vanishing angular momentum "to be a real-world prediction of a fully relational theory", as reported by [Gomes & Gryb, arXiv:2011.01693](https://arxiv.org/abs/2011.01693).

### Inferences
- **Best matching = horizontality.** Mercati's finite d_BM, Montgomery's quotient metric (eq. 23) and Barbour–Bertotti best matching between two configurations all minimise a (mass-)weighted Euclidean distance over the Euclidean group. That is ordinary Procrustes superimposition without scaling, i.e. the size-and-shape metric. The first-order condition for the optimal rotation, Σ m_a Δr_a^i × Ω Δr_a^f = 0, becomes Σ m_a r_a × dr_a = 0, i.e. J = 0, as the step shrinks. So infinitesimal best matching is the horizontality condition of the mechanical connection (Montgomery Lemma 1; Barbour's orthogonal "strut").
- **Chained vs one-shot = holonomy.** Let a path of shapes γ(t) run from shape s0 to shape s1.
  - Chained infinitesimal best matching gives the horizontal lift of γ.
  - A one-shot Procrustes fit between the two endpoint configurations gives the horizontal lift of the shape-space **geodesic** from s0 to s1. Montgomery's Thm 3 says geodesics are projections of horizontal lines, and the optimal Procrustes pair lies on such a line.
  - The two results differ by a rotation equal to the holonomy of the closed loop "γ, then the geodesic back".
  - For planar triangles, that rotation is Δθ = ½·Ω(loop) (mod 2π, sign by orientation), with Ω the solid angle on the unit sphere enclosed by the loop. Equivalently, Δθ = 2·Area on Kendall's S²(1/2). For a closed shape loop (s1 = s0) the one-shot fit is the identity and the chained fit rotates by Ω/2.
- **Numerical check (my computation this session; not from the literature).** I built triangles from w through the inverse Hopf map, for unit masses and for masses (1, 2, 5), in mass-weighted Jacobi coordinates. I then chained exact weighted Procrustes rotation fits around circles on the shape sphere. Pseudocode:
  ```
  for each k: c = Σ m_i q_k,i · conj(ref_i)
              ref ← e^{−i·arg c} · q_k
  ```
  - For a cap of polar radius α the net rotation converged to π(1 − cos α) = Ω/2. For example, α = 1.0 gave 1.41630 at N = 12 steps, 1.44379 at N = 100, and 1.44418 at N = 2000, against Ω/2 = 1.44418.
  - The sign flipped with loop orientation. The result was the same for loops centred on a collinear (equatorial) shape and for unequal masses.
  - The one-shot fit of last against first returned rotation ≈ 0 (about 10⁻¹⁷).
  - The finite-N chained value equalled, to all printed digits, half the solid angle of the **geodesic polygon** through the sampled shapes. So a discretely chained Procrustes fit is exact parallel transport along a geodesic polygon, and its error relative to the smooth loop is purely the polygon-vs-curve area difference.
- **Eckart frame vs zero-angular-momentum frame.** In molecular physics this is the same distinction as in the bullet above. The Eckart frame (one-shot Procrustes to a fixed reference, per Kudin & Dymarsky) is a gauge choice, a section of the bundle. The Guichardet frame (zero angular momentum) is parallel transport. They differ by the holonomy, which Littlejohn & Reinsch study as gauge-field curvature.
- **Barbour's Jacobi factor does not change the holonomy.** Barbour's best matching weights the metric by (E − V), but a conformal factor rescales lengths, not orthogonality. The horizontal space, and hence the holonomy, is the same as for the plain mass metric.

### Gaps
- I did not see the body of Littlejohn & Reinsch (paywalled), so I cannot quote their three-body monopole formula or their Eckart-frame section directly.
- I did not see the full text of Montgomery's 1993 falling-cat chapter.
- The sign convention in "Δθ = ±Ω/2" depends on the orientation of the shape sphere and of the loop. Montgomery states the planar geometric phase as the normalised area form with total 2π, i.e. Ω/2 up to orientation.

---

## 5. Does any literature explicitly call Barbour's best matching "Procrustes"?

### Takeaway
I found only one explicit statement: Anderson (2015) says best matching "bears much conceptual similarity to the Procrustes procedure of Shape Geometry", while stressing that they extremise different objects. The core shape-dynamics texts I checked (Barbour, Mercati, Gryb, Gomes–Gryb) define best matching in a way that is mathematically a Procrustes fit, but none of them uses the word. Kendall's 1984 title ("Procrustean metrics") and Kudin & Dymarsky's Eckart–RMSD result are the closest links from the statistics and chemistry side.

### Cited Findings
- Anderson: "Matching is a procedure that bears much conceptual similarity to the Procrustes procedure of Shape Geometry [10]; however these each involve extremizing a different object" — [Anderson, arXiv:1307.1923 / Stud. Hist. Phil. Mod. Phys. 51 (2015) 1–8](https://arxiv.org/abs/1307.1923)
- Mercati's tutorial defines the "best-matched distance" d_BM = inf over Ω, θ of ‖r^i − Ωr^f − θ‖. This is the ordinary Procrustes criterion, but the word "Procrustes" does not appear in the text — [arXiv:1409.0105](https://arxiv.org/abs/1409.0105). I checked this by full-text search this session.
- Full-text searches this session found no occurrence of "Procrustes" in:
  - Barbour, arXiv:1105.0183, arXiv:1007.3368 and arXiv:gr-qc/0211021;
  - Gryb, arXiv:1003.1973;
  - Gomes & Gryb, arXiv:2011.01693;
  - Anderson, arXiv:1111.1472, arXiv:1503.01507, arXiv:1711.10054 and arXiv:1711.11492;
  - "The Physics and Metaphysics of Pure Shape Dynamics", arXiv:2203.02773.

  Sources: [arXiv:1105.0183](https://arxiv.org/abs/1105.0183), [arXiv:1007.3368](https://arxiv.org/abs/1007.3368), [arXiv:gr-qc/0211021](https://arxiv.org/abs/gr-qc/0211021), [arXiv:1003.1973](https://arxiv.org/abs/1003.1973), [arXiv:2011.01693](https://arxiv.org/abs/2011.01693), [arXiv:1111.1472](https://arxiv.org/abs/1111.1472), [arXiv:1503.01507](https://arxiv.org/abs/1503.01507), [arXiv:1711.10054](https://arxiv.org/abs/1711.10054), [arXiv:1711.11492](https://arxiv.org/abs/1711.11492), [arXiv:2203.02773](https://arxiv.org/abs/2203.02773)
- A web search for "best matching" together with "Procrustes" in relational mechanics surfaced only Anderson's paper as an explicit link — [arXiv HTML of 1307.1923](https://arxiv.org/html/1307.1923)
- Montgomery's shape-space metric (eq. 23, the infimum over rigid motions) is the Procrustes construction, but he does not use the term — [arXiv:1402.0841](https://arxiv.org/abs/1402.0841)
- Kendall's 1984 paper is titled "Shape manifolds, Procrustean metrics, and complex projective spaces" — [DOI 10.1112/blms/16.2.81](https://doi.org/10.1112/blms/16.2.81)
- Kudin & Dymarsky: the RMSD-minimising rotation (Kabsch) satisfies the Eckart conditions — [DOI 10.1063/1.1929739](https://doi.org/10.1063/1.1929739); [Kabsch 1976, DOI 10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873)

### Inferences
- Calling Barbour–Bertotti's two-configuration best matching "Procrustes superimposition (without scaling, mass-weighted)" is mathematically exact. Calling its infinitesimal, path-chained use "Procrustes" hides the holonomy discussed in §4, which may be why Anderson says the two extremise different objects.

### Gaps
- I could not search the full texts of Anderson's The Problem of Time (2017), Mercati's Shape Dynamics: Relativity and Relationalism (OUP 2018), Barbour's The End of Time (1999) or The Janus Point (2020) for "Procrustes".
- I did not run a systematic philosophy-database search (PhilPapers, PhilSci-Archive). A 2021–2026 philosophy paper using the term could exist without being found here.

---

## Source catalogue

"verified" = yes means I saw a bibliographic record this session (arXiv API, Crossref API, Math-Net.ru, Numdam/EuDML search results, or the author PDF). "partial" means the record was seen but some fields come from secondary sources or memory. "no" means the entry is from memory or secondary citation only.

| key | authors | year | title | venue/volume/pages | DOI | arXiv ID | open-access PDF URL | verified |
|---|---|---|---|---|---|---|---|---|
| kendall1977diffusion | Kendall, David G. | 1977 | The diffusion of shape | Advances in Applied Probability 9(3):428–430 | 10.2307/1426091 | | | yes |
| kendall1984shape | Kendall, David G. | 1984 | Shape manifolds, Procrustean metrics, and complex projective spaces | Bulletin of the London Mathematical Society 16(2):81–121 | 10.1112/blms/16.2.81 | | | yes |
| kendall1985exact | Kendall, David G. | 1985 | Exact distributions for shapes of random triangles in convex sets | Advances in Applied Probability 17(2):308–329 | 10.2307/1427143 | | | yes |
| kendall1989survey | Kendall, David G. | 1989 | A survey of the statistical theory of shape | Statistical Science 4(2):87–99 | 10.1214/ss/1177012582 | | | yes |
| le1993riemannian | Le, Huiling; Kendall, David G. | 1993 | The Riemannian structure of Euclidean shape spaces: a novel environment for statistics | Annals of Statistics 21(3):1225–1271 | 10.1214/aos/1176349259 | | | yes |
| kendall1999shape | Kendall, David G.; Barden, Dennis; Carne, T. K.; Le, Huiling | 1999 | Shape and Shape Theory | Wiley Series in Probability and Statistics, Wiley, Chichester | 10.1002/9780470317006 | | | partial (Crossref record lacks author list; authors from memory) |
| goodall1991procrustes | Goodall, Colin | 1991 | Procrustes methods in the statistical analysis of shape | Journal of the Royal Statistical Society Series B 53(2):285–321 | 10.1111/j.2517-6161.1991.tb01825.x | | | yes |
| dryden2016statistical | Dryden, Ian L.; Mardia, Kanti V. | 2016 | Statistical Shape Analysis, with Applications in R (2nd ed.) | Wiley Series in Probability and Statistics, Wiley | 10.1002/9781119072492 | | | yes |
| kabsch1976solution | Kabsch, Wolfgang | 1976 | A solution for the best rotation to relate two sets of vectors | Acta Crystallographica Section A 32(5):922–923 | 10.1107/S0567739476001873 | | | yes |
| kudin2005eckart | Kudin, Konstantin N.; Dymarsky, Anatoly Y. | 2005 | Eckart axis conditions and the minimization of the root-mean-square deviation: two closely related problems | Journal of Chemical Physics 122(22):224105 | 10.1063/1.1929739 | | | yes |
| montgomery2015shape | Montgomery, Richard | 2015 | The three-body problem and the shape sphere | American Mathematical Monthly 122(4):299–321 | 10.4169/amer.math.monthly.122.04.299 | 1402.0841 | https://arxiv.org/pdf/1402.0841 | yes (end page 321 from a citing source; Crossref gives first page only) |
| montgomery1996geometric | Montgomery, Richard | 1996 | The geometric phase of the three-body problem | Nonlinearity 9(5):1341–1360 | 10.1088/0951-7715/9/5/014 | dg-ga/9510005 | https://arxiv.org/pdf/dg-ga/9510005 | yes |
| montgomery2002syzygies | Montgomery, Richard | 2002 | Infinitely many syzygies | Archive for Rational Mechanics and Analysis 164(4):311–340 | 10.1007/s00205-002-0211-z | | https://people.ucsc.edu/~rmont/papers/InfinitelyManySyzygies.pdf | yes |
| montgomery2005pants | Montgomery, Richard | 2005 | Fitting hyperbolic pants to a three-body problem | Ergodic Theory and Dynamical Systems 25(3):921–947 | 10.1017/S0143385704000653 | math/0405014 | https://arxiv.org/pdf/math/0405014 | yes |
| montgomery2007zero | Montgomery, Richard | 2007 | The zero angular momentum, three-body problem: all but one solution has syzygies | Ergodic Theory and Dynamical Systems 27(6):1933–1946 | 10.1017/S0143385707000338 | | | yes |
| chenciner2000remarkable | Chenciner, Alain; Montgomery, Richard | 2000 | A remarkable periodic solution of the three-body problem in the case of equal masses | Annals of Mathematics (2) 152(3):881–901 | 10.2307/2661357 | math/0011268 | https://arxiv.org/pdf/math/0011268 | yes |
| montgomery1998nbody | Montgomery, Richard | 1998 | The N-body problem, the braid group, and action-minimizing periodic solutions | Nonlinearity 11(2):363–376 | 10.1088/0951-7715/11/2/011 | | | yes |
| moore1993braids | Moore, Cristopher | 1993 | Braids in classical dynamics | Physical Review Letters 70(24):3675–3679 | 10.1103/PhysRevLett.70.3675 | | | yes (Montgomery's preprint cites it as "Braids in Classical Gravity"; Crossref title used) |
| moeckel2013symmetric | Moeckel, Richard; Montgomery, Richard | 2013 | Symmetric regularization, reduction and blow-up of the planar three-body problem | Pacific Journal of Mathematics 262(1):129–189 | 10.2140/pjm.2013.262.129 | 1202.0972 | https://arxiv.org/pdf/1202.0972 | yes |
| moeckel2014realizing | Moeckel, Richard; Montgomery, Richard | 2014 | Realizing all free homotopy classes for the Newtonian three-body problem | arXiv preprint (journal version not checked) | | 1412.2263 | https://arxiv.org/pdf/1412.2263 | yes (arXiv record only) |
| moeckel2012brake | Moeckel, Richard; Montgomery, Richard; Venturelli, Andrea | 2012 | From brake to syzygy | Archive for Rational Mechanics and Analysis 204(3):1009–1060 | 10.1007/s00205-012-0502-y | 1109.1089 | | partial (arXiv ID from a citing bibliography) |
| moeckel1988qualitative | Moeckel, Richard | 1988 | Some qualitative features of the three-body problem | Contemporary Mathematics 81 (Hamiltonian Dynamical Systems, eds. K. Meyer, D. Saari), AMS, 1–22 | | | | no (from Montgomery 1996 bibliography) |
| moeckel2015central | Moeckel, Richard | 2015 | Central configurations | In: Llibre, J.; Moeckel, R.; Simó, C., Central Configurations, Periodic Orbits, and Hamiltonian Systems, Advanced Courses in Mathematics – CRM Barcelona, Birkhäuser, 105–167 | 10.1007/978-3-0348-0933-7_2 | | | yes |
| hsiang2007kinematic | Hsiang, Wu-Yi; Straume, Eldar | 2007 | Kinematic geometry of triangles and the study of the three-body problem | Lobachevskii Journal of Mathematics 25:9–130 (MR2312347; Zbl 1118.53011) | | math-ph/0608060 | https://arxiv.org/pdf/math-ph/0608060 | yes |
| lemaitre1954regularisation | Lemaître, Georges | 1954 | Régularisation dans le problème des trois corps | Bulletin de la Classe des Sciences, Académie Royale de Belgique 40:759–767 | | | | no (from Moeckel–Montgomery 2013 bibliography) |
| lemaitre1964three | Lemaître, Georges | 1964 | The three body problem | NASA Contractor Report CR-110 | | | https://ntrs.nasa.gov/ (per Moeckel–Montgomery) | no |
| murnaghan1936symmetric | Murnaghan, F. D. | 1936 | A symmetric reduction of the planar three-body problem | American Journal of Mathematics 58(4):829–832 | 10.2307/2371252 | | | yes (end page from Moeckel–Montgomery bibliography) |
| smale1970topology1 | Smale, Stephen | 1970 | Topology and mechanics. I | Inventiones Mathematicae 10(4):305–331 | 10.1007/BF01418778 | | | yes |
| smale1970topology2 | Smale, Stephen | 1970 | Topology and mechanics. II (The planar n-body problem) | Inventiones Mathematicae 11(1):45–64 | 10.1007/BF01389805 | | | yes |
| guichardet1984rotation | Guichardet, Alain | 1984 | On rotation and vibration motions of molecules | Annales de l'I.H.P. Physique théorique 40(3):329–342 | | | https://www.numdam.org/item/AIHPA_1984__40_3_329_0.pdf | yes |
| iwai1987geometric | Iwai, Toshihiro | 1987 | A geometric setting for classical molecular dynamics | Annales de l'I.H.P. Physique théorique 47(2):199–219 | | | http://www.numdam.org/item/AIHPA_1987__47_2_199_0/ | yes |
| iwai1987gauge | Iwai, Toshihiro | 1987 | A gauge theory for the quantum planar three-body problem | Journal of Mathematical Physics 28(4):964–974 | 10.1063/1.527588 | | | yes |
| iwai1987internal | Iwai, Toshihiro | 1987 | A geometric setting for internal motions of the quantum three-body system | Journal of Mathematical Physics 28(6):1315–1326 | 10.1063/1.527534 | | | yes |
| shapere1987self | Shapere, Alfred; Wilczek, Frank | 1987 | Self-propulsion at low Reynolds number | Physical Review Letters 58(20):2051–2054 | 10.1103/PhysRevLett.58.2051 | | | yes |
| shapere1989geometry | Shapere, Alfred; Wilczek, Frank | 1989 | Geometry of self-propulsion at low Reynolds number | Journal of Fluid Mechanics 198:557–585 | 10.1017/S002211208900025X | | | yes |
| shapere1989gauge | Shapere, Alfred; Wilczek, Frank | 1989 | Gauge kinematics of deformable bodies | American Journal of Physics 57(6):514–518 | 10.1119/1.15986 | | | yes |
| shapere1989geometric | Shapere, Alfred; Wilczek, Frank (eds.) | 1989 | Geometric Phases in Physics | World Scientific, Singapore | | | | no |
| littlejohn1995internal | Littlejohn, Robert G.; Reinsch, Matthias | 1995 | Internal or shape coordinates in the n-body problem | Physical Review A 52(3):2035–2051 | 10.1103/PhysRevA.52.2035 | | | yes |
| littlejohn1997gauge | Littlejohn, Robert G.; Reinsch, Matthias | 1997 | Gauge fields in the separation of rotations and internal motions in the n-body problem | Reviews of Modern Physics 69(1):213–276 | 10.1103/RevModPhys.69.213 | | | yes |
| montgomery1990isoholonomic | Montgomery, Richard | 1990 | Isoholonomic problems and some applications | Communications in Mathematical Physics 128(3):565–592 | 10.1007/BF02096874 | | | yes |
| montgomery1991rigid | Montgomery, Richard | 1991 | How much does the rigid body rotate? A Berry's phase from the 18th century | American Journal of Physics 59(5):394–398 | 10.1119/1.16514 | | | yes |
| montgomery1993falling | Montgomery, Richard | 1993 | Gauge theory of the falling cat | In: Enos, M. J. (ed.), Dynamics and Control of Mechanical Systems: The Falling Cat and Related Problems, Fields Institute Communications 1, American Mathematical Society, 193–218 | | | | partial (secondary sources only) |
| marsden1990reduction | Marsden, Jerrold E.; Montgomery, Richard; Ratiu, Tudor | 1990 | Reduction, symmetry, and phases in mechanics | Memoirs of the American Mathematical Society 88(436) | 10.1090/memo/0436 | | | yes |
| batterman2003falling | Batterman, Robert W. | 2003 | Falling cats, parallel parking, and polarized light | Studies in History and Philosophy of Modern Physics 34(4):527–557 | 10.1016/S1355-2198(03)00062-5 | | | yes |
| barbour1982mach | Barbour, Julian B.; Bertotti, Bruno | 1982 | Mach's principle and the structure of dynamical theories | Proceedings of the Royal Society of London A 382(1783):295–306 | 10.1098/rspa.1982.0102 | | | yes |
| barbour2003scale | Barbour, Julian | 2003 | Scale-invariant gravity: particle dynamics | Classical and Quantum Gravity 20(8):1543–1570 | 10.1088/0264-9381/20/8/310 | gr-qc/0211021 | https://arxiv.org/pdf/gr-qc/0211021 | yes |
| barbour2010definition | Barbour, Julian | 2010 | The definition of Mach's principle | Foundations of Physics 40:1263–1284 | 10.1007/s10701-010-9490-7 | 1007.3368 | https://arxiv.org/pdf/1007.3368 | yes |
| barbour2012shape | Barbour, Julian | 2012 | Shape dynamics. An introduction | In: Finster, F. et al. (eds.), Quantum Field Theory and Gravity, Birkhäuser, 257–297 | 10.1007/978-3-0348-0043-3_13 | 1105.0183 | https://arxiv.org/pdf/1105.0183 | yes |
| mercati2014tutorial | Mercati, Flavio | 2014 | A shape dynamics tutorial | arXiv preprint | | 1409.0105 | https://arxiv.org/pdf/1409.0105 | yes |
| mercati2018shape | Mercati, Flavio | 2018 | Shape Dynamics: Relativity and Relationalism | Oxford University Press | 10.1093/oso/9780198789475.001.0001 | | | yes |
| gryb2010definition | Gryb, Sean | 2010 | A definition of background independence | Classical and Quantum Gravity 27(21):215018 | 10.1088/0264-9381/27/21/215018 | 1003.1973 | https://arxiv.org/pdf/1003.1973 | yes |
| gomes2021angular | Gomes, Henrique; Gryb, Sean | 2021 | Angular momentum without rotation: turbocharging relationalism | Studies in History and Philosophy of Science 88:138–155 | 10.1016/j.shpsa.2021.05.006 | 2011.01693 | https://arxiv.org/pdf/2011.01693 | yes |
| pooley2002relationalism | Pooley, Oliver; Brown, Harvey R. | 2002 | Relationalism rehabilitated? I: Classical mechanics | British Journal for the Philosophy of Science 53(2):183–204 | 10.1093/bjps/53.2.183 | | | yes |
| anderson2008foundations | Anderson, Edward | 2008 | Foundations of relational particle dynamics | Classical and Quantum Gravity 25(2):025003 | 10.1088/0264-9381/25/2/025003 | 0706.3934 | https://arxiv.org/pdf/0706.3934 | yes |
| anderson2009triangleland1 | Anderson, Edward | 2009 | Triangleland. I. Classical dynamics with exchange of relative angular momentum | Classical and Quantum Gravity 26(13):135020 | 10.1088/0264-9381/26/13/135020 | 0809.1168 | https://arxiv.org/pdf/0809.1168 | yes |
| anderson2009triangleland2 | Anderson, Edward | 2009 | Triangleland. II. Quantum mechanics of pure shape | Classical and Quantum Gravity 26(13):135021 | 10.1088/0264-9381/26/13/135021 | 0809.3523 | https://arxiv.org/pdf/0809.3523 | yes |
| anderson2011problem | Anderson, Edward | 2011 | The problem of time and quantum cosmology in the relational particle mechanics arena | arXiv preprint | | 1111.1472 | https://arxiv.org/pdf/1111.1472 | yes |
| anderson2014quadrilateralland1 | Anderson, Edward | 2014 | Relational quadrilateralland I: the classical theory | International Journal of Modern Physics D 23(2):1450014 | 10.1142/S021827181450014X | 1202.4186 | https://arxiv.org/pdf/1202.4186 | yes |
| anderson2014quadrilateralland2 | Anderson, Edward; Kneller, Sophie | 2014 | Relational quadrilateralland II: the quantum theory | International Journal of Modern Physics D 23(5):1450052 | 10.1142/S0218271814500527 | 1303.5645 | https://arxiv.org/pdf/1303.5645 | yes (arXiv journal-ref lists Paper I's DOI; Crossref DOI used) |
| anderson2015kendall | Anderson, Edward | 2015 | Kendall's shape statistics as a classical realization of Barbour-type timeless records theory approach to quantum gravity | Studies in History and Philosophy of Modern Physics 51:1–8 | 10.1016/j.shpsb.2015.05.002 | 1307.1923 | https://arxiv.org/pdf/1307.1923 | yes |
| anderson2015configuration | Anderson, Edward | 2015 | Configuration spaces in fundamental physics | arXiv preprint | | 1503.01507 | https://arxiv.org/pdf/1503.01507 | yes |
| anderson2016six | Anderson, Edward | 2016 | Six new mechanics corresponding to further shape theories | International Journal of Modern Physics D 25:1650044 | 10.1142/S0218271816500449 | 1505.00488 | https://arxiv.org/pdf/1505.00488 | yes |
| anderson2017problem | Anderson, Edward | 2017 | The Problem of Time: Quantum Mechanics Versus General Relativity | Fundamental Theories of Physics 190, Springer, Cham | 10.1007/978-3-319-58848-3 | | | partial (Crossref title "The Problem of Time" and series seen; subtitle and vol. 190 from memory or brief) |
| anderson2017smallest1 | Anderson, Edward | 2017 | The smallest shape spaces. I. Shape theory posed, with example of 3 points on the line | arXiv preprint | | 1711.10054 | https://arxiv.org/pdf/1711.10054 | yes |
| anderson2017smallest2 | Anderson, Edward | 2017 | The smallest shape spaces. II. 4 points in 1-d suffices to have a complex background-independent theory of inhomogeneity | arXiv preprint | | 1711.10073 | https://arxiv.org/pdf/1711.10073 | yes |
| anderson2017smallest3 | Anderson, Edward | 2017 | The smallest shape spaces. III. Triangles in 2- and 3-d | arXiv preprint | | 1711.10115 | https://arxiv.org/pdf/1711.10115 | yes |
| anderson2017alice | Anderson, Edward | 2017 | Alice in Triangleland: Lewis Carroll's pillow problem and variants solved on shape space of triangles | arXiv preprint | | 1711.11492 | https://arxiv.org/pdf/1711.11492 | yes |
| anderson2018monopoles | Anderson, Edward | 2018 | Monopoles of twelve types in 3-body problems | arXiv preprint | | 1802.03465 | https://arxiv.org/pdf/1802.03465 | yes |
| anderson2018topological | Anderson, Edward | 2018 | Topological shape theory | arXiv preprint | | 1803.11126 | https://arxiv.org/pdf/1803.11126 | yes |
| kuang2021shape | Kuang, Wentian | 2021 | On shape sphere and rotation of three-body motion | arXiv preprint | | 2106.15106 | https://arxiv.org/pdf/2106.15106 | yes |
| krishnaswami2016curvature | Krishnaswami, Govind S.; Senapati, Himalaya | 2016 | Curvature and geodesic instabilities in a geometrical approach to the planar three-body problem | Journal of Mathematical Physics 57:102901 | 10.1063/1.4964340 | 1606.05091 | https://arxiv.org/pdf/1606.05091 | yes |
| garaud2026geometric | Garaud, Julien; Niemi, Antti J. | 2026 | Geometric qubits in programmable atomic trimers | arXiv preprint | | 2607.16846 | https://arxiv.org/pdf/2607.16846 | yes |
| vassallo2022physics | Vassallo, Antonio; Naranjo, Pedro; Koslowski, Tim | 2022 | The physics and metaphysics of Pure Shape Dynamics | arXiv preprint (cited only as a negative full-text search) | | 2203.02773 | https://arxiv.org/pdf/2203.02773 | yes (arXiv record) |
