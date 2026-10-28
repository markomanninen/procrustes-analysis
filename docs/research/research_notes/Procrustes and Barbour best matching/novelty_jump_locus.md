# Prior-art check: the jump locus and the rank rule for unlabelled, mirror-allowed triangle Procrustes

These notes check two claims made about the triangle page (`js/procrustes.js`), which fits DEF onto ABC by trying 6 vertex pairings × {proper, mirrored} = 12 candidates and keeping the largest |s|, s = Σ z̄_i u_π(i):

- **Claim A (jump locus).** With mirrors allowed, the best candidate changes exactly when the target ABC passes through an isosceles or collinear shape, whatever DEF is. When DEF passes through such a shape, only the labels change and the fitted triangle moves continuously. With mirrors off, the switching set depends on both triangles.
- **Claim B (rank rule).** With mirrors allowed, the best candidate pairs vertices by the rank of the opposite side length (equivalently by angle), and mirrors exactly when ABC and the relabelled DEF have opposite signed area. The caller reports that this matched exhaustive search in 40,000 of 40,000 random trials, with and without scaling. Ranking alone matched about 60% of the time with mirrors off. These numbers come from the caller's session and were not re-run here.
- **Proposed mechanism.** For a finite reflection group W and x, y in one closed chamber, ⟨x, w·y⟩ ≤ ⟨x, y⟩ for every w ∈ W.

**Search status legend used below.** "Read" means I saw the full text or the relevant passage this session. "Record only" means I saw a bibliographic record, abstract or search-engine extract, not the text. "Not searched" means no search was made. About 30 tool calls were used: web searches, fetches of ADS, arXiv and author PDFs, and local text extraction from downloaded PDFs. Sibling notes in this folder (`shape_sphere_bridge.md`, `correspondence_soft_alignment.md`) already cover the shape-sphere geometry, Rubix, and molecular RMSD. They are referred to, not repeated.

---

## Q1. Is Claim A or Claim B stated in Kendall's shape-theory literature (spherical blackboard), Small (1996), or Dryden & Mardia ch. 14?

### Takeaway
The geometry behind both claims is classical. It includes the 12-tile division of the triangle shape sphere by the collinear equator and the three isosceles great circles, and the reading of one tile as "unlabelled, unoriented triangles" (Kendall's spherical blackboard). In the material I could read, I found no statement that ranking by side length plus an orientation check gives the least-squares Procrustes optimum, and none locating the jumps of the hard fit. Kendall's primary texts, Small (1996), Kendall–Barden–Carne–Le (1999) and Dryden & Mardia ch. 14 were **not readable** this session. So for this literature the status is "searched, not found in accessible material", not "absent".

### Cited Findings
- **Blackboard as a tile of the shape sphere (secondary sources).** Anderson calls the smallest triangle quotient "the Leibniz space demilune of width π/3: Kendall's spherical blackboard", "the 1/6th hemisphere bounded by arcs of isoscelesness and collinearity" and "the Leibs(3,2) 1/12-sphere". These quotes were verified in the sibling note — [Anderson, arXiv:1711.10115](https://arxiv.org/abs/1711.10115); [Anderson, arXiv:1711.11492](https://arxiv.org/abs/1711.11492); see `shape_sphere_bridge.md`.
- **Behrend's tutorial on the blackboard.** A search-engine extract quotes it: "The 6 ways of labelling a given triangle will usually produce 6 different shapes. The exceptions are: equilateral triangle (2 shapes); and collinear triangle with one point midway between the other two (3 shapes)". It also says the dividing planes contain "collinear shapes – separating – clockwise and anticlockwise labelled hemispheres of triangles, isosceles shapes – separating left- and right-leaning hemispheres of scalene triangles". A direct fetch returned HTTP 403, so the text is unverified beyond the extract — [Behrend, "Kendall's spherical blackboard"](https://cantab.net/users/michael.behrend/ley_stats/spher_black/pages/spher_black.html)
- **Mirror structure on the sphere (equal masses).** The isosceles triangles form "three great circles on the shape sphere, each circle passing through the two Lagrange points". A transposition acts as a half-twist about an axis through an Euler point and a binary-collision point. Composed with the planar reflection, it acts as "reflection about the plane ISOSC₂". Verified in the sibling notes — [Montgomery, arXiv:1402.0841](https://arxiv.org/abs/1402.0841)
- **Kendall 1984 (Bull. LMS)** has sections on "the manifold carrying the shapes of triangles" and "the use of the 'spherical blackboard'" (search-engine summary). The PDF mirror returned no readable text — [scispace/typeset mirror](https://typeset.io/pdf/shape-manifolds-procrustean-metrics-and-complex-projective-1sjgkp2ynn.pdf)
- **Kendall 1985** (Adv. Appl. Prob. 17:308–329), the paper where the blackboard is used for exact shape distributions of random triangles. Record only (Crossref, sibling note); text not read — [DOI 10.2307/1427143](https://doi.org/10.2307/1427143)
- **Kendall 1989 survey** (Statist. Sci. 4(2):87–99) is open access on Project Euclid, but the fetch returned an empty page. Content not checked — [Project Euclid](https://projecteuclid.org/euclid.ss/1177012582)
- **Shape and Shape Theory** (Kendall, Barden, Carne & Le, Wiley 1999). Record only, ISBN 978-0-471-96823-8. Chapter content not seen — [Wiley listing](https://gov.wiley.com/WileyCDA/WileyTitle/productCd-0471968234.html)
- **Dryden & Mardia (2016), ch. 14 "Unlabelled size-and-shape and shape analysis"**, pp. 341–352. The table of contents (sibling note) shows sections on the Green–Mardia model, a Procrustes model and related methods. These are general unlabelled point-set matching by MCMC or optimisation. I saw nothing about an analytic triangle rule — [DOI 10.1002/9781119072492.ch14](https://doi.org/10.1002/9781119072492.ch14)
- **A closely related tiling in a different metric (read).** Bowden, Haynes, Shonkwiler & Shukert (Geom. Dedicata 2019) model triangles by the Hausmann–Knutson/Grassmannian sphere, not Kendall's.
  - "The fundamental domain of this action is a region … whose boundary corresponds to those n-gons which are unchanged by some group element".
  - For triangles: "Its boundary is precisely the set of isosceles and degenerate triangles and its interior is the set of scalene triangles".
  - The tile with |z| ≤ y ≤ x satisfies a ≤ b ≤ c, "so D parametrizes triangles with side lengths written in ascending order".
  - In that metric the sphere is cut into "48 congruent 45–60–90 triangles" (group B₃), not Kendall's 12 tiles of 60–90–90.
  - They use the tile to find the "least symmetric triangle", the point furthest from the walls. They do not use it for matching two triangles.
  — [arXiv:1708.01559](https://arxiv.org/abs/1708.01559); [Geom. Dedicata 198:19–34, DOI 10.1007/s10711-018-0327-4](https://doi.org/10.1007/s10711-018-0327-4)

### Inferences
- **The chamber is determined by side order and orientation.** On Kendall's sphere each isosceles great circle is the locus where two labelled side lengths are equal (Montgomery). So the six lunes between the three circles are the six strict orderings of (|AB|, |BC|, |CA|), and the equator splits each lune by the sign of the area. Naming the chamber of a labelled triangle is therefore the same as naming its side-length order and its orientation. The rank rule (Claim B) is the folding map written in these coordinates. Bowden et al. state the side-order half of this dictionary for their own metric. I did not find it stated for Kendall's metric in accessible text, but it is elementary.
- **The blackboard picture already implies Claim A in quotient form.** If the blackboard is a fundamental chamber of a reflection group acting by isometries, then the quotient distance between two unlabelled, unoriented triangles is the spherical distance between their blackboard representatives. "Fold both into the blackboard and measure there" is the quotient-metric statement of Claims A and B. It would be unsurprising if Kendall, Small or Le stated it. I could not confirm or exclude this.
- Bowden et al.'s 48-tile picture shows that the wall set depends on the metric. In the Hausmann–Knutson metric, degenerate (collinear) triangles bound the domain through a different group action. Claim A is specific to the Procrustes (Kendall) metric, in which |s| is the quantity being maximised.

### Gaps
- Kendall 1984, Kendall 1985, Kendall & Kendall 1980 (flat triangles, standing stones), Kendall et al. 1999 and Small 1996 were **searched but not read**: no open-access text was reachable, and the Project Euclid fetch of the 1989 survey was empty. A statement equivalent to Claim A or B in these sources cannot be ruled out. The most likely places are Kendall 1985 §§ on the blackboard and the triangle chapter of Kendall et al. 1999. A library check is needed.
- Mardia, Edwards & Puri (1977, central-place theory) and Mardia (1989, "Shape analysis of triangles through directional techniques") handle unlabelled triangles by ordering angles. **Not searched** this session.
- Dryden & Mardia ch. 14: only the table of contents was seen, in the sibling note.

---

## Q2. Do astronomical triangle pattern-matching papers state the rank rule, and do they claim it is least-squares optimal?

### Takeaway
This is the strongest partial prior art. Groth (1986) labels triangle vertices by side rank, stores the orientation sense, and explicitly counts the same 12 candidates ("six ways to match the vertices and two ways to match the orientations"). He then uses the rank labelling to compare one canonical pairing per triangle pair, plus a same/opposite-sense flag. That is Claim B's rule as an algorithm. Pál & Bakos (2006) observe that an orientation-preserving, side-ordered triangle space is discontinuous at isosceles triangles, the invariant-space counterpart of Claim A. None of these papers states or needs least-squares optimality for non-similar triangles. The labelling is an invariant for matching nearly similar triangles within a tolerance.

### Cited Findings
- **Groth 1986, labelling (read, ADS full text).** Each triangle stores "the three points forming the vertices arranged so that the side between vertices 1 and 2 is the shortest side, the side between vertices 2 and 3 is the intermediate side, and the side between vertices 3 and 1 is the longest side". It also stores "the orientation—whether vertices 1,2, and 3 are traversed in a clockwise or counterclockwise sense", the ratio of longest to shortest side, and "the cosine of the angle at vertex 1" — [Groth, AJ 91:1244 (1986)](https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G)
- **Groth, the 12 candidates.** He counts "12n_t² potential matches (there are six ways to match the vertices and two ways to match the orientations of a pair of triangles)". He states: "Triangles may be matched in the same or opposite senses; that is, if both triangles of a matched pair have clockwise or both have counterclockwise orientations, the pair have the same sense, otherwise they have the opposite sense." — [Groth 1986](https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G)
- **Groth, matching criterion (not least squares).** Triangles match when (R_A − R_B)² < τ_RA² + τ_RB² and likewise for the cosines. "In the event that more than one triangle in list B matches a triangle in list A … the closest is chosen", meaning closest in (R, C) invariant space. The global transformation is a later "polynomial fit of data coordinates versus catalog coordinates". Triangles with ratio above a limit (typically 10) are discarded. The extracted text contains no discussion of isosceles or near-isosceles triangles — [Groth 1986](https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G)
- **Tabur 2007 (read, arXiv full text).** "Each set of 3 stars (triplet) may be matched in 6 different ways with a triplet from the other list. Using an optimization introduced by Groth (1986), the number of candidates is reduced by a factor of 6 by assigning the vertices of the triangle such that vertices A and B define the shortest side, B and C the longest side, and A and C define the intermediate length side." Also: "Valdes et al. (1995) used the ratios of two sides, (b/a, c/a) where a, b and c are the side lengths in decreasing order". Matches are declared within a 2% tolerance — [Tabur, PASA 24 (2007), arXiv:0710.3618](https://arxiv.org/abs/0710.3618)
- **Pál & Bakos 2006 (read via ar5iv, §III.1.1 "Triangle spaces").**
  - Their "mixed" space (attributed to Valdes et al. 1995) sorts the sides, so "the orientation information is lost", and "using mixed triangle space coordinates can yield false triangle matchings".
  - Their "chiral" space takes the sides "in counter-clockwise order" with the longest side first. "The shortcoming of T(chir) is that it is not continuous": "a small perturbation of an isosceles triangle can result in a new coordinate".
  - Their proposed "continuous" space is "a continuous function of the sides for all non-singular triangles, and also preserves chirality information".
  - The paper attributes triangle matching to Groth 1986, Stetson 1989 and Valdes et al. 1995.
  — [Pál & Bakos, PASP 118:1474 (2006), arXiv:astro-ph/0609658](https://arxiv.org/abs/astro-ph/0609658)
- **Astrometry.net (Lang et al. 2010)** hashes "asterisms (sets of four or five stars)", not triangles. Abstract only; how it breaks label symmetry inside a quad was not checked — [arXiv:0910.2233](https://arxiv.org/abs/0910.2233)
- **Valdes et al. 1995 (FOCAS) and Stetson (1989/1990)**: record only, through the citations in Tabur and Pál & Bakos. Not read.

### Inferences
- **Groth's procedure is Claim B's rule.** Matching vertex k of one triangle to vertex k of the other after side-rank labelling, with a same/opposite-sense flag, is exactly "pair by opposite-side rank, mirror iff orientations differ". So the rule as an algorithm dates from at least 1986, and the count "6 × 2 = 12" is in Groth's text.
- **What Groth does not claim.** For exactly similar triangles, any invariant labelling gives the correct correspondence, so optimality is trivial. Groth works in that regime (tolerance ε, near-similar triangles), and his closeness measure is in (R, C) invariant space, not a Procrustes residual. The content of Claim B beyond Groth is that the same labelling is the **exact least-squares Procrustes optimum for arbitrary, possibly very dissimilar triangles**, with or without scaling. That is a theorem about the residual, which none of these papers states.
- **Pál & Bakos is Claim A seen in invariant coordinates.** Their chiral triangle space is a coordinate chart on Kendall's 12-tile structure modulo labels. The discontinuity they report at isosceles triangles is the label flip at a wall. Their mixed space quotients the orientation out, which is the collinear wall. They treat it as a nuisance for hashing. They do not connect it to where a least-squares fit jumps, nor to the fact that only the target's crossing moves the fitted triangle.
- **Practical corollary for astrometry (my inference).** Claim B implies that rank-plus-sense labelling remains the Procrustes-optimal correspondence under arbitrarily large noise. When noise pushes a triangle across a wall, the rank labelling changes exactly when the least-squares optimum does, so the two cannot disagree. Ambiguity arises only from near-ties at walls.

### Gaps
- Star-tracker triangle algorithms (Liebe 1992; Cole & Crassidis 2006; Mortari's pyramid) were **not searched**. They typically use inter-star angular distances on the celestial sphere, and I would not expect an optimality claim, but this is unchecked.
- Valdes et al. 1995 and Stetson were not read directly. I relied on Tabur's and Pál & Bakos's descriptions.
- A direct full-text search of ADS for phrases like "least squares" together with "isosceles" in triangle-matching papers was not run. ADS full-text search needs an API token.

---

## Q3. Do molecular symmetry-corrected RMSD, bilateral-symmetry Procrustes, or computer-vision triangle methods contain the rule?

### Takeaway
I found no closed-form rank rule here. Molecular tools solve relabelling by enumeration, Hungarian assignment or graph isomorphism. The bilateral-symmetry literature (Mardia, Bookstein & Moreton 2000; Kent & Mardia 2001) uses the order-2 "relabel + reflect" map, a single mirror of the kind that D3h has three of. That is the same mechanism in its simplest form, applied to symmetry assessment, not to choosing a correspondence.

### Cited Findings
- **spyrmsd** (symmetry-corrected RMSD, J. Cheminform. 2020) offers two approaches.
  - It treats symmetry as a per-element linear assignment problem solved by the Hungarian method, which "can give artificially low values". In their docking example it reported 1.00 Å where the correct value was 2.46 Å.
  - It also offers graph isomorphism: "All possible graph isomorphisms are computed using the VF2 algorithm and the lowest RMSD among all isomorphisms is retained". Minimised RMSD uses QCP.
  — [spyrmsd paper](https://link.springer.com/doi/10.1186/s13321-020-00455-2); [spyrmsd docs](https://spyrmsd.readthedocs.io/en/latest/api/spyrmsd.rmsd.html)
- **Kent & Mardia (2001)**, "Shape, Procrustes tangent projections and bilateral symmetry", Biometrika 88(2):469–485. The abstract gives "an analytical description of bilateral symmetry of objects". Mardia, Bookstein & Moreton's earlier bilateral-symmetry paper is in Biometrika 87:285–300 (2000). Record only — [DOI 10.1093/biomet/88.2.469](https://doi.org/10.1093/biomet/88.2.469); [scilit record](https://scilit.com/publications/1b2e262fe83d107095aa0f6fa2a90968)
- **Rubix (Oct 2026)** solves "the equally weighted planar problem globally under squared Euclidean loss" via the "permutation polygon" of complex correlations, with "a sharp bound of n(n−1) vertices". Its abstract does not mention triangles, reflections or symmetric configurations (abstract read) — [arXiv:2610.10408](https://arxiv.org/abs/2610.10408)
- The sibling note `correspondence_soft_alignment.md` (Q5) covers molecular cases further, including Hungarian versus automorphism mappings and cyclic or dihedral symmetry. None there gives a sorting rule for three equivalent atoms.

### Inferences
- **3-D triangles.** For three points in R³, a proper rotation in SO(3) can flip the triangle's plane, so the planar mirror comes free. Superposing two 3-D triangles over all relabellings is the "mirrors allowed" case. By Claim B the optimal correspondence is opposite-side-rank pairing, with no orientation test. This applies to three interchangeable points superposed alone. In a molecule the other atoms fix the rotation, so the rule does not transfer directly to symmetry-corrected RMSD of, for example, methyl hydrogens.
- **The rule's place in the bilateral-symmetry picture.** In shape space, each "relabel + reflect" map (transposition ∘ mirror) is a reflection. A bilaterally symmetric shape is a point on its mirror. Kent & Mardia's setting has one such mirror. The triangle with mirrors allowed has three, plus the collinear equator, and these generate D3h.

### Gaps
- Computer-vision triangle invariants and geometric hashing (Lamdan & Wolfson; triangle-based affine or similarity hashing) were **not searched** beyond the astronomy papers. I would expect the same Groth-style canonical labelling without an optimality claim, but this is unchecked.
- Kent & Mardia 2001 and Mardia et al. 2000 were not read in full.

---

## Q4. Is the proposed mechanism (chamber lemma, nearest orbit point by folding) published, and has it been applied to nearest-orbit-point problems?

### Takeaway
Yes, the mechanism is textbook and is used for exactly this kind of nearest-orbit-point problem in coding theory. It appears in four places:
- Kuznetsov & Polishchuk state it as Corollary 2.4 for Weyl groups.
- Eaton triples take ⟨x, gy⟩ ≤ ⟨x, y⟩ on a cone as an axiom, and Eaton & Perlman (1977) extend majorization (the rearrangement inequality) to reflection groups.
- Slepian's permutation-modulation decoder (1965) is the type-A, sorting case.
- Reflection-group codes (Mittelholzer & Lahtonen 1996; Fossorier, Nation & Peterson 2010; Kim, Nation & Shepler) decode by folding the received vector into the fundamental chamber.

D3h = W(A₂) × W(A₁) is a Weyl group, so the lemma applies to the shape sphere without change.

### Cited Findings
- **Chamber lemma, verbatim (read).** "Corollary 2.4. For a pair of weights λ and µ the maximum (resp., minimum) of the scalar product (wλ, µ) when w runs through the Weyl group W is achieved when wλ and µ lie in the same Weyl chamber (resp., opposite Weyl chambers)." The proof reduces to "if λ is also dominant then (wλ, µ) ≤ (λ, µ) for any w ∈ W", via "Lemma 2.3. If λ is dominant then for any w ∈ W_G one has λ − wλ ∈ Q⁺_G" (non-negative combination of positive roots) — [Kuznetsov & Polishchuk, arXiv:1110.5607](https://arxiv.org/abs/1110.5607)
- **Eaton triples.** In an Eaton triple (V, G, F), G is a closed subgroup of the orthogonal group and F a closed convex cone with ⟨x, Ay⟩ ≤ ⟨x, y⟩ for all x, y ∈ F and A ∈ G. The orbit of each x meets F in exactly one point (search-engine extracts of these sources; not read in full) — [Results Math. 2023 article](https://link.springer.com/article/10.1007/s00025-023-01845-2); [Gowda, FTvN systems](https://userpages.umbc.edu/~gowda/papers/ftvn.pdf)
- **Eaton & Perlman (1977)**, Ann. Probab. 5(6):829–860. "The classical results of Schur (1923), Ostrowski (1952), Rado (1952), and Hardy, Littlewood and Polya (1952) concerning majorization and Schur functions are generalized to reflection groups" (abstract, record only) — [Project Euclid](https://www.projecteuclid.org/euclid.aop/1176995655)
- **Slepian's sorting decoder (type A).** Slepian's "Permutation modulation" (Proc. IEEE 53(3):228–236, 1965) decodes by sorting the received vector. For full symmetric-group permutation codes this simple decoder is maximum likelihood (search-engine extract citing Peterson 1997 and Biglieri, Karlof & Viterbo 1999; not read) — [Biglieri–Karlof–Viterbo PDF](https://ecse.monash.edu/staff/eviterbo/papers/itsep99r.pdf); [Nordio–Viterbo PDF](https://ecse.monash.edu/staff/eviterbo/papers/ict03.pdf)
- **Reflection-group codes (read, ar5iv full text).** Kim, Nation & Shepler define the fundamental region as "one vector closest to x₀ from each H-orbit". They then state: "If G is a Coxeter group, then the fundamental region is just a fundamental chamber that contains x₀", and "the decoding region of g in G is just the chamber containing g⁻¹x₀". Decoding folds the received vector: "d_m⋯d₁r is in the fundamental region of G and r decodes as d_m⋯d₁". They credit Mittelholzer & Lahtonen with solving the initial-vector problem for reflection groups, and Fossorier, Nation & Peterson with a coset-based decoder — [arXiv:1311.7038](https://arxiv.org/abs/1311.7038)
- **Mittelholzer & Lahtonen (1996)**, "Group codes generated by finite reflection groups", IEEE Trans. Inf. Theory 42(2). They give "fast maximum-likelihood decoding algorithms … relying on Slepian's decoding algorithms for permutation modulation" (search-engine extract; record only) — [Bohrium record](https://www.bohrium.com/paper-details/group-codes-generated-by-finite-reflection-groups/812461814884335618-2561)

### Inferences
- **The triangle problem is this setting exactly.** Map unit-size, centred labelled triangles to unit vectors x(z) ∈ R³ (Kendall's sphere, scaled to radius 1). The standard Hopf identity gives ⟨x(z), x(u)⟩ = 2|Σ z̄_i u_i|² − 1 = 2|s|² − 1, which is monotone in |s|. The page maximises |s| in every mode, and the size-and-shape residual Σ|z|² + Σ|u|² − 2|s| is also monotone in |s| because relabelling preserves Σ|u|². So the best of the 12 candidates is the orbit point of DEF nearest ABC under D3h acting on R³. D3h is a Weyl group W(A₂) × W(A₁): S₃ acts on the (x₁, x₂) plane as W(A₂), and the mirror acts as W(A₁) on the axis x₃ ∝ signed area. Kuznetsov–Polishchuk's proof uses only linearity and the positivity of dominant pairings, so it holds for real vectors. Claims A and B follow:
  - The Voronoi cells of DEF's orbit are the chambers (Kim–Nation–Shepler, Remark 5). The argmax therefore changes only when ABC crosses a wall, which gives Claim A's first half.
  - The selected point is DEF folded into ABC's chamber. Folding is continuous in DEF and only relabels when DEF crosses a wall, which gives Claim A's second half.
  - Reading the chamber as (side order, orientation), per Q1, gives Claim B.
  - In Slepian's terms, the page decodes a reflection-group code of length 12, with DEF as the initial vector and ABC as the received vector.
- **Why mirrors-off fails.** Without mirrors the group is D3 (three 120° turns and three half-turns), generated by rotations. It contains no reflections of the sphere, so there are no walls and no chambers. The Voronoi cells of a 6-point D3-orbit depend on DEF's position, which is the second half of Claim A, and ranking alone can be wrong. That fits the caller's report that ranking alone matched about 60% of the time.
- **A closed form for the gap across a wall.** For the reflection r_n(y) = y − 2⟨y,n⟩n in a wall with unit normal n, ⟨x, y⟩ − ⟨x, r_n y⟩ = 2⟨x,n⟩⟨y,n⟩. With the identity above, the two candidates on either side of a wall differ by Δ(|s|²) = ⟨x,n⟩⟨y,n⟩, the **product** of ABC's and folded DEF's (sine-)distances to that wall. So the jump vanishes if either triangle is on the wall, and the size of the fitted-triangle jump scales with DEF's distance from the wall. This matches the numerical jump sizes in `correspondence_soft_alignment.md`. I derived this; I did not find it in a source.

### Gaps
- Humphreys (1990), the canonical reference for chambers and fundamental domains of finite reflection groups, was **not fetched**. Kuznetsov–Polishchuk Cor. 2.4 serves as the verified written statement.
- Eaton & Perlman 1977, Mittelholzer & Lahtonen 1996 and Slepian 1965 were not read in full. Their content is from abstracts and secondary descriptions.
- I found no source in the Coxeter or coding literature that applies the lemma to Kendall's shape sphere or to Procrustes matching of landmark configurations.

---

## Q5. For each hit: optimality for the least-squares fit, or only an invariant labelling heuristic?

### Takeaway
No source found states either claim as a least-squares optimality result for triangle Procrustes. The pieces that exist fall into three kinds: invariant labelling heuristics (astronomy); quotient geometry used for probability densities or symmetry measures (Kendall, Bowden et al., Montgomery, Anderson); and general nearest-orbit theorems with no shape-analysis application (Coxeter, Eaton, coding).

### Cited Findings
| Source | What it contains | Least-squares optimality stated? | Relation to Claims A/B |
|---|---|---|---|
| [Groth 1986](https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G) (read) | Side-rank vertex labelling, orientation sense, "six ways … two ways", same/opposite-sense matches | No. Invariant matching within tolerance; global fit is a later polynomial fit | Claim B's rule as an algorithm, without the optimality theorem |
| [Tabur 2007](https://arxiv.org/abs/0710.3618) (read) | Groth labelling described as reducing candidates "by a factor of 6" | No | Same as Groth |
| [Pál & Bakos 2006](https://arxiv.org/abs/astro-ph/0609658) (read) | Chiral side-ordered triangle space "is not continuous" at isosceles triangles; mixed space loses orientation | No | Invariant-space version of the Claim A walls; nothing on fit jumps or target/source asymmetry |
| Valdes et al. 1995 (record via Tabur) | Sides sorted in decreasing order; orientation dropped | No | Heuristic labelling |
| [Lang et al. 2010](https://arxiv.org/abs/0910.2233) (abstract) | Quads/quints, not triangles | No | Not directly relevant |
| Kendall blackboard ([Kendall 1985](https://doi.org/10.2307/1427143), not read; [Behrend](https://cantab.net/users/michael.behrend/ley_stats/spher_black/pages/spher_black.html), extract; [Anderson](https://arxiv.org/abs/1711.10115)) | 12-tile structure; blackboard = unlabelled, unoriented triangles; isosceles/collinear boundaries | Not found in accessible text | Geometry underlying both claims; quotient-metric version plausible but unconfirmed |
| [Montgomery](https://arxiv.org/abs/1402.0841) (sibling note) | Isosceles great circles; transposition ∘ mirror = reflection in isosceles plane | No (dynamics, not matching) | Identifies D3h's mirrors |
| [Bowden et al. 2019](https://arxiv.org/abs/1708.01559) (read) | Fundamental domain bounded by isosceles and degenerate triangles; tile = ascending side lengths (Hausmann–Knutson metric, B₃) | No (least-symmetric triangle) | Side-order ↔ chamber dictionary, different metric |
| [Dryden & Mardia 2016 ch. 14](https://doi.org/10.1002/9781119072492.ch14) (TOC only) | Unlabelled matching by Bayesian/Procrustes models | Not seen | General framework, numerical |
| [spyrmsd](https://link.springer.com/doi/10.1186/s13321-020-00455-2) | Hungarian or graph-isomorphism enumeration | Optimal by enumeration only | No closed form |
| [Kent & Mardia 2001](https://doi.org/10.1093/biomet/88.2.469) (record) | Bilateral symmetry via relabel + reflect | Not applicable | One-mirror analogue |
| [Kuznetsov & Polishchuk](https://arxiv.org/abs/1110.5607) Cor. 2.4 (read) | max over w of (wλ, µ) attained in the same Weyl chamber | Yes, as an inner-product theorem | The mechanism, with no shape application |
| [Kim, Nation & Shepler](https://arxiv.org/abs/1311.7038) (read); Mittelholzer & Lahtonen 1996; Slepian 1965 | Nearest orbit point = fold into chamber; sorting decoder is ML | Yes, for decoding (Euclidean nearest codeword) | Same mathematics as Claim B; no shape application |
| [Rubix 2026](https://arxiv.org/abs/2610.10408) (abstract) | Global planar unlabelled Procrustes via the permutation polygon | Yes, for general n, by geometry and assignment queries | General-n algorithm; n = 3 rule not singled out |

### Inferences
- The only sources that state optimality are general theorems (Coxeter/Eaton/coding) or general-n algorithms (Rubix, enumeration). The triangle-specific sources (Groth, Pál & Bakos, Kendall, Bowden) state no optimality for matching.

### Gaps
- Classifications marked "record", "TOC" or "not read" are provisional.

---

## Q6. Verdict: novel, known, or known in pieces?

### Takeaway
**Known in pieces.** Each ingredient is published:
- the 12-tile geometry and the blackboard (Kendall, 1980s);
- the side-rank-plus-orientation labelling and the 6 × 2 = 12 count (Groth 1986);
- the discontinuity of a side-ordered chiral labelling at isosceles triangles (Pál & Bakos 2006);
- the chamber lemma and nearest-orbit-point decoding by folding (Coxeter and Weyl theory; Eaton & Perlman 1977; Slepian 1965; Mittelholzer & Lahtonen 1996).

What I did **not** find stated anywhere I could read is the assembled statement for triangle Procrustes:
- (i) rank + orientation is the exact least-squares optimum among the 12 candidates for arbitrary, non-similar triangles, with or without scaling;
- (ii) the hard fit jumps only when the target crosses an isosceles or collinear shape, independently of the source, while crossings of the source only relabel;
- (iii) both fail with mirrors off, because D3 is not a reflection group.

These follow in a few lines from published ingredients. They are best presented as a clean observation or corollary that connects Kendall's blackboard, Groth's labelling and reflection-group decoding, not as a new theorem. Because Kendall 1985, Kendall et al. 1999 and Small 1996 were not read, the quotient-metric form of (i) and (ii) may already be in print there.

### Cited Findings
- Groth's rule and 12-count — [Groth 1986](https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G)
- Isosceles discontinuity of chiral triangle space — [Pál & Bakos 2006](https://arxiv.org/abs/astro-ph/0609658)
- Chamber lemma — [Kuznetsov & Polishchuk, Cor. 2.4](https://arxiv.org/abs/1110.5607); decoding regions are chambers — [Kim, Nation & Shepler, Remark 5](https://arxiv.org/abs/1311.7038)
- Blackboard geometry — [Anderson](https://arxiv.org/abs/1711.10115); [Montgomery](https://arxiv.org/abs/1402.0841); [Kendall 1985 record](https://doi.org/10.2307/1427143)

### Inferences
- **Suggested wording.** "The optimal 12-candidate fit is the folding of DEF into ABC's Kendall chamber (a standard reflection-group fact). In side-length language this is Groth's (1986) labelling, which is therefore not only an invariant but the least-squares optimum." Avoid "we discovered"; "we observe" or "we note that this classical labelling is in fact optimal" is defensible.
- **The most defensible new pieces** are the target/source asymmetry of the jump locus, the explicit failure with mirrors off, and the product formula Δ(|s|²) = ⟨x,n⟩⟨y,n⟩ for the gap across a wall. I found these nowhere, but each is a short derivation.

### Gaps
- A library check of Kendall 1985 (§ on the blackboard), Kendall et al. 1999 (triangle chapters), Small 1996 (ch. on triangles) and Mardia, Edwards & Puri 1977 is needed before any claim to priority.
- Google Scholar forward-citation searches on Groth 1986 for "Procrustes" or "least squares optimal" were not run. Groth has hundreds of citing papers.

---

## Q7. What is the result useful for, and what open questions does it raise?

### Takeaway
Computationally it saves little, since 12 candidates are cheap. Its value is explanatory and diagnostic:
- it predicts where the drawn fit jumps and why only the target matters;
- it certifies Groth-type labelling as least-squares optimal under any noise;
- it gives canonical "blackboard" coordinates for unlabelled triangle statistics;
- it gives a closed form for the residual gap near a wall, which is what controls the soft fit.

The interesting open questions are what survives when the group stops being a reflection group: n ≥ 4, unequal masses, and non-isotropic weights.

### Cited Findings
- Rubix proves a sharp n(n−1) bound on vertices of the planar permutation polygon and solves the planar unlabelled problem globally. That is the natural general-n reference point — [arXiv:2610.10408](https://arxiv.org/abs/2610.10408)
- In a different (Hausmann–Knutson) metric, fundamental domains of the hyperoctahedral group on G₂(Rⁿ) parametrise unordered n-gons, and their walls are n-gons with a symmetry. The authors say they "don't know the general answer" to the furthest-from-walls question for n > 3 — [Bowden et al., arXiv:1708.01559](https://arxiv.org/abs/1708.01559)
- Montgomery's mirror description assumes equal masses. The poles are equilateral "only when the masses are equal" (sibling note) — [arXiv:1402.0841](https://arxiv.org/abs/1402.0841); see `shape_sphere_bridge.md`

### Inferences
- **Uses.**
  1. **Closed-form unlabelled triangle Procrustes.** Sort sides, compare signed areas, then do one complex Procrustes. It is exact, branch-light, and gives the same answer as enumeration.
  2. **A diagnostic for jumps.** A hard-fit jump means the target crossed a wall. Its size is proportional to the source's distance from that wall.
  3. **Astrometry and hashing.** The rule justifies Groth/Valdes labelling as residual-optimal, not just invariant. Near-isosceles or near-flat triangles are exactly where the labelling is fragile, which Pál & Bakos observed empirically.
  4. **Statistics on unlabelled triangles.** Folding all data into one chamber realises the quotient metric. Means or Fréchet means computed there are valid only if the data cloud stays away from walls; otherwise the quotient is an orbifold and folding biases the mean. This is a standing caution, not a sourced result.
  5. **3-D triangles.** Rotation-only superposition of 3-D triangles automatically includes the mirror. The optimal correspondence is plain opposite-side-rank pairing, and the walls are only the isosceles targets. Collinear triangles become the boundary of the hemisphere-shaped 3-D shape space, not a wall.
- **Open question 1: n ≥ 4.** Planar shape space is CP^{n−2}. Relabellings act holomorphically. Each "relabel ∘ mirror" involution is antiholomorphic, with a fixed set that is a real form of real dimension n−2, so codimension n−2 in a space of real dimension 2n−4. Only for n = 3 is this codimension 1, a wall. For n ≥ 4 the fixed sets of S_n × Z₂ do not separate the space, there are no chambers, and no folding rule of Claim-B type can exist in general. The Voronoi cells of an orbit depend on both configurations, as in the mirrors-off triangle case. Rubix's permutation polygon is the right replacement. A natural question is which special families recover a reflection-group structure. For collinear (1-D) configurations the rearrangement inequality gives sort-and-match with an optional flip, the type-A case. Configurations with a symmetry axis are another candidate.
- **Open question 2: unequal masses (Barbour).** With mass-weighted Jacobi coordinates, only permutations of equal-mass bodies are isometries.
  - With m₁ = m₂ ≠ m₃ the group is {id, (12)} × {proper, mirror} ≅ A₁ × A₁, still a reflection group. The rule should become: swap 1 and 2 according to the order of r₁₃ versus r₂₃, and mirror by orientation.
  - With all masses distinct, only the mirror remains. The rule reduces to "mirror iff the signed areas differ in sign", which holds by the wall-gap identity with n along the area axis.
  - I have not checked whether the r₁₃ = r₂₃ locus is still a great circle in mass-weighted coordinates. It should be, since (12) acts linearly, but this is unverified. The equilateral pole and the "rank" picture no longer coincide for unequal masses.
- **Open question 3: weighted or anisotropic Procrustes.** Unequal landmark weights or anisotropic noise break the S₃ isometry in the same way as unequal masses. Does a reflection subgroup survive, and what replaces "rank"?
- **Open question 4: soft and probabilistic fits near walls.** With weights e^{−Δrms²/σ²}, the weight ratio of the two candidates straddling a wall depends on the product ⟨x,n⟩⟨y,n⟩ (Q4). So the blend's transition width in ABC-space scales as σ² divided by DEF's distance from the wall. Exactly on a wall the 50/50 blend projects onto the wall: a flat fit for collinear targets, an isosceles fit for isosceles targets, as found numerically in `correspondence_soft_alignment.md`. Under isotropic Gaussian noise with a flat prior over the 12 labellings, the MAP labelling is the rank rule. The posterior over labellings is a softmax over 12 orbit points, which could be written in closed form with the wall-gap identity. Whether this matches Green–Mardia-type posteriors for n = 3 is open.
- **Open question 5: 3-D configurations with n ≥ 4.** SO(3) relabelling problems are not reflection-group problems either. Do any closed-form rules exist beyond n = 3?

### Gaps
- None of the open-question statements above was checked against literature this session. They are derivations or conjectures and should be labelled as such in any write-up.
- The great-circle status of partial isosceles loci for unequal masses was not verified.

---

## Source catalogue

"verified" = yes if I saw the bibliographic record or text during this session. "Sibling" means it was verified in a sibling note in this folder, not re-checked here. Fields marked "(memory)" were filled from memory and not checked. OA URLs are legal sources only.

| key | authors | year | title | venue/volume/pages | DOI | arXiv ID | open-access PDF URL | verified |
|---|---|---|---|---|---|---|---|---|
| groth1986pattern | Groth, Edward J. | 1986 | A pattern-matching algorithm for two-dimensional coordinate lists | Astronomical Journal 91(5):1244–1248 | 10.1086/114099 | — | https://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G | yes (full text read) |
| pal2006astrometry | Pál, András; Bakos, Gáspár Á. | 2006 | Astrometry in Wide-Field Surveys | PASP 118(848):1474–1483 | 10.1086/508573 | astro-ph/0609658 | https://arxiv.org/pdf/astro-ph/0609658 | yes (§III.1 read) |
| tabur2007fast | Tabur, Vladimir | 2007 | Fast algorithms for matching CCD images to a stellar catalogue | Publ. Astron. Soc. Australia 24:189–198 (pages memory) | 10.1071/AS07028 (memory) | 0710.3618 | https://arxiv.org/pdf/0710.3618 | yes (text read; DOI/pages no) |
| valdes1995focas | Valdes, F. G.; Campusano, L. E.; Velásquez, J. D.; Stetson, P. B. | 1995 | FOCAS Automatic Catalog Matching Algorithms (title memory) | PASP 107:1119 | 10.1086/133667 (memory) | — | — | no (cited in Tabur) |
| lang2010astrometry | Lang, Dustin; Hogg, David W.; Mierle, Keir; Blanton, Michael; Roweis, Sam | 2010 | Astrometry.net: Blind astrometric calibration of arbitrary astronomical images | Astronomical Journal 139:1782 | 10.1088/0004-6256/139/5/1782 (memory) | 0910.2233 | https://arxiv.org/pdf/0910.2233 | yes (abstract) |
| kendall1984shape | Kendall, David G. | 1984 | Shape manifolds, Procrustean metrics, and complex projective spaces | Bull. London Math. Soc. 16(2):81–121 (memory) | 10.1112/blms/16.2.81 (memory) | — | — | no (text not readable) |
| kendall1985exact | Kendall, David G. | 1985 | Exact distributions for shapes of random triangles in convex sets | Adv. Appl. Probab. 17(2):308–329 | 10.2307/1427143 | — | — | sibling |
| kendall1980alignments | Kendall, David G.; Kendall, Wilfrid S. | 1980 | Alignments in two-dimensional random sets of points | Adv. Appl. Probab. 12(2):380–424 (memory) | 10.2307/1426603 (memory) | — | — | no (not searched) |
| kendall1989survey | Kendall, David G. | 1989 | A Survey of the Statistical Theory of Shape | Statistical Science 4(2):87–99 | 10.1214/ss/1177012582 | — | https://projecteuclid.org/euclid.ss/1177012582 (OA landing page) | yes (record; text fetch empty) |
| kendall1999shape | Kendall, D. G.; Barden, D.; Carne, T. K.; Le, H. | 1999 | Shape and Shape Theory | Wiley, Chichester; ISBN 978-0-471-96823-8 | 10.1002/9780470317006 (memory) | — | — | yes (listing); content no |
| small1996statistical | Small, Christopher G. | 1996 | The Statistical Theory of Shape | Springer Series in Statistics | 10.1007/978-1-4612-4032-7 (memory) | — | — | no (not searched) |
| dryden2016unlabelled | Dryden, Ian L.; Mardia, Kanti V. | 2016 | Unlabelled size-and-shape and shape analysis (ch. 14) | Statistical Shape Analysis, with Applications in R, 2nd ed., Wiley, pp. 341–352 | 10.1002/9781119072492.ch14 | — | — | sibling (TOC only) |
| behrend_blackboard | Behrend, Michael | n.d. | Kendall's "spherical blackboard" (tutorial web pages) | cantab.net personal pages | — | — | https://cantab.net/users/michael.behrend/ley_stats/spher_black/pages/spher_black.html | no (HTTP 403; extract only) |
| montgomery2015shape | Montgomery, Richard | 2015 | The three-body problem and the shape sphere | Amer. Math. Monthly 122(4):299–321 (memory) | 10.4169/amer.math.monthly.122.04.299 (memory) | 1402.0841 | https://arxiv.org/pdf/1402.0841 | sibling |
| anderson2017smallest3 | Anderson, Edward | 2017 | The Smallest Shape Spaces III. Triangles in the Plane and in 3-d | arXiv preprint | — | 1711.10115 | https://arxiv.org/pdf/1711.10115 | sibling |
| anderson2017alice | Anderson, Edward | 2017 | Alice in Triangleland | arXiv preprint | — | 1711.11492 | https://arxiv.org/pdf/1711.11492 | sibling |
| bowden2019least | Bowden, Laney; Haynes, Andrea; Shonkwiler, Clayton; Shukert, Aaron | 2019 | Spherical geometry and the least symmetric triangle | Geometriae Dedicata 198:19–34 | 10.1007/s10711-018-0327-4 | 1708.01559 | https://arxiv.org/pdf/1708.01559 | yes (text read) |
| kuznetsov2016exceptional | Kuznetsov, Alexander; Polishchuk, Alexander | 2016 | Exceptional collections on isotropic Grassmannians | J. Eur. Math. Soc. (venue memory) | — | 1110.5607 | https://arxiv.org/pdf/1110.5607 | yes (arXiv text; Cor. 2.4 read) |
| humphreys1990reflection | Humphreys, James E. | 1990 | Reflection Groups and Coxeter Groups | Cambridge Univ. Press, Cambridge Studies in Adv. Math. 29 (memory) | 10.1017/CBO9780511623646 (memory) | — | — | no (not fetched) |
| eaton1977reflection | Eaton, Morris L.; Perlman, Michael D. | 1977 | Reflection Groups, Generalized Schur Functions, and the Geometry of Majorization | Annals of Probability 5(6):829–860 | 10.1214/aop/1176995655 | — | https://www.projecteuclid.org/euclid.aop/1176995655 (OA landing page) | yes (record) |
| slepian1965permutation | Slepian, David | 1965 | Permutation modulation | Proc. IEEE 53(3):228–236 | 10.1109/PROC.1965.3680 (memory) | — | — | no (secondary only) |
| mittelholzer1996group | Mittelholzer, Thomas; Lahtonen, Jyrki | 1996 | Group codes generated by finite reflection groups | IEEE Trans. Inf. Theory 42(2):519–528 (pages memory) | 10.1109/18.485721 (memory) | — | — | no (search record only) |
| fossorier2010reflection | Fossorier, Marc P. C.; Nation, J. B.; Peterson, W. W. | 2010 | Reflection group codes and their decoding | IEEE Trans. Inf. Theory 56(12) (memory) | — | — | — | no (search record only) |
| kim2013group | Kim, Hye Jung; Nation, J. B.; Shepler, Anne V. | 2013 | Group Coding with Complex Isometries | arXiv preprint (journal version not checked) | — | 1311.7038 | https://arxiv.org/pdf/1311.7038 | yes (text read) |
| biglieri1999permutation | Biglieri, E.; Karlof, J.; Viterbo, E. (authors per search extract) | 1999 | Representing group codes as permutation codes (title per search; unverified) | IEEE Trans. Inf. Theory (memory) | — | — | https://ecse.monash.edu/staff/eviterbo/papers/itsep99r.pdf | no |
| meli2020spyrmsd | Meli, Rocco; Biggin, Philip C. (authors memory) | 2020 | spyrmsd: symmetry-corrected RMSD calculations in Python | J. Cheminformatics 12:49 (memory) | 10.1186/s13321-020-00455-2 | — | https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00455-2 | yes (search record) |
| kent2001shape | Kent, John T.; Mardia, Kanti V. | 2001 | Shape, Procrustes tangent projections and bilateral symmetry | Biometrika 88(2):469–485 | 10.1093/biomet/88.2.469 | — | — | yes (record; authors via secondary) |
| mardia2000bilateral | Mardia, K. V.; Bookstein, F. L.; Moreton, I. J. | 2000 | Statistical assessment of bilateral symmetry of shapes | Biometrika 87(2):285–300 | 10.1093/biomet/87.2.285 (memory) | — | — | no (search mention only) |
| bhattacharjee2026rubix | Bhattacharjee, Subhransu S.; Campbell, Dylan; Shome, Rahul | 2026 | Rubix: Global Correspondence-Free Point Set Alignment through Assignment Geometry | arXiv preprint (7 Oct 2026) | — | 2610.10408 | https://arxiv.org/pdf/2610.10408 | yes (abstract) |
| mardia1977central | Mardia, K. V.; Edwards, R.; Puri, M. L. | 1977 | Analysis of central place theory | Bull. Int. Statist. Inst. 47:93–110 (memory) | — | — | — | no (not searched) |
