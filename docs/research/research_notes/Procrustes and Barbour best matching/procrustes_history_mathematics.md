# Procrustes analysis: history, mathematics, variants, and the link to Kendall's shape spaces

Conventions used throughout these notes (they follow Dryden & Mardia, and match this archive's `js/procrustes.js`):

- A planar configuration of k landmarks is a complex vector in C^k. `*` is the conjugate transpose, so z\*w = Σ z̄ᵢ wᵢ. `ᵀ` is the plain transpose, so wᵀz = Σ wᵢ zᵢ, with no conjugation.
- z is the **target** and w is the configuration **fitted onto** it. In the code these are `z` (ABC) and `u` (DEF). The code's sum `s = Σ z̄ᵢ u_perm(i)` is z\*u, its `phi` is the rotation θ, and its `rho` is the **scale** β. In Dryden & Mardia, ρ is the **Riemannian shape distance**. The two ρ's are different things.
- For real m-dimensional data, X is a k×m configuration matrix, H is the (k−1)×k sub-Helmert matrix, and Z = HX/‖HX‖ is the pre-shape.

---

## 1. Originators: who did what, with exact citations

### Takeaway
The problem started in factor analysis. Mosier (1939) fitted one factor matrix to a target. Green (1952) solved the orthogonal least-squares version for full-rank matrices. Hurley & Cattell (1962) named it "Procrustes". Schönemann (1966) gave the general SVD solution, and Schönemann & Carroll (1970) added translation and a central dilation. Gower (1975) generalised the method to many configurations (GPA), and Ten Berge (1977) corrected its scaling step and studied its optimality. The same mathematics was found again, independently, in satellite attitude (Wahba 1965), crystallography (Kabsch 1976/1978), metrology (Hanson & Norris 1981), and computer vision and photogrammetry (Arun et al. 1987; Horn 1987; Umeyama 1991). Kendall (1977, 1984) and the shape-statistics school (Bookstein, Goodall, Kent, Dryden & Mardia) then turned Procrustes distances into the metric geometry of shape spaces.

### Cited Findings

**Psychometric origins**
- **Mosier (1939)**, Psychometrika 4(2):149–162. The problem: "Given a primary trait matrix for n tests and r₁ traits, and a matrix for the same n tests and r₂ reference axes, to discover the transformation which will transform the second matrix into the first." Mosier gave "a rigorous and an approximate solution", applied them to hypothetical tests and to two batteries sharing 8 tests, and presented the method as a way of "comparing directly and quantitatively the results of two different factor studies." The transformation is unconstrained, which makes this the "oblique" Procrustes problem. — [Crossref record/abstract, DOI 10.1007/BF02288493](https://doi.org/10.1007/BF02288493)
- **Green (1952)**, Psychometrika 17(4):429–440. "A procedure is derived for obtaining an orthogonal transformation which most nearly transforms one given matrix into another given matrix, according to some least-squares criterion of fit." From it Green derived three analytic methods for approximating an oblique factor matrix by an orthogonal one. — [DOI 10.1007/BF02288918](https://doi.org/10.1007/BF02288918)
- **Hurley & Cattell (1962)**, "The Procrustes program: Producing direct rotation to test a hypothesized factor structure", Behavioral Science 7(2):258–262. — [DOI 10.1002/bs.3830070216](https://doi.org/10.1002/bs.3830070216). Crossref gives an online date of 2007 (a digitisation); volume 7 is 1962.
- Hurley & Cattell are credited with naming the method. One indexed abstract says it was "named Procrustes by Hurley & Cattell (1962)", and the R package EFA.dimensions documents that its *oblique* Procrustes rotation "is based on Hurley and Cattell (1962)" while its *orthogonal* rotation is based on Schönemann (1966). — [Kansai Univ. repository record](https://kansai-u.repo.nii.ac.jp/records/17221); [EFA.dimensions PROCRUSTES docs](https://rdrr.io/cran/EFA.dimensions/man/PROCRUSTES.html). Both are secondary sources.
- **Schönemann (1966)**, "A generalized solution of the orthogonal Procrustes problem", Psychometrika 31(1):1–10. "A solution T of the least-squares problem AT = B + E, given A and B so that trace(E′E) = minimum and T′T = I is presented. It is compared with a less general solution of the same problem which was given by Green. The present solution, in contrast to Green's, is applicable to matrices A and B which are of less than full column rank." — [DOI 10.1007/BF02289451](https://doi.org/10.1007/BF02289451)
- **Cliff (1966)**, "Orthogonal rotation to congruence", Psychometrika 31(1):33–42. Cliff treated two problems: rotating two factor solutions orthogonally "to a position where corresponding factors are as similar as possible", and rotating one factor matrix orthogonally to a specified target. — [DOI 10.1007/BF02289455](https://doi.org/10.1007/BF02289455)
- **Browne (1967)**, "On oblique Procrustes rotation", Psychometrika 32(2):125–132. — [DOI 10.1007/BF02289420](https://doi.org/10.1007/BF02289420). **Ten Berge & Nevels (1977)**, "A general solution to Mosier's oblique Procrustes problem", Psychometrika 42(4):593–600. — [DOI 10.1007/BF02295981](https://doi.org/10.1007/BF02295981)
- **Schönemann & Carroll (1970)**, "Fitting one matrix to another under choice of a central dilation and a rigid motion", Psychometrika 35(2):245–255: "A least squares method is presented for fitting a given matrix A to another given matrix B under choice of an unknown rotation, an unknown translation, and an unknown central dilation." The stated motivation was comparing nonmetric-scaling results across samples. This is the similarity (scaled) Procrustes problem. — [DOI 10.1007/BF02291266](https://doi.org/10.1007/BF02291266)
- **Gower (1975)**, "Generalized Procrustes analysis", Psychometrika 40(1):33–51. The problem is "translating, rotating, reflecting and scaling the m configurations to minimize the goodness-of-fit criterion Σᵢ Σⱼ Δ²(Pⱼ⁽ⁱ⁾, Gⱼ), where Gⱼ is the centroid". The centroid configuration acts as a "consensus". Results can be "summarized in analysis of variance form". "The special case m = 2 corresponds to Classical Procrustes analysis but the choice of criterion that fits each configuration to the common centroid configuration avoids difficulties that arise when one set is fitted to the other, regarded as fixed." — [DOI 10.1007/BF02291478](https://doi.org/10.1007/BF02291478)
- **Ten Berge (1977)**, "Orthogonal Procrustes rotation for two or more matrices", Psychometrika 42(2):267–276. He gave "a more straightforward proof" of a theorem by Fischer and Roppert that solves the two-matrix case. He showed that a sufficient condition for a best least-squares fit of more than two matrices "is ... not necessary", and that the "necessary conditions suggested by Kristof and Wingersky are shown to be not sufficient". He proposed an alternative rotation procedure, derived "upper bounds ... for determining the extent to which the procedure falls short of attaining the best least-squares fit", treated "scaling matrices to maximal agreement", and suggested "modifications of Gower's method of generalized Procrustes analysis." — [DOI 10.1007/BF02294053](https://doi.org/10.1007/BF02294053)
- **Sibson (1978)**, "Studies in the robustness of multidimensional scaling: Procrustes statistics", JRSS-B 40(2):234–238. — [DOI 10.1111/j.2517-6161.1978.tb01669.x](https://doi.org/10.1111/j.2517-6161.1978.tb01669.x)

**Engineering, crystallography and vision: the same problem found again**
- **Wahba (1965)**, "A least squares estimate of satellite attitude", SIAM Review 7(3):409. Posed as a SIAM problem. **Farrell, Stuelpnagel, Wessner, Velman & Brook (1966)** published solutions in SIAM Review 8(3):384–386. — [DOI 10.1137/1007077](https://doi.org/10.1137/1007077); [DOI 10.1137/1008080](https://doi.org/10.1137/1008080)
- **Kabsch (1976)**, "A solution for the best rotation to relate two sets of vectors", Acta Cryst. A32(5):922–923. **Kabsch (1978)**, "A discussion of the solution for the best rotation to relate two sets of vectors", Acta Cryst. A34(5):827–828. — [DOI 10.1107/S0567739476001873](https://doi.org/10.1107/S0567739476001873); [DOI 10.1107/S0567739478001680](https://doi.org/10.1107/S0567739478001680)
- Lawrence, Bernal & Witzgall (2019, NIST) describe the "Kabsch-Umeyama algorithm" as "the algorithm of choice" for the constrained (rotation-only) orthogonal Procrustes problem. It "is effectively no more than the computation of the singular value decomposition of a particular matrix". They note that "its justification, as presented separately by Kabsch and Umeyama, is not totally algebraic since it is based on solving the minimization problem via Lagrange multipliers." They cite Gower & Dijksterhuis (2004) and Schönemann (1966) for the unconstrained problem, and Kabsch 1976, Kabsch 1978 and Umeyama 1991 for the constrained one. — [Lawrence et al. 2019, J. Res. NIST 124:124028 (open access)](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **Hanson & Norris (1981)**, "Analysis of measurements based on the singular value decomposition", SIAM J. Sci. Stat. Comput. 2(3):363–373. Quality control: matching measured points to drawing points in a different coordinate system, assuming "a certain rigid transformation, an explicit least squares solution is obtained. This solution requires the singular value decomposition of a related matrix." — [DOI 10.1137/0902029](https://doi.org/10.1137/0902029)
- **Arun, Huang & Blostein (1987)**, "Least-squares fitting of two 3-D point sets", IEEE TPAMI PAMI-9(5):698–700. — [DOI 10.1109/TPAMI.1987.4767965](https://doi.org/10.1109/TPAMI.1987.4767965)
- **Horn (1987)**, "Closed-form solution of absolute orientation using unit quaternions", JOSA A 4(4):629–642. In his words: "The best translational offset is the difference between the centroid of the coordinates in one system and the rotated and scaled centroid of the coordinates in the other system. The best scale is equal to the ratio of the root-mean-square deviations of the coordinates in the two systems from their respective centroids ... The unit quaternion representing the best rotation is the eigenvector associated with the most positive eigenvalue of a symmetric 4 × 4 matrix." Horn credits Schut and Thompson with laying "the groundwork for the application of quaternions in photogrammetry". — [Horn 1987, author-hosted PDF](https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf); [DOI 10.1364/JOSAA.4.000629](https://doi.org/10.1364/JOSAA.4.000629)
- **Umeyama (1991)**, "Least-squares estimation of transformation parameters between two point patterns", IEEE TPAMI 13(4):376–380. His abstract states that the earlier solutions of Arun et al. and Horn et al. can return a reflection instead of a rotation "when the data is severely corrupted", and presents a solution that always gives correct parameters. — [DOI 10.1109/34.88573](https://doi.org/10.1109/34.88573); abstract seen via an [aggregator mirror](https://manuscript.isc.ac/Inventory/49/1408853.htm), not the IEEE page.

**Shape statistics**
- **D. G. Kendall (1977)**, "The diffusion of shape", Adv. Appl. Prob. 9(3):428–430. — [DOI 10.2307/1426091](https://doi.org/10.2307/1426091). Kendall's 1989 survey describes "the 'theory of shape' introduced by the author in 1977". — [Kendall 1989, DOI 10.1214/ss/1177012582](https://doi.org/10.1214/ss/1177012582)
- **Kendall (1984)**, "Shape manifolds, Procrustean metrics, and complex projective spaces", Bull. London Math. Soc. 16(2):81–121. — [DOI 10.1112/blms/16.2.81](https://doi.org/10.1112/blms/16.2.81). The term "pre-shape" "was coined by Kendall (1984)". — [Dryden, STAT718A notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Bookstein (1986)**, "Size and shape spaces for landmark data in two dimensions", Statistical Science 1(2):181–222. — [DOI 10.1214/ss/1177013696](https://doi.org/10.1214/ss/1177013696). **Bookstein (1991)**, *Morphometric Tools for Landmark Data*, CUP. — [DOI 10.1017/CBO9780511573064](https://doi.org/10.1017/CBO9780511573064). Crossref gives 1992 for this record; the book is conventionally cited as 1991, as in Dryden's notes.
- **Goodall (1991)**, "Procrustes methods in the statistical analysis of shape", JRSS-B 53(2):285–321 per Crossref. The paper presents "a model-based Procrustes approach to analysing sets of shapes". It develops understanding "through the familiar QR and singular value decompositions". It defines "a Gaussian model for landmark data" for one population, generalised to two-sample, ANOVA and regression models. "Maximum likelihood estimation is by least squares superimposition of the figures". It also generalises Procrustes "to allow non-isotropic errors at and between landmarks." — [DOI 10.1111/j.2517-6161.1991.tb01825.x](https://doi.org/10.1111/j.2517-6161.1991.tb01825.x)
- **Kent (1994)**, "The complex Bingham distribution and shape analysis", JRSS-B 56(2):285–299. The distribution is defined on the complex unit sphere in C^k, and its "symmetries ... make it a natural tool for the analysis of the shape of landmark data in two dimensions". The paper also discusses links with Bookstein coordinates. — [DOI 10.1111/j.2517-6161.1994.tb01978.x](https://doi.org/10.1111/j.2517-6161.1994.tb01978.x)
- **Kendall, Barden, Carne & Le (1999)**, *Shape and Shape Theory*, Wiley. — [DOI 10.1002/9780470317006](https://doi.org/10.1002/9780470317006)
- **Gower & Dijksterhuis (2004)**, *Procrustes Problems*, Oxford Statistical Science Series, OUP, ISBN 978-0-19-851058-1. — [DOI 10.1093/acprof:oso/9780198510581.001.0001](https://doi.org/10.1093/acprof:oso/9780198510581.001.0001)
- **Dryden & Mardia**, *Statistical Shape Analysis* (1st ed., Wiley 1998, xvii+347 pp., ISBN 0-471-95816-6), per [Rudemo's review in Stat. Med., DOI 10.1002/1097-0258(20001015)19:19<2716::aid-sim590>3.0.co;2-o](https://doi.org/10.1002/1097-0258(20001015)19:19%3C2716::aid-sim590%3E3.0.co;2-o). The 2nd edition is *Statistical Shape Analysis, with Applications in R* (Wiley 2016, ISBN 978-0-470-69962-1). — [DOI 10.1002/9781119072492](https://doi.org/10.1002/9781119072492)

### Inferences
- The algorithm was discovered at least four separate times: psychometrics (Green/Schönemann), attitude estimation (Wahba), crystallography (Kabsch), and metrology and vision (Hanson & Norris, Arun, Horn, Umeyama). Each community kept its own name: "orthogonal Procrustes", "Wahba's problem", "Kabsch algorithm", "absolute orientation". The psychometric line allowed any orthogonal T (T′T = I, so reflections are allowed). The physical-science lines needed proper rotations, which is why the det = +1 correction carries Kabsch's and Umeyama's names.
- Gower's 1975 change was more than going from two matrices to many. Fitting every configuration to a shared centroid makes the criterion symmetric in the configurations, and that symmetry is what later let Procrustes sums of squares serve as shape distances.

### Gaps
- I could not see the full text of Hurley & Cattell 1962. The naming claim rests on secondary sources, and the cited abstract and R documentation agree on it. The myth itself (Procrustes stretching or cutting his guests to fit his bed) is the usual etymology; I did not trace it to a primary citation.
- I could not see the text of Kabsch 1976 or 1978. That the 1978 note handles the improper-rotation (det = −1) case comes from secondary sources: [MathWorld "Kabsch Algorithm"](https://mathworld.wolfram.com/KabschAlgorithm.html) and the [BTK RMSD theory notes](https://btk.sourceforge.net/html_docs/0.8.1/rmsd_theory.html). Treat the exact 1978 wording as unverified.
- I did not see Kristof & Wingersky (1971, Proc. 79th Annual Convention APA 6:89–90) or Fischer & Roppert. They are known here only as they are cited in Ten Berge's 1977 abstract.
- I did not open the content of Gower & Dijksterhuis (2004) or Kendall et al. (1999). Only their bibliographic records were verified.

---

## 2. Solving the orthogonal Procrustes problem by SVD, the det = +1 constraint, and how scale, weights and translation enter

### Takeaway
Remove translation by centring. The rotation then maximises a trace, tr(ΓᵀA) with A the cross-covariance. If A = VΛUᵀ, the unconstrained orthogonal optimum is Γ = UVᵀ (in the right orientation of the product). To force det = +1, multiply through by diag(1,…,1, det(UVᵀ)) (the Kabsch–Umeyama correction). This flips the sign of the smallest singular direction, and the fit pays for it with 2σ_min of the trace. The optimal scale is the trace of the corrected singular values divided by the source sum of squares. Weights that act per landmark keep the closed form. Non-isotropic weights generally do not.

### Cited Findings
- **Reduction to a trace maximisation.** For d × n matrices P, Q, ‖UQ − P‖²_F = tr(QᵀQ) + tr(PᵀP) − 2 tr(PᵀUQ). It is therefore enough to find a rotation U maximising tr(UM) with M = QPᵀ. "It is well known that one such U can be computed from the singular value decomposition of M." — [Lawrence et al. 2019](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **The Kabsch–Umeyama algorithm, verbatim from Lawrence et al.** "Compute d × d matrix M = QPᵀ. Compute SVD of M, i.e., identify d × d matrices V, S, W, so that M = VSWᵀ ... Set s₁ = … = s_{d−1} = 1. If det(VW) > 0, then set s_d = 1, else set s_d = −1. Set S̃ = diag{s₁, …, s_d}. Return d × d rotation matrix U = WS̃Vᵀ." — [Lawrence et al. 2019](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **Why the correction is optimal.** Proposition 1 of Lawrence et al.: if D = diag(σ₁,…,σ_d) with σⱼ ≥ 0 and W is orthogonal, then (1) tr(WD) ≤ Σσⱼ, and (3) if det(W) = −1 and σ_d ≤ σⱼ for all j, then tr(WD) ≤ Σ_{j<d} σⱼ − σ_d. The authors note that a justification "that requires the conclusion in (3) but lacks a proof for it, is not exactly complete." — [Lawrence et al. 2019](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **Rigid motion reduces to the centred problem.** Translating both sets to their centroids reduces "the orientation-preserving rigid motion problem" to the constrained orthogonal Procrustes problem. — [Lawrence et al. 2019](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **Similarity version (Umeyama), as implemented in Eigen with Umeyama's equation numbers.** Minimise (1/n)Σ‖yᵢ − (cRxᵢ + t)‖². Then:
  - Eq. (38): Σ_xy = (1/n)Σ(yᵢ − μ_y)(xᵢ − μ_x)ᵀ, with SVD Σ_xy = UDVᵀ.
  - Eq. (39): S = I, except S_mm = −1 if det Σ_xy < 0.
  - Eqs. (40)/(43): R = USVᵀ.
  - Eq. (42): c = (1/σ_x²) tr(DS), where σ_x² is the source variance.
  - Eq. (41): t = μ_y − cRμ_x.

  — [Eigen `Umeyama.h`](https://eigen.googlesource.com/mirror/+/17af8c763d5aa58c40fc2c40ff6e8d6fc28a8d63/Eigen/src/Geometry/Umeyama.h). MathWorld states the same correction as D = diag(1, …, 1, det(UVᵀ)), citing Kabsch 1976/1978 and Lawrence et al. 2019. — [MathWorld](https://mathworld.wolfram.com/KabschAlgorithm.html)
- **The shape-analysis form (Dryden & Mardia, lecture-note version).** For centred X₁, X₂ (k × m), minimise D²_OPA = ‖X₂ − βX₁Γ − 1_kγᵀ‖² over Γ ∈ SO(m), β > 0, γ ∈ R^m:
  - γ̂ = 0.
  - With X₂ᵀX₁ = ‖X₁‖‖X₂‖ VΛUᵀ, U, V ∈ SO(m), the rotation is Γ̂ = UVᵀ, where λ₁ ≥ … ≥ λ_{m−1} ≥ |λ_m| and "λ_m is the negative square root iff det(Z₁ᵀZ₂) < 0."
  - β̂ = tr(X₂ᵀX₁Γ̂)/tr(X₁ᵀX₁) = (‖X₂‖/‖X₁‖) cos ρ(X₁, X₂).
  - The minimised sum is OSS(X₁, X₂) = ‖X₂‖² sin²ρ(X₁, X₂).

  — [Dryden STAT718A notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf), [ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **When reflection matters.** "λ_m will be negative in the cases where an orthogonal transformation (reflection and rotation) would produce a smaller sum of squares than just a rotation. In practice, for fairly close shapes λ_m will usually be positive – in which case the solution is the same as minimizing over the orthogonal matrices O(m) instead of SO(m)." — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **OPA is asymmetric.** Reversing roles gives γ̂ᴿ = −γ̂ and Γ̂ᴿ = Γ̂ᵀ, but β̂ᴿ ≠ 1/β̂ in general, and OSS(X₂, X₁) ≠ OSS(X₁, X₂) unless the figures have the same size. Normalising both to unit size gives √OSS = sin ρ = d_F, a proper distance. — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **Scale: least squares versus symmetric.** Horn shows that the asymmetric least-squares scale depends on which set is transformed. He prefers a symmetric scale, equal to "the ratio of the root-mean-square deviations of the coordinates in the two systems from their respective centroids", so that the inverse transformation gives the inverse scale. — [Horn 1987](https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf)
- **Quaternion alternative (3D).** The optimal rotation quaternion is the eigenvector for the most positive eigenvalue of a symmetric 4 × 4 matrix N built from sums of products of corresponding coordinates. — [Horn 1987](https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf). Coutsias, Seok & Dill (2004) use the same idea for RMSD. — [DOI 10.1002/jcc.20110](https://doi.org/10.1002/jcc.20110)
- **The reflection failure Umeyama fixed.** Arun et al. and Horn et al. could return a reflection "when the data is severely corrupted". Umeyama's correction keys on the sign of det(Σ_xy). One secondary note warns that ad hoc fixes, such as flipping a column of U after the fact, do not always give the least-squares answer. — [aggregator abstract](https://manuscript.isc.ac/Inventory/49/1408853.htm); [Heidenreich note](https://hunterheidenreich.pages.dev/notes/biology/computational-biology/umeyama-similarity-transformation/), a secondary blog source.
- **Non-isotropic weighting.** Goodall (1991) generalised Procrustes "to allow non-isotropic errors at and between landmarks", and observed that "the superimposition metric used for fitting, and the model metric, or covariance, used for testing, may not coincide." — [Goodall 1991](https://doi.org/10.1111/j.2517-6161.1991.tb01825.x)

### Inferences
- **Notation.** Sources write the same algorithm with different SVD labels. Lawrence writes M = QPᵀ = VSWᵀ and U = WS̃Vᵀ. Umeyama writes Σ_xy = UDVᵀ and R = USVᵀ. Dryden writes X₂ᵀX₁ ∝ VΛUᵀ and Γ̂ = UVᵀ (right-multiplication, since configurations are k × m rows). In every case the rule is the same: rotation = (left/right factors of the cross-covariance arranged so the trace is maximised) × diag(1,…,1, sign det). The diag factor should be computed from the determinant of the cross-covariance, or of the product of the singular-vector matrices, and never by flipping signs afterwards.
- **Translation always separates.** For isotropic, unweighted least squares the optimum puts the transformed source centroid on the target centroid: t = μ_y − cRμ_x (Umeyama eq. 41; Horn's abstract). Centring first is therefore without loss of generality.
- **Per-landmark weights wᵢ keep the closed form.** Use weighted centroids and the weighted cross-covariance Σwᵢ(yᵢ − μ_y^w)(xᵢ − μ_x^w)ᵀ, then the same SVD. Weights per coordinate axis, or a full landmark covariance as in Goodall's non-isotropic models, break the trace structure and in general need iteration. I did not verify this claim against a primary weighted-Procrustes paper.
- **The planar case needs no SVD or determinant test** (see §4). Complex multiplication by βe^{iθ} is always a proper similarity. "Allowing reflection" just means also trying the conjugate w̄ and keeping the better fit.

### Gaps
- I did not see the original Wahba (1965) problem text. The usual statement is: find the orthogonal matrix with determinant +1 that best brings one set of unit vectors into least-squares coincidence with another. The weighted form L(A) = ½Σaᵢ‖bᵢ − Arᵢ‖² is the later standard formulation. I did not verify that it appears in the 1965 note.
- I could not directly confirm the rank-deficient branch of Umeyama (eq. 43: when rank Σ_xy = m − 1, choose S by det(U)det(V)). Eigen's code has a separate rank-deficient branch citing "(40) and (43)", and its details here are from memory.
- Primary sources for weighted orthogonal Procrustes in psychometrics (for example Lissitz, Schönemann & Lingoes 1976; Koschat & Swayne 1991) were not retrieved.

---

## 3. Variants: ordinary, generalized (GPA), partial, full, weighted, extended; the GPA iteration and its convergence

### Takeaway
The adjectives sit on two independent axes:

- **How many configurations are fitted.** Ordinary (OPA) fits one onto another; generalized (GPA) fits n to a common mean.
- **Which transformations are allowed.** Partial uses translation and rotation; full adds isotropic scale; reflection Procrustes uses O(m) instead of SO(m). Weighted and oblique are further variants, and "extended" usually means rotation + translation + scale.

GPA alternates OPA steps against a leave-one-out mean, with an optional Ten Berge rescaling step. The loss never increases, but the algorithm is not guaranteed to reach a global optimum. For planar data the full-Procrustes mean has a closed form (a complex eigenvector), so no iteration is needed.

### Cited Findings
- **OPA versus GPA.** "When several objects are fitted using Procrustes superimposition the method has been called generalized Procrustes analysis (GPA) (Gower, 1975), whereas when a single object is fitted to one other ... ordinary Procrustes analysis (OPA). Note that OPA is not symmetrical in the ordering of the objects, whereas GPA is invariant under re-orderings of the objects." — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf)
- **Full, partial, reflection.** "Partial Procrustes involves just superimposition by translation and rotation (not scaling) as opposed to full Procrustes which involves the full set of similarity transformations. Another alternative could be matching using orthogonal matrices instead of rotation matrices (reflection Procrustes)." In the R package `shapes`, `procGPA(x, scale=TRUE)` "carries out full Procrustes analysis (the default)", `scale=FALSE` gives partial Procrustes, and "The option reflection=TRUE uses reflection invariance whereas the default value does not have reflection invariance." — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **"Full" means the full similarity group.** "The term full is used because the full set of Euclidean similarity transformations is estimated in the matching (translation, rotation and scale), but note that y and w are pre-scaled to unit size." Also, "the optimal rotation is the same whether or not scaling is in the minimization." — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf), [ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **The GPA algorithm, as written in Dryden's notes:**
  1. *Translations.* Centre each configuration; set Xᵢᴾ = Xᵢ.
  2. *Rotations.* For each i let X̄₍ᵢ₎ = (1/(n−1))Σ_{j≠i} Xⱼᴾ, and replace Xᵢᴾ by its OPA rotation-only fit onto X̄₍ᵢ₎. Rotate the n figures in turn, and repeat "until the Procrustes sum of squares ... cannot be reduced further".
  3. *Scaling.* "Let Φ be the n × n correlation matrix of the vec(Xᵢᴾ) ... with eigenvector φ ... corresponding to the largest eigenvalue. Then from Ten Berge (1977) take β̂ᵢ = (Σₖ‖Xₖᴾ‖² / ‖Xᵢᴾ‖²)^{1/2} φᵢ."
  4. Repeat steps 2–3 until the Procrustes sum of squares cannot be reduced further. "The algorithm usually converges quickly."

  For m = 3 and higher "the full Procrustes mean shape has to be found iteratively ... but an explicit eigenvector solution is available for two dimensional data." — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **Optimality of multi-matrix rotation.** Ten Berge (1977) showed that the Kristof–Wingersky necessary conditions are not sufficient for a best fit with more than two matrices. He also gave upper bounds on how far a rotation procedure can fall short of the best fit. — [Ten Berge 1977](https://doi.org/10.1007/BF02294053)
- **Gower's design and output.** Each configuration is fitted to the common centroid, which "avoids difficulties that arise when one set is fitted to the other, regarded as fixed". Results "can be summarized in analysis of variance form". — [Gower 1975](https://doi.org/10.1007/BF02291478)
- **Convergence and consistency literature (records verified):**
  - Groisser (2005), "On the convergence of some Procrustean averaging algorithms", Stochastics 77(1):31–60. — [DOI 10.1080/17442500512331341059](https://doi.org/10.1080/17442500512331341059)
  - Kent & Mardia (1997): "Lele has shown that the Procrustes estimator of form is inconsistent". The Procrustes estimator of *shape* "is consistent under the assumption of an isotropic error distribution and ... consistency breaks down if the assumption of isotropy is relaxed." — [DOI 10.1111/1467-9868.00069](https://doi.org/10.1111/1467-9868.00069)
  - Le (1998): under isotropy, "an increasing sequence of procrustean mean shapes defined in either of the two ways used in practice will tend almost surely to the shape of the means." — [DOI 10.1239/aap/1035227991](https://doi.org/10.1239/aap/1035227991)
- **Extended (similarity) Procrustes.** Schönemann & Carroll (1970) solve rotation + translation + "central dilation". — [DOI 10.1007/BF02291266](https://doi.org/10.1007/BF02291266)
- **Oblique Procrustes** (the transformation is not orthogonal): Mosier 1939; Hurley & Cattell 1962; Browne 1967; Ten Berge & Nevels 1977, "A general solution to Mosier's oblique Procrustes problem". — [Browne](https://doi.org/10.1007/BF02289420); [Ten Berge & Nevels](https://doi.org/10.1007/BF02295981)
- **Robust and resistant variants (records verified, content not seen):** Siegel & Benson (1982), "A robust comparison of biological shapes", Biometrics 38(2):341 ff. — [DOI 10.2307/2530448](https://doi.org/10.2307/2530448). Rohlf & Slice (1990), "Extensions of the Procrustes method for the optimal superimposition of landmarks", Syst. Zool. 39(1):40 ff. — [DOI 10.2307/2992207](https://doi.org/10.2307/2992207)

### Inferences
- **GPA cannot increase the loss.** Each rotation step replaces Xᵢ by its exact least-squares fit to the mean of the others, which cannot increase the total GPA sum of squares. Each scaling step solves its subproblem exactly. The loss is bounded below by 0, so the sequence of loss values converges. Whether the iterates converge to a global minimiser is a separate matter: Ten Berge's (1977) results show that stationary points need not be global optima for n > 2. Groisser (2005) is the rigorous reference, but I did not see its content.
- **Ten Berge's scaling preserves total size.** With ‖φ‖ = 1, Σᵢ β̂ᵢ²‖Xᵢᴾ‖² = Σₖ‖Xₖᴾ‖² · Σᵢφᵢ² = Σₖ‖Xₖᴾ‖². This constraint stops GPA from shrinking every configuration towards zero, which would otherwise trivially minimise the sum of squares.
- **In 2D, full GPA is one eigenproblem.** The full Procrustes mean is the leading eigenvector of S = Σ zᵢzᵢ\* (Kent 1994; see §4), and each fit is then a single complex multiplication. The leave-one-out iteration is only needed for m ≥ 3, or when a different loss is used (partial, weighted, resistant).
- **Terminology varies.** "Extended orthogonal Procrustes" is used, especially in geodesy and photogrammetry, for the 7-parameter similarity fit (rotation + translation + one scale), which is Schönemann & Carroll's problem. I did not retrieve a primary source fixing the term, so the label should be used with care.

### Gaps
- I did not retrieve the abstract or content of Groisser (2005). Its precise theorem (local convergence conditions for Procrustean averaging) is not stated here.
- I did not verify the content of Gower & Dijksterhuis's (2004) chapter taxonomy (two-sided, projection, oblique, weighted, missing-value Procrustes).

---

## 4. The planar case with complex numbers: optimal rotation and scale, and the full, partial and Riemannian distances

### Takeaway
For centred z, w ∈ C^k, the optimal similarity fit of w onto z is a complex regression with no intercept: β̂e^{iθ̂} = w\*z / w\*w. The residual gives the full Procrustes distance, d_F² = 1 − |w\*z|²/(‖w‖²‖z‖²). Write cos ρ = |w\*z|/(‖w‖‖z‖). Then d_F = sin ρ, d_P = 2 sin(ρ/2), and ρ = arccos(·) is the Riemannian distance in Kendall's shape space CP^{k−2}. This archive's code already implements these formulas exactly.

### Cited Findings
- **The model** (Dryden & Mardia, planar chapter): y = (a+ib)1_k + βe^{iθ}w + ε, a complex regression with design X_D = [1_k, w] and parameters A = (a+ib, βe^{iθ})ᵀ. The full Procrustes fit minimises D²(y, w) = ‖y − wβe^{iθ} − (a+ib)1_k‖². — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf)
- **Result 3.1** (y and w centred). The fit has â + ib̂ = 0 (3.2), θ̂ = arg(w\*y) = −arg(y\*w) (3.3), and β̂ = (w\*y y\*w)^{1/2}/(w\*w) = |y\*w|/(w\*w) (3.4). Proof: write y\*w = γe^{iφ}; the cross term is 2βγ cos(θ + φ), maximised at θ = −φ; then ∂D²/∂β = 0 gives β̂ = γ/(w\*w). — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf)
- **The fit as ordinary least squares.** Â = (X_D\*X_D)⁻¹X_D\*y gives Â₁ = 0 and Â₂ = w\*y/(w\*w) (3.6). The full Procrustes fit is wᴾ = β̂e^{iθ̂}w = (w\*y) w/(w\*w). The residual sum of squares is r\*r = y\*y − (y\*w w\*y)/(w\*w) (3.7). "The full Procrustes fit of w onto y is actually obtained by complex linear regression of y on w." — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf)
- **Full Procrustes distance** (Definition 3.2, eq. 3.8): d_F(w, y) = inf_{β,θ,a,b} ‖y/‖y‖ − (w/‖w‖)βe^{iθ} − 1_k(a+ib)‖ = (1 − (y\*w w\*y)/(w\*w y\*y))^{1/2}. — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf)
- **The three distances**, from Dryden & Mardia, Table 1 and eqs. 4.11–4.14, general m; here λᵢ are the signed singular values of Z₁ᵀZ₂ for pre-shapes:

  | Distance | Formula | Range |
  |---|---|---|
  | Full Procrustes, d_F | {1 − (Σλᵢ)²}^{1/2} | 0 ≤ d_F ≤ 1 |
  | Partial Procrustes, d_P | √2 (1 − Σλᵢ)^{1/2} | 0 ≤ d_P ≤ √2 |
  | Procrustes (Riemannian), ρ | arccos(Σλᵢ) | 0 ≤ ρ ≤ π/2 |

  ρ is "the closest great circle distance between Z₁ and Z₂ on the pre-shape sphere", minimising over rotations; ρ = 2 arcsin(d_P/2). "Indeed d_F(X₁,X₂) = sin ρ, d_P(X₁,X₂) = 2 sin(ρ/2)." For close shapes, d_P = d_F + O(d_F³) and ρ = d_F + O(d_F³). "ρ can be considered as the smallest angle between the complex vectors Z₁ and Z₂ over rotations of Z₁ and Z₂." — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Partial distance** (Definition 4.13): d_P = inf_{Γ∈SO(m)} ‖Z₂ − Z₁Γ‖ on pre-shapes, which is "keeping β = 1 fixed". — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Planar full Procrustes mean** (Kent 1994, Result 3.2). [µ̂] is the eigenvector for the largest eigenvalue of the complex sum-of-squares-and-products matrix S = Σᵢ wᵢwᵢ\*/(wᵢ\*wᵢ) = Σ zᵢzᵢ\*, with zᵢ = wᵢ/‖wᵢ‖. Proof: Σd_F²(wᵢ, µ) = n − µ\*Sµ/(µ\*µ). The full Procrustes coordinates are wᵢᴾ = (wᵢ\*µ̂) wᵢ/(wᵢ\*wᵢ) (3.12). The solution coincides with the maximum-likelihood mode under the complex Bingham model. — [Dryden notes ch. 3](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf); [Kent 1994](https://doi.org/10.1111/j.2517-6161.1994.tb01978.x)
- **Similarity group in complex form.** The Euclidean similarity transformations of z° are {ηz° + 1_kξ : η = βe^{iθ} ∈ C, ξ ∈ C}. — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Software.** `shapes::procdist(x, y, type=…)` returns "full", "partial" or "Riemannian" (and size-and-shape variants), and accepts 2D configurations. — [CRAN `procdist` reference](https://search.r-project.org/CRAN/refmans/shapes/html/procdist.html)
- **Canonical book references:** Kendall (1984), [DOI 10.1112/blms/16.2.81](https://doi.org/10.1112/blms/16.2.81); Dryden & Mardia (2016), [DOI 10.1002/9781119072492](https://doi.org/10.1002/9781119072492); Mardia & Dryden (1989), "The statistical analysis of shape data", Biometrika 76(2):271–281, [DOI 10.1093/biomet/76.2.271](https://doi.org/10.1093/biomet/76.2.271); Kent (1994).

### Inferences
These were derived by me and checked numerically on 2,000 random configurations for each identity (script run this session).

- **Formulas in the question's notation** (target z, fitted w, both centred):
  - β̂ e^{iθ̂} = w\*z / (w\*w), so β̂ = |w\*z|/‖w‖² and θ̂ = arg(w\*z) = −arg(z\*w).
  - min‖z − βe^{iθ}w‖² = ‖z‖² − |w\*z|²/‖w‖² = ‖z‖² d_F².
  - cos ρ = |w\*z|/(‖w‖‖z‖).
  - d_F² = 1 − cos²ρ = sin²ρ.
  - d_P² = 2(1 − cos ρ) = 4 sin²(ρ/2).
  - Uncentred data: ξ̂ = z̄_c − β̂e^{iθ̂} w̄_c, where z̄_c and w̄_c are the centroids.
- **Mapping to `js/procrustes.js`.** `s = Σ conj(z_i)·u_perm(i)` is z\*u. `phi = −arg(s)` is θ̂. `rho = |s|/Σ|u|²` (when `scaling`) is β̂. So `candidate()` is exactly Dryden & Mardia's Result 3.1 applied to each vertex pairing.
  - The comment in the code is correct that GeoGebra's θ = atan(y/x) is wrong when Re s < 0. The derivation needs cos(θ + φ) = +1. An angle determined only modulo π can land on cos = −1, which is the *worst* rotation. A full-quadrant `atan2` (arg) is required.
- **Scale versus distance.** The code's `rho` (scale) and Dryden's ρ (distance) are related by β̂ = (‖z‖/‖w‖) cos ρ, the 2D case of Dryden's eq. (4.7).
- **The 2D singular values are the two complex inner products.** For the real 2×2 cross-product A = Zᵀ W, with configurations as k×2 real matrices:
  - |w\*z| = σ₁ + sign(det A) σ₂ and |wᵀz| = σ₁ − sign(det A) σ₂.
  - det A = (|w\*z|² − |wᵀz|²)/4.
  - So the Kabsch correction (det A < 0) fires exactly when the mirrored fit beats the proper one, |wᵀz| > |w\*z|. That is the complex-number form of Dryden's "λ_m negative iff det < 0".
- **The code's shape-sphere identity is cos²ρ.** In `js/procrustes.js`, |s|²/(‖z‖²‖u‖²) = (1 + n·n′)/2. The left side is cos²ρ. So n·n′ = cos 2ρ, and on the unit sphere the angle between the shape points is 2ρ. Equivalently, ρ is the great-circle distance on the sphere of radius ½. It also follows that d_F = sin ρ = ½‖n − n′‖, half the chord on the unit sphere, which is the chord on S²(½).

### Gaps
- I did not see the published book pages, only Dryden's own lecture-note version of the text. Equation numbers (3.2–3.8, 4.10–4.14) follow those notes and may differ from the printed 1998 and 2016 editions.
- I did not see Kendall 1984's own statement of these formulas. The attribution of ρ as the "Procrustean metric" follows Dryden's notes and the paper's title.

---

## 5. Reflections in shape analysis: reflection shape, size-and-reflection shape, and when to allow them

### Takeaway
Shape analysis distinguishes four spaces: shape (quotient by translation, scale, SO(m)), size-and-shape (no scaling), and their *reflection* versions (quotient by O(m)). The default in morphometric shape analysis is **no reflection**: Kendall's spaces quotient by rotations only, and `procGPA` has `reflection=FALSE` by default. The psychometric and MDS tradition (Green, Schönemann, Gower's GPA) allows reflections naturally. So do analyses built on inter-landmark distances, because distances cannot see handedness. Robotics, crystallography and attitude estimation must exclude reflections (det = +1).

### Cited Findings
- **Definition 4.10.** "The reflection shape of a configuration matrix X is all the geometrical information that is invariant under the similarity transformations and reflection ... [X]_R = {ZR : R ∈ O(m)}", where O(m) is the set of orthogonal matrices with |R| = ±1 and Z is the pre-shape.
- **Definition 4.11.** "The reflection size-and-shape ... is invariant under translation, rotation and reflection ... [X]_RS = {X_H R : R ∈ O(m)}", where X_H are the Helmertized coordinates.

  Dryden's hierarchy figure (after Goodall & Mardia 1992) runs from configuration, through Helmertized/centred and pre-shape or size-and-shape, to shape and reflection shape, with "remove reflection" arrows. — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Defaults in software.** In `shapes`, "The option reflection=TRUE uses reflection invariance whereas the default value does not have reflection invariance." — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **When the distinction matters.** "For fairly close shapes λ_m will usually be positive – in which case the solution is the same as minimizing over the orthogonal matrices O(m) instead of SO(m)". The distinction only matters for dissimilar or nearly degenerate shapes. — [Dryden notes ch. 5](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap5.pdf)
- **Psychometric GPA includes reflection.** Gower's criterion includes "translating, rotating, reflecting and scaling". — [Gower 1975](https://doi.org/10.1007/BF02291478). Schönemann's constraint is T′T = I, with no determinant condition. — [Schönemann 1966](https://doi.org/10.1007/BF02289451)
- **Physical problems constrain to rotations.** The *constrained* orthogonal Procrustes problem requires "U ... to be a rotation matrix, i.e., an orthogonal matrix of determinant 1"; the *unconstrained* problem allows any orthogonal U. — [Lawrence et al. 2019](https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf)
- **Reflection shape in 3D statistics (records verified, content not seen):**
  - Preston & Wood (2011), "Bootstrap inference for mean reflection shape and size-and-shape with three-dimensional landmark data", Biometrika 98(1):49–63. — [DOI 10.1093/biomet/asq065](https://doi.org/10.1093/biomet/asq065)
  - Dryden, Kume, Le & Wood (2008), "A multi-dimensional scaling approach to shape analysis", Biometrika 95(4):779–798. — [DOI 10.1093/biomet/asn050](https://doi.org/10.1093/biomet/asn050)
- **Bilateral symmetry.** Kent & Mardia (2001), "Shape, Procrustes tangent projections and bilateral symmetry", Biometrika 88(2):469–485. — [DOI 10.1093/biomet/88.2.469](https://doi.org/10.1093/biomet/88.2.469). The Crossref record lists only J. T. Kent as author.

### Inferences
- **Complex form of reflection.** A reflection of a planar configuration is complex conjugation, w ↦ w̄, composed with a similarity. The mirror fit is z ≈ a w̄ with â = wᵀz/‖w‖²; note that (w̄)\*z = wᵀz, with no conjugation. The reflection-shape full Procrustes distance is therefore:

  d_{F,R}² = 1 − max(|w\*z|², |wᵀz|²)/(‖w‖²‖z‖²)

  This matches the code's `mirror` branch, where `conj` is applied to `u` and the sum becomes conj(Σ zᵢuᵢ), whose modulus is |wᵀz|. Checked numerically.
- **On the triangle shape sphere, reflection flips the height.** Reflection fixes the first two coordinates and negates the third: (x, y, h) ↦ (x, y, −h). Checked numerically. The equator (collinear triangles) is fixed, and the two equilateral poles (the two orientations) swap. The reflection shape space of triangles is therefore a closed hemisphere, i.e. a disc. Allowing reflection can make the best fit jump between a proper and a mirrored candidate where the two images of DEF are equidistant from ABC; that is what the code's "mirror" option and hysteresis manage.
- **When to allow reflection.**
  - *Standard to allow it:* when only distances or similarity structure matter and handedness carries no meaning, as in MDS and factor-analysis configurations, sensory-panel GPA, and analyses built from distance matrices. Distances are reflection-invariant.
  - *Not standard (excluded):* biological landmarks with chirality or left/right labelling, molecules (enantiomers), physical rigid-body pose (robotics, attitude, crystallography), and any study of handedness or asymmetry. There, reflection is either physically impossible or exactly the signal.

### Gaps
- I did not retrieve explicit guidance text (for example from Dryden & Mardia 2016 or Gower & Dijksterhuis 2004) stating when reflection invariance is recommended. The "when to allow" guidance above is my synthesis from the defaults, the definitions and the problem domains.
- The claim that Kendall's (1984) spaces quotient only by SO(m) follows from Dryden's definitions; Kendall's own text was not seen.

---

## 6. Procrustes and Kendall's shape spaces: pre-shape sphere, CP^{k−2}, triangles on S²(½)

### Takeaway
The steps from configuration to shape are:

- Centre (Helmert) and scale to unit size. This gives the **pre-shape**, a point on the unit sphere S^{(k−1)m−1}.
- Quotient by SO(m). This gives Kendall's shape space Σ^k_m.

The Riemannian (Procrustean) distance ρ is the smallest great-circle distance between rotation fibres, which is what Procrustes rotation optimises. For m = 2, Σ^k_2 is CP^{k−2} with the Fubini–Study metric. For triangles, CP¹ is a 2-sphere of radius ½, so ρ is ordinary great-circle distance on S²(½). Full and partial Procrustes distances are monotone functions of ρ: d_F = sin ρ and d_P = 2 sin(ρ/2).

### Cited Findings
- **Pre-shape** (Definitions 4.4–4.5): Z = HX/‖HX‖. "The pre-shape space S^k_m ≡ S^{(k−1)m−1} is a hypersphere of unit radius in (k − 1)m real dimensions, since ‖Z‖ = 1. The term 'pre-shape' signifies that we are one step away from shape – rotation still has to be removed. The term was coined by Kendall (1984)." Shape is obtained by identifying "all rotated versions of the pre-shape with each other". — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Fibres and distances.** "The shapes of configurations are represented by fibres on the pre-shape sphere", and shape distance is "the closest distance between the fibres". ρ is the closest great-circle distance and d_P the closest chordal distance, and their relationships are illustrated on a cross-section of the pre-shape sphere. — [Dryden notes ch. 4](https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap4.pdf)
- **Kendall 1984, planar case.** The paper's abstract (as indexed) says that for m = 2 "the shape-space can be identified with a version of CP^{k−2}". — [zendy.io index of Bull. LMS 16(2):81–121](https://zendy.io/title/10.1112/blms/16.2.81), a secondary index of the abstract. An arXiv paper restates it: the planar shape space "is simply a complex projective space ... with the Fubini-Study metric defining the geodesic distances". — [Anderson, arXiv:1307.1923](https://arxiv.org/abs/1307.1923)
- **Triangles.** For triangles the shape space "is isometric to the 2-sphere of radius 1/2 [K1984]", and the distance between two triangles is "the arc-length of the arc joining the two projected points on the sphere". — [Geomstats tutorial notebook](https://nbviewer.org/github/geomstats/geomstats/blob/sam/notebooks/16_real_world_applications__visualizations_in_kendall_shape_spaces.ipynb), a secondary source. Teaching material on "Kendall's sphere of triangles" is also at [J. Feydy, ENS](https://www.math.ens.fr/~feydy/Teaching/GeomData_2018/Kendall_1.html). One morphometrics mailing-list post describes Kendall's theorem as a "three-dimensional sphere"; this conflicts with the standard 2-sphere statement and appears to be an error.
- **General m.** "The Riemannian metric structure of the shape space Σ^k_m for k labelled points in R^m was given by Kendall for the atypically simple situations in which m = 1 or 2 and k ≥ 2. Here we deal with the general case ... by using the properties of Riemannian submersions and warped products ... size-and-shape space ... is the warped product of the shape space and the half-line R₊ (carrying size), the warping function being equal to the square of the size." — [Le & Kendall 1993, Ann. Statist. 21(3):1225–1271](https://doi.org/10.1214/aos/1176349259)
- **Kendall's 1989 survey.** Shape is defined "for a set of k points in m dimensions". "The first task is to identify the shape spaces in which such objects naturally live, and then to examine the probability structures induced". The survey also outlines "size-and-shape", shapes in Riemannian spaces, and random Delaunay tessellations. Applications include archaeology, astronomy, geography and physical chemistry. — [Kendall 1989, Stat. Sci. 4(2):87–99](https://doi.org/10.1214/ss/1177012582)
- **Goodall's Procrustes view of shape.** "Two geometrical figures ... have the same shape if they differ by at most a rotation, a translation and isotropic scaling." "With few exceptions, the metric geometry of shape spaces is quite complicated." The strategy is "to work directly with the N × K co-ordinate matrix, while allowing for an arbitrary similarity transformation at all stages". — [Goodall 1991](https://doi.org/10.1111/j.2517-6161.1991.tb01825.x)
- **Complex Bingham on the pre-shape sphere.** A natural model for planar shape, invariant under the rotation fibres. — [Kent 1994](https://doi.org/10.1111/j.2517-6161.1994.tb01978.x)
- **Bookstein's coordinate approach** for 2D landmark size and shape: [Bookstein 1986](https://doi.org/10.1214/ss/1177013696); [Bookstein 1991](https://doi.org/10.1017/CBO9780511573064). Exact distributions of random-triangle shapes on Kendall's sphere: [Kendall 1985, Adv. Appl. Prob. 17(2):308–329](https://doi.org/10.2307/1427143).
- **Overviews:** Small, *The Statistical Theory of Shape* (Springer 1996), [DOI 10.1007/978-1-4612-4032-7](https://doi.org/10.1007/978-1-4612-4032-7); a 2026 review, [Choi et al., arXiv:2606.17022](https://arxiv.org/abs/2606.17022), which notes that "Quotienting by rigid motions and scalings yields Kendall's classical shape spaces."
- **Link to the Barbour topic of this folder.** Anderson argues that Kendall's shape statistics is a classical realisation of a "Barbour-type timeless records theory". — [arXiv:1307.1923](https://arxiv.org/abs/1307.1923)

### Inferences
- **Procrustes is the computational form of Kendall's quotient geometry.** OPA's rotation step finds the closest point on one fibre to a given point on another. Partial OPA on unit-size configurations returns the chord d_P. Full OPA, by also optimising scale, returns the perpendicular distance sin ρ = d_F. The geodesic ρ = arccos|z\*w| (unit z, w) is the intrinsic Fubini–Study-type distance on CP^{k−2}. GPA's full Procrustes mean is an extrinsic mean (Kent's eigenvector). An intrinsic Fréchet mean under ρ would instead require iteration (Le 1998 discusses uniqueness of Fréchet means on planar shape space).
- **Triangles in this codebase.** `shapePoint()` uses sub-Helmert coordinates w₁ = (p₂ − p₁)/√2 and w₂ = (2p₃ − p₁ − p₂)/√6. It then maps (w₁, w₂) ∈ C² to (|w₁|² − |w₂|², 2 Re(w̄₁w₂), 2 Im(w̄₁w₂))/(|w₁|² + |w₂|²) on the unit sphere, which is the Hopf map CP¹ → S².
  - The unit sphere is Kendall's S²(½) scaled by 2, so angles on it are 2ρ. Equilateral triangles map to the poles (checked numerically: (0, 0, ±1)), and collinear triangles map to the equator (h = 0).
  - The best vertex pairing (one of 3! relabellings, times 2 with mirroring) is the image of DEF nearest to ABC in ρ. Ties lie on the great circles equidistant from two images.
- Dryden's dimension formula for Σ^k_m is (k−1)m − 1 − m(m−1)/2. This gives real dimension 2k − 4 for m = 2 (matching CP^{k−2}) and 2 for triangles. I did not re-verify this formula against the printed text this session.

### Gaps
- Kendall 1984 itself (and the Kendall et al. 1999 book) were not readable this session. The CP^{k−2} and S²(½) statements are supported by the indexed abstract and by secondary and teaching sources, not by the primary text. The Fubini–Study curvature normalisation (holomorphic sectional curvature 4, matching the radius-½ sphere for CP¹) is from memory and unverified.
- I did not open Bookstein 1986 or the discussion contributions (including Kendall's), so the specifics of Bookstein's shape coordinates and their relation to Kendall's coordinates are not cited here beyond Kent's (1994) mention of "links and contrasts with Bookstein co-ordinates".

---

## Source catalogue

"verified = yes" means the bibliographic record (authors, title, venue, volume, issue, pages, DOI) was seen this session via the Crossref API, the arXiv API, or the publisher or author PDF. Notes in brackets flag discrepancies. "no" means the entry is from memory.

| key | authors | year | title | venue/volume/pages | DOI | arXiv ID | open-access PDF URL | verified |
|---|---|---|---|---|---|---|---|---|
| Mosier1939 | Mosier, Charles I. | 1939 | Determining a Simple Structure When Loadings for Certain Tests are Known | Psychometrika 4(2):149–162 | 10.1007/BF02288493 | | | yes |
| Green1952 | Green, Bert F. | 1952 | The Orthogonal Approximation of an Oblique Structure in Factor Analysis | Psychometrika 17(4):429–440 | 10.1007/BF02288918 | | | yes |
| HurleyCattell1962 | Hurley, John R.; Cattell, Raymond B. | 1962 | The Procrustes program: Producing direct rotation to test a hypothesized factor structure | Behavioral Science 7(2):258–262 | 10.1002/bs.3830070216 | | | yes [Crossref "issued" 2007 = online digitisation; vol. 7 = 1962] |
| Wahba1965 | Wahba, Grace | 1965 | A Least Squares Estimate of Satellite Attitude | SIAM Review 7(3):409 | 10.1137/1007077 | | | yes |
| FarrellStuelpnagel1966 | Farrell, J. L.; Stuelpnagel, J. C.; Wessner, R. H.; Velman, J. R.; Brook, J. E. | 1966 | A Least Squares Estimate of Satellite Attitude (Grace Wahba) [problem solutions] | SIAM Review 8(3):384–386 | 10.1137/1008080 | | | yes |
| Schonemann1966 | Schönemann, Peter H. | 1966 | A Generalized Solution of the Orthogonal Procrustes Problem | Psychometrika 31(1):1–10 | 10.1007/BF02289451 | | | yes |
| Cliff1966 | Cliff, Norman | 1966 | Orthogonal Rotation to Congruence | Psychometrika 31(1):33–42 | 10.1007/BF02289455 | | | yes |
| Browne1967 | Browne, Michael W. | 1967 | On Oblique Procrustes Rotation | Psychometrika 32(2):125–132 | 10.1007/BF02289420 | | | yes |
| SchonemannCarroll1970 | Schönemann, Peter H.; Carroll, Robert M. | 1970 | Fitting One Matrix to Another Under Choice of a Central Dilation and a Rigid Motion | Psychometrika 35(2):245–255 | 10.1007/BF02291266 | | | yes |
| KristofWingersky1971 | Kristof, W.; Wingersky, B. | 1971 | Generalization of the orthogonal Procrustes rotation procedure to more than two matrices | Proceedings of the 79th Annual Convention of the American Psychological Association 6:89–90 | | | | no |
| Gower1975 | Gower, J. C. | 1975 | Generalized Procrustes Analysis | Psychometrika 40(1):33–51 | 10.1007/BF02291478 | | | yes |
| Kabsch1976 | Kabsch, W. | 1976 | A solution for the best rotation to relate two sets of vectors | Acta Crystallographica Section A 32(5):922–923 | 10.1107/S0567739476001873 | | | yes |
| TenBerge1977 | Ten Berge, Jos M. F. | 1977 | Orthogonal Procrustes Rotation for Two or More Matrices | Psychometrika 42(2):267–276 | 10.1007/BF02294053 | | | yes |
| TenBergeNevels1977 | ten Berge, Jos M. F.; Nevels, Klaas | 1977 | A General Solution to Mosier's Oblique Procrustes Problem | Psychometrika 42(4):593–600 | 10.1007/BF02295981 | | | yes |
| Kendall1977 | Kendall, D. G. | 1977 | The diffusion of shape | Advances in Applied Probability 9(3):428–430 | 10.2307/1426091 | | | yes |
| Kabsch1978 | Kabsch, W. | 1978 | A discussion of the solution for the best rotation to relate two sets of vectors | Acta Crystallographica Section A 34(5):827–828 | 10.1107/S0567739478001680 | | | yes |
| Sibson1978 | Sibson, Robin | 1978 | Studies in the Robustness of Multidimensional Scaling: Procrustes Statistics | Journal of the Royal Statistical Society Series B 40(2):234–238 | 10.1111/j.2517-6161.1978.tb01669.x | | | yes |
| HansonNorris1981 | Hanson, Richard J.; Norris, Michael J. | 1981 | Analysis of Measurements Based on the Singular Value Decomposition | SIAM Journal on Scientific and Statistical Computing 2(3):363–373 | 10.1137/0902029 | | | yes |
| SiegelBenson1982 | Siegel, Andrew F.; Benson, Richard H. | 1982 | A Robust Comparison of Biological Shapes | Biometrics 38(2):341–350 | 10.2307/2530448 | | | yes [Crossref gives start page 341 only; end page from memory] |
| Kendall1984 | Kendall, David G. | 1984 | Shape Manifolds, Procrustean Metrics, and Complex Projective Spaces | Bulletin of the London Mathematical Society 16(2):81–121 | 10.1112/blms/16.2.81 | | | yes |
| Kendall1985 | Kendall, D. G. | 1985 | Exact distributions for shapes of random triangles in convex sets | Advances in Applied Probability 17(2):308–329 | 10.2307/1427143 | | | yes |
| Bookstein1986 | Bookstein, Fred L. | 1986 | Size and Shape Spaces for Landmark Data in Two Dimensions | Statistical Science 1(2):181–222 | 10.1214/ss/1177013696 | | | yes |
| Arun1987 | Arun, K. S.; Huang, T. S.; Blostein, S. D. | 1987 | Least-Squares Fitting of Two 3-D Point Sets | IEEE Transactions on Pattern Analysis and Machine Intelligence PAMI-9(5):698–700 | 10.1109/TPAMI.1987.4767965 | | | yes |
| Horn1987 | Horn, Berthold K. P. | 1987 | Closed-form solution of absolute orientation using unit quaternions | Journal of the Optical Society of America A 4(4):629–642 | 10.1364/JOSAA.4.000629 | | https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf | yes [Crossref gives start page 629; end page 642 from memory; author PDF opened] |
| Kendall1989 | Kendall, David G. | 1989 | A Survey of the Statistical Theory of Shape | Statistical Science 4(2):87–99 | 10.1214/ss/1177012582 | | | yes |
| MardiaDryden1989 | Mardia, K. V.; Dryden, I. L. | 1989 | The statistical analysis of shape data | Biometrika 76(2):271–281 | 10.1093/biomet/76.2.271 | | | yes |
| RohlfSlice1990 | Rohlf, F. James; Slice, Dennis | 1990 | Extensions of the Procrustes Method for the Optimal Superimposition of Landmarks | Systematic Zoology 39(1):40–59 | 10.2307/2992207 | | | yes [Crossref gives start page 40; end page from memory] |
| Goodall1991 | Goodall, Colin | 1991 | Procrustes Methods in the Statistical Analysis of Shape | Journal of the Royal Statistical Society Series B 53(2):285–321 | 10.1111/j.2517-6161.1991.tb01825.x | | | yes [Crossref pages 285–321; often cited as 285–339 including discussion] |
| Umeyama1991 | Umeyama, S. | 1991 | Least-squares estimation of transformation parameters between two point patterns | IEEE Transactions on Pattern Analysis and Machine Intelligence 13(4):376–380 | 10.1109/34.88573 | | | yes |
| Bookstein1991 | Bookstein, Fred L. | 1991 | Morphometric Tools for Landmark Data: Geometry and Biology | Cambridge University Press | 10.1017/CBO9780511573064 | | | yes [Crossref lists 1992; subtitle from memory] |
| LeKendall1993 | Le, Huiling; Kendall, David G. | 1993 | The Riemannian Structure of Euclidean Shape Spaces: A Novel Environment for Statistics | The Annals of Statistics 21(3):1225–1271 | 10.1214/aos/1176349259 | | | yes |
| Lele1993 | Lele, Subhash | 1993 | Euclidean Distance Matrix Analysis (EDMA): Estimation of mean form and mean form difference | Mathematical Geology 25(5):573–602 | 10.1007/bf00890247 | | | yes |
| Kent1994 | Kent, John T. | 1994 | The Complex Bingham Distribution and Shape Analysis | Journal of the Royal Statistical Society Series B 56(2):285–299 | 10.1111/j.2517-6161.1994.tb01978.x | | | yes |
| Small1996 | Small, Christopher G. | 1996 | The Statistical Theory of Shape | Springer Series in Statistics, Springer New York | 10.1007/978-1-4612-4032-7 | | | yes |
| KentMardia1997 | Kent, John T.; Mardia, Kanti V. | 1997 | Consistency of Procrustes Estimators | Journal of the Royal Statistical Society Series B 59(1):281–290 | 10.1111/1467-9868.00069 | | | yes |
| Le1998 | Le, Huiling | 1998 | On the consistency of procrustean mean shapes | Advances in Applied Probability 30(1):53–63 | 10.1239/aap/1035227991 | | | yes |
| DrydenMardia1998 | Dryden, Ian L.; Mardia, Kanti V. | 1998 | Statistical Shape Analysis | Wiley, Chichester (xvii+347 pp.; ISBN 0-471-95816-6) | | | | yes [details via Rudemo 2000 review record, 10.1002/1097-0258(20001015)19:19<2716::aid-sim590>3.0.co;2-o] |
| KendallBardenCarneLe1999 | Kendall, D. G.; Barden, D.; Carne, T. K.; Le, H. | 1999 | Shape and Shape Theory | Wiley Series in Probability and Statistics, John Wiley & Sons (ISBN 978-0-471-96823-8) | 10.1002/9780470317006 | | | yes |
| KentMardia2001 | Kent, John T.; Mardia, Kanti V. | 2001 | Shape, Procrustes tangent projections and bilateral symmetry | Biometrika 88(2):469–485 | 10.1093/biomet/88.2.469 | | | yes [Crossref lists only Kent, J. T.; Mardia as co-author from memory] |
| GowerDijksterhuis2004 | Gower, John C.; Dijksterhuis, Garmt B. | 2004 | Procrustes Problems | Oxford Statistical Science Series 30, Oxford University Press (ISBN 978-0-19-851058-1) | 10.1093/acprof:oso/9780198510581.001.0001 | | | yes [series number 30 from memory] |
| Coutsias2004 | Coutsias, Evangelos A.; Seok, Chaok; Dill, Ken A. | 2004 | Using quaternions to calculate RMSD | Journal of Computational Chemistry 25(15):1849–1857 | 10.1002/jcc.20110 | | | yes |
| Groisser2005 | Groisser, David | 2005 | On the convergence of some Procrustean averaging algorithms | Stochastics 77(1):31–60 | 10.1080/17442500512331341059 | | | yes |
| DrydenKumeLeWood2008 | Dryden, I. L.; Kume, A.; Le, H.; Wood, A. T. A. | 2008 | A multi-dimensional scaling approach to shape analysis | Biometrika 95(4):779–798 | 10.1093/biomet/asn050 | | | yes |
| MyronenkoSong2009 | Myronenko, Andriy; Song, Xubo | 2009 | On the closed-form solution of the rotation matrix arising in computer vision problems | arXiv preprint | | 0904.1613 | https://arxiv.org/pdf/0904.1613 | yes |
| PrestonWood2011 | Preston, S. P.; Wood, A. T. A. | 2011 | Bootstrap inference for mean reflection shape and size-and-shape with three-dimensional landmark data | Biometrika 98(1):49–63 | 10.1093/biomet/asq065 | | | yes |
| Anderson2013 | Anderson, Edward | 2013 | Kendall's Shape Statistics as a Classical Realization of Barbour-type Timeless Records Theory approach to Quantum Gravity | arXiv preprint (v5) | | 1307.1923 | https://arxiv.org/pdf/1307.1923 | yes [journal version, if any, not checked] |
| DrydenMardia2016 | Dryden, Ian L.; Mardia, Kanti V. | 2016 | Statistical Shape Analysis, with Applications in R (2nd ed.) | Wiley Series in Probability and Statistics, Wiley (ISBN 978-0-470-69962-1) | 10.1002/9781119072492 | | | yes |
| Lawrence2019 | Lawrence, Jim; Bernal, Javier; Witzgall, Christoph | 2019 | A Purely Algebraic Justification of the Kabsch-Umeyama Algorithm | Journal of Research of the National Institute of Standards and Technology 124:124028 | 10.6028/jres.124.028 | | https://nvlpubs.nist.gov/nistpubs/jres/124/jres.124.028.pdf | yes |
| Brussel2024 | Brussel, Eric; Goertz, Madeleine; Guptill, Elijah; Lyle, Kelly | 2024 | The Stack of Similarity Classes of Triangles | arXiv preprint (v3; v2 titled "The Stack of Triangles up to Similarity") | | 2408.07792 | https://arxiv.org/pdf/2408.07792 | yes |
| Choi2026 | Choi, Gary P. T.; Dao Duc, Khanh; Faigenbaum-Golovin, Shira; Habermann, Karen; Hartman, Emmanuel; von Tycowicz, Christoph; Zhang, Chi; Zhao, Wenjun; Zhou, Felix | 2026 | Learning the Geometry of Data: A Mathematical Review of Shape Space Analysis | arXiv preprint (v1, 15 Jun 2026) | | 2606.17022 | https://arxiv.org/pdf/2606.17022 | yes |
| DrydenNotesSTAT718A | Dryden, Ian L. | n.d. | STAT 718A lecture notes, chapters 3 (Planar Procrustes Analysis), 4 (Shape Space and Distances), 5 (Procrustes Analysis) | Course notes, University of South Carolina (lecture version of Dryden & Mardia) | | | https://people.stat.sc.edu/dryden/STAT718A/notes/shape-chap3.pdf ; …/shape-chap4.pdf ; …/shape-chap5.pdf | yes [year not stated on the page] |
| ShapesProcdist | Dryden, Ian L. (package author) | n.d. | procdist: Procrustes distance (R package `shapes`) | CRAN reference manual | | | https://search.r-project.org/CRAN/refmans/shapes/html/procdist.html | yes [page reachable; version not recorded] |
| EigenUmeyama | Heibel, Hauke (Eigen contributors) | n.d. | Eigen/src/Geometry/Umeyama.h | Eigen C++ library source (implements Umeyama 1991, eqs. 38–43) | | | https://eigen.googlesource.com/mirror/+/17af8c763d5aa58c40fc2c40ff6e8d6fc28a8d63/Eigen/src/Geometry/Umeyama.h | yes |
