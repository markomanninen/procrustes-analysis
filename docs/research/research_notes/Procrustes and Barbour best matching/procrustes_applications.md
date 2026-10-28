# Applications of Procrustes analysis across scientific fields

Scope and method note (read first). These notes were compiled on 2026-10-28. Bibliographic records were checked this session against Crossref (`api.crossref.org/works/{DOI}`), doi.org content negotiation (for the Hystrix DOIs, which are not in Crossref), the arXiv API, and OpenAlex. OpenAlex supplied abstracts, open-access locations and citation counts. Content claims were checked against abstracts (via OpenAlex), full texts (PMC, arXiv PDFs) or official documentation (CRAN, SciPy, vegan, MDAnalysis, MorphoJ, rdrr.io). Citation counts are given as "Crossref cited-by / OpenAlex cited-by, as of 2026-10-28". They are only a rough adoption proxy: the two databases disagree, and neither covers all citations. Anything recalled from background knowledge and not re-verified this session is labelled as such, or is placed under Gaps.

Terminology used below:
- **Ordinary/orthogonal Procrustes (OPA)**: fit one configuration to another.
- **Generalized Procrustes analysis (GPA)**: fit many configurations iteratively to their mean.
- **Full vs. partial Procrustes**: with or without optimal scaling.
- **"Proper" rotation**: det = +1, so reflections are excluded.
- **Absolute orientation / Kabsch / Wahba**: the same SVD/quaternion core under the names used by different fields.

---

## Q1. Geometric morphometrics (GM): what is the standard Procrustes workflow, which variant is used, what software, and which works are landmarks?

### Takeaway
Geometric morphometrics is probably the largest and most institutionalised user of Procrustes methods. The field describes itself as built on a "Procrustes paradigm". The standard pipeline is:
1. Digitise landmarks and semilandmarks.
2. Run GPA (translation, unit centroid size, rotation; semilandmarks slid by bending energy or Procrustes distance).
3. Take the Procrustes shape coordinates, projected to a tangent space.
4. Run PCA, then Procrustes ANOVA or regression with permutation (RRPP), allometry and symmetry decompositions, and phylogenetic comparative analyses.

The dominant tools are the R packages geomorph and RRPP, MorphoJ, and Rohlf's tps series, now increasingly fed by 3D pipelines (SlicerMorph, ALPACA).

### Cited Findings
**Founding methods and theory**
- Rohlf & Slice (1990) reviewed superimposition methods. In these methods, shape differences are the residuals left after rotating, translating and scaling specimens into best alignment. They generalised Siegel & Benson's two-object resistant-fit (theta-rho) method to many objects at once (the generalized resistant fit, GRF), and extended both least-squares and resistant fits to affine (uniform) transformations. Their test data were 18 landmarks on the wings of 127 species of North American mosquitoes. — [Rohlf & Slice 1990](https://doi.org/10.2307/2992207). Citations: 3135 Crossref / 4287 OpenAlex.
- Siegel & Benson (1982) is the original "robust comparison of biological shapes" that the resistant fit builds on. — [Siegel & Benson 1982](https://doi.org/10.2307/2530448)
- Gower (1975) defined GPA as translating, rotating, **reflecting** and scaling m configurations to minimise the summed squared distances to their centroid configuration, which is read as a consensus. — [Gower 1975](https://doi.org/10.1007/BF02291478). Citations: 2885 / 3265.
- Goodall (1991): two figures of N landmarks in R^K have the same shape if they differ at most by rotation, translation and isotropic scaling. He presented a model-based Procrustes approach using QR and SVD decompositions. — [Goodall 1991](https://doi.org/10.1111/j.2517-6161.1991.tb01825.x)
- Kendall (1984) introduced the shape space Σ^k_m of k-ads in R^m as a quotient space with the quotient (Procrustean) metric. For planar shapes (m = 2) it can be identified with CP^(k−2). — [Kendall 1984](https://doi.org/10.1112/blms/16.2.81)
- After size, position and orientation are removed, the Procrustes shape coordinates are analysed statistically. Kendall shape space is a metric space that "can be approximated locally by a Euclidean tangent space", which is why GM statistics are run in tangent coordinates. — [Mitteroecker & Gunz 2009](https://doi.org/10.1007/s11692-009-9055-x); see also [Rohlf 1999, "Shape statistics: Procrustes superimpositions and tangent spaces"](https://doi.org/10.1007/s003579900054) (title verified; full text not read this session).
- Bookstein's thin-plate spline / principal warps paper is the standard source for visualising Procrustes shape differences as deformation grids and for bending energy. — [Bookstein 1989](https://doi.org/10.1109/34.24792) (3973 / 4923 cites); [Bookstein 1991 book](https://doi.org/10.1017/CBO9780511573064) (Crossref dates the online monograph 1992).

**History and reviews**
- The field dates its "revolution" to the shift towards landmark coordinates. Rohlf & Marcus (1993) proclaimed it. — [Rohlf & Marcus 1993](https://doi.org/10.1016/0169-5347(93)90024-J)
- Adams, Rohlf & Slice (2004) describe a fundamental change in the 1980s towards analysing landmark coordinates, which also allowed multivariate results to be visualised as shapes. — [Adams et al. 2004](https://doi.org/10.1080/11250000409356545)
- Adams, Rohlf & Slice (2013): GM "has matured into a rich and cohesive discipline", and "the development of the field is identified with the Procrustes paradigm". — [Adams et al. 2013, Hystrix](https://doi.org/10.4404/hystrix-24.1-6283)
- Mitteroecker & Schaefer (2022), the state-of-the-art review, list the building blocks of GM:
  - landmarks and semilandmarks;
  - shape or form variables computed by superimposition;
  - visualisation of statistical results as actual shapes;
  - decomposition into symmetric and asymmetric components and across spatial scales;
  - the geometry of shape and form spaces;
  - models linking shape to environmental, genetic or behavioural data.

  They also flag current challenges "arising from the increasing number of landmarks and semilandmarks". — [Mitteroecker & Schaefer 2022](https://doi.org/10.1002/ajpa.24531) (open access, [PMC9545184](https://pmc.ncbi.nlm.nih.gov/articles/PMC9545184/))

**Workflow components**
- Semilandmarks: sliding semilandmarks quantify homologous curves and surfaces together with traditional landmarks, and the same sliding algorithm can estimate missing data. — [Gunz & Mitteroecker 2013](https://doi.org/10.4404/hystrix-24.1-6292). A practical guide to sliding and surface semilandmarks was motivated by high-resolution CT and surface-scan datasets. — [Bardua et al. 2019](https://doi.org/10.1093/iob/obz016)
- geomorph's `gpagen` performs "Generalized Procrustes analysis of points, curves, and surfaces" with fixed landmarks plus sliding semilandmarks on curves and surfaces. Its options include:
  - `ProcD`: slide by Procrustes distance, otherwise by bending energy;
  - an approximate thin-plate-spline mode for bending energy, which is faster and allows more iterations;
  - `max.iter`: GPA iterations.

  — [geomorph `gpagen` documentation](https://rdrr.io/cran/geomorph/man/gpagen.html). geomorph also provides `gm.prcomp` ("Principal and phylogenetically-aligned components analysis") and `bilat.symmetry`. — [same page, function index](https://rdrr.io/cran/geomorph/man/gpagen.html)
- Procrustes ANOVA and linear models: RRPP ("residual randomization in permutation procedures") builds empirical sampling distributions for ANOVA statistics on high-dimensional data. It supports OLS or GLS, analysis of data or dissimilarity matrices, and type I/II/III sums of squares. — [Collyer & Adams 2018](https://doi.org/10.1111/2041-210X.13029). The phenotypic-change method for high-dimensional data is in [Collyer, Sekora & Adams 2015](https://doi.org/10.1038/hdy.2014.75) (title verified).
- Allometry: Klingenberg's (2016) systematic review contrasts two schools. The Gould–Mosimann school defines allometry as the covariation of shape with size, and the review examines how each concept is implemented in Procrustes shape space. — [Klingenberg 2016](https://doi.org/10.1007/s00427-016-0539-2)
- Symmetry and asymmetry: for structures with object symmetry, the configuration and a **reflected copy** are Procrustes-superimposed together. Variation is then partitioned into symmetric variation among individuals and asymmetry, in a way compatible with two-factor ANOVA (the basis of "Procrustes ANOVA" for fluctuating asymmetry). — [Klingenberg, Barluenga & Meyer 2002](https://doi.org/10.1111/j.0014-3820.2002.tb00117.x); [Klingenberg & McIntyre 1998](https://doi.org/10.1111/j.1558-5646.1998.tb02018.x)

**Software (versions current as of 2026-10-28)**
- **geomorph** (R), introduced by Adams & Otárola-Castillo (2013). It provides routines for "all stages of landmark-based geometric morphometric analyses in two and three-dimensions", including Procrustes superimposition of landmarks and of semilandmarks from curves and surfaces. — [Adams & Otárola-Castillo 2013](https://doi.org/10.1111/2041-210X.12035). Citations: 2069 / 2177.
  - geomorph v4.0 added gmShiny, a GUI, and positions geomorph "to be the primary tool for GM analyses, particularly those employing a phylogenetic comparative approach". — [Baken et al. 2021](https://doi.org/10.1111/2041-210X.13723) (680 / 666 cites)
  - Current CRAN version: **4.1.1, published 2026-07-05**. — [CRAN geomorph](https://cran.r-project.org/package=geomorph)
- **RRPP** CRAN version **2.2.0 (2026-07-03)**. — [CRAN RRPP](https://cran.r-project.org/package=RRPP)
- **MorphoJ** (Java, Apache 2.0 licence). It extracts shape information "by Procrustes superimposition" and offers PCA, discriminant analysis, multivariate regression and more. — [Klingenberg 2011](https://doi.org/10.1111/j.1755-0998.2010.02924.x) (3392 / 3776 cites, the most-cited GM software paper in this set)
  - Feature list: Procrustes fit for 2D/3D data with and without object symmetry; outlier identification; pooled within-group regression "for allometry correction in multiple groups"; two-block PLS; modularity tests; CVA/DFA; phylogenetic mapping; quantitative-genetics analyses.
  - Current version **1.08.02**. — [MorphoJ page](https://morphometrics.uk/MorphoJ_page.html)
- **tps series** (Rohlf; Windows): "used in hundreds of studies in mammals and other organisms". — [Rohlf 2015](https://doi.org/10.4404/hystrix-26.1-11264) (OpenAlex 961 cites)
- **Morpho / Rvcg** (R) for 3D meshes and sliding semilandmarks. — [Schlager 2017](https://doi.org/10.1016/B978-0-12-810493-4.00011-0). CRAN Morpho **2.13 (2025-09-04)**. — [CRAN Morpho](https://cran.r-project.org/package=Morpho)
- **shapes** (R, Dryden). CRAN **1.2.8 (2025-10-24)**. — [CRAN shapes](https://cran.r-project.org/package=shapes). The companion textbook is [Dryden & Mardia 2016](https://doi.org/10.1002/9781119072492). It cites uses such as brain shape in schizophrenia and protein molecules.
- **3D data and automation (state of the art)**: SlicerMorph (3D Slicer extension) responds to projects such as #ScanAllFishes and oVert, which generate microCT scans "by the thousands". MorphoSource already hosts "tens of thousands of 3D scans". — [Rolfe et al. 2021](https://doi.org/10.1111/2041-210X.13669). ALPACA automates 3D landmarking because 3D data "are still mainly gathered by trained experts". — [Porto, Rolfe & Maga 2021](https://doi.org/10.1111/2041-210X.13689)
- The standard textbook is [Zelditch, Swiderski & Sheets, *Geometric Morphometrics for Biologists: A Primer*, 2nd ed. 2012](https://doi.org/10.1016/b978-0-12-386903-6.00014-9) (a chapter DOI from the 2012 edition). Crossref also lists 2026-dated chapters by the same authors ("Allometry", "Morphometrics and systematics") under the ISBN prefix 978-0-443-28965-1, which suggests a new edition in 2026. — [Crossref chapter record](https://doi.org/10.1016/b978-0-443-28965-1.00003-8)

**Quantitative indicators**
- Typical data sizes in the verified sources:
  - 18 landmarks × 127 species ([Rohlf & Slice 1990](https://doi.org/10.2307/2992207));
  - simulation figures of 6–20 landmarks ([Walker 2000](https://doi.org/10.1080/106351500750049770));
  - modern 3D projects with thousands of scans ([Rolfe et al. 2021](https://doi.org/10.1111/2041-210X.13669)).
- A dedicated study of sample size and sampling error in GM exists: [Cardini & Elton 2007](https://doi.org/10.1007/s00435-007-0036-2). Its title was verified; the abstract was not retrievable this session.

### Inferences
- The "standard workflow" in the question (landmarks/semilandmarks → GPA → Procrustes coordinates → tangent projection → PCA / Procrustes ANOVA / allometry) is consistent with every GM source checked. The modern additions are:
  - semilandmark sliding inside GPA (`gpagen`);
  - permutation-based Procrustes ANOVA (RRPP) rather than parametric F-tests;
  - phylogenetic comparative extensions (geomorph v4);
  - symmetry decomposition via reflected copies (MorphoJ, `bilat.symmetry`).
- GM generally uses **full GPA with scaling to unit centroid size and no reflections**. Reflection enters only deliberately, in the object-symmetry analysis. This contrasts with Gower's original GPA and with ecology's PROTEST, which allow reflection.
- Citation counts suggest that MorphoJ and geomorph together account for the bulk of software citations in GM (about 3.4k and 2.1k Crossref cites for the founding papers). The tps series is comparably cited (about 1k via its 2015 paper, but older tps versions were cited by URL and are undercounted).

### Gaps
- I found no systematic survey of typical sample sizes or landmark counts in published GM studies. Cardini & Elton 2007 and Cardini 2020 (on "high-density" morphometrics, [DOI](https://doi.org/10.1007/s00435-020-00499-w)) are the relevant sources, but their abstracts and full texts were not read this session.
- Page ranges for Adams et al. 2013 (7–14), Gunz & Mitteroecker 2013 (103–109) and Rohlf 2015 (9–12) are from memory. The doi.org metadata did not include pages.
- The geomorph function `procD.lm` (Procrustes ANOVA / regression) is standard but was not seen in this session's documentation snippet. Treat the function name as background knowledge.
- The exact edition number and date of the 2026 Zelditch et al. volume were not confirmed.

---

## Q2. Structural biology and chemistry: how are Kabsch, quaternion and QCP superposition used, and how is RMSD used in the PDB and in structure comparison?

### Takeaway
In structural biology the Procrustes problem is "optimal rigid superposition minimising RMSD": proper rotation plus translation, no scaling.
- **Algorithms**: Kabsch's SVD-based solution (1976/1978), the quaternion eigenproblem (McLachlan, Kearsley, Horn, Coutsias et al. 2004), and Theobald's QCP (2005), which is now the default in high-throughput tools such as MDAnalysis.
- **How RMSD is reported**: in Å, together with the number of aligned atoms or residues. It is increasingly complemented or replaced by size-normalised or superposition-free scores (TM-score, GDT/LGA, lDDT). RCSB PDB's alignment tool reports RMSD (Å) alongside TM-score.
- **Ligands**: "pocket-aligned ligand RMSD < 2 Å" is the de facto success criterion (AlphaFold 3, PoseBusters).
- **Refinements**: maximum-likelihood superposition (THESEUS) and symmetry-aware RMSD address known pitfalls.

### Cited Findings
**Algorithms**
- Kabsch (1976) derived "a best rotation of a given vector set into a second vector set by minimizing the weighted sum of squared deviations". — [Kabsch 1976](https://doi.org/10.1107/S0567739476001873) (2938 / 3205 cites). The follow-up discussion paper is [Kabsch 1978](https://doi.org/10.1107/S0567739478001680) (994 / 1145).
- Quaternion formulations:
  - Kearsley (1989) found the optimal rotation "as an eigenvalue problem in quaternion parameters". — [Kearsley 1989](https://doi.org/10.1107/S0108767388010128)
  - McLachlan's "Rapid comparison of protein structures" is an earlier crystallographic solution. — [McLachlan 1982](https://doi.org/10.1107/S0567739482001806)
  - Horn (1987) gave the closed-form unit-quaternion solution of "absolute orientation", a "classic photogrammetric task" also used in robotics. — [Horn 1987](https://doi.org/10.1364/JOSAA.4.000629)
- Coutsias, Seok & Dill (2004) derive the quaternion method and prove it is **equivalent to Kabsch's formula**. They enumerate the special cases via the eigenvalues of a traceless 4×4 symmetric matrix and give the gradient of RMSD with respect to model parameters, which is useful for elastic-band reaction paths and for fitting. — [Coutsias et al. 2004](https://doi.org/10.1002/jcc.20110)
- Theobald (2005), the QCP method: the rotations are often not needed, and common algorithms require eigendecomposition or inversion of 3×3 or 4×4 matrices. QCP "rapidly and stably determines RMSDs by circumventing the decomposition and inversion problems". — [Theobald 2005](https://doi.org/10.1107/S0108767305015266)
  - Liu, Agrafiotis & Theobald (2010) added the rotation matrix via Newton–Raphson plus an adjoint matrix, "at least an order of magnitude more efficient than conventional inversion/decomposition methods". — [Liu et al. 2010](https://doi.org/10.1002/jcc.21439)
- MDAnalysis's `rms` module "uses the fast QCP algorithm [Theobald2005]" for superposition and RMSD, and asks users to cite Theobald 2005. — [MDAnalysis rms docs](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html); package paper [Michaud-Agrawal et al. 2011](https://doi.org/10.1002/jcc.21787) (4914 / 4541 cites).

**Beyond least squares**
- Theobald & Wuttke (2006): the ordinary unweighted least-squares Procrustes criterion "can give erroneous solutions when the errors have heterogeneous variances (heteroscedasticity) or the errors are correlated". They solve the non-isotropic maximum-likelihood Procrustes problem under a matrix Gaussian model. — [Theobald & Wuttke 2006, PNAS](https://doi.org/10.1073/pnas.0508445103)
  - Their program THESEUS "down-weight[s] variable structural regions", corrects for correlations among atoms, and is "insensitive to the specific atoms included", so it avoids subjective pruning. — [THESEUS, Bioinformatics 2006](https://doi.org/10.1093/bioinformatics/btl332)
- Symmetry: Coutsias & Wester (2019) treat alternative alignments due to indistinguishable atoms, and degenerate superpositions, for which a necessary condition is optimal superimposability to one's mirror image. Their software is frmsd. — [Coutsias & Wester 2019](https://doi.org/10.1002/jcc.25802)
  - For docking, assuming direct atomic correspondence for symmetric ligands "can lead to an artificial inflation" of RMSD. DockRMSD recasts symmetry correction as graph isomorphism. — [Bell & Zhang 2019](https://doi.org/10.1186/s13321-019-0362-7)
- Molecular dynamics: essential dynamics (protein PCA) separates configurational space into a few "essential" degrees of freedom once "overall translation and rotation are eliminated". This is the same superimpose-then-PCA pattern as in GM. — [Amadei, Linssen & Berendsen 1993](https://doi.org/10.1002/prot.340170408) (3324 / 3504 cites)

**How RMSD is used and reported**
- RMSD between equivalent atom pairs is "usually" the similarity measure, but it "depends on the dimension of the proteins, that is, on the number of equivalent atom pairs". Hence normalised RMSD (RMSD100). — [Carugo & Pongor 2001](https://doi.org/10.1110/ps.690101). The statistical significance of RMSD for globular proteins is analysed in [Maiorov & Crippen 1994](https://doi.org/10.1006/jmbi.1994.1017) (title verified).
- TM-score extends GDT and MaxSub with a protein-size-dependent scale "to eliminate the inherent protein size dependence" and to account for random structure pairs. — [Zhang & Skolnick 2004](https://doi.org/10.1002/prot.20264) (2409 cites). TM-align: [Zhang & Skolnick 2005](https://doi.org/10.1093/nar/gki524) (3859 cites).
- LGA (Local-Global Alignment) compares structures in sequence-dependent and sequence-independent modes, and its output can be used to rank similarity and classify structures. — [Zemla 2003](https://doi.org/10.1093/nar/gkg571)
- lDDT: "conventional similarity measures based on a global superposition of carbon α atoms are strongly influenced by domain motions and do not assess the accuracy of local atomic details". lDDT is superposition-free. — [Mariani et al. 2013](https://doi.org/10.1093/bioinformatics/btt473)
- **RCSB PDB**: the pairwise structure alignment tool (rcsb.org/alignment; API alignment.rcsb.org) offers FATCAT, CE, TM-align and Smith–Waterman 3D. Its output panel reports "Root-Mean-Square-Deviation or RMSD in Ångstrom Units; TM-score, range 0.00–1.00", plus % sequence identity and the number of aligned residues. — [Bittrich et al. 2024](https://doi.org/10.1093/bioinformatics/btae370) (full text [PMC11212067](https://pmc.ncbi.nlm.nih.gov/articles/PMC11212067/)). The underlying aligners are [CE](https://doi.org/10.1093/protein/11.9.739) and [FATCAT](https://doi.org/10.1093/bioinformatics/btg1086).
- **AlphaFold 3 and PoseBusters**: protein–ligand accuracy is reported as the "percentage of pocket-aligned ligand r.m.s.d. < 2 Å".
  - The pocket "is used to align the predicted structure to the ground-truth structure with least-squares rigid alignment and then the r.m.s.d. is computed on all heavy atoms of the ligand".
  - The PoseBusters set has 428 protein–ligand structures released to the PDB in 2021 or later. AF3's general "recent PDB" evaluation set has 8,856 complexes.

  — [Abramson et al. 2024](https://doi.org/10.1038/s41586-024-07487-w) ([PMC11168924](https://pmc.ncbi.nlm.nih.gov/articles/PMC11168924/))
- PoseBusters warns that state-of-the-art claims "in terms of crystallographic RMSD" hide physically implausible poses. Evaluation "solely by RMSD to a native binding mode" is insufficient, so steric and energetic checks are added. — [Buttenschoen, Morris & Deane 2024](https://doi.org/10.1039/D3SC04185A)
- Crystallography: COMPACK compares crystal structures through interatomic distances of molecular packing environments, avoiding space-group and cell information. — [Chisholm & Motherwell 2005](https://doi.org/10.1107/S0021889804027074)
- Reviews of structure-comparison methods: [Kufareva & Abagyan 2012](https://doi.org/10.1007/978-1-61779-588-6_10) ([PMC4321859](https://pmc.ncbi.nlm.nih.gov/articles/PMC4321859/)). The classic sequence–structure divergence relation is [Chothia & Lesk 1986](https://doi.org/10.1002/j.1460-2075.1986.tb04288.x).

### Inferences
- Structural biology uses the **two-set, proper-rotation, no-scaling** Procrustes problem almost exclusively. It is the one field where computational speed of the core solver mattered enough to produce a dedicated algorithmic line (Kabsch → quaternion → QCP).
- Reporting has shifted. "Cα RMSD over N aligned residues" remains universal in the PDB and in MD, but headline benchmarks (CASP, AF3) now emphasise TM-score, GDT and lDDT for proteins. RMSD survives mainly as a thresholded success rate (< 2 Å) for ligands, always after a defined alignment frame such as pocket alignment. Which atoms define the superposition is thus part of the metric.
- ML-Procrustes (THESEUS) is the structural-biology analogue of GM's resistant fit and Pinocchio-effect concerns: both respond to least squares spreading localised variation across all points.

### Gaps
- I did not verify the specific CASP/AF2 RMSD figures (e.g., AF2's median backbone accuracy) this session. The AF2 paper record ([Jumper et al. 2021](https://doi.org/10.1038/s41586-021-03819-2)) was verified bibliographically only.
- Background knowledge, not verified this session:
  - Kabsch 1978 is the paper that handles the improper-rotation (reflection) case.
  - AF3's diffusion training loss uses a weighted rigid (Kabsch) alignment of prediction to ground truth.
- No source was found this session that quantifies how many PDB-era papers report RMSD. Adoption is therefore indicated only by citation counts.

---

## Q3. Computer vision, robotics and aerospace: Wahba's problem, ICP, Active Shape Models, SLAM and point-cloud registration

### Takeaway
The same closed-form solution appears under different names in different fields:
- "Wahba's problem" in spacecraft attitude estimation (1965; QUEST 1981);
- "absolute orientation" in photogrammetry and vision (Arun 1987 SVD; Horn 1987 quaternions; Umeyama 1991 with scale and reflection-safe);
- the inner step of ICP (Besl & McKay 1992; Chen & Medioni 1992) and of robust registration (TEASER 2021).

Statistical shape models (Active Shape Models, 1995) use GPA to align training shapes before PCA. In SLAM and pose estimation, Procrustes alignment is also built into evaluation metrics: absolute trajectory error after alignment, and Procrustes-aligned MPJPE ("Protocol #2") in 3D human pose.

### Cited Findings
**Aerospace (Wahba's problem)**
- Wahba posed "A least squares estimate of satellite attitude" as SIAM Review Problem 65-1. — [Wahba 1965](https://doi.org/10.1137/1007077) (943 / 1116 cites). Published solutions followed. — [Farrell et al. 1966](https://doi.org/10.1137/1008080)
- Shuster & Oh (1981) present TRIAD, a deterministic non-optimal solution from two vector observations, and QUEST, an optimal algorithm for the "best weighted overlap of an arbitrary number of reference and observation vectors", with covariance expressions. — [Shuster & Oh 1981](https://doi.org/10.2514/3.19717) (1136 / 1399 cites). TRIAD origin: [Black 1964](https://doi.org/10.2514/3.2555). Textbook: [Markley & Crassidis 2014](https://doi.org/10.1007/978-1-4939-0802-8).
- SciPy's `Rotation.align_vectors(a, b, weights=None, return_sensitivity=False)`:
  - estimates the rotation "with Kabsch algorithm, and solves what is known as the 'pointing problem', or 'Wahba's problem'";
  - supports per-vector weights, including an infinite weight that makes "primary" vectors align exactly;
  - optionally returns a sensitivity matrix, and returns `rssd`.

  — [SciPy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html)

**Absolute orientation**
- Arun, Huang & Blostein (1987): least-squares R and T from the SVD of a 3×3 matrix. — [Arun et al. 1987](https://doi.org/10.1109/TPAMI.1987.4767965) (3409 / 3961 cites)
- Horn's orthonormal-matrix variant: [Horn, Hilden & Negahdaripour 1988](https://doi.org/10.1364/JOSAA.5.001127)
- Umeyama (1991) generalises to m-dimensional similarity transforms (rotation, translation, scaling). He notes that the Arun and Horn solutions "sometimes fail[] to give a correct rotation matrix and give[] a reflection instead when the data is severely corrupted", and gives a strict solution. — [Umeyama 1991](https://doi.org/10.1109/34.88573) (2165 / 2436 cites)
- A comparison of four closed-form rigid-transform algorithms: [Eggert, Lorusso & Fisher 1997](https://doi.org/10.1007/s001380050048) (title verified).

**ICP and registration**
- Besl & McKay (1992): ICP "always converges monotonically to the nearest local minimum of a mean-square distance metric". Global minimisation requires testing "an adequate set of initial rotations and translations". The headline application was registering sensed data to an ideal model before inspection. — [Besl & McKay 1992](https://doi.org/10.1109/34.121791) (14609 / 18600 cites, by far the most-cited paper in this catalogue)
- Parallel and variant ICP lines:
  - [Chen & Medioni 1992](https://doi.org/10.1016/0262-8856(92)90066-C) (2315 / 3164);
  - [Zhang 1994](https://doi.org/10.1007/BF01427149);
  - [Rusinkiewicz & Levoy 2001, "Efficient variants of the ICP algorithm"](https://doi.org/10.1109/IM.2001.924423);
  - [Segal, Haehnel & Thrun 2009, Generalized-ICP](https://doi.org/10.15607/RSS.2009.V.021).
- Robotics review: registration "associate[s] sets of data into a common coordinate system", is used in reconstruction, inspection, medicine and mobile-robot localisation, and the "plethora of solutions" can be organised by a few elements. — [Pomerleau, Colas & Siegwart 2015](https://doi.org/10.1561/2300000035)
- State of the art in robust and certifiable registration: TEASER reformulates registration with a truncated-least-squares cost that is "insensitive to a large fraction of spurious correspondences". It decouples scale, rotation and translation, and certifies optimality. — [Yang, Shi & Carlone 2021](https://doi.org/10.1109/TRO.2020.3033695) (arXiv 2001.07715)
- Learned registration: [Deep Closest Point (Wang & Solomon 2019)](https://doi.org/10.1109/ICCV.2019.00362); [Deep Global Registration (Choy, Dong & Koltun 2020)](https://doi.org/10.1109/CVPR42600.2020.00259) (titles and records verified).

**SLAM and odometry evaluation**
- TUM RGB-D benchmark: 39 Kinect sequences with ground truth from an 8-camera, 100 Hz motion-capture system, plus automatic evaluation tools. — [Sturm et al. 2012](https://doi.org/10.1109/IROS.2012.6385773) (3521 cites)
  - Its tools page defines the absolute trajectory error (ATE), "well-suited for measuring the performance of visual SLAM systems". Estimated poses are first associated with ground truth by timestamp. — [TUM tools page](https://cvg.cit.tum.de/data/datasets/rgbd-dataset/tools)
- Zhang & Scaramuzza (2018) show "how to determine the transformation type to use in trajectory alignment based on the specific sensing modality (i.e., monocular, stereo and visual-inertial)" before computing ATE or relative error. — [Zhang & Scaramuzza 2018](https://doi.org/10.1109/IROS.2018.8593941)
- The SLAM system record [ORB-SLAM, Mur-Artal et al. 2015](https://doi.org/10.1109/TRO.2015.2463671) was verified (6733 cites).

**Statistical shape models in vision**
- Active Shape Models: [Cootes et al. 1995](https://doi.org/10.1006/cviu.1995.1004) (5734 / 7190 cites). Active Appearance Models: [Cootes, Edwards & Taylor 2001](https://doi.org/10.1109/34.927467).
- Davies et al. (2002) identify establishing "dense correspondences between all members of a set of training shapes" as the main drawback. Manual landmarking is "time consuming and subjective in two dimensions and almost impossible in three dimensions", so they pose correspondence as minimum description length. — [Davies et al. 2002](https://doi.org/10.1109/TMI.2002.1009388)

**3D human pose (PA-MPJPE)**
- Human3.6M has 3.6 million 3D poses from 11 subjects (5 female, 6 male) in 4 viewpoints. — [Ionescu et al. 2014](https://doi.org/10.1109/TPAMI.2013.248)
- Martinez et al. (2017): "the prediction has been further aligned with the ground truth via a rigid transformation… We call this post-processing protocol #2". On HumanEva "the error is always computed after a rigid transformation". Their Table 3 (Protocol #2) gives a 52.1 mm average for their single-action model on Stacked-Hourglass detections. — [Martinez et al. 2017](https://doi.org/10.1109/ICCV.2017.288) ([arXiv 1705.03098](https://arxiv.org/abs/1705.03098))
- The same aligned-error convention is used in mesh recovery: [HMR, Kanazawa et al. 2018](https://doi.org/10.1109/CVPR.2018.00744) (record verified).

### Inferences
- These fields solve the same orthogonal-Procrustes core under different group constraints: SO(3) for attitude, SE(3) for rigid registration and stereo/VIO trajectories, Sim(3) for monocular SLAM and Umeyama scaling. They differ mainly in two things:
  1. **Where correspondences come from**: known vectors (Wahba), nearest neighbours (ICP), learned features (DCP/DGR), or hypothesis-and-verify / robust costs (RANSAC, TEASER).
  2. **How outliers are handled.**
- In evaluation, Procrustes alignment acts as a *gauge-fixing* step: it removes unobservable degrees of freedom (global frame, and scale for monocular) before measuring error. That makes the choice of transform group a reporting decision that changes the numbers. Protocol #1 vs. #2 differ in this way, as does SE(3) vs. Sim(3) ATE.

### Gaps
- Background knowledge, not verified this session:
  - Chen & Medioni (1992) introduced point-to-plane ICP.
  - ORB-SLAM computes Sim(3) loop-closure transforms with Horn's method inside RANSAC.
  - The TUM ATE tool aligns trajectories with Horn's method (the tools page text retrieved did not show this sentence).
  - Zhang & Scaramuzza recommend Sim(3) for monocular, SE(3) for stereo and 4-DoF (position + yaw) for VIO.
  - DCP uses a differentiable SVD (Procrustes) head.
  - In ASMs, training shapes are aligned by generalised Procrustes analysis before PCA builds the point distribution model. The Cootes 1995 open-access PDF could not be fetched this session.
- No verified numbers on how many point-cloud registration benchmarks report Procrustes-type errors.

---

## Q4. Sensory science and psychometrics: GPA for panel consensus, and target (Procrustes) rotation of factors

### Takeaway
Sensory science adopted Gower's GPA in the 1980s to reconcile individual assessors' configurations, especially in free-choice profiling. Each assessor's configuration is translated, rotated or reflected, and isotropically scaled to a consensus. The main methodological issue since then has been testing whether the consensus is real, using Monte Carlo and permutation tests. GPA remains in use in 2025 for cross-cultural panel comparison, alongside alternatives such as MFA and napping.

Psychometrics is the historical origin of the name "Procrustes" (Hurley & Cattell 1962). There, orthogonal and oblique target rotation (Green 1952; Schönemann 1966; Browne 1967; ten Berge 1977) remains the standard way to test factor-structure replication across cultures, e.g. the Big Five.

### Cited Findings
**Sensory science**
- Williams & Langron (1984) introduced free-choice profiling. Each assessor uses his or her own terms, the individual configurations "are rationalised by generalised Procrustes statistics", and the result is "a consensus configuration". The demonstration used ten assessors (expert and non-expert) profiling commercial ports. — [Williams & Langron 1984](https://doi.org/10.1002/jsfa.2740350513) (321 / 356 cites)
- GPA is defined by Gower (1975) as above (translation, rotation, reflection, scaling; centroid = consensus). — [Gower 1975](https://doi.org/10.1007/BF02291478). An interpretive guide for sensory GPA: [Dijksterhuis & Gower 1991](https://doi.org/10.1016/0950-3293(91)90027-C) (record verified). The standard monograph: [Gower & Dijksterhuis 2004, *Procrustes Problems*, OUP](https://doi.org/10.1093/acprof:oso/9780198510581.001.0001) (618 cites).
- Consensus testing:
  - King & Arents (1991): GPA "always produce[s] a matrix of maximal agreement", so whether consensus plots show "true consensus among panelists, as opposed to being merely artifacts of the analysis" needs a test. They proposed a Monte Carlo test (fruit-flavour data). — [King & Arents 1991](https://doi.org/10.1111/j.1745-459X.1991.tb00500.x)
  - Wakeling, Raats & MacFie (1992) permute data rows while retaining the assessor configurations. This gives "a more powerful statistical test", shown on free-choice profiling of eight milk products. — [Wakeling et al. 1992](https://doi.org/10.1111/j.1745-459X.1992.tb00526.x)
  - Later permutation tests: [Xiong et al. 2008](https://doi.org/10.1016/j.foodqual.2007.03.003) (record verified).
- Current use (2025): "Understanding cultural differences across global descriptive sensory panels using Generalized Procrustes Analysis", *Food Quality and Preference* 133:105633. — [Yu et al. 2025](https://doi.org/10.1016/j.foodqual.2025.105633) (record verified; abstract not available)
- Context and alternatives: the review of novel rapid sensory methods [Varela & Ares 2012](https://doi.org/10.1016/j.foodres.2012.06.037) (706 Crossref cites), and napping analysed by multiple factor analysis [Pagès 2005](https://doi.org/10.1016/j.foodqual.2005.01.006). FactoMineR is the common R package for these multivariate analyses. — [Lê, Josse & Husson 2008](https://doi.org/10.18637/jss.v025.i01) (9201 Crossref cites)

**Psychometrics**
- Pre-history of target rotation:
  - Mosier (1939) rotated to a known simple structure. — [Mosier 1939](https://doi.org/10.1007/BF02288493)
  - Green (1952) found the orthogonal approximation of an oblique structure. — [Green 1952](https://doi.org/10.1007/BF02288918)
- Hurley & Cattell, "The Procrustes program: producing direct rotation to test a hypothesized factor structure", *Behavioral Science* 7(2), 1962 (Crossref shows a 2007 online date). — [Hurley & Cattell 1962](https://doi.org/10.1002/bs.3830070216)
- Schönemann (1966) solved AT = B + E with T′T = I, minimising trace(E′E), including matrices of less than full column rank, generalising Green. — [Schönemann 1966](https://doi.org/10.1007/BF02289451) (1601 / 2045 cites). This is the reference SciPy cites for `orthogonal_procrustes` ([SciPy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html)).
- Extensions:
  - dilation and rigid motion ([Schönemann & Carroll 1970](https://doi.org/10.1007/BF02291266));
  - oblique Procrustes ([Browne 1967](https://doi.org/10.1007/BF02289420));
  - rotation of two or more matrices to maximal agreement, i.e. the psychometric GPA ([ten Berge 1977](https://doi.org/10.1007/BF02294053)).
- Applied landmark (cross-cultural personality):
  - McCrae et al. (1996): CFA fit indices for the NEO-PI-R were not high (n = 229). An "alternative method that uses targeted rotation… clearly replicated the 5-factor structure" in Chinese and Japanese versions. — [McCrae et al. 1996](https://doi.org/10.1037/0022-3514.70.3.552)
  - McCrae & Costa (1997): German, Portuguese, Hebrew, Chinese, Korean and Japanese samples (**N = 7,134**) were compared with the American structure. "When targeted rotations were used, the American factor structure was closely reproduced, even at the level of secondary loadings." — [McCrae & Costa 1997](https://doi.org/10.1037/0003-066X.52.5.509) (2896 Crossref cites)
- Software and tutorials: gradient-projection rotation algorithms (the R package GPArotation lineage) are in [Bernaards & Jennrich 2005](https://doi.org/10.1177/0013164404272507). A tutorial on three major classes of cross-cultural invariance testing in R: [Fischer & Karl 2019](https://doi.org/10.3389/fpsyg.2019.01507).

### Inferences
- Sensory GPA is the variant closest to Gower's original. It allows **reflection and per-assessor isotropic scaling**: assessors use scales differently and vocabularies are arbitrary. The scaling step handles "scale usage" differences between assessors, and the consensus is a mean configuration, not a reference object.
- The psychometric use is OPA (one-to-target) rather than GPA. The quantity reported is typically a congruence coefficient between rotated and target loadings, not a residual sum of squares. The "Procrustes" stigma (forcing data onto a hypothesis) originated here: Hurley & Cattell's choice of name.

### Gaps
- Not verified this session:
  - FactoMineR and SensoMineR include a `GPA` function, and XLSTAT offers GPA (common in industry).
  - Fischer & Karl 2019 include Procrustes/target rotation among their three technique classes (the abstract says only "three major classes").
  - Tucker's congruence coefficient is the standard reported statistic.
- No quantitative adoption data for GPA in sensory journals was found. Arnold & Williams (1986) and Dijksterhuis's 1996 chapter (early sensory-GPA reviews) were not located with DOIs.

---

## Q5. Ecology (and adjacent population genetics, microbiome and spatial omics): PROTEST and vegan::procrustes

### Takeaway
In community ecology, Procrustes analysis compares two ordinations or data tables on the same sites. Examples: species vs. environment, two taxonomic groups, two sampling methods, two distance metrics. Significance comes from Jackson's PROTEST permutation test (1995). Peres-Neto & Jackson (2001) argued it is preferable to the Mantel test. The default implementation is `vegan::procrustes()` / `protest()`, which reports a symmetric Procrustes statistic and a correlation-like r = √(1 − ss) with permutation p-values.

The same "Procrustes concordance" logic spread to:
- microbiome studies (16S vs. metagenome ordinations);
- human population genetics (genetic PCA maps vs. geography; projecting samples into reference PCA spaces);
- spatial transcriptomics (PASTE aligns tissue slices by weighted Procrustes).

### Cited Findings
- **PROTEST** (Jackson 1995): one matrix is subjected to "reflection, rigid rotation, translation, and dilation to minimize the sum of the squared residual deviations" to the target. A randomization test asks whether the residual sum is smaller than expected by chance. It works on raw data or on ordination summaries (direct or indirect gradient analysis). — [Jackson 1995](https://doi.org/10.1080/11956860.1995.11682297) (421 / 479 cites)
- Peres-Neto & Jackson (2001), "How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test", *Oecologia* 129:169–178. — [Peres-Neto & Jackson 2001](https://doi.org/10.1007/s004420100720) (915 / 985 cites; abstract not retrievable)
- Lisboa et al. (2014) promote the **Procrustean association metric**, the vector of per-observation Procrustes residuals, to estimate *local* levels of association between tables. In the Mantel framework, converting to pairwise distance matrices loses "the ability to contrast homologous data points across dimensions". — [Lisboa et al. 2014](https://doi.org/10.1371/journal.pone.0101238)
- **vegan**: `procrustes(X, Y, scale = TRUE, symmetric = FALSE, ...)` "rotates a configuration to maximum similarity with another configuration". It is "typically used in comparison of ordination results".
  - With `symmetric = TRUE` both configurations are scaled to equal dispersions.
  - `protest(X, Y, permutations = how(nperm = 999))` "uses a correlation-like statistic derived from the symmetric Procrustes sum of squares ss as r = √(1−ss), and also prints the sum of squares of the symmetric analysis, sometimes called m12²".

  — [vegan reference](https://vegandevs.github.io/vegan/reference/procrustes.html). Current CRAN vegan **2.7-6 (2026-08-31)**; the pkgdown site shows a 2.8-0 development version. — [CRAN vegan](https://cran.r-project.org/package=vegan)
- Microbiome:
  - Muegge et al. (2011) combined shotgun metagenomes and 16S rRNA data from 33 mammalian species and 18 humans. — [Muegge et al. 2011](https://doi.org/10.1126/science.1198719)
  - EMPeror was built to visualise ordinations for the Earth Microbiome Project, profiling "tens of thousands of samples". — [Vázquez-Baeza et al. 2013](https://doi.org/10.1186/2047-217X-2-16)
  - QIIME 2: [Bolyen et al. 2019](https://doi.org/10.1038/s41587-019-0209-9)
- Population genetics:
  - Wang et al. (2010) describe "a Procrustes analysis approach for quantitatively assessing the similarity of population-genetic and geographic maps". They confirm high concordance for SNP data from Europe and worldwide, PCA vs. MDS concordance, and CNV vs. SNP maps. — [Wang et al. 2010](https://doi.org/10.2202/1544-6115.1493) ([PMC2861313](https://pmc.ncbi.nlm.nih.gov/articles/PMC2861313/))
  - Wang, Zöllner & Rosenberg (2012) use genome-wide SNPs from **128 populations** and apply Procrustes analysis per geographic region to compare PCA maps with geography. — [Wang et al. 2012](https://doi.org/10.1371/journal.pgen.1002886)
  - The motivating observation is [Novembre et al. 2008, "Genes mirror geography within Europe"](https://doi.org/10.1038/nature07331). Procrustes projection of samples into a reference ancestry space is used for sequencing data in [Wang et al. 2014 (LASER)](https://doi.org/10.1038/ng.2924) (record verified).
- **Spatial transcriptomics (2022)**: PASTE stacks adjacent tissue slices. It uses optimal-transport mappings Π^(k) "to solve a generalized weighted Procrustes problem", finding R ∈ R^{2×2} with RᵀR = I and v minimising Σ_ij π_ij ‖z_i^(k) − R z_j^(k+1) − v‖², solved by SVD. — [Zeira et al. 2022](https://doi.org/10.1038/s41592-022-01459-6) ([PMC9334025](https://pmc.ncbi.nlm.nih.gov/articles/PMC9334025/))

### Inferences
- Ecology uses Procrustes not to analyse shape but as a **matrix-association statistic** (a multivariate correlation between two configurations of the same sites) plus residual diagnostics. The usual variant is OPA with reflection allowed and symmetric scaling. Sign and orientation of ordination axes are arbitrary, so reflection must be allowed.
- PASTE shows a modern pattern: soft correspondences (OT couplings) become weights in a weighted Procrustes fit. This is the same structure as ICP with soft assignments and as Wasserstein-Procrustes in NLP (Q6).

### Gaps
- Background knowledge, not verified this session:
  - Muegge et al. specifically used Procrustes to compare 16S- and metagenome-based PCoA.
  - EMPeror and QIIME 2 provide Procrustes plots (q2-diversity `procrustes-analysis`).
  - Wang et al.'s similarity statistic t0 and their permutation-test details.
- Peres-Neto & Jackson 2001's specific power comparisons with the Mantel test could not be read (paywalled; abstract not in OpenAlex).

---

## Q6. Machine learning and neuroscience: cross-lingual embedding alignment, representation comparison, hyperalignment

### Takeaway
Orthogonal Procrustes became a workhorse in NLP after 2015:
- Xing et al. (2015) constrained the Mikolov et al. (2013) linear mapping between embedding spaces to be orthogonal.
- Smith et al. (2017) argued that a self-consistent map should be orthogonal and solved it by SVD.
- MUSE (Conneau et al. 2018) used adversarial initialisation followed by iterated Procrustes refinement with CSLS retrieval.
- VecMap (Artetxe et al. 2016–2018) used Procrustes inside a robust self-learning loop.

Its key assumption, approximate isomorphism of embedding spaces, fails for distant languages and mismatched domains (Søgaard et al. 2018). In representation-similarity research, the Procrustes distance was rehabilitated:
- Williams et al. (2021) cast neural-representation comparison as statistical shape analysis with proper metrics.
- Ding et al. (2021) found orthogonal Procrustes distance "balances gracefully" between CKA and CCA.

In fMRI, hyperalignment (Haxby et al. 2011) is an iterative GPA-like orthogonal Procrustes alignment of subjects' voxel spaces into a common model space. It raised between-subject classification from 32% to 70.6%.

### Cited Findings
**Cross-lingual embeddings**
- Mikolov, Le & Sutskever (2013), "Exploiting similarities among languages for machine translation", introduced linear mapping between monolingual embedding spaces. — [arXiv 1309.4168](https://arxiv.org/abs/1309.4168)
- Xing et al. (2015) found inconsistencies between the embedding objective, the transform objective and the distance measure. They normalise word vectors onto a hypersphere and constrain the linear transform to be orthogonal, improving word similarity and English→Spanish translation. — [Xing et al. 2015](https://doi.org/10.3115/v1/N15-1104)
- Smith et al. (2017):
  - "a self-consistent linear transformation between vector spaces should be orthogonal… found using the singular value decomposition";
  - training dictionary of 5,000 pairs;
  - the inverted softmax raises precision@1 of Mikolov's mapping "from 34% to 43%" for English→Italian;
  - a pseudo-dictionary of identical character strings gives 40%;
  - sentence retrieval from 200k Italian sentences reaches P@1 68%.

  — [Smith et al. 2017, ICLR](https://arxiv.org/abs/1702.03859)
- Conneau et al. (2018), MUSE:
  - With an orthogonality constraint the mapping problem "boils down to the Procrustes problem" with closed form W* = UVᵀ, where UΣVᵀ = SVD(YXᵀ).
  - Pipeline: adversarial initialisation, then a synthetic dictionary of mutual nearest neighbours among frequent words, then Procrustes refinement, applied iteratively.
  - Further iterations usually improve word translation by "below 1%".
  - CSLS reduces hubness.

  — [Conneau et al. 2018, arXiv 1710.04087](https://arxiv.org/abs/1710.04087). Hubness background: [Dinu, Lazaridou & Baroni 2015](https://arxiv.org/abs/1412.6568).
- Artetxe, Labaka & Agirre:
  - orthogonal mappings preserving monolingual invariance ([EMNLP 2016](https://doi.org/10.18653/v1/D16-1250));
  - self-learning from almost no bilingual data ([ACL 2017](https://doi.org/10.18653/v1/P17-1042));
  - a fully unsupervised initialisation exploiting structural similarity plus robust self-learning, which succeeds where adversarial methods "often fail in more realistic scenarios" ([ACL 2018](https://doi.org/10.18653/v1/P18-1073)).
- Successors relaxing or replacing Procrustes:
  - RCSLS ([Joulin et al. 2018](https://doi.org/10.18653/v1/D18-1330));
  - Gromov-Wasserstein, which casts correspondence as optimal transport ([Alvarez-Melis & Jaakkola 2018](https://doi.org/10.18653/v1/D18-1214));
  - Wasserstein Procrustes ([Grave, Joulin & Berthet, arXiv 1805.11222](https://arxiv.org/abs/1805.11222)).
- Limits and evaluation:
  - Unsupervised alignment "performs much worse on morphologically rich languages that are not dependent marking, when monolingual corpora from different domains or different embedding algorithms are used". Identical words provide a robust weak signal. — [Søgaard, Ruder & Vulić 2018](https://doi.org/10.18653/v1/P18-1072)
  - Performance "largely depends on the task at hand", and optimising for bilingual lexicon induction "can result in deteriorated downstream performance". — [Glavaš et al. 2019](https://doi.org/10.18653/v1/P19-1070)
  - Survey: [Ruder, Vulić & Søgaard 2019, JAIR 65](https://doi.org/10.1613/jair.1.11640)
- Diachronic semantics: six historical corpora in four languages over two centuries; laws of conformity and innovation. — [Hamilton, Leskovec & Jurafsky 2016](https://doi.org/10.18653/v1/P16-1141)
- 2025 frontier: unsupervised translation between text-embedding spaces of different models, "Harnessing the universal geometry of embeddings". — [Jha et al. 2025, arXiv 2505.12540](https://arxiv.org/abs/2505.12540) (record verified; content not read)

**Representation comparison**
- Williams et al. (2021), "Generalized shape metrics on neural representations" (NeurIPS 2021):
  - They draw on statistical shape analysis to define dissimilarities that are proper metrics (symmetric, satisfying the triangle inequality), enabling k-NN and hierarchical clustering over "thousands of networks".
  - Their d₁(X_i, X_j) = min_{Q∈O} ‖X_i^φ − X_j^φ Q‖ "is equivalent to the Procrustes size-and-shape distance with reflections".
  - The angular version "is closely related to the Riemannian distance on Kendall's shape space".
  - Both are solved as an orthogonal Procrustes problem in O(n³).

  — [Williams et al. 2021, arXiv 2110.14739](https://arxiv.org/abs/2110.14739)
- Ding, Denain & Steinhardt (2021):
  - They compare CKA, three CCA-based measures and orthogonal Procrustes.
  - d_Proc(A, B) = ‖A‖²_F + ‖B‖²_F − 2‖AᵀB‖_* (nuclear norm). It is not normalised, but lies in [0, 2] for normalised inputs.
  - "The Orthogonal Procrustes distance, a classical but often overlooked dissimilarity measure, balances gracefully between CKA and CCA and consistently performs well."

  — [Ding et al. 2021, arXiv 2108.01661](https://arxiv.org/abs/2108.01661)
- Related works:
  - CKA: [Kornblith et al. 2019](https://arxiv.org/abs/1905.00414)
  - stochastic-network metrics: [Duong et al. 2023, ICLR](https://arxiv.org/abs/2211.11665)
  - decodability interpretation: [Harvey, Lipshutz & Williams 2024](https://arxiv.org/abs/2411.08197)
  - dynamical similarity analysis: [Ostrow et al. 2023](https://arxiv.org/abs/2306.10168)
  - survey: [Klabunde et al. 2025, ACM Computing Surveys 57(9):242](https://doi.org/10.1145/3728458)

**Neuroscience (hyperalignment)**
- Haxby et al. (2011):
  - "The non-reference subjects were iteratively chosen and their time-point vectors were aligned to the time-point vectors of the current reference using the Procrustean transformation (procrustes as implemented in Matlab)." The reference is updated by averaging, which makes this a GPA-like scheme.
  - Between-subject classification (BSC) of movie time segments reached **70.6% (SE 2.6%, chance < 1%) after hyperalignment vs. 32.0% (SE 2.5%) after anatomical alignment**.
  - A 35-dimensional common model matched the full 1000 dimensions: 68.3% vs. 70.6% for movie segments, 64.8% vs. 63.9% for faces and objects, 67.6% vs. 68.0% for animal species.

  — [Haxby et al. 2011](https://doi.org/10.1016/j.neuron.2011.08.026) ([PMC3201764](https://pmc.ncbi.nlm.nih.gov/articles/PMC3201764/))
- Whole-cortex searchlight hyperalignment: [Guntupalli et al. 2016](https://doi.org/10.1093/cercor/bhw068). Conceptual framework: individual transformation matrices project into a common information space "preserving the geometry of pairwise dissimilarities". — [Haxby et al. 2020, eLife](https://doi.org/10.7554/eLife.56601). Toolbox: [PyMVPA, Hanke et al. 2009](https://doi.org/10.1007/s12021-008-9041-y).

**Software**
- `scipy.linalg.orthogonal_procrustes(A, B)` minimises ‖AR − B‖_F subject to RᴴR = I. It returns R and `scale` (the sum of singular values of AᴴB). It uses "only orthogonal transformations like rotations and reflections, and it does not use scaling or translation". It cites Schönemann 1966 and has experimental Array-API (CuPy/PyTorch/JAX) support. — [SciPy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html)
- `scipy.spatial.procrustes(data1, data2)`:
  - standardises both inputs (centred, tr(AAᵀ) = 1), then fits data2 to data1 "including scaling/dilation, rotations, and reflections";
  - returns `disparity` M² = Σ(data1 − data2)²;
  - is not designed for different numbers of rows; for different numbers of columns, zero-pad;
  - "the disparity scales as the number of points per input matrix".

  — [SciPy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html); [Virtanen et al. 2020](https://doi.org/10.1038/s41592-019-0686-2)

### Inferences
- ML predominantly uses OPA in O(d) without scaling, after centring and length-normalising. Reflections are usually permitted because orientation is meaningless in embedding spaces.
- In BLI the Procrustes fit is rarely the bottleneck. Two things are: (a) obtaining correspondences (seed dictionary, adversarial, OT or self-learning), and (b) retrieval in the mapped space (hubness, CSLS). The same "correspondence is the hard part" lesson appears in ICP (Q3) and in SSMs (Q7).
- In representation comparison, Procrustes became attractive because the shape-analysis view yields **metrics**, and because empirically it sits between CKA's insensitivity and CCA's oversensitivity. These are arguably the two main 2021–2025 developments.

### Gaps
- Background knowledge, not verified this session:
  - Hamilton et al. (2016) align decade-wise embeddings with orthogonal Procrustes (the abstract retrieved did not mention it).
  - Ostrow et al. (2023) use a Procrustes-over-vector-fields formulation.
  - Grave et al. (2019) appeared at AISTATS 2019 (PMLR 89).
- Typical BLI numbers (e.g., MUSE en–es P@1 with Procrustes + CSLS) were not extracted from tables this session. Only the Smith et al. 34% → 43% figure is verified.

---

## Q7. Other notable uses: geography and cartography, geodesy, archaeology, medical-imaging shape models, gait, face and biomechanics

### Takeaway
- **Geography**: Tobler's map-comparison work (1965) and bidimensional regression (1994) are Procrustes and regression analogues for comparing maps, including old vs. modern maps and cognitive maps.
- **Cartography**: MapAnalyst-style distortion analysis of historical maps.
- **Geodesy**: the 7-parameter Helmert similarity transformation between datums is exactly a Procrustes problem, with closed-form Procrustes solutions (and 9-parameter and weighted extensions) that need no initial values or iteration.
- **Archaeology and statistics**: Kendall's shape theory grew partly out of testing whether standing-stone and ley-line alignments exceed chance, i.e. the distribution of "flat" triangles.
- **Medical imaging**: statistical shape models (ASM, MDL, GPMM, ShapeWorks) depend on Procrustes alignment plus correspondence.
- **Biometrics and biomechanics**: Procrustes mean shapes for gait recognition, SVD rigid-body fits for marker-based motion capture, and morphable face models.

### Cited Findings
**Geography and cartography**
- Tobler (1965), "Computation of the correspondence of geographical patterns". — [Tobler 1965](https://doi.org/10.1007/BF01947869)
- Tobler (1994): bidimensional regression allows "geographic maps [to] be compared with each other", including "ancient and modern maps", "mental maps" and GIS rubber-sheeting. It also applies to shapes of "leaves, fish, faces, or skulls after the manner of D'Arcy W. Thompson". The linear case yields "a Pearsonian-like correlation coefficient". — [Tobler 1994](https://doi.org/10.1111/j.1538-4632.1994.tb00320.x)
- Statistical inference for bidimensional regression: [Nakaya 1997](https://doi.org/10.1111/j.1538-4632.1997.tb00954.x). Psychology adoption for cognitive maps, with a new distortion index and Euclidean vs. non-Euclidean geometries: [Friedman & Kohler 2003](https://doi.org/10.1037/1082-989X.8.4.468).
- Historical-map accuracy: MapAnalyst ([Jenny, Weber & Hurni 2007](https://doi.org/10.3138/carto-v42-1-089)); distortion analysis and visualisation ([Jenny & Hurni 2011](https://doi.org/10.1016/j.cag.2011.01.005)) (records verified).

**Geodesy**
- "The Procrustes method is a very effective method for determining the Helmert's datum transformation parameters since it requires neither initial starting values nor iteration." Extending ABC-Procrustes to 3D affine transformations with anisotropic scale works only for mild anisotropy. The example was AGD84 → GDA94 (Australia). — [Paláncz et al. 2010](https://doi.org/10.5047/eps.2010.10.004)
- Further work:
  - 9-parameter Procrustes: [Awange, Bae & Claessens 2008](https://doi.org/10.1186/BF03353115)
  - nonlinear analysis of the conformal-group C7(3) datum transformation: [Grafarend & Awange 2003](https://doi.org/10.1007/s00190-002-0299-9)
  - Procrustes for transformation parameters between ITRS realisations in Turkey: [Yetkin & Ansari 2017](https://doi.org/10.1515/jag-2016-0048)
- Ligas (2020) gives point-wise weighted Procrustes solutions for similarity transformations under Gauss–Markov (errors in one system, closed form) and Gauss–Helmert / errors-in-variables (errors in both systems; generally iterative) models. — [Ligas 2020](https://doi.org/10.1007/s12518-020-00306-7)
- Reference-frame context (state of the art): ITRF2020 combines DORIS, GNSS, SLR and VLBI solutions spanning 28, 27, 38 and 41 years of observations. — [Altamimi et al. 2023](https://doi.org/10.1007/s00190-023-01738-w)

**Archaeology and the origin of shape theory**
- Kendall & Kendall (1980) count "ε-blunt" triads (triangles with largest angle ≥ π − ε) among n random points. They introduce the pontogram and simulation tests "by random lateral perturbations of the data" to test for too many alignments. — [Kendall & Kendall 1980](https://doi.org/10.2307/1426603)
- The problem was posed by Broadbent's "Simulating the ley hunter". — [Broadbent 1980](https://doi.org/10.2307/2981985)
- Kendall (1989) lists applications of shape theory "to archeology, astronomy, geography and physical chemistry". — [Kendall 1989](https://doi.org/10.1214/ss/1177012582)
- Secondary sources describe Kendall's Land's End (Cornwall) standing-stones analysis as testing near-collinear ("flat") triangles on the shape sphere within a convex region bounded by the coastline. They also connect this to Barbour-type relational mechanics. — [Anderson, arXiv 1307.1923](https://arxiv.org/abs/1307.1923); [Anderson, arXiv 1505.00488](https://arxiv.org/abs/1505.00488) (seen via search summaries only)

**Medical imaging**
- Statistical shape models for 3D segmentation: review by [Heimann & Meinzer 2009](https://doi.org/10.1016/j.media.2009.05.004) (1248 cites).
- The automatic-correspondence problem: [Davies et al. 2002](https://doi.org/10.1109/TMI.2002.1009388), see Q3.
- Gaussian process morphable models: [Lüthi et al. 2018](https://doi.org/10.1109/TPAMI.2017.2739743)
- Overview chapter: [Ambellan et al. 2019](https://doi.org/10.1007/978-3-030-19385-0_5)
- Particle-based correspondence software ShapeWorks: [Cates, Elhabian & Whitaker 2017](https://doi.org/10.1016/B978-0-12-810493-4.00012-2)

**Face, gait and biomechanics**
- Face: the original 3D morphable face model is [Blanz & Vetter 1999](https://doi.org/10.1145/311535.311556) (3511 / 4710 cites). The 20-year survey is [Egger et al. 2020](https://doi.org/10.1145/3395208).
- Gait: silhouettes over a gait cycle are represented as complex configurations, "analyzed using the Procrustes shape analysis method to obtain mean shape as gait signature". Recognition uses "the full Procrustes distance measure". — [Wang et al. 2003](https://doi.org/10.1109/TIP.2003.815251)
- Biomechanics and motion capture: [Söderkvist & Wedin 1993](https://doi.org/10.1016/0021-9290(93)90098-Y) (621 / 734 cites) and [Challis 1995](https://doi.org/10.1016/0021-9290(94)00116-L) are the standard rigid-body-from-markers references (records verified).

### Inferences
- Geodesy and cartography use OPA with a similarity transform (Sim(2)/Sim(3)), and treat residuals as distortion or error fields to be mapped. This is the same as Lisboa et al.'s Procrustean residual vectors in ecology and as GM's landmark-wise residual visualisations.
- Kendall's archaeological motivation is directly relevant to this repository's triangle and shape-sphere work: flat (collinear) triangles form a distinguished locus on the shape sphere, and their statistics motivated the theory.

### Gaps
- Background knowledge, not verified this session:
  - ITRF2020 transformations use 14-parameter (7 + rates) Helmert similarity transformations.
  - MapAnalyst uses Helmert and affine transformations.
  - Söderkvist & Wedin and Challis use the SVD (Procrustes) solution (only the titles were verified).
  - Face-recognition pipelines (e.g., ArcFace) align faces with a 5-point similarity (Umeyama) transform.
- The Anderson arXiv sources were not bibliographically verified. The original Kendall Land's End data analysis (likely in Kendall & Kendall 1980 or a later Kendall paper) was not read this session.

---

## Q8. Cross-field pitfalls: the Pinocchio effect, outliers and landmark choice, reflections, correspondence assumptions, and metric reporting

### Takeaway
The pitfalls recur across fields because they come from the least-squares Procrustes criterion itself.
1. **Localised variation is smeared over all points** (the "Pinocchio effect" in GM; heteroscedasticity in structural biology). Remedies: resistant fits, ML superposition (THESEUS), robust/TLS costs (TEASER).
2. **Results depend on which points define the fit.** Examples: landmark or semilandmark choice, the atoms or pocket used for alignment, the trajectory alignment group.
3. **Reflection handling must match the domain**: forbidden for molecules and rigid bodies, allowed for ordinations, embeddings and sensory configurations.
4. **Known correspondence is assumed.** Where it is absent, it is the hard part: ICP, MDL shape models, symmetric ligands, BLI seed dictionaries.
5. **Raw Procrustes residuals or RMSD are size- and sample-dependent and can be produced by chance.** Hence normalised scores (TM-score, nRMSD), permutation tests (PROTEST, sensory consensus tests, RRPP), and caution about global-superposition metrics.

### Cited Findings
**Pinocchio effect and outliers**
- Using GPA to assess landmark repeatability suffers from "the 'Pinocchio effect'", in which variation at one landmark is spread to others by the superimposition. Alternatives (centroid-distance based) are also potentially misleading. — [von Cramon-Taubadel, Frazier & Lahr 2007](https://doi.org/10.1002/ajpa.20616). Measurement-error review: [Fruciano 2016](https://doi.org/10.1007/s00427-016-0537-4).
- Errors in estimating the "nuisance" parameters (translation, rotation, scale) "will be reflected in the covariance structure of the coordinates". Generalized resistant fit vs. relative-warp performance depends on landmarks with excessive variance (simulations with 6–20 landmarks). — [Walker 2000](https://doi.org/10.1080/106351500750049770)
- Resistant (robust) superimposition, generalised to many specimens: [Rohlf & Slice 1990](https://doi.org/10.2307/2992207); [Siegel & Benson 1982](https://doi.org/10.2307/2530448).
- In structural biology, LS superposition "can give erroneous solutions" under heteroscedastic or correlated errors. ML superposition down-weights variable regions and is "insensitive to the specific atoms included". — [Theobald & Wuttke 2006](https://doi.org/10.1073/pnas.0508445103); [THESEUS](https://doi.org/10.1093/bioinformatics/btl332)
- Global Cα-superposition scores are "strongly influenced by domain motions". — [Mariani et al. 2013](https://doi.org/10.1093/bioinformatics/btt473)
- Robust registration with large outlier fractions requires non-LS costs, e.g. truncated least squares. — [Yang et al. 2021](https://doi.org/10.1109/TRO.2020.3033695)

**Landmark choice and density**
- The increasing number of landmarks and semilandmarks is a current methodological challenge. — [Mitteroecker & Schaefer 2022](https://doi.org/10.1002/ajpa.24531); critique of "high-density" morphometrics in [Cardini 2020](https://doi.org/10.1007/s00435-020-00499-w) (title verified)
- geomorph lets users choose which points drive the rotation ("fixed", "sliders" or a custom vector). — [gpagen docs](https://rdrr.io/cran/geomorph/man/gpagen.html)

**Reflections**
- Umeyama (1991): earlier SVD and quaternion solutions can return a reflection instead of a rotation when data are corrupted; his strict solution always returns a proper rotation. — [Umeyama 1991](https://doi.org/10.1109/34.88573)
- Defaults differ by toolkit and field:
  - Gower's GPA explicitly includes reflection ([Gower 1975](https://doi.org/10.1007/BF02291478)), as does PROTEST ([Jackson 1995](https://doi.org/10.1080/11956860.1995.11682297)).
  - `scipy.spatial.procrustes` includes reflections ([docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html)), as does `scipy.linalg.orthogonal_procrustes` ([docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html)).
  - Williams et al.'s representational metric is the size-and-shape distance "with reflections" ([arXiv 2110.14739](https://arxiv.org/abs/2110.14739)).
  - `Rotation.align_vectors` (Kabsch) returns rotations ([docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html)).
- Symmetric molecules can be optimally superimposable on their mirror image, giving degenerate superpositions. — [Coutsias & Wester 2019](https://doi.org/10.1002/jcc.25802)
- GM's object-symmetry analysis deliberately superimposes a reflected copy. — [Klingenberg et al. 2002](https://doi.org/10.1111/j.0014-3820.2002.tb00117.x)

**Correspondence**
- ICP converges only to the nearest local minimum and needs adequate initial poses. — [Besl & McKay 1992](https://doi.org/10.1109/34.121791)
- Dense correspondence for shape models is "almost impossible in three dimensions" by hand. — [Davies et al. 2002](https://doi.org/10.1109/TMI.2002.1009388)
- Symmetric ligands: naive correspondence inflates RMSD. — [Bell & Zhang 2019](https://doi.org/10.1186/s13321-019-0362-7)
- Embedding alignment assumes near-isomorphic spaces, which breaks for typologically distant languages, different domains or different embedding algorithms. — [Søgaard et al. 2018](https://doi.org/10.18653/v1/P18-1072)
- `scipy.spatial.procrustes` requires equal row counts, i.e. known one-to-one correspondence. — [docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html)

**Metric reporting and chance agreement**
- RMSD depends on the number of aligned atoms, which motivated nRMSD. — [Carugo & Pongor 2001](https://doi.org/10.1110/ps.690101). TM-score removes size dependence. — [Zhang & Skolnick 2004](https://doi.org/10.1002/prot.20264)
- RMSD alone is insufficient for docking evaluation. — [Buttenschoen et al. 2024](https://doi.org/10.1039/D3SC04185A)
- SciPy disparity "scales as the number of points". — [docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html). The orthogonal Procrustes distance is unnormalised unless inputs are normalised. — [Ding et al. 2021](https://arxiv.org/abs/2108.01661)
- GPA always produces a consensus, so it must be tested: [King & Arents 1991](https://doi.org/10.1111/j.1745-459X.1991.tb00500.x); [Wakeling et al. 1992](https://doi.org/10.1111/j.1745-459X.1992.tb00526.x).
- Ecological concordance is tested by permutation (PROTEST, vegan default 999 permutations). — [vegan docs](https://vegandevs.github.io/vegan/reference/procrustes.html)
- In GM, Procrustes ANOVA uses residual randomisation. — [Collyer & Adams 2018](https://doi.org/10.1111/2041-210X.13029)

**Synthesis table: variant by field** (each cell is supported by the sources cited in Q1–Q7; entries marked † are background knowledge not verified this session)

| Field | Typical Procrustes variant | Transform allowed | Reflection | Standard software | What is reported |
|---|---|---|---|---|---|
| Geometric morphometrics | GPA (iterative, to mean), sliding semilandmarks | translation + rotation + scaling to unit centroid size | no (except deliberate reflected copies for symmetry) | geomorph/RRPP, MorphoJ, tps, Morpho, shapes | Procrustes coordinates/distances, PCA, Procrustes ANOVA (permutation) |
| Structural biology / MD | two-set rigid superposition (Kabsch, quaternion, QCP); ML variant (THESEUS) | rotation + translation | no | MDAnalysis (QCP), THESEUS, TM-align/CE/FATCAT (RCSB) | RMSD in Å with N aligned atoms or residues; TM-score; ligand RMSD < 2 Å success rate |
| Aerospace attitude | Wahba's problem (weighted unit vectors) | rotation only | no | QUEST/TRIAD; SciPy `align_vectors` | attitude estimate, covariance; rssd |
| Vision / robotics | absolute orientation inside ICP / RANSAC / TLS | SE(3) or Sim(3) (Umeyama) | no | TEASER; ICP variants; benchmark tools (TUM) | ATE after alignment; PA-MPJPE (Protocol #2) |
| Sensory science | GPA of assessor configurations | translation + rotation + isotropic scaling | yes | GPA implementations† (FactoMineR, XLSTAT) | consensus configuration; residual and consensus tests (Monte Carlo, permutation) |
| Psychometrics | target (Procrustes) rotation, orthogonal or oblique | rotation (or oblique transform) | yes (orthogonal) | GPArotation lineage | replicated loadings; congruence coefficients† |
| Ecology | OPA plus PROTEST | translation + rotation + dilation | yes | vegan `procrustes`/`protest` | m12², r = √(1−ss), permutation p |
| Population genetics / spatial omics | OPA of PCA maps vs. geography; weighted Procrustes with OT weights (PASTE) | similarity / rigid | — | PASTE; LASER/TRACE† | Procrustes similarity; aligned 3D tissue stacks |
| NLP / representation similarity | orthogonal Procrustes | orthogonal only (after centring and normalising) | usually yes | `scipy.linalg.orthogonal_procrustes`, MUSE, VecMap | BLI precision@k; Procrustes or angular shape distance |
| fMRI | iterative Procrustes to an evolving group mean (hyperalignment) | orthogonal | (Matlab `procrustes` default) | PyMVPA | between-subject classification accuracy |
| Geodesy / cartography | Helmert similarity, Procrustes closed form; bidimensional regression | Sim(3) / Sim(2) (also affine extensions) | no | MapAnalyst† | transformation parameters; residual and distortion fields; bidimensional r |

### Inferences
- The fields fall into two families:
  1. **Gauge-removal users** (structural biology, robotics, aerospace, geodesy, evaluation metrics). They want the rigid or similarity transform, or the residual after removing it, and treat reflection as an error.
  2. **Configuration-comparison users** (GM, sensory science, ecology, psychometrics, ML representation analysis). They want a consensus or a similarity statistic, and often allow reflection because axes are arbitrary.

  Most practical pitfalls come from using one family's defaults (e.g., SciPy's reflection-allowing `procrustes`) in the other family's problem.
- The modern trend across fields is the same move away from plain least squares, in three forms:
  - robust or probabilistic weighting (THESEUS, TLS/TEASER, resistant fit);
  - soft correspondences fed into weighted Procrustes (PASTE, Wasserstein-Procrustes, ICP variants);
  - normalised, permutation-tested or superposition-free summaries (TM-score, lDDT, RRPP, PROTEST).

### Gaps
- The origin of the term "Pinocchio effect" (usually attributed to Chapman 1990, Proceedings of the Michigan Morphometrics Workshop) could not be verified. No DOI exists for it, and the attribution here rests on secondary use.
- No cross-field quantitative survey of how often reflections are mishandled was found.
- Default reflection behaviour in MATLAB's `procrustes` (used by Haxby 2011) was not checked this session.

---

## Source catalogue

Key conventions:
- **verified = yes**: the bibliographic record (Crossref, doi.org metadata, arXiv API or CRAN) was seen this session. Any field filled from memory is noted in the "verified" cell.
- **verified = no**: the record is from memory or search snippets only.
- Open-access URLs are legal sources only: arXiv, PMC, ACL Anthology, institutional repositories, or publisher OA copies listed by OpenAlex.
- Long author lists are truncated with "et al." and the Crossref author count in parentheses.

| key | authors | year | title | venue/volume/pages | DOI | arXiv ID | open-access PDF URL | verified |
|---|---|---|---|---|---|---|---|---|
| rohlf1990extensions | Rohlf, F. James; Slice, Dennis | 1990 | Extensions of the Procrustes method for the optimal superimposition of landmarks | Systematic Zoology 39(1):40–59 | 10.2307/2992207 | | | yes (end page 59 from memory; Crossref gives first page only) |
| siegel1982robust | Siegel, Andrew F.; Benson, Richard H. | 1982 | A robust comparison of biological shapes | Biometrics 38(2):341–350 | 10.2307/2530448 | | | yes (end page from memory) |
| gower1975generalized | Gower, J. C. | 1975 | Generalized Procrustes analysis | Psychometrika 40(1):33–51 | 10.1007/BF02291478 | | | yes |
| goodall1991procrustes | Goodall, Colin | 1991 | Procrustes methods in the statistical analysis of shape | Journal of the Royal Statistical Society: Series B 53(2):285–321 | 10.1111/j.2517-6161.1991.tb01825.x | | | yes |
| kendall1984shape | Kendall, David G. | 1984 | Shape manifolds, Procrustean metrics, and complex projective spaces | Bulletin of the London Mathematical Society 16(2):81–121 | 10.1112/blms/16.2.81 | | | yes |
| kendall1989survey | Kendall, David G. | 1989 | A survey of the statistical theory of shape | Statistical Science 4(2):87–99 | 10.1214/ss/1177012582 | | https://projecteuclid.org/journals/statistical-science/volume-4/issue-2/A-Survey-of-the-Statistical-Theory-of-Shape/10.1214/ss/1177012582.pdf | yes |
| kendall1980alignments | Kendall, David G.; Kendall, Wilfrid S. | 1980 | Alignments in two-dimensional random sets of points | Advances in Applied Probability 12(2):380–424 | 10.2307/1426603 | | | yes (OpenAlex also lists DOI 10.1017/s0001867800050230) |
| broadbent1980simulating | Broadbent, Simon | 1980 | Simulating the ley hunter | Journal of the Royal Statistical Society: Series A (General) 143(2):109–140 | 10.2307/2981985 | | | yes (end page from memory) |
| bookstein1989principal | Bookstein, F. L. | 1989 | Principal warps: thin-plate splines and the decomposition of deformations | IEEE Transactions on Pattern Analysis and Machine Intelligence 11(6):567–585 | 10.1109/34.24792 | | | yes |
| bookstein1991morphometric | Bookstein, Fred L. | 1991 | Morphometric Tools for Landmark Data: Geometry and Biology | Cambridge University Press | 10.1017/CBO9780511573064 | | | yes (Crossref online date 1992; subtitle from memory) |
| rohlf1993revolution | Rohlf, F. James; Marcus, Leslie F. | 1993 | A revolution in morphometrics | Trends in Ecology & Evolution 8(4):129–132 | 10.1016/0169-5347(93)90024-J | | | yes (Crossref title string garbled as "A revolution morphometrics") |
| rohlf1999shape | Rohlf, F. James | 1999 | Shape statistics: Procrustes superimpositions and tangent spaces | Journal of Classification 16(2):197–223 | 10.1007/s003579900054 | | | yes |
| adams2004geometric | Adams, Dean C.; Rohlf, F. James; Slice, Dennis E. | 2004 | Geometric morphometrics: ten years of progress following the 'revolution' | Italian Journal of Zoology 71(1):5–16 | 10.1080/11250000409356545 | | https://www.tandfonline.com/doi/pdf/10.1080/11250000409356545 | yes |
| adams2013field | Adams, Dean C.; Rohlf, F. James; Slice, Dennis E. | 2013 | A field comes of age: geometric morphometrics in the 21st century | Hystrix, the Italian Journal of Mammalogy 24(1):7–14 | 10.4404/hystrix-24.1-6283 | | https://lib.dr.iastate.edu/eeob_ag_pubs/121 (repository landing page) | yes (doi.org metadata; pages from memory) |
| zelditch2012geometric | Zelditch, Miriam L.; Swiderski, Donald L.; Sheets, H. David | 2012 | Geometric Morphometrics for Biologists: A Primer (2nd ed.) | Academic Press (Elsevier), ISBN 978-0-12-386903-6 | | | | no (book-level record not found; ISBN inferred from 2012 chapter DOI 10.1016/b978-0-12-386903-6.00014-9) |
| klingenberg2011morphoj | Klingenberg, Christian Peter | 2011 | MorphoJ: an integrated software package for geometric morphometrics | Molecular Ecology Resources 11(2):353–357 | 10.1111/j.1755-0998.2010.02924.x | | | yes |
| adams2013geomorph | Adams, Dean C.; Otárola-Castillo, Erik | 2013 | geomorph: an R package for the collection and analysis of geometric morphometric shape data | Methods in Ecology and Evolution 4(4):393–399 | 10.1111/2041-210X.12035 | | https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.12035 | yes |
| baken2021geomorph | Baken, Erica K.; Collyer, Michael L.; Kaliontzopoulou, Antigoni; Adams, Dean C. | 2021 | geomorph v4.0 and gmShiny: Enhanced analytics and a new graphical interface for a comprehensive morphometric experience | Methods in Ecology and Evolution 12(12):2355–2363 | 10.1111/2041-210X.13723 | | https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.13723 | yes |
| collyer2018rrpp | Collyer, Michael L.; Adams, Dean C. | 2018 | RRPP: An R package for fitting linear models to high-dimensional data using residual randomization | Methods in Ecology and Evolution 9(7):1772–1779 | 10.1111/2041-210X.13029 | | | yes |
| collyer2015method | Collyer, M. L.; Sekora, D. J.; Adams, D. C. | 2015 | A method for analysis of phenotypic change for phenotypes described by high-dimensional data | Heredity 115(4):357–365 | 10.1038/hdy.2014.75 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4815463/ | yes |
| rohlf2015tps | Rohlf, F. James | 2015 | The tps series of software | Hystrix, the Italian Journal of Mammalogy 26(1):9–12 | 10.4404/hystrix-26.1-11264 | | | yes (doi.org metadata; pages from memory) |
| gunz2013semilandmarks | Gunz, Philipp; Mitteroecker, Philipp | 2013 | Semilandmarks: a method for quantifying curves and surfaces | Hystrix, the Italian Journal of Mammalogy 24(1):103–109 | 10.4404/hystrix-24.1-6292 | | | yes (doi.org metadata; pages from memory) |
| mitteroecker2009advances | Mitteroecker, Philipp; Gunz, Philipp | 2009 | Advances in geometric morphometrics | Evolutionary Biology 36(2):235–247 | 10.1007/s11692-009-9055-x | | https://link.springer.com/content/pdf/10.1007/s11692-009-9055-x.pdf | yes |
| klingenberg2016size | Klingenberg, Christian Peter | 2016 | Size, shape, and form: concepts of allometry in geometric morphometrics | Development Genes and Evolution 226(3):113–137 | 10.1007/s00427-016-0539-2 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4896994/ | yes |
| klingenberg2002shape | Klingenberg, Christian Peter; Barluenga, Marta; Meyer, Axel | 2002 | Shape analysis of symmetric structures: quantifying variation among individuals and asymmetry | Evolution 56(10):1909–1920 | 10.1111/j.0014-3820.2002.tb00117.x | | https://academic.oup.com/evolut/article-pdf/56/10/1909/47652105/evolut1909.pdf | yes |
| klingenberg1998geometric | Klingenberg, Christian Peter; McIntyre, Grant S. | 1998 | Geometric morphometrics of developmental instability: analyzing patterns of fluctuating asymmetry with Procrustes methods | Evolution 52(5):1363–1375 | 10.1111/j.1558-5646.1998.tb02018.x | | | yes |
| dryden2016statistical | Dryden, Ian L.; Mardia, Kanti V. | 2016 | Statistical Shape Analysis, with Applications in R (2nd ed.) | Wiley Series in Probability and Statistics, Wiley | 10.1002/9781119072492 | | | yes |
| walker2000ability | Walker, Jeffrey A. | 2000 | Ability of geometric morphometric methods to estimate a known covariance matrix | Systematic Biology 49(4):686–696 | 10.1080/106351500750049770 | | | yes |
| voncramontaubadel2007problem | von Cramon-Taubadel, Noreen; Frazier, Brenda C.; Lahr, Marta Mirazón | 2007 | The problem of assessing landmark error in geometric morphometrics: theory, methods, and modifications | American Journal of Physical Anthropology 134(1):24–35 | 10.1002/ajpa.20616 | | | yes |
| fruciano2016measurement | Fruciano, Carmelo | 2016 | Measurement error in geometric morphometrics | Development Genes and Evolution 226(3):139–158 | 10.1007/s00427-016-0537-4 | | | yes |
| cardini2007sample | Cardini, Andrea; Elton, Sarah | 2007 | Sample size and sampling error in geometric morphometric studies of size and shape | Zoomorphology 126(2):121–134 | 10.1007/s00435-007-0036-2 | | | yes |
| cardini2020less | Cardini, Andrea | 2020 | Less tautology, more biology? A comment on "high-density" morphometrics | Zoomorphology 139(4):513–529 | 10.1007/s00435-020-00499-w | | | yes |
| mitteroecker2022thirty | Mitteroecker, Philipp; Schaefer, Katrin | 2022 | Thirty years of geometric morphometrics: Achievements, challenges, and the ongoing quest for biological meaningfulness | American Journal of Biological Anthropology 178(S74):181–210 | 10.1002/ajpa.24531 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC9545184/ | yes |
| bardua2019practical | Bardua, C.; Felice, R. N.; Watanabe, A.; Fabre, A.-C.; Goswami, A. | 2019 | A practical guide to sliding and surface semilandmarks in morphometric analyses | Integrative Organismal Biology 1(1):obz016 | 10.1093/iob/obz016 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC7780474/ | yes |
| schlager2017morpho | Schlager, Stefan | 2017 | Morpho and Rvcg – Shape analysis in R | In: Statistical Shape and Deformation Analysis (Academic Press), pp. 217–256 | 10.1016/B978-0-12-810493-4.00011-0 | | | yes (book editors not retrieved) |
| porto2021alpaca | Porto, Arthur; Rolfe, Sara; Maga, A. Murat | 2021 | ALPACA: A fast and accurate computer vision approach for automated landmarking of three-dimensional biological structures | Methods in Ecology and Evolution 12(11):2129–2144 | 10.1111/2041-210X.13689 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC9291522/ | yes |
| rolfe2021slicermorph | Rolfe, Sara; Pieper, Steve; Porto, Arthur; Diamond, Kelly; Winchester, Julie; Shan, Shan; et al. (10 authors) | 2021 | SlicerMorph: An open and extensible platform to retrieve, visualize and analyse 3D morphology | Methods in Ecology and Evolution 12(10):1816–1825 | 10.1111/2041-210X.13669 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC12094517/ | yes |
| adams2026geomorphcran | Adams, Dean C.; et al. | 2026 | geomorph: Geometric Morphometric Analyses of 2D and 3D Landmark Data, R package version 4.1.1 | CRAN (published 2026-07-05) | | | https://cran.r-project.org/package=geomorph | yes (CRAN page; full author list and exact title not captured) |
| collyer2026rrppcran | Collyer, Michael L.; et al. | 2026 | RRPP: Linear Model Evaluation with Randomized Residuals in a Permutation Procedure, R package version 2.2.0 | CRAN (published 2026-07-03) | | | https://cran.r-project.org/package=RRPP | yes (version/date; title from memory) |
| schlager2025morphocran | Schlager, Stefan; et al. | 2025 | Morpho: Calculations and Visualisations Related to Geometric Morphometrics, R package version 2.13 | CRAN (published 2025-09-04) | | | https://cran.r-project.org/package=Morpho | yes (version/date; title from memory) |
| dryden2025shapescran | Dryden, Ian L. | 2025 | shapes: Statistical Shape Analysis, R package version 1.2.8 | CRAN (published 2025-10-24) | | | https://cran.r-project.org/package=shapes | yes (version/date) |
| klingenberg2026morphojsoftware | Klingenberg, Christian Peter | 2026 | MorphoJ, version 1.08.02 (software) | https://morphometrics.uk/MorphoJ_page.html | | | | yes (web page seen 2026-10-28) |
| kabsch1976solution | Kabsch, W. | 1976 | A solution for the best rotation to relate two sets of vectors | Acta Crystallographica Section A 32(5):922–923 | 10.1107/S0567739476001873 | | | yes |
| kabsch1978discussion | Kabsch, W. | 1978 | A discussion of the solution for the best rotation to relate two sets of vectors | Acta Crystallographica Section A 34(5):827–828 | 10.1107/S0567739478001680 | | | yes |
| mclachlan1982rapid | McLachlan, A. D. | 1982 | Rapid comparison of protein structures | Acta Crystallographica Section A 38(6):871–873 | 10.1107/S0567739482001806 | | | yes |
| kearsley1989orthogonal | Kearsley, S. K. | 1989 | On the orthogonal transformation used for structural comparisons | Acta Crystallographica Section A 45(2):208–210 | 10.1107/S0108767388010128 | | | yes |
| horn1987closed | Horn, Berthold K. P. | 1987 | Closed-form solution of absolute orientation using unit quaternions | Journal of the Optical Society of America A 4(4):629–642 | 10.1364/JOSAA.4.000629 | | | yes (end page from memory) |
| horn1988closed | Horn, Berthold K. P.; Hilden, Hugh M.; Negahdaripour, Shahriar | 1988 | Closed-form solution of absolute orientation using orthonormal matrices | Journal of the Optical Society of America A 5(7):1127–1135 | 10.1364/JOSAA.5.001127 | | | yes (end page from memory) |
| coutsias2004using | Coutsias, Evangelos A.; Seok, Chaok; Dill, Ken A. | 2004 | Using quaternions to calculate RMSD | Journal of Computational Chemistry 25(15):1849–1857 | 10.1002/jcc.20110 | | | yes |
| theobald2005rapid | Theobald, Douglas L. | 2005 | Rapid calculation of RMSDs using a quaternion-based characteristic polynomial | Acta Crystallographica Section A 61(4):478–480 | 10.1107/S0108767305015266 | | https://journals.iucr.org/a/issues/2005/04/00/sh5029/sh5029.pdf | yes |
| liu2010fast | Liu, Pu; Agrafiotis, Dimitris K.; Theobald, Douglas L. | 2010 | Fast determination of the optimal rotational matrix for macromolecular superpositions | Journal of Computational Chemistry 31(7):1561–1563 | 10.1002/jcc.21439 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/ | yes |
| theobald2006empirical | Theobald, Douglas L.; Wuttke, Deborah S. | 2006 | Empirical Bayes hierarchical models for regularizing maximum likelihood estimation in the matrix Gaussian Procrustes problem | Proceedings of the National Academy of Sciences 103(49):18521–18527 | 10.1073/pnas.0508445103 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC1664551/ | yes |
| theobald2006theseus | Theobald, Douglas L.; Wuttke, Deborah S. | 2006 | THESEUS: maximum likelihood superpositioning and analysis of macromolecular structures | Bioinformatics 22(17):2171–2172 | 10.1093/bioinformatics/btl332 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2584349/ | yes |
| coutsias2019rmsd | Coutsias, Evangelos A.; Wester, Michael J. | 2019 | RMSD and symmetry | Journal of Computational Chemistry 40(15):1496–1508 | 10.1002/jcc.25802 | | | yes |
| michaudagrawal2011mdanalysis | Michaud-Agrawal, Naveen; Denning, Elizabeth J.; Woolf, Thomas B.; Beckstein, Oliver | 2011 | MDAnalysis: A toolkit for the analysis of molecular dynamics simulations | Journal of Computational Chemistry 32(10):2319–2327 | 10.1002/jcc.21787 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC3144279/ | yes |
| amadei1993essential | Amadei, Andrea; Linssen, Antonius B. M.; Berendsen, Herman J. C. | 1993 | Essential dynamics of proteins | Proteins: Structure, Function, and Bioinformatics 17(4):412–425 | 10.1002/prot.340170408 | | | yes |
| maiorov1994significance | Maiorov, Vladimir N.; Crippen, Gordon M. | 1994 | Significance of root-mean-square deviation in comparing three-dimensional structures of globular proteins | Journal of Molecular Biology 235(2):625–634 | 10.1006/jmbi.1994.1017 | | | yes |
| carugo2001normalized | Carugo, Oliviero; Pongor, Sándor | 2001 | A normalized root-mean-square distance for comparing protein three-dimensional structures | Protein Science 10(7):1470–1473 | 10.1110/ps.690101 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2374114/ | yes (Crossref title has typo "spuare") |
| chothia1986relation | Chothia, C.; Lesk, A. M. | 1986 | The relation between the divergence of sequence and structure in proteins | The EMBO Journal 5(4):823–826 | 10.1002/j.1460-2075.1986.tb04288.x | | https://pmc.ncbi.nlm.nih.gov/articles/PMC1166865/ | yes |
| zemla2003lga | Zemla, A. | 2003 | LGA: a method for finding 3D similarities in protein structures | Nucleic Acids Research 31(13):3370–3374 | 10.1093/nar/gkg571 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC168977/ | yes |
| zhang2004scoring | Zhang, Yang; Skolnick, Jeffrey | 2004 | Scoring function for automated assessment of protein structure template quality | Proteins: Structure, Function, and Bioinformatics 57(4):702–710 | 10.1002/prot.20264 | | | yes |
| zhang2005tmalign | Zhang, Yang; Skolnick, Jeffrey | 2005 | TM-align: a protein structure alignment algorithm based on the TM-score | Nucleic Acids Research 33(7):2302–2309 | 10.1093/nar/gki524 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC1084323/ | yes (Crossref lists only "Zhang, Y."; Skolnick per RCSB 2024 reference list) |
| mariani2013lddt | Mariani, Valerio; Biasini, Marco; Barbato, Alessandro; Schwede, Torsten | 2013 | lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests | Bioinformatics 29(21):2722–2728 | 10.1093/bioinformatics/btt473 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC3799472/ | yes |
| shindyalov1998protein | Shindyalov, I. N.; Bourne, P. E. | 1998 | Protein structure alignment by incremental combinatorial extension (CE) of the optimal path | Protein Engineering 11(9):739–747 | 10.1093/protein/11.9.739 | | | yes |
| ye2003flexible | Ye, Yuzhen; Godzik, Adam | 2003 | Flexible structure alignment by chaining aligned fragment pairs allowing twists | Bioinformatics 19(suppl_2):ii246–ii255 | 10.1093/bioinformatics/btg1086 | | | yes |
| bittrich2024rcsb | Bittrich, Sebastian; Segura, Joan; Duarte, Jose M.; Burley, Stephen K.; Rose, Yana | 2024 | RCSB Protein Data Bank: exploring protein 3D similarities via comprehensive structural alignments | Bioinformatics 40(6):btae370 | 10.1093/bioinformatics/btae370 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC11212067/ | yes |
| kufareva2012methods | Kufareva, Irina; Abagyan, Ruben | 2012 | Methods of protein structure comparison | Methods in Molecular Biology 857:231–257 | 10.1007/978-1-61779-588-6_10 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4321859/ | yes (Crossref online year 2011; volume 857 from memory) |
| jumper2021highly | Jumper, John; Evans, Richard; Pritzel, Alexander; Green, Tim; Figurnov, Michael; Ronneberger, Olaf; et al. (34 authors) | 2021 | Highly accurate protein structure prediction with AlphaFold | Nature 596(7873):583–589 | 10.1038/s41586-021-03819-2 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC8371605/ | yes |
| abramson2024accurate | Abramson, Josh; Adler, Jonas; Dunger, Jack; Evans, Richard; Green, Tim; Pritzel, Alexander; et al. (48 authors) | 2024 | Accurate structure prediction of biomolecular interactions with AlphaFold 3 | Nature 630(8016):493–500 | 10.1038/s41586-024-07487-w | | https://pmc.ncbi.nlm.nih.gov/articles/PMC11168924/ | yes |
| buttenschoen2024posebusters | Buttenschoen, Martin; Morris, Garrett M.; Deane, Charlotte M. | 2024 | PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences | Chemical Science 15(9):3130–3139 | 10.1039/D3SC04185A | | https://pmc.ncbi.nlm.nih.gov/articles/PMC10901501/ | yes |
| bell2019dockrmsd | Bell, Eric W.; Zhang, Yang | 2019 | DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism | Journal of Cheminformatics 11(1):40 | 10.1186/s13321-019-0362-7 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC6556049/ | yes |
| chisholm2005compack | Chisholm, James A.; Motherwell, Sam | 2005 | COMPACK: a program for identifying crystal structure similarity using distances | Journal of Applied Crystallography 38(1):228–231 | 10.1107/S0021889804027074 | | | yes |
| wahba1965least | Wahba, Grace | 1965 | A least squares estimate of satellite attitude (Problem 65-1) | SIAM Review 7(3):409 | 10.1137/1007077 | | | yes |
| farrell1966least | Farrell, J. L.; Stuelpnagel, J. C.; Wessner, R. H.; Velman, J. R.; Brook, J. E. | 1966 | A least squares estimate of satellite attitude (Grace Wahba) [problem solutions] | SIAM Review 8(3):384–386 | 10.1137/1008080 | | | yes |
| black1964passive | Black, Harold D. | 1964 | A passive system for determining the attitude of a satellite | AIAA Journal 2(7):1350–1351 | 10.2514/3.2555 | | | yes |
| shuster1981three | Shuster, M. D.; Oh, S. D. | 1981 | Three-axis attitude determination from vector observations | Journal of Guidance and Control 4(1):70–77 | 10.2514/3.19717 | | | yes |
| markley2014fundamentals | Markley, F. Landis; Crassidis, John L. | 2014 | Fundamentals of Spacecraft Attitude Determination and Control | Springer (Space Technology Library) | 10.1007/978-1-4939-0802-8 | | | yes (series name from memory) |
| arun1987least | Arun, K. S.; Huang, T. S.; Blostein, S. D. | 1987 | Least-squares fitting of two 3-D point sets | IEEE Transactions on Pattern Analysis and Machine Intelligence PAMI-9(5):698–700 | 10.1109/TPAMI.1987.4767965 | | | yes |
| umeyama1991least | Umeyama, S. | 1991 | Least-squares estimation of transformation parameters between two point patterns | IEEE Transactions on Pattern Analysis and Machine Intelligence 13(4):376–380 | 10.1109/34.88573 | | | yes |
| eggert1997estimating | Eggert, D. W.; Lorusso, A.; Fisher, R. B. | 1997 | Estimating 3-D rigid body transformations: a comparison of four major algorithms | Machine Vision and Applications 9(5–6):272–290 | 10.1007/s001380050048 | | | yes |
| besl1992method | Besl, P. J.; McKay, Neil D. | 1992 | A method for registration of 3-D shapes | IEEE Transactions on Pattern Analysis and Machine Intelligence 14(2):239–256 | 10.1109/34.121791 | | | yes |
| chen1992object | Chen, Yang; Medioni, Gérard | 1992 | Object modelling by registration of multiple range images | Image and Vision Computing 10(3):145–155 | 10.1016/0262-8856(92)90066-C | | | yes |
| zhang1994iterative | Zhang, Zhengyou | 1994 | Iterative point matching for registration of free-form curves and surfaces | International Journal of Computer Vision 13(2):119–152 | 10.1007/BF01427149 | | | yes |
| rusinkiewicz2001efficient | Rusinkiewicz, S.; Levoy, M. | 2001 | Efficient variants of the ICP algorithm | Proceedings Third International Conference on 3-D Digital Imaging and Modeling (3DIM), pp. 145–152 | 10.1109/IM.2001.924423 | | | yes (year not in Crossref record; 2001 from conference) |
| segal2009generalized | Segal, A.; Haehnel, D.; Thrun, S. | 2009 | Generalized-ICP | Robotics: Science and Systems V | 10.15607/RSS.2009.V.021 | | https://doi.org/10.15607/rss.2009.v.021 | yes |
| pomerleau2015review | Pomerleau, François; Colas, Francis; Siegwart, Roland | 2015 | A review of point cloud registration algorithms for mobile robotics | Foundations and Trends in Robotics 4(1):1–104 | 10.1561/2300000035 | | https://hal.science/hal-01178661/document | yes |
| yang2021teaser | Yang, Heng; Shi, Jingnan; Carlone, Luca | 2021 | TEASER: Fast and certifiable point cloud registration | IEEE Transactions on Robotics 37(2):314–333 | 10.1109/TRO.2020.3033695 | 2001.07715 | https://arxiv.org/pdf/2001.07715 | yes |
| wang2019deep | Wang, Yue; Solomon, Justin M. | 2019 | Deep Closest Point: Learning representations for point cloud registration | 2019 IEEE/CVF International Conference on Computer Vision (ICCV), pp. 3522–3531 | 10.1109/ICCV.2019.00362 | 1905.03304 | https://arxiv.org/pdf/1905.03304 | yes |
| choy2020deep | Choy, Christopher; Dong, Wei; Koltun, Vladlen | 2020 | Deep Global Registration | 2020 IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), pp. 2511–2520 | 10.1109/CVPR42600.2020.00259 | 2004.11540 | https://arxiv.org/pdf/2004.11540 | yes |
| sturm2012benchmark | Sturm, Jürgen; Engelhard, Nikolas; Endres, Felix; Burgard, Wolfram; Cremers, Daniel | 2012 | A benchmark for the evaluation of RGB-D SLAM systems | 2012 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), pp. 573–580 | 10.1109/IROS.2012.6385773 | | | yes |
| zhang2018tutorial | Zhang, Zichao; Scaramuzza, Davide | 2018 | A tutorial on quantitative trajectory evaluation for visual(-inertial) odometry | 2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), pp. 7244–7251 | 10.1109/IROS.2018.8593941 | | https://www.zora.uzh.ch/id/eprint/175991/1/IROS18_Zhang.pdf | yes |
| murartal2015orb | Mur-Artal, Raúl; Montiel, J. M. M.; Tardós, Juan D. | 2015 | ORB-SLAM: A versatile and accurate monocular SLAM system | IEEE Transactions on Robotics 31(5):1147–1163 | 10.1109/TRO.2015.2463671 | 1502.00956 | https://arxiv.org/pdf/1502.00956 | yes |
| cootes1995active | Cootes, T. F.; Taylor, C. J.; Cooper, D. H.; Graham, J. | 1995 | Active shape models—their training and application | Computer Vision and Image Understanding 61(1):38–59 | 10.1006/cviu.1995.1004 | | https://pure.manchester.ac.uk/ws/files/30124448/POST-PEER-REVIEW-PUBLISHERS.PDF (listed by OpenAlex; not retrievable this session) | yes |
| cootes2001active | Cootes, T. F.; Edwards, G. J.; Taylor, C. J. | 2001 | Active appearance models | IEEE Transactions on Pattern Analysis and Machine Intelligence 23(6):681–685 | 10.1109/34.927467 | | | yes |
| heimann2009statistical | Heimann, Tobias; Meinzer, Hans-Peter | 2009 | Statistical shape models for 3D medical image segmentation: A review | Medical Image Analysis 13(4):543–563 | 10.1016/j.media.2009.05.004 | | | yes |
| davies2002minimum | Davies, R. H.; Twining, C. J.; Cootes, T. F.; Waterton, J. C.; Taylor, C. J. | 2002 | A minimum description length approach to statistical shape modeling | IEEE Transactions on Medical Imaging 21(5):525–537 | 10.1109/TMI.2002.1009388 | | | yes |
| luthi2018gaussian | Lüthi, Marcel; Gerig, Thomas; Jud, Christoph; Vetter, Thomas | 2018 | Gaussian process morphable models | IEEE Transactions on Pattern Analysis and Machine Intelligence 40(8):1860–1873 | 10.1109/TPAMI.2017.2739743 | 1603.07254 | https://arxiv.org/pdf/1603.07254 | yes |
| ambellan2019statistical | Ambellan, Felix; Lamecker, Hans; von Tycowicz, Christoph; Zachow, Stefan | 2019 | Statistical shape models: understanding and mastering variation in anatomy | Advances in Experimental Medicine and Biology 1156:67–84 | 10.1007/978-3-030-19385-0_5 | | | yes (volume 1156 from memory) |
| cates2017shapeworks | Cates, Joshua; Elhabian, Shireen; Whitaker, Ross | 2017 | ShapeWorks: particle-based shape correspondence and visualization software | In: Statistical Shape and Deformation Analysis (Academic Press), pp. 257–298 | 10.1016/B978-0-12-810493-4.00012-2 | | | yes (Crossref title is "ShapeWorks"; subtitle from memory) |
| ionescu2014human36m | Ionescu, Catalin; Papava, Dragos; Olaru, Vlad; Sminchisescu, Cristian | 2014 | Human3.6M: Large scale datasets and predictive methods for 3D human sensing in natural environments | IEEE Transactions on Pattern Analysis and Machine Intelligence 36(7):1325–1339 | 10.1109/TPAMI.2013.248 | | | yes |
| martinez2017simple | Martinez, Julieta; Hossain, Rayat; Romero, Javier; Little, James J. | 2017 | A simple yet effective baseline for 3d human pose estimation | 2017 IEEE International Conference on Computer Vision (ICCV), pp. 2659–2668 | 10.1109/ICCV.2017.288 | 1705.03098 | https://arxiv.org/pdf/1705.03098 | yes |
| kanazawa2018end | Kanazawa, Angjoo; Black, Michael J.; Jacobs, David W.; Malik, Jitendra | 2018 | End-to-end recovery of human shape and pose | 2018 IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), pp. 7122–7131 | 10.1109/CVPR.2018.00744 | 1712.06584 | https://arxiv.org/pdf/1712.06584 | yes |
| blanz1999morphable | Blanz, Volker; Vetter, Thomas | 1999 | A morphable model for the synthesis of 3D faces | Proceedings of SIGGRAPH '99, pp. 187–194 | 10.1145/311535.311556 | | https://pure.mpg.de/pubman/item/item_1793809_3/component/file_3271006/SIGGRAPH-1999-Blanz.pdf | yes |
| egger20203d | Egger, Bernhard; Smith, William A. P.; Tewari, Ayush; Wuhrer, Stefanie; Zollhoefer, Michael; Beeler, Thabo; et al. (13 authors) | 2020 | 3D morphable face models—past, present, and future | ACM Transactions on Graphics 39(5):157 (pp. 1–38) | 10.1145/3395208 | | | yes (article number 157 from memory) |
| wang2003automatic | Wang, Liang; Tan, Tieniu; Hu, Weiming; Ning, Huazhong | 2003 | Automatic gait recognition based on statistical shape analysis | IEEE Transactions on Image Processing 12(9):1120–1131 | 10.1109/TIP.2003.815251 | | | yes |
| soderkvist1993determining | Söderkvist, Inge; Wedin, Per-Åke | 1993 | Determining the movements of the skeleton using well-configured markers | Journal of Biomechanics 26(12):1473–1477 | 10.1016/0021-9290(93)90098-Y | | | yes |
| challis1995procedure | Challis, John H. | 1995 | A procedure for determining rigid body transformation parameters | Journal of Biomechanics 28(6):733–737 | 10.1016/0021-9290(94)00116-L | | | yes |
| williams1984use | Williams, Anthony A.; Langron, Steven P. | 1984 | The use of free-choice profiling for the evaluation of commercial ports | Journal of the Science of Food and Agriculture 35(5):558–568 | 10.1002/jsfa.2740350513 | | | yes |
| dijksterhuis1991interpretation | Dijksterhuis, G. B.; Gower, J. C. | 1991 | The interpretation of generalized Procrustes analysis and allied methods | Food Quality and Preference 3(2):67–87 | 10.1016/0950-3293(91)90027-C | | | yes |
| gower2004procrustes | Gower, John C.; Dijksterhuis, Garmt B. | 2004 | Procrustes Problems | Oxford University Press (Oxford Statistical Science Series) | 10.1093/acprof:oso/9780198510581.001.0001 | | | yes (series name from memory) |
| king1991statistical | King, Bonnie M.; Arents, Paul | 1991 | A statistical test of consensus obtained from generalized Procrustes analysis of sensory data | Journal of Sensory Studies 6(1):37–48 | 10.1111/j.1745-459X.1991.tb00500.x | | | yes |
| wakeling1992new | Wakeling, Ian N.; Raats, Monique M.; MacFie, Halliday J. H. | 1992 | A new significance test for consensus in generalized Procrustes analysis | Journal of Sensory Studies 7(2):91–96 | 10.1111/j.1745-459X.1992.tb00526.x | | | yes |
| xiong2008permutation | Xiong, R.; Blot, K.; Meullenet, J. F.; Dessirier, J. M. | 2008 | Permutation tests for generalized Procrustes analysis | Food Quality and Preference 19(2):146–155 | 10.1016/j.foodqual.2007.03.003 | | | yes |
| yu2025understanding | Yu, [given name not retrieved]; Logan, [—]; Xie, [—]; King, [—] | 2025 | Understanding cultural differences across global descriptive sensory panels using Generalized Procrustes Analysis | Food Quality and Preference 133:105633 | 10.1016/j.foodqual.2025.105633 | | | yes (family names only) |
| varela2012sensory | Varela, Paula; Ares, Gastón | 2012 | Sensory profiling, the blurred line between sensory and consumer science. A review of novel methods for product characterization | Food Research International 48(2):893–908 | 10.1016/j.foodres.2012.06.037 | | | yes |
| pages2005collection | Pagès, Jérôme | 2005 | Collection and analysis of perceived product inter-distances using multiple factor analysis: Application to the study of 10 white wines from the Loire Valley | Food Quality and Preference 16(7):642–649 | 10.1016/j.foodqual.2005.01.006 | | | yes |
| le2008factominer | Lê, Sébastien; Josse, Julie; Husson, François | 2008 | FactoMineR: An R package for multivariate analysis | Journal of Statistical Software 25(1):1–18 | 10.18637/jss.v025.i01 | | | yes (pages from memory) |
| mosier1939determining | Mosier, Charles I. | 1939 | Determining a simple structure when loadings for certain tests are known | Psychometrika 4(2):149–162 | 10.1007/BF02288493 | | | yes |
| green1952orthogonal | Green, Bert F. | 1952 | The orthogonal approximation of an oblique structure in factor analysis | Psychometrika 17(4):429–440 | 10.1007/BF02288918 | | | yes |
| hurley1962procrustes | Hurley, John R.; Cattell, Raymond B. | 1962 | The Procrustes program: producing direct rotation to test a hypothesized factor structure | Behavioral Science 7(2):258–262 | 10.1002/bs.3830070216 | | | yes (Crossref shows 2007 online date; print 1962) |
| schonemann1966generalized | Schönemann, Peter H. | 1966 | A generalized solution of the orthogonal Procrustes problem | Psychometrika 31(1):1–10 | 10.1007/BF02289451 | | | yes |
| browne1967oblique | Browne, Michael W. | 1967 | On oblique Procrustes rotation | Psychometrika 32(2):125–132 | 10.1007/BF02289420 | | | yes |
| schonemann1970fitting | Schönemann, Peter H.; Carroll, Robert M. | 1970 | Fitting one matrix to another under choice of a central dilation and a rigid motion | Psychometrika 35(2):245–255 | 10.1007/BF02291266 | | | yes |
| tenberge1977orthogonal | ten Berge, Jos M. F. | 1977 | Orthogonal Procrustes rotation for two or more matrices | Psychometrika 42(2):267–276 | 10.1007/BF02294053 | | | yes |
| mccrae1996evaluating | McCrae, Robert R.; Zonderman, Alan B.; Costa, Paul T.; Bond, Michael H.; Paunonen, Sampo V. | 1996 | Evaluating replicability of factors in the Revised NEO Personality Inventory: Confirmatory factor analysis versus Procrustes rotation | Journal of Personality and Social Psychology 70(3):552–566 | 10.1037/0022-3514.70.3.552 | | | yes |
| mccrae1997personality | McCrae, Robert R.; Costa, Paul T. | 1997 | Personality trait structure as a human universal | American Psychologist 52(5):509–516 | 10.1037/0003-066X.52.5.509 | | | yes |
| bernaards2005gradient | Bernaards, Coen A.; Jennrich, Robert I. | 2005 | Gradient projection algorithms and software for arbitrary rotation criteria in factor analysis | Educational and Psychological Measurement 65(5):676–696 | 10.1177/0013164404272507 | | | yes |
| fischer2019primer | Fischer, Ronald; Karl, Johannes A. | 2019 | A primer to (cross-cultural) multi-group invariance testing possibilities in R | Frontiers in Psychology 10:1507 | 10.3389/fpsyg.2019.01507 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC6657455/ | yes |
| jackson1995protest | Jackson, Donald A. | 1995 | PROTEST: A PROcrustean Randomization TEST of community environment concordance | Écoscience 2(3):297–303 | 10.1080/11956860.1995.11682297 | | | yes |
| peresneto2001how | Peres-Neto, Pedro R.; Jackson, Donald A. | 2001 | How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test | Oecologia 129(2):169–178 | 10.1007/s004420100720 | | | yes |
| lisboa2014much | Lisboa, Francy Junio Gonçalves; Peres-Neto, Pedro R.; Chaer, Guilherme Montandon; Jesus, Ederson da Conceição; Mitchell, Ruth Joy; Chapman, Stephen James; Berbara, Ricardo Luis Louro | 2014 | Much beyond Mantel: Bringing Procrustes association metric to the plant and soil ecologist's toolbox | PLoS ONE 9(6):e101238 | 10.1371/journal.pone.0101238 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4074130/ | yes (Berbara's given names from memory) |
| oksanen2026vegan | Oksanen, Jari; et al. | 2026 | vegan: Community Ecology Package, R package version 2.7-6 | CRAN (published 2026-08-31) | | | https://cran.r-project.org/package=vegan | yes (CRAN page; full author list not captured) |
| muegge2011diet | Muegge, Brian D.; Kuczynski, Justin; Knights, Dan; Clemente, Jose C.; González, Antonio; Fontana, Luigi; et al. (9 authors) | 2011 | Diet drives convergence in gut microbiome functions across mammalian phylogeny and within humans | Science 332(6032):970–974 | 10.1126/science.1198719 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC3303602/ | yes |
| vazquezbaeza2013emperor | Vázquez-Baeza, Yoshiki; Pirrung, Meg; Gonzalez, Antonio; Knight, Rob | 2013 | EMPeror: a tool for visualizing high-throughput microbial community data | GigaScience 2(1):16 | 10.1186/2047-217X-2-16 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4076506/ | yes |
| bolyen2019qiime2 | Bolyen, Evan; Rideout, Jai Ram; Dillon, Matthew R.; Bokulich, Nicholas A.; Abnet, Christian C.; Al-Ghalith, Gabriel A.; et al. (112 authors) | 2019 | Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2 | Nature Biotechnology 37(8):852–857 | 10.1038/s41587-019-0209-9 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC7015180/ | yes |
| wang2010comparing | Wang, Chaolong; Szpiech, Zachary A.; Degnan, James H.; Jakobsson, Mattias; Pemberton, Trevor J.; Hardy, John A.; Singleton, Andrew B.; Rosenberg, Noah A. | 2010 | Comparing spatial maps of human population-genetic variation using Procrustes analysis | Statistical Applications in Genetics and Molecular Biology 9(1): Article 13 | 10.2202/1544-6115.1493 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2861313/ | yes (article number 13 from memory) |
| wang2012quantitative | Wang, Chaolong; Zöllner, Sebastian; Rosenberg, Noah A. | 2012 | A quantitative comparison of the similarity between genes and geography in worldwide human populations | PLoS Genetics 8(8):e1002886 | 10.1371/journal.pgen.1002886 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC3426559/ | yes |
| wang2014ancestry | Wang, Chaolong; Zhan, Xiaowei; Bragg-Gresham, Jennifer; Kang, Hyun Min; Stambolian, Dwight; et al. | 2014 | Ancestry estimation and control of population stratification for sequence-based association studies | Nature Genetics 46(4):409–415 | 10.1038/ng.2924 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4084909/ | yes (Crossref lists consortium "The FUSION Study" as an author) |
| novembre2008genes | Novembre, John; Johnson, Toby; Bryc, Katarzyna; Kutalik, Zoltán; Boyko, Adam R.; Auton, Adam; et al. (12 authors) | 2008 | Genes mirror geography within Europe | Nature 456(7218):98–101 | 10.1038/nature07331 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2735096/ | yes |
| zeira2022alignment | Zeira, Ron; Land, Max; Strzalkowski, Alexander; Raphael, Benjamin J. | 2022 | Alignment and integration of spatial transcriptomics data | Nature Methods 19(5):567–575 | 10.1038/s41592-022-01459-6 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC9334025/ | yes |
| mikolov2013exploiting | Mikolov, Tomas; Le, Quoc V.; Sutskever, Ilya | 2013 | Exploiting similarities among languages for machine translation | arXiv preprint | | 1309.4168 | https://arxiv.org/pdf/1309.4168 | yes |
| xing2015normalized | Xing, Chao; Wang, Dong; Liu, Chao; Lin, Yiye | 2015 | Normalized word embedding and orthogonal transform for bilingual word translation | Proceedings of NAACL-HLT 2015, pp. 1006–1011 | 10.3115/v1/N15-1104 | | https://aclanthology.org/N15-1104.pdf | yes |
| artetxe2016learning | Artetxe, Mikel; Labaka, Gorka; Agirre, Eneko | 2016 | Learning principled bilingual mappings of word embeddings while preserving monolingual invariance | Proceedings of EMNLP 2016, pp. 2289–2294 | 10.18653/v1/D16-1250 | | https://aclanthology.org/D16-1250.pdf | yes |
| artetxe2017learning | Artetxe, Mikel; Labaka, Gorka; Agirre, Eneko | 2017 | Learning bilingual word embeddings with (almost) no bilingual data | Proceedings of ACL 2017 (Vol. 1: Long Papers), pp. 451–462 | 10.18653/v1/P17-1042 | | https://aclanthology.org/P17-1042.pdf | yes |
| artetxe2018robust | Artetxe, Mikel; Labaka, Gorka; Agirre, Eneko | 2018 | A robust self-learning method for fully unsupervised cross-lingual mappings of word embeddings | Proceedings of ACL 2018 (Vol. 1: Long Papers), pp. 789–798 | 10.18653/v1/P18-1073 | 1805.06297 | https://aclanthology.org/P18-1073.pdf | yes |
| smith2017offline | Smith, Samuel L.; Turban, David H. P.; Hamblin, Steven; Hammerla, Nils Y. | 2017 | Offline bilingual word vectors, orthogonal transformations and the inverted softmax | International Conference on Learning Representations (ICLR 2017) | | 1702.03859 | https://arxiv.org/pdf/1702.03859 | yes |
| conneau2018word | Conneau, Alexis; Lample, Guillaume; Ranzato, Marc'Aurelio; Denoyer, Ludovic; Jégou, Hervé | 2018 | Word translation without parallel data | International Conference on Learning Representations (ICLR 2018) | | 1710.04087 | https://arxiv.org/pdf/1710.04087 | yes |
| dinu2015improving | Dinu, Georgiana; Lazaridou, Angeliki; Baroni, Marco | 2015 | Improving zero-shot learning by mitigating the hubness problem | ICLR 2015 Workshop track | | 1412.6568 | https://arxiv.org/pdf/1412.6568 | yes (venue from memory) |
| joulin2018loss | Joulin, Armand; Bojanowski, Piotr; Mikolov, Tomas; Jégou, Hervé; Grave, Edouard | 2018 | Loss in translation: Learning bilingual word mapping with a retrieval criterion | Proceedings of EMNLP 2018, pp. 2979–2984 | 10.18653/v1/D18-1330 | 1804.07745 | https://aclanthology.org/D18-1330.pdf | yes |
| grave2019unsupervised | Grave, Edouard; Joulin, Armand; Berthet, Quentin | 2019 | Unsupervised alignment of embeddings with Wasserstein Procrustes | Proceedings of AISTATS 2019 (PMLR 89) | | 1805.11222 | https://arxiv.org/pdf/1805.11222 | yes (arXiv record; venue from memory) |
| alvarezmelis2018gromov | Alvarez-Melis, David; Jaakkola, Tommi | 2018 | Gromov-Wasserstein alignment of word embedding spaces | Proceedings of EMNLP 2018, pp. 1881–1890 | 10.18653/v1/D18-1214 | | https://aclanthology.org/D18-1214.pdf | yes |
| ruder2019survey | Ruder, Sebastian; Vulić, Ivan; Søgaard, Anders | 2019 | A survey of cross-lingual word embedding models | Journal of Artificial Intelligence Research 65:569–631 | 10.1613/jair.1.11640 | 1706.04902 | https://arxiv.org/pdf/1706.04902 | yes |
| sogaard2018limitations | Søgaard, Anders; Ruder, Sebastian; Vulić, Ivan | 2018 | On the limitations of unsupervised bilingual dictionary induction | Proceedings of ACL 2018 (Vol. 1: Long Papers), pp. 778–788 | 10.18653/v1/P18-1072 | | https://aclanthology.org/P18-1072.pdf | yes |
| glavas2019properly | Glavaš, Goran; Litschko, Robert; Ruder, Sebastian; Vulić, Ivan | 2019 | How to (properly) evaluate cross-lingual word embeddings: On strong baselines, comparative analyses, and some misconceptions | Proceedings of ACL 2019, pp. 710–721 | 10.18653/v1/P19-1070 | | | yes |
| hamilton2016diachronic | Hamilton, William L.; Leskovec, Jure; Jurafsky, Dan | 2016 | Diachronic word embeddings reveal statistical laws of semantic change | Proceedings of ACL 2016 (Vol. 1: Long Papers), pp. 1489–1501 | 10.18653/v1/P16-1141 | | https://aclanthology.org/P16-1141.pdf | yes |
| jha2025harnessing | Jha, Rishi; Zhang, Collin; Shmatikov, Vitaly; Morris, John X. | 2025 | Harnessing the universal geometry of embeddings | arXiv preprint | | 2505.12540 | https://arxiv.org/pdf/2505.12540 | yes |
| williams2021generalized | Williams, Alex H.; Kunz, Erin; Kornblith, Simon; Linderman, Scott W. | 2021 | Generalized shape metrics on neural representations | Advances in Neural Information Processing Systems 34 (NeurIPS 2021) | | 2110.14739 | https://arxiv.org/pdf/2110.14739 | yes |
| ding2021grounding | Ding, Frances; Denain, Jean-Stanislas; Steinhardt, Jacob | 2021 | Grounding representation similarity with statistical testing | Advances in Neural Information Processing Systems 34 (NeurIPS 2021) | | 2108.01661 | https://arxiv.org/pdf/2108.01661 | yes |
| kornblith2019similarity | Kornblith, Simon; Norouzi, Mohammad; Lee, Honglak; Hinton, Geoffrey | 2019 | Similarity of neural network representations revisited | Proceedings of ICML 2019 (PMLR 97) | | 1905.00414 | https://arxiv.org/pdf/1905.00414 | yes (PMLR volume from memory) |
| klabunde2025similarity | Klabunde, Max; Schumacher, Tobias; Strohmaier, Markus; Lemmerich, Florian | 2025 | Similarity of neural network models: A survey of functional and representational measures | ACM Computing Surveys 57(9):242 (52 pp.) | 10.1145/3728458 | 2305.06329 | https://dl.acm.org/doi/pdf/10.1145/3728458 | yes |
| duong2023representational | Duong, Lyndon R.; Zhou, Jingyang; Nassar, Josue; Berman, Jules; Olieslagers, Jeroen; Williams, Alex H. | 2023 | Representational dissimilarity metric spaces for stochastic neural networks | International Conference on Learning Representations (ICLR 2023) | | 2211.11665 | https://arxiv.org/pdf/2211.11665 | yes |
| harvey2024what | Harvey, Sarah E.; Lipshutz, David; Williams, Alex H. | 2024 | What representational similarity measures imply about decodable information | arXiv preprint | | 2411.08197 | https://arxiv.org/pdf/2411.08197 | yes |
| ostrow2023beyond | Ostrow, Mitchell; Eisen, Adam; Kozachkov, Leo; Fiete, Ila | 2023 | Beyond geometry: Comparing the temporal structure of computation in neural circuits with dynamical similarity analysis | Advances in Neural Information Processing Systems 36 (NeurIPS 2023) | | 2306.10168 | https://arxiv.org/pdf/2306.10168 | yes (arXiv record; venue from memory) |
| haxby2011common | Haxby, James V.; Guntupalli, J. Swaroop; Connolly, Andrew C.; Halchenko, Yaroslav O.; Conroy, Bryan R.; Gobbini, M. Ida; Hanke, Michael; Ramadge, Peter J. | 2011 | A common, high-dimensional model of the representational space in human ventral temporal cortex | Neuron 72(2):404–416 | 10.1016/j.neuron.2011.08.026 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC3201764/ | yes |
| haxby2020hyperalignment | Haxby, James V.; Guntupalli, J. Swaroop; Nastase, Samuel A.; Feilong, Ma | 2020 | Hyperalignment: Modeling shared information encoded in idiosyncratic cortical topographies | eLife 9:e56601 | 10.7554/eLife.56601 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC7266639/ | yes |
| guntupalli2016model | Guntupalli, J. Swaroop; Hanke, Michael; Halchenko, Yaroslav O.; Connolly, Andrew C.; Ramadge, Peter J.; Haxby, James V. | 2016 | A model of representational spaces in human cortex | Cerebral Cortex 26(6):2919–2934 | 10.1093/cercor/bhw068 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC4869822/ | yes |
| hanke2009pymvpa | Hanke, Michael; Halchenko, Yaroslav O.; Sederberg, Per B.; Hanson, Stephen José; Haxby, James V.; Pollmann, Stefan | 2009 | PyMVPA: A Python toolbox for multivariate pattern analysis of fMRI data | Neuroinformatics 7(1):37–53 | 10.1007/s12021-008-9041-y | | https://pmc.ncbi.nlm.nih.gov/articles/PMC2664559/ | yes |
| virtanen2020scipy | Virtanen, Pauli; Gommers, Ralf; Oliphant, Travis E.; Haberland, Matt; Reddy, Tyler; Cournapeau, David; et al. (112 authors) | 2020 | SciPy 1.0: fundamental algorithms for scientific computing in Python | Nature Methods 17(3):261–272 | 10.1038/s41592-019-0686-2 | | https://pmc.ncbi.nlm.nih.gov/articles/PMC7056644/ | yes |
| tobler1965computation | Tobler, Waldo R. | 1965 | Computation of the correspondence of geographical patterns | Papers of the Regional Science Association 15(1):131–139 | 10.1007/BF01947869 | | | yes |
| tobler1994bidimensional | Tobler, W. R. | 1994 | Bidimensional regression | Geographical Analysis 26(3):187–212 | 10.1111/j.1538-4632.1994.tb00320.x | | | yes |
| nakaya1997statistical | Nakaya, Tomoki | 1997 | Statistical inferences in bidimensional regression models | Geographical Analysis 29(2):169–186 | 10.1111/j.1538-4632.1997.tb00954.x | | | yes |
| friedman2003bidimensional | Friedman, Alinda; Kohler, Bernd | 2003 | Bidimensional regression: Assessing the configural similarity and accuracy of cognitive maps and other two-dimensional data sets | Psychological Methods 8(4):468–491 | 10.1037/1082-989X.8.4.468 | | | yes |
| jenny2007visualizing | Jenny, Bernhard; Weber, Adrian; Hurni, Lorenz | 2007 | Visualizing the planimetric accuracy of historical maps with MapAnalyst | Cartographica 42(1):89–94 | 10.3138/carto-v42-1-089 | | | yes |
| jenny2011studying | Jenny, Bernhard; Hurni, Lorenz | 2011 | Studying cartographic heritage: Analysis and visualization of geometric distortions | Computers & Graphics 35(2):402–411 | 10.1016/j.cag.2011.01.005 | | | yes |
| grafarend2003nonlinear | Grafarend, E. W.; Awange, J. L. | 2003 | Nonlinear analysis of the three-dimensional datum transformation [conformal group C7(3)] | Journal of Geodesy 77(1–2):66–76 | 10.1007/s00190-002-0299-9 | | | yes (group symbol garbled in Crossref title) |
| awange2008procrustean | Awange, J. L.; Bae, K.-H.; Claessens, S. J. | 2008 | Procrustean solution of the 9-parameter transformation problem | Earth, Planets and Space 60(6):529–537 | 10.1186/BF03353115 | | | yes |
| palancz2010extension | Paláncz, B.; Zaletnyik, P.; Awange, J. L.; Heck, B. | 2010 | Extension of the ABC-Procrustes algorithm for 3D affine coordinate transformation | Earth, Planets and Space 62(11):857–862 | 10.5047/eps.2010.10.004 | | | yes |
| ligas2020point | Ligas, Marcin | 2020 | Point-wise weighted solution for the similarity transformation parameters | Applied Geomatics 12(3):371–377 | 10.1007/s12518-020-00306-7 | | | yes |
| yetkin2017determination | Yetkin, Mevlut; Ansari, Kutubuddin | 2017 | On the determination of transformation parameters between different ITRS realizations using Procrustes approach in Turkey | Journal of Applied Geodesy 11(3):123–129 | 10.1515/jag-2016-0048 | | | yes |
| altamimi2023itrf2020 | Altamimi, Zuheir; Rebischung, Paul; Collilieux, Xavier; Métivier, Laurent; Chanard, Kristel | 2023 | ITRF2020: an augmented reference frame refining the modeling of nonlinear station motions | Journal of Geodesy 97(5):47 | 10.1007/s00190-023-01738-w | | https://link.springer.com/content/pdf/10.1007/s00190-023-01738-w.pdf | yes |
| anderson2013kendall | Anderson, Edward | 2013 | Kendall's shape statistics as a classical realization of Barbour-type timeless records theory approach to quantum gravity | arXiv preprint | | 1307.1923 | https://arxiv.org/pdf/1307.1923 | no (title and ID seen in search results; author from memory) |
| anderson2015six | Anderson, Edward | 2015 | Six new mechanics corresponding to further shape theories | arXiv preprint | | 1505.00488 | https://arxiv.org/pdf/1505.00488 | no (title and ID seen in search results; author from memory) |
