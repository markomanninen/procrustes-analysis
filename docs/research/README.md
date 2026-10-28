# Research archive: Procrustes analysis and Barbour's best matching

October 2026

This folder holds the research behind the triangle page: what Procrustes fitting is, where it
comes from, and how it relates to Julian Barbour's best matching. It contains eight sets of notes,
one report that draws them together, a merged bibliography, a list of open-access papers with a
script to fetch them, and the scripts that checked the numerical claims.

The main report is
[reports/Procrustes and Barbour best matching.md](reports/Procrustes%20and%20Barbour%20best%20matching.md).

## What it covers

- **Procrustes analysis and its uses.** Fitting one configuration onto another by translation,
  rotation and possibly scale, from the factor-analysis rotations of 1939–1966 through Kabsch,
  Wahba, Horn and Umeyama to Kendall's shape spaces. In the plane the whole fit is one complex
  number, w\*z/(w\*w). It is used in morphometrics, protein superposition, spacecraft attitude,
  point-cloud registration, sensory science, ecology, population genetics, word embeddings and
  geodesy.
- **Barbour's best matching and shape dynamics.** Best matching compares two configurations
  after shifting and turning one of them (and from 2003 also scaling it) to minimise a
  mass-weighted distance. That is a Procrustes fit, applied between neighbouring instants of a
  motion. In the twelve Barbour texts searched, Kendall is cited once and the word Procrustes never
  appears. The notes follow
  the programme to shape dynamics for gravity, the Janus point and work up to October 2026.
- **Kendall's shape sphere, Anderson and holonomy.** Triangles up to similarity form a sphere.
  Edward Anderson identified Kendall's shape spaces with Barbour's reduced configuration spaces.
  Chained best matching round a closed loop of shapes turns the triangle by half the solid angle
  the loop encloses, Ω/2. A single Procrustes fit of the last shape onto the first does not turn
  it at all. This is the "falling cat" holonomy of the gauge theory of deformable bodies
  (Shapere and Wilczek, Montgomery, Littlejohn and Reinsch).
- **Labelling and covariance in general relativity.** Coordinates are labels. Comparing fields at
  neighbouring points needs the Christoffel connection, and the hole argument makes the labels
  gauge. The ADM shift is the best-matching variable for spatial relabellings, and superspace
  (metrics up to diffeomorphism) plays the part of shape space.
- **Fitting with unknown correspondence.** When the vertex pairing is unknown, the best of
  several fits is chosen, and that choice jumps. ICP, softassign, Coherent Point Drift and optimal
  transport replace the hard choice by weights on point pairs. The page's soft fit instead
  weights each of the 12 whole candidate fits by e^(−Δrms²/σ²), so it never jumps. The notes
  compare the two.
- **Novelty of the page's jump locus and side-length rule.** With mirrors allowed, the best of
  the 12 fits pairs vertices by the rank of the opposite side length, and mirrors when the two
  triangles have opposite orientation. The fit jumps only when the target triangle passes through
  an isosceles or flat shape. The geometry is classical: Kendall's "spherical blackboard" and the
  chambers of a reflection group. The labelling is Groth's (1986) star-matching rule. The notes
  found no earlier statement that this labelling is the least-squares optimum, or of where the fit
  jumps. Kendall (1985), Kendall et al. (1999) and Small (1996) could not be read, so the result
  is best stated as an observation, not a new theorem.

## The notes

The notes are in
[research_notes/Procrustes and Barbour best matching/](research_notes/Procrustes%20and%20Barbour%20best%20matching/).
Each ends with a Source catalogue table.

| file | contents |
|---|---|
| [procrustes_history_mathematics.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/procrustes_history_mathematics.md) | History of Procrustes problems from Mosier (1939) on; the planar solution; full, partial and Riemannian distances; reflections; Kendall's shape spaces; how `js/procrustes.js` maps onto all of it. |
| [procrustes_applications.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/procrustes_applications.md) | Where Procrustes analysis is used, field by field, with the software that implements it. |
| [barbour_best_matching.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/barbour_best_matching.md) | Barbour's best matching from 1974 to the 2003 scale-invariant theory, its critics, and how it differs from statistical Procrustes. |
| [shape_dynamics_janus_point.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/shape_dynamics_janus_point.md) | Shape dynamics for fields and gravity, shape complexity, the Janus point and arrows of time, and the debate up to October 2026. |
| [shape_sphere_bridge.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/shape_sphere_bridge.md) | The triangle shape sphere, its use in the three-body problem, Anderson's link from Kendall to Barbour, and the holonomy Ω/2. |
| [correspondence_soft_alignment.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/correspondence_soft_alignment.md) | Alignment with unknown pairing: unlabelled shape analysis, ICP, softassign, CPD, optimal transport, and why the hard fit jumps. |
| [labelling_covariance_gr.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/labelling_covariance_gr.md) | Labels, the Christoffel connection, the hole argument, the ADM shift as best matching, and superspace as the field version of shape space. |
| [novelty_jump_locus.md](research_notes/Procrustes%20and%20Barbour%20best%20matching/novelty_jump_locus.md) | Prior-art search for the page's two claims: where the best fit jumps, and the side-length rule. |

## Other files

- [bibliography.bib](bibliography.bib): all 481 sources as BibTeX, merged from the eight
  catalogues (604 rows). Duplicates were merged by DOI, then arXiv ID, then title. Each entry uses
  a catalogue key. The other keys of the same work are listed in `ids`, which biblatex accepts as
  aliases. BibTeX and biber process it without errors.
- [sources.md](sources.md): the same sources as tables, grouped by topic and sorted by year,
  with links, verification status and the notes that cite each one.
- [papers-manifest.tsv](papers-manifest.tsv): 240 open-access PDFs for 238 sources, from arXiv,
  PMC, author and institution pages, and publishers' open-access copies. For each it gives the
  target file name, URL, host, licence where one is stated, size from an HTTP HEAD request, and
  status code.
- [papers-core.tsv](papers-core.tsv): the core set, 33 files for 31 key sources (39 MB), in the
  same format.
- [fetch-papers.sh](fetch-papers.sh): downloads the files in a manifest.
- [scripts/](scripts/): the scripts that checked the numerical claims, with their output.

## Fetching the PDFs

The PDFs are not in the repository. To download the core set, or every listed paper, run

    ./fetch-papers.sh papers-core.tsv
    ./fetch-papers.sh

from this folder, or by its path from anywhere. The core set covers the Procrustes solutions
(Horn, Dryden's notes, Green and Mardia, Coherent Point Drift, Williams et al., Smith et al.),
Barbour's best matching and shape dynamics (Barbour–Foster–Ó Murchadha, Barbour 2003, 2009,
2010, 2012 and 2024, Mercati's tutorial, Gomes–Gryb–Koslowski), the Janus-point debate
(Barbour–Koslowski–Mercati 2014 and 2016, Zeh, Lazarovici–Reichert), the shape sphere and
holonomy (Montgomery 1996 and 2015, Anderson 2008 and 2015), labelling in GR (ADM, Gryb,
Giulini, Carroll), and the prior art for the side-length rule (Groth, Pál–Bakos,
Kuznetsov–Polishchuk, Bowden et al., Rubix). It was downloaded on 2026-10-28; all 33 files
arrived. The files go to `papers/`, which is git-ignored
(`docs/research/papers/` in the repository's `.gitignore`). They are kept out of git because not
every paper's licence allows redistribution. Of the 240 files, 117 carry arXiv's non-exclusive
distribution licence, which grants the right to distribute to arXiv alone, and 71 state no
licence.

The known sizes add up to 422 MB, for 208 of the files. The script skips files that are already
there, waits 3 seconds between arXiv requests (set `ARXIV_DELAY` to change this), keeps a
download only if it is a PDF, and prints a summary at the end. 209 of the URLs answered the HEAD
request with a PDF. The rest refuse scripted downloads (publisher bot checks, PMC's browser
check) or failed the HEAD request. The script lists the files it could not fetch at the end, so
they can be fetched by browser.

## Verification convention

Every source in the catalogues is marked **verified**, **partial** or **no**.

- **yes**: the bibliographic record (authors, title, venue, volume, pages, DOI) was seen during
  the research session of 27–28 October 2026, in Crossref, the arXiv API, or a publisher or author
  page. A note in the catalogue flags any field that was still filled from memory.
- **partial**: the record was seen, but some fields come from secondary sources or memory.
- **no**: the entry is from memory, or was seen only as a citation inside another source.

Of the 481 sources, 435 are verified, 5 partial and 41 no. In bibliography.bib, partial and no
entries carry `note = {unverified: ...}` with the reason. Verified entries with a field from
memory carry the caveat in `annote`. Eleven DOIs that the notes gave from memory were looked up in
Crossref while the archive was built. All of them resolve to a record with the expected title.
