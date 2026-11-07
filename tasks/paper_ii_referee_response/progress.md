# Paper II response: progress and evidence

## 2026-10-24 — framework initialized

- Registered the active clean source as
  `papers/II/submission_a&a_v2/lancet-II-v2.tex`.
- Set the revision macro to render changed text in red, with a final clean
  rendering switch. No scientific manuscript text was changed at setup.
- Preserved the reserved bold and diff sources. Their setup hashes are in
  `task.toml` for subsequent checks.
- Mapped all 13 numbered referee items and the final editorial work to
  manuscript edits, isolated diagnostics, or maintained task/script updates.
- Recorded the existing reply as an initial v1-based draft that needs
  reconciliation with v2 before reuse.
- The current v2 abstract already contains updated mass fractions; earlier
  v1 numerical concerns must be rechecked against current task outputs rather
  than copied into the revision as unresolved findings.
- Identified existing task settings for region masks and opacity conversion;
  these are candidate provenance sources for R1/R3, not newly verified results.
- Verified the task TOML, all 14 work-package entries and their existing
  task/script references. Checked the frozen-file hashes after setup.
- Compiled an isolated revision-marker example and checked its PDF text
  colors: red during preparation and black with marking disabled, with no
  bracketed revision label. This checks the macro, not the full manuscript.
- Verified that the clean-source edits are confined to the revision-marker
  configuration and a final newline. The manuscript body is unchanged.

Scientific subtasks remain pending. No new scientific test or pipeline rerun
has been performed as part of initialization.

## 2026-10-24 — VS Code response build configuration

- Reproduced `spawn latexmk ENOENT`: no `latexmk` executable is installed.
- Configured the `papers/II/` VS Code workspace to use the available Tectonic
  executable through LaTeX Workshop, preserving existing file exclusions.
- Ran the configured recipe against `response-draft.tex`; it produced the
  10-page PDF and SyncTeX file successfully. No manuscript text was changed.
- Follow-up investigation in VS Code found a file-only window with Explorer
  reporting **No Folder Opened**; it did not load the `papers/II/` settings.
  After the `II` folder workspace was opened, invoked the actual LaTeX
  Workshop build button. The extension log confirmed the `Tectonic` recipe,
  the configured executable and arguments, and `Successfully built` followed
  by `BUILD_DONE`. This verifies editor integration, not only a terminal run.

## 2026-10-25 — comment dependency review

- Added `dependencies.md` with start-versus-completion dependencies for every
  referee item, conditional rerun routes and a proposed revision sequence.
- Recommended starting with R5 PB/selection provenance and R3 boundary
  definitions. R1 terminology and R2a narrative can proceed early; final
  mass fractions depend on selection and physical-systematic checks.
- Checked current task inputs and runner routing. Identified distinct region
  mask variants and fixed fragment/cloud mass inputs as issues to account for
  when designing a consistent boundary-sensitivity test, not established
  scientific errors.
- Confirmed existing configured CASA image and measurement-set paths, so
  product inspection should precede deciding whether new imaging is needed.
- No scientific pipeline run, manuscript/reply edit or scientific status
  completion is claimed. Reserved bold/diff sources remain frozen.

## 2026-10-25 — R5/R3 product and selection audit

Subtask: R5 and the boundary-definition part of R3, with R6 implications.

Question / scope: Establish the delivered-product relationships, selection
gains, regional masks and catalogue coverage before changing the analysis.

Mode: diagnostic. No adopted scientific or manuscript/reply update.

Inputs and configuration: Existing `prepare_continuum_images`,
`source_catalogue`, `column_density`, `core_formation` and
`plot_allscale_continuum` configurations; delivered FITS products; retained
triple-getsf images/settings/components/catalogue; six configured regional
rectangles; published Paper-I 65-fragment catalogue and display ellipses.
Paths, sizes, mtimes, header hashes and units are in the diagnostic inventory.

Command / verification: Existing `py312` environment, with `PYTHONPATH=src`,
running
`results/paper_ii_referee_response/diagnostics/R5_R3/2026-10-25/audit.py`.
The command and numerical tolerances are in that directory's `report.md`.
Full-footprint comparisons, frame-aware catalogue sampling and mask checks
are recorded in `summary.json`. An initial plotting attempt failed because
native PB maps are already clipped at 0.2 and contain no numeric crossing;
the corrected overlay explicitly shows their retained-mask edges. This
display fix did not alter scientific products.

Evidence and outcome:

- Both adopted extraction inputs have exactly the same finite mask as the
  adopted G≥0.2 selection, containing 4,546,473 pixels. Their getsf observed
  component images match the current extraction inputs on this footprint.
- The 12m7m input matches the unmultiplied MIRIAD restored image after unit
  conversion. The TP input matches the **dust-only** combined product.
  Neither matches multiplication by the adopted selection gain. These file
  identities do not by themselves establish that the restored MIRIAD image
  has no residual attenuation.
- The delivered `pbant` product matches restored × the **12m-only**
  reprojected PB. That PB differs from the selection gain configured under
  `7MTM2TM1/cont/casa/primary_beam_rp2miriad.fits`.
- All 231 adopted source flux rows equal their raw getsf columns; the
  maintained catalogue code applies no later gain division. This conflicts
  with the current manuscript processing account, but does not authorize an
  additional flux correction before the residual response is established.
- All 231 core centroids pass a candidate intersection of ACA, native 12m
  and adopted-selection G≥0.2 footprints. Only 55 of 65 published fragment
  centroids pass, split 11/31/13 by their published regional labels. Adopted
  selection gain alone retains 58. No mass or completeness correction is
  inferred from these centroid counts.
- Configured original boxes are 300″ × 400″; shrink boxes are approximately
  288″ × 400″. These are distinct from published Paper-I emitting-region
  dimensions. Unshrunk boxes double-count four adopted cores when regional
  memberships are naively summed; shrink boxes retain all 231 once.
- The current core/fragment calculation is a regional summed-mass ratio
  with fixed fragment/cloud denominators, not a parent-associated ratio.
  No revised mass fractions have been calculated.

Owning task/script changes: None. The self-contained audit lives under the
required diagnostic root and is not an adopted scientific pipeline.

Affected downstream products: None changed. Response, catalogue selection,
mass fractions and environmental measurements may need changes depending on
the unresolved gain/combination provenance.

Clean-manuscript location / red revision: None. Existing `cat:extract` and
`cat:measure` descriptions need reconciliation after the actual response
and photometry chain is established; no unsupported replacement was made.

Response location: R5 and R3 remain provisional in the response draft.
The new evidence is author-internal in the diagnostic `report.md`; no claim
of a completed scientific revision was added to the reply.

Remaining work / status: R3 and R5 are `in_progress`. Need the actual J-comb
and extraction-preparation script/notebook, the MIRIAD residual response and
noise convention, Paper-I aperture provenance, fragment-boundary handling,
matched mass accounting, boundary sensitivity and maintained Fig. 2/3
selection overlays. Asked the user for the preparation-script location.
Other scientific subtasks remain pending; frozen variants remain untouched.

## 2026-10-25 — provenance availability and request-scope discussion

- User reports J-comb was performed by Sihan and the original script is not
  available. The getsf preparation was performed on a remote server that is
  presently inaccessible. Historical provenance remains unresolved; do not
  treat finding those scripts locally as the next assumed step.
- User requested a breakdown of the referee's requests before discussing
  how to address them. Set further scientific execution aside for this
  discussion; no scientific status was marked complete or blocked.
- Added `referee_requests.md`, grounded in the preserved original comment
  boxes, separating direct requests, suggested methods, explicitly optional
  TP joint deconvolution, and the additional checks proposed in our draft.
- In particular, extensive injection/recovery, parent matching, a prescribed
  boundary-perturbation suite and environmental covariance reanalysis were
  not explicitly prescribed by the referee. Their necessity depends on the
  claims and evidence retained in the revision. Existing draft promises are
  proposals, not a record of completed work or an agreed immutable scope.
- No manuscript/reply source or scientific product changed in this work unit.

## 2026-10-25 — R5 reply, manuscript and selection-boundary figures

Subtask: R5(a–b). Mode: manuscript clarification and maintained figure update.

Author confirmation: The user states that PB-corrected 12m7m and TP were
combined first; both PB-corrected 12m7m and 12m7mTP were then multiplied by
the same combined 12m7m PB for getsf extraction, followed by division of
source fluxes by that gain. This is author-confirmed processing history,
not a newly recovered original script or an independent noise-map test.
Earlier file-identity tests do not establish the PB state by themselves.

Inputs and configuration: Both maintained figure tasks now read
`data/delivered/alma/7MTM2TM1/cont/casa/primary_beam_rp2miriad.fits` and
draw its numerical 0.2 contour in its own celestial WCS. Continuum/surface
density data, display scales and free–free contour levels are retained.

Owning task/script changes: Updated `plot_allscale_continuum` and
`column_density` task settings/runners and the shared continuum plotting
helper. Replaced the former DS9 0.1 display region for the continuum view.
Added a WCS-alignment regression test. Removed duplicate overlay-axis tick
labels, moved the resolution label clear of contours, and clipped the
displayed free–free contour paths to the existing viewport for figure QA.

Command / verification: Ran `plot_allscale_continuum.py --stage combined
--overwrite` and `column_density.py --stage plot --overwrite` in `py312`.
All 11 tests in `test_continuum_plotting.py` and `test_column_density.py`
passed. Both final figure PDFs pass rendered collision checks and the
5-point glyph-floor check; their minimum detected glyph size is 7 points.
Visually inspected both final renders. Package-consistency checker returned
no findings. The desktop compiler reports success for the existing response
source and the clean v2 manuscript.

Adopted figure copies:

- Fig. 2 source: `results/plot_allscale_continuum/figures/allscale_continuum_dust_freefree.pdf`.
  Paper copy: `papers/II/submission_a&a_v2/figures/ImageCombine_Continuum_Dust_FF_TP7MTM2TM1_DualColor.pdf`.
- Fig. 3 source: `results/column_density/figures/surface_density_allscale_freefree.pdf`.
  Paper copy: `papers/II/submission_a&a_v2/figures/ImageCombine_Sigma_DualColor_OverlayFF.pdf`.

Full source/copy hashes, backups, diffs, build logs and figure QA are under
`results/paper_ii_referee_response/evidence/R5/2026-10-25/`.

Clean-manuscript location / red revision: Replaced the combination sentence
in `combine:interfTP`, the extraction PB explanation in `cat:extract`, and
the flux-correction sentence in `cat:measure`. Updated captions of
`fig:allscale` and `fig:sigma`. New/replaced text uses `\revise{...}` and
`\markrevisionstrue` remains enabled. Removed the old Fig. 2 caption's
0.1/resolution-transition sentence because the new 0.2 selection boundary
does not itself define a transition to 18″ resolution.

Response location: Replaced only R5's reply and its superseded author note
in `referee-response-v1/response-draft.tex`. The reply explains the confirmed
combination/attenuation/photometry order and reports the actual Fig. 2/3
updates. It describes a more nearly uniform rms without claiming a new
quantitative noise validation.

Affected downstream scientific products: None. No imaging, extraction,
catalogue, mass or environmental-analysis rerun was performed. Plot stages
regenerated their associated figure products; only the requested two PDFs
were copied into the paper.

Remaining work / status: R5 is `ready_for_review`, based on author-confirmed
methods and verified text/figure updates. Original-script provenance and
independent response/noise diagnostics remain distinct from this completed
clarification. R3 remains in progress. The clean manuscript is not finalized;
reserved bold/diff sources remain frozen.

## 2026-10-25 — R5 caption wording follow-up

- At the user's request, removed “, defining the footprint adopted for
  source extraction” from the Fig. 2 caption and the corresponding caption
  update statement from the R5 reply. The figures and Fig. 3 caption were
  unchanged.
- Both saved LaTeX sources compile successfully with the desktop compiler;
  the package-consistency check returns no findings.

## 2026-10-25 — visible response status markers

- Added an explicit `Solved` or `To be solved` argument to every referee
  comment box in the working response, displayed beside its title. R5 is
  `Solved`; the other 13 boxes, including the introductory marked/diff
  request, are `To be solved`. Editorial requirements are also marked
  `To be solved`.
- Markers follow the existing `\ifauthornotes` switch, so finalization can
  hide the working-status labels together with internal author notes.
- Set the R5 tracker to `complete` for the delivered clarification and
  Fig. 2/3 updates. This does not claim original-script recovery or a new
  noise-uniformity test. All other subtask states are retained.
- Verified that referee wording is unchanged and all 14 comment boxes have
  markers. The saved response compiles successfully with the desktop compiler.

## 2026-10-25 — R3 Paper-I classification and clump-extension rationale

- Author confirms that the evolutionary classification, segmentation and
  characterization were established from the original lower-resolution
  data in Paper I. The segmentation extends around three prominent massive
  star-forming clumps associated with cluster formation to enclose their
  surrounding parental-cloud material.
- Read Paper I's regional overview, Table 1 and cloud-property measurement
  procedure. Replaced the brief target overview in the clean v2 Introduction
  with red-marked explanations of the independent infrared/radio signposts,
  clump L/M values 0.7/26/156 on 0.2–0.4 pc scales, clump-extension
  segmentation, 18″ column-density selection threshold and cloud masses.
- Updated only R3's reply and author note. The reply lists the regional
  signposts and provides a table of Paper-I clump/cloud L/M, cloud masses,
  mean temperatures and moment-based cloud dimensions, with reported
  uncertainties. It explicitly distinguishes cloud measurements from
  segmentation-box dimensions and from evolutionary-age uncertainties.
- Removed the draft's unagreed future promises to add a diagnostic table
  and perform boundary perturbations. The table is now present in the reply;
  no new scientific calculation or boundary-sensitivity result is claimed.
- Both saved documents compile successfully with the desktop compiler.
  Table measurements were checked against Paper I; the package-consistency
  check reports no findings. Reserved bold/diff sources remain unchanged.
- R3 remains `To be solved` / `in_progress` while the exact mask/boundary
  provenance and remaining uncertainty wording are discussed. No data,
  catalogue, mass sums or figure products changed. Evidence and source
  snapshots are under `results/paper_ii_referee_response/evidence/R3/2026-10-25/`.

## 2026-10-25 — restore deleted opening author note

- Recovered the opening working-draft author note verbatim from the saved
  pre-R5 response snapshot and restored it before the Referee 1 section.
  Other notes and the current response wording were preserved.
- The note remains controlled by `\authornotestrue` / `\authornotesfalse`.
  The existing response compiles successfully with the desktop compiler.

## 2026-10-27 — polish author's conceptual target overview

- Located the author's conceptual edit at line 146 of the clean manuscript
  (the open response draft's line 146 contained only an environment ending).
- Polished only that target-overview passage. Preserved the sequence from
  LANCET motivation to Paper-I classification, three principal clumps and
  the extension of regional boundaries into the parental clouds. Split the
  classification and cloud-selection material into two red-marked spans.
- All numbers and surrounding author edits are preserved. The response and
  neighboring manuscript paragraphs were not edited. The clean manuscript
  compiles successfully; the package-consistency check has no findings.
- Pre-edit source and a before/after record are in
  `results/paper_ii_referee_response/evidence/R3/2026-10-27/`.

## 2026-10-27 — preserve the regional cloud-scale definition

- Corrected the earlier polish's loss of the author's key distinction.
  The manuscript now defines the subregions as clouds on approximately
  5 pc scales, each hosting a representative clump associated with a
  different stage of cluster formation. Only that sentence was edited;
  other current author edits are preserved.
- The clean manuscript compiles successfully. A pre-edit snapshot and the
  exact sentence change are recorded in the 2026-10-27 R3 evidence folder.

## 2026-10-27 — separate complete Introduction draft

- Created `papers/II/intro.tex` at the user's request as a removable,
  standalone temporary draft. It is not included in the working manuscript.
- Reorganized the Introduction into seven paragraphs covering the scientific
  questions, spatial requirements, evolutionary comparison, LANCET/G316.8,
  regional cloud definitions, quantitative classification and this paper's
  approach. Retained the author's current approximately 4 pc cloud-region
  definition, representative-clump distinction and quantitative properties.
- The draft positions combination as an established observational method
  and emphasizes extensive coverage at high resolution. It frames the
  analyses as instantaneous mass fractions and core/environment relations.
- Reused existing bibliography entries, checked all citation keys and
  compiled the draft successfully with the desktop compiler. Verified that
  the working manuscript and response remain byte-identical to their state
  before drafting. Hashes and drafting checks are recorded in
  `results/paper_ii_referee_response/evidence/R3/2026-10-27/intro_draft_guard.json`.
- This draft is for discussion only. No revision was adopted into the main
  source and no referee status was changed.

## 2026-10-27 — retain multiscale context in the temporary Introduction

- Expanded the opening of `papers/II/intro.tex` following the author's
  co-author feedback. Restored the hierarchy from GMFs and GMCs through
  filamentary clouds, clumps and dense cores, with characteristic scales
  and their roles as gas reservoirs for stars and clusters.
- Kept the research questions in a separate following paragraph and
  preserved all later draft paragraphs, including the approximately 4 pc
  cloud-region definition. Reused the existing manuscript bibliography.
- The temporary file compiles successfully. The working manuscript and
  response were not edited. A pre-edit snapshot and scoped change record
  are stored in the 2026-10-27 R3 evidence directory.

## 2026-10-27 — temporary Introduction removed

- The author finished the Introduction and adopted their selected wording
  in the main manuscript. Removed only `papers/II/intro.tex`, as requested,
  retaining an archival copy in the R3 evidence directory.
- Did not overwrite the author's adopted Introduction or change revision
  markup. The recorded rule in `papers/II/AGENTS.md` still requires red
  `\revise{...}` spans during preparation and leaves the reserved variants
  frozen until finalization. Finishing this Introduction is not treated as
  finalizing the whole manuscript.

## 2026-10-27 — mark the adopted Introduction revisions

- Compared the author-adopted Introduction with the submitted v1 baseline
  and added seven `\revise{...}` spans: the substantively replaced sentence
  within the retained opening hierarchy paragraph, plus separate spans for
  the six new/replaced following paragraphs. Unchanged opening material
  remains unmarked.
- Verified that removing only these new wrappers exactly reproduces the
  pre-edit Introduction. No wording, citations, values or other manuscript
  sections were altered. `\markrevisionstrue` remains enabled.
- The saved manuscript compiles successfully and the package-consistency
  check reports no findings. Reserved bold/diff sources are unchanged.
  Pre-edit source and the markup record are in the 2026-10-27 R3 evidence
  directory. No latexdiff or reserved variant was generated.

## 2026-10-27 — R3 closed and reply shortened

- The author accepted the evolutionary-sequence clarification as solved.
  Shortened the R3 reply to two sentences pointing to the revised Introduction,
  its quantitative diagnostics and the fuller description in Paper I.
- Marked the response comment Solved and the R3 tracker complete. Replaced
  the obsolete pending author note with a completion note and updated the
  dependency review so proposed boundary perturbations do not block progress.
- The saved response compiles successfully with the desktop compiler.
  The pre-edit response is archived in the R3 evidence directory as
  `2026-10-27/response-draft.before-R3-closure.tex`.
- Next discussion item: R4a, quantitative 12 m/7 m overlap consistency.
  No imaging or scientific pipeline was run for this wording update.

## 2026-10-27 — R4a maintained combination diagnostics

- At the author's request, renamed `compare_interferometer_arrays` to
  `combination_diagnostic` across the task, runner, result root and active
  references. Moved existing results; retained an old-name result symlink
  for historical references. Added task/result symlinks from the referee area.
- Inspected the supplied read-only continuum MSs, combination script and TP
  model. Rebuilt the uv map from unflagged cross-correlations and redesigned
  it as a horizontal full/short-baseline pair with consistent array colours.
- Added maintained amplitude, mean/summed weight and row-count diagnostics,
  full/zoom views, per-field/per-DDID tables, metadata and an averaging
  consistency control. Added a separate windowed TP Fourier/beam diagnostic.
- Measured ACA/12 m radial-support intersection near 8.8–31.6 kλ. Raw profiles
  are not a common-PB flux-scale assessment; R4a is in progress. Supplied-script
  settings differ from adopted image filenames and are documented for later
  imaging work rather than silently treated as identical.
- Nine relevant tests pass. Final rendered alignment, PDF text and collision
  checks pass for all four figures. The source validator's config-driven
  export false positives are recorded with actual artifact verification.
- Evidence: `results/combination_diagnostic/reports/R4a_initial_visibility_diagnostics.md`
  and its linked referee-result directory. No reduction input, manuscript,
  paper figure or response text was changed. Figure 1 placement is deferred.

## 2026-10-27 — R4a noise, matched-sky and PB tests

- Added the explicitly requested controlled tests to the maintained
  `combination_diagnostic` task, with reusable numerical helpers and tests.
  No scientific pipeline is duplicated in the referee coordination area.
- Empirically estimated component noise from adjacent same-baseline,
  same-field integrations and checked independent-hand cross-power.
  ACA broad-window mean amplitude decreases 18.6% after coherent averaging
  in 26–32 kλ, with negligible predicted fringe smearing. Narrow windows
  contribute a smaller additional scalar bias to the initial all-window plot.
- Matched 19 pointing pairs within 5″ and common sideband/radius/orientation
  cells. Noise-corrected ACA/12 m RMS ratios remain 1.50–1.82.
- Forward-modelled two restored sky maps through native versus common
  Gaussian PBs from the supplied script. Much of the matched excess is
  reproduced by the wider ACA response. Conditional common-PB ratios are
  0.98–1.17; a stricter three-pair check reaches 1.32, so precise calibration
  agreement is not established and R4a remains in progress.
- Fifteen tests pass. All three final diagnostic figures pass rendered
  alignment, text and collision checks; numerical sampling changes are
  below 0.91% on the checked pointing pair.
- Evidence: `results/combination_diagnostic/reports/R4a_noise_sky_pb_tests.md`
  and `results/combination_diagnostic/overlap_tests/`, accessible through the
  existing referee-result symlink. No reduction input, paper figure,
  manuscript or reply text is changed. No flux rescaling is adopted.

## 2026-10-27 — simple ACA PB-adjusted amplitude curve

- Added the author's requested before/after PB comparison as a single plot,
  reusing the completed matched-cell, noise-corrected RMS analysis. The
  original mean scalar statistic is not silently substituted or relabelled.
- Orange dashed/solid curves show native/common-PB ACA; the blue reference
  is the measured 12 m RMS translated to the same pair-centred PB. A shaded
  band shows the range of adjustments from two existing sky models. No
  fitted normalization forces agreement; residual ratios remain 0.98–1.17.
- The maintained `pb-curve` stage exports PDF/SVG/PNG and CSV without
  rereading visibilities or rerunning model FFTs. Future full overlap-test
  runs generate the same additional plot.
- The six relevant numerical tests pass. Final PDF text and collision checks
  pass, and the exported CSV ratios reproduce the completed test to numerical
  precision. The final plot was visually inspected.
- No reduction data, adopted calibration, paper figure or manuscript text
  is changed. R4a remains in progress. Evidence is in the owning task's
  `overlap_tests/aca_pb_amplitude_curve.*` files and curve manifest.

## 2026-10-27 — comprehensive Figure 1 candidate

- Assembled the author's requested four-panel figure in the maintained
  combination task: (a) pointings/matched pairs, (b) native/common-PB ACA
  and common-PB 12 m RMS amplitudes, (c) full uv coverage, (d) short-baseline
  coverage plus TP zero spacing. All panel titles and legends are internal.
- Tried TP in panel b as an explicitly labelled beam-smoothed model
  reference, through the common PB and lowest-bin matched population.
  Its native 18.1″ beam is retained, and plotting stops at 10% beam response
  below native pixel Nyquist. It is not a measured visibility/calibration curve.
- Panel d restores TP context with a zero-spacing point and model-derived
  50%/10% beam-response radii, 6.70/12.21 m at 225 GHz, rather than arbitrary
  hard support radii. These are model guides, not actual TP uv sampling.
- Exported PDF/SVG/600-dpi PNG, TP source CSV, proposed caption, figure
  contract and provenance manifest under `results/combination_diagnostic/figures/`.
  The maintained runner has an explicit `figure1` stage.
- Verified equal plot areas, PDF text size and source/numerical consistency.
  Collision audit has no FAIL; two raster-edge WARNs were visually reviewed
  and documented as acceptable internal annotation overlays. Six numerical
  helper tests pass; no new scientific calibration is adopted.
- This is a reviewable Figure 1 draft. The working manuscript, paper figures,
  old Appendix-A figure/caption, response text and frozen variants are unchanged.
  R4a remains in progress.

## 2026-10-27 — correct Figure 1 sampling scope and match TP to the 12 m PB

- The author correctly identified that the first amplitude panel used only
  overlap-selected bins, concealing the distinct sampling ranges. Rebuilt
  profiles over all available radii in the same 19 paired pointings: ACA
  5.065–31.505 kλ and 12 m 8.825–956.332 kλ. Actual baseline-range bars and
  shaded overlap now show complementarity directly.
- Retained common occupied pointing/sideband/orientation cells inside
  overlap and each array's own cells outside it. The extended signal RMS
  uses independent-hand cross-power throughout; negative/weak long-baseline
  estimates produce conditional limits rather than fabricated detections.
- Evaluated ACA model PB factors over its whole sampled range, targeting
  the actual matched 12 m pointing. Corrected overlap ratios remain
  0.981–1.172, without a fitted flux normalization.
- At the author's follow-up request, TP now has that exact target 12 m PB
  multiplied into its brightness model before FFT, using the shortest-ACA-
  bin sky weights. Its native 18.1″ restoring beam remains explicit. An
  analytic uniform-sky test verifies PB solid-angle normalization and sign.
- Added maintained `full-sampling` products and two regression tests. All
  17 relevant tests pass. Four-panel alignment and PDF text checks pass;
  collision audit has no FAIL and two reviewed raster-edge WARNs.
- Updated the Figure 1 draft, caption, task README, manifest and report.
  Evidence: `results/combination_diagnostic/full_sampling/report.md` and
  `results/combination_diagnostic/figures/Figure1_combination_diagnostics.*`.
  No observed data, adopted calibration, manuscript or paper figure is
  changed. R4a remains in progress.

## 2026-10-27 — full log uv-distance scale and no matched-pair boxes

- At the author's request, panel b now uses a fully logarithmic positive
  uv-distance axis. The zero-frequency TP value remains in the source CSV
  and the zero-spacing point remains in panel d.
- Removed panel-a matched-pair squares and their legend entry. The numerical
  analysis cohort and all amplitude/source-cell products are unchanged.
- Adjusted internal label placement for the log view and regenerated the
  PDF/SVG/600-dpi PNG. Updated the proposed caption and task documentation.
- Final alignment and PDF glyph checks pass; collision audit has no FAIL
  and the two existing panel-c raster-edge WARNs are visually reviewed.
  Verified amplitude and cell-statistics files against their saved hashes.
  No new scientific pipeline, manuscript edit or paper-figure replacement.

## 2026-10-27 — simplify Figure 1 and display observation sampling

- Kept only (a)–(d) as panel labels. Panel b now spans 1–1000 kλ on its
  log axis and omits the native-ACA curve and both shaded regions.
- Added eleven sampling tracks below the amplitude curves: eight ACA,
  one C2 and two C5 observation blocks. Dot area represents valid-row
  counts in 45 logarithmic bins using a single global area scale. This
  shows observation-level density without mixing raw scalar amplitudes
  with noise-corrected, PB-translated RMS. Counts use the same 19 pointing
  pairs and broad windows before common-cell selection, without thinning
  or imaging weights. All 279,245 rows are accounted for in the source CSV.
- Regenerated PDF/SVG/600-dpi PNG and updated the caption, contract, task
  README and manifest. Alignment and PDF text checks pass; collision audit
  has no FAIL and two visually reviewed panel-c raster-edge WARNs. Static
  preflight's helper/configuration limitations are recorded in QA notes.
- Verified unchanged amplitude/cell-statistics hashes. No scientific
  recalibration, manuscript edit or paper-figure replacement; R4a remains
  in progress. Evidence: `results/combination_diagnostic/figures/Figure1_*`
  and `results/combination_diagnostic/qa/Figure1_combination_diagnostics.*`.

## 2026-10-27 — separate C2/C5, A&A text size and TP beam-response test

- Figure 1 now shows separate C2/C5 cross-power RMS profiles from the existing
  row caches. In each radial bin the arrays present share common occupied
  cells where possible, retaining unique modes elsewhere. The displayed
  ACA profile uses the corresponding common population. Original pooled
  CSV/NPZ products are unchanged. Matched C5/C2 has median ratio 1.008 in
  14 bins (range 0.896–1.217); long-baseline limits belong to C5.
- Renamed the lower strip "Visibility counts" and removed its area key
  from the panel; the caption explains dot areas. Export is exactly 184 mm
  wide with 10 pt base/legend/axis text and 11 pt panel letters, following
  the local A&A class. Math superscripts are 7 pt. Canvas-preserving export
  retains these sizes when inserted at text width.
- The author chose matched-beam testing before any vertical normalization.
  Added the maintained `tp-beam-test` stage and `src/lancet/tp_beam_test.py`.
  Sky models are convolved to the native 18.1″ TP beam before the actual
  paired 12 m PB is applied. TP is sampled at the same measured ACA uv
  coordinates and the same weighted cells. No flux factor is fitted.
- Beam/PB-matched ACA/TP is 2.11–2.58 with the all-scale model and 0.63–0.81
  with the 12 m-only model. The strong model dependence does not establish
  a unique smooth join; Figure 1 remains on its physical amplitude scale.
  The TP interval has 4,670 selected ACA rows, 4,253 eligible rows in 489
  cells. C2/C5 have 42/27 rows and no cells meeting the four-row threshold.
  Evidence: `results/combination_diagnostic/tp_beam_test/report.md`.
- Nineteen relevant tests pass, including new common-cell/limit and analytic
  beam-kernel flux/covariance tests. Figure 1 alignment passes with zero
  collision FAIL and one reviewed raster-edge WARN; the TP diagnostic has
  zero FAIL/WARN. Both PDFs pass text checks. Input/source-product hashes
  are verified; static save-helper recognition limitations remain documented.
- Updated the task README, proposed caption, contracts, reports and QA notes.
  No observed visibility, calibration, manuscript, response claim or adopted
  paper figure is changed. R4a remains in progress.

## 2026-10-27 — pointing-resolved ACA/C2/C5/TP amplitude diagnostic

- At the author's request, added a separate plot for each of the original
  19 matched pointing pairs: orange ACA, blue C2, purple C5 and grey TP.
  Exported a large-format 19-panel overview plus a four-page A4-text-width
  atlas with six panels per page, 10 pt body text and common log-axis limits.
  Every panel includes its pair number and native field IDs; coordinates
  and matching separations are in `pointings.csv`.
- Added the maintained `pointing-profiles` stage and reusable
  `src/lancet/pointing_profiles.py`. Used the original full-range row caches
  and cross-power estimator, calculating each field separately. Retained
  minimum four-row/common-cell rules and refined the grid at 9.165 kλ.
- TP uses the exact measured ACA uv coordinates, row weights and selected
  cells within each pointing/bin; all 76 TP points have identical population
  counts and mean radii to ACA. TP retains its 18.1″ beam, and the ACA PB
  translation remains conditional on the all-scale model. No flux shift or
  restoring-beam matching is introduced by this diagnostic.
- Accounted for input, eligible and selected rows for every array in the
  result README/manifest. No pointing is selected by brightness or omitted.
  First-bin PB-translated ACA/TP ratios span 1.16–2.70 across the 19 fields;
  this is descriptive, with beam/model limitations, not a calibration fit.
- Five focused tests pass, including a new field-isolation/TP-population
  regression. All comparable layouts pass alignment; the final single-panel
  page is NOT APPLICABLE. Atlas/overview collisions have zero FAIL/WARN and
  PDF text checks pass. Inspected all panels and verified source hashes and
  exactly-once P01–P19 inclusion. Static helper-recognition limits are recorded.
- Evidence: `results/combination_diagnostic/pointing_profiles/`. No Figure 1,
  observed input, calibration, manuscript, response claim or adopted result
  changes. R4a remains in progress.

## 2026-10-28 — overlay all pointings and extend TP to the log-axis left edge

- Clarified that the earlier TP curves started near 5 kλ because they were
  evaluated only at actual ACA samples, not because TP lacks shorter modes.
- Added the maintained `pointing-range` stage, overlaying all 19 curves per
  array on a single set of log axes, with consistent colours and faint
  upper-limit arrows. The spread is the actual pointing curves; no fitted
  normalization, mean curve or confidence interval is added.
- Evaluated continuous TP model curves through each actual matched 12 m PB
  from zero to 9.165 kλ, using each pointing's fixed first-bin sideband/
  orientation-cell population. Positive 1–9.165 kλ values reach the plot's
  left edge; zero values remain in the CSV. Grey diamonds preserve exact
  ACA-sampled comparisons, whose cell population varies with bin.
- Exported `pointing_profiles_overlay.{pdf,svg,png}` at 184 × 135 mm with
  10 pt text, plus `tp_continuous_profiles.csv` and a dedicated manifest.
  Verified all 19 × 97 TP samples and 92 reference cells. The reference-cell-
  weighted ensemble reproduces the prior TP overview within 0.075% after
  radial interpolation. Source hashes match.
- Collision audit passes with zero FAIL/WARN and PDF text passes (7 pt
  superscript floor); single-panel alignment is NOT APPLICABLE. Visually
  inspected the complete figure. Static source limitations remain documented.
- Evidence: `results/combination_diagnostic/pointing_profiles/overlay.qa.md`
  and accompanying plot/source files. Previous atlases, scientific profiles,
  Figure 1, manuscript, response claims and adopted calibration are unchanged.
  R4a remains in progress.

## 2026-10-28 — adopt the compact amplitude overlay as manuscript Figure 1

- The author requested a TP upward shift, short y-axis label, tighter layout
  and first-figure placement. Applied a common ×2.7 display factor to TP
  curves/diamonds only, rounded from the median of 74 detected matched
  ACA/TP bins (unrounded 2.710615348473197). The caption discloses the factor
  and states that it is not a flux-calibration correction. No measured
  visibility, physical model/source amplitude or adopted calibration changes.
- The maintained `pointing_range.tp_display_scale` configuration and
  `pointing-range` plotting stage save a separate `_display_scaled` version,
  preserving the unscaled overlay and the existing scientific source tables.
  Y label is "Amplitude (mJy)"; its RMS definition remains in the caption.
  The compact canvas is 184 × 110 mm with 10 pt body text.
- Final graphic passes collision/text QA (zero FAIL/WARN, 7 pt math-script
  floor); one-panel alignment is NOT APPLICABLE. Visually inspected all
  curves, upper limits, labels and margins. Source-preflight save-helper
  recognition limitations are preserved in the QA report.
- Copied the verified PDF to
  `papers/II/submission_a&a_v2/figures/Visibility_amplitude_pointing_ranges.pdf`
  and inserted it before the original first figure, labelled
  `fig:visibility_amplitudes`. Added a short introduction and caption wrapped
  in `\revise`; retained `\markrevisionstrue`. Existing references renumber
  automatically. Updated only the R5 author's response references to Figures
  3/4; original referee comments retain their submitted numbering.
- Built-in compilation succeeds for the existing working manuscript and
  response. Figure-copy hashes, pre-edit backups and source diffs are in
  `results/paper_ii_referee_response/evidence/R4a/2026-10-28/`.
  Reserved bold/diff manuscript variants are untouched. R4a remains in
  progress: display normalization does not establish absolute flux agreement.

## 2026-10-28 — correct Figure 1 adoption to the four-panel version

- The author clarified that the four-panel combination figure, rather than
  the separate 19-pointing overlay, is intended as the manuscript's first
  figure. Replaced the erroneous insertion with (a) continuum pointings,
  (b) amplitude profiles/counts, (c) full uv coverage and (d) short baselines/TP.
- Carried the requested display edits into panel b: "Amplitude (mJy)",
  separate C2/C5 curves and TP ×2.7 for display only, disclosed in the
  caption. The physical TP CSV, visibilities and adopted calibration are
  unchanged. Tightened margins/gutters, retaining 184 mm width and 10 pt
  text with equal square plot areas.
- Copied the verified maintained PDF to
  `papers/II/submission_a&a_v2/figures/Continuum_combination_diagnostics.pdf`.
  Rewrote the new caption/intro for all four panels with `\revise`, labelled
  `fig:combination_diagnostics`. Removed the unused single-axis paper copy
  after preserving it in evidence; the diagnostic outputs remain available.
  Figure numbering and the R5 response's Figure 3/4 references are unchanged.
- Alignment and text checks pass. Collision audit has zero FAIL and one
  visually reviewed panel-c raster-edge WARN. Built-in compilation of the
  working manuscript succeeds. Corrected source/copy hashes and superseded
  adoption metadata are retained in the R4a evidence directory. The frozen
  variants remain untouched. R4a scientific assessment remains in progress.

## 2026-10-28 — separate three-panel horizontal layout candidate

- The author requested preserving the current four-panel version and trying
  (a), (b), (c) in one row, with the former (d) as an upper-right inset of
  (c). Added the maintained `figure1-three-panel` stage, reusing the same
  data/plotting pipeline and preserving the existing figure outputs.
- The inset spans −25 to +25 m in both u and v, highlighting ACA sampling,
  TP zero spacing and the 18.1″ beam-response guides. It has no (d) label
  or separate legend. The combined legend and caption are in/under (c).
  Parent uv axes use km and the inset uses m. Full C2/C5 sampling remains
  in the parent coverage. The original four-panel view remains available.
- Three equal square main plot areas fit one 184 mm-wide row at 10 pt body
  text, with a 60.454 mm canvas height. Shortened repeated legend labels
  and retained all eleven visibility-count rows, identified by array colour.
  Amplitudes, selection, source counts and the disclosed TP ×2.7 display
  factor are unchanged.
- Exported `Figure1_combination_diagnostics_three_panel.{pdf,svg,png}` and
  separate caption/contract/manifest/source tables. The inset is excluded
  from the equal-main-panel alignment group with a recorded reason and no
  tolerance relaxation. Main alignment passes with equal 126.811 pt widths.
  PDF text passes (7 pt superscript floor); collision audit has zero FAIL
  and three visually reviewed intentional parent-raster/inset WARNs.
- Verified saved hashes: current four-panel bundle, adopted paper PDF and
  manuscript are byte-identical. Scientific input hashes also agree. No
  TeX source or referee response was edited, and no alternative was adopted.
  Evidence: `results/combination_diagnostic/figures/` and matching QA notes.

## 2026-10-28 — apply two PDF comments to the three-panel candidate

- Reduced panel-b curve markers from 3 to 1.5 pt and line widths from
  1.4 to 0.9 pt in the three-panel version only; reduced corresponding
  upper-limit caps/markers. No data points or source values are changed.
- Removed the white inset backing and made its axes background transparent.
  Moved the inset flush to the upper/right parent edges and added exactly
  two solid black connectors from the central zoom square. Insets retain
  the ±25 m range, TP response disks and zero-spacing marker. Horizontal
  tick labels are inside the box to avoid the connectors, with 10 pt text.
- Preserved the annotated pre-review PDF/manifest in the owning QA directory,
  updated the alternative's caption/contract/manifest, and regenerated its
  PDF/SVG/PNG. Main alignment and glyph checks pass; collision audit has
  zero FAIL and three visually reviewed raster-overlay WARNs. Verified
  exactly two black connector paths and 1.5 pt markers in output metadata.
- Original four-panel bundle, adopted paper PDF, working manuscript and
  scientific inputs remain byte-identical. No manuscript or response edit
  and no scientific recalibration. Evidence: the three-panel figure files
  and `qa/Figure1_combination_diagnostics_three_panel.notes.md`.

## 2026-10-28 — remove three-panel inset labels

- Removed axis and numeric tick labels from the inset in (c), retaining tick
  marks, zoom connectors, transparent background and the ±25 m caption extent.
- Regenerated the separate three-panel exports. Alignment passes; PDF text
  meets the 7 pt minimum; collision QA reports zero FAIL and one visually
  reviewed panel-label/raster-edge WARN.
- The adopted four-panel files and manuscript remain unchanged. No scientific
  values or normalization were changed.

## 2026-10-28 — white inset background without ticks

- Set the separate three-panel candidate's inset background to white and
  removed all inset tick marks; axis labels and numeric labels remain hidden.
  The border, zoom connectors and ±25 m limits are retained.
- This changes only presentation; the adopted four-panel figure and manuscript
  are preserved.

## 2026-10-28 — adopt reordered three-panel Figure 1

- Swapped the former (b)/(c): (a) pointings, (b) full uv with the white,
  tick-free ACA/TP inset, (c) amplitudes and visibility counts. The rightmost
  x-axis tick label is aligned inward to prevent clipping at the page edge.
- Regenerated and visually checked PDF/SVG/PNG. Alignment passes; PDF glyph
  minimum is 7 pt; collision QA has zero FAIL and one reviewed warning for
  the readable (b) label over the uv raster boundary. Static export/alignment
  false positives remain recorded; export and alignment are checked directly.
- Copied the verified PDF into the v2 Figure 1 path, reordered the caption
  under `\revise`, and compiled the working manuscript successfully with
  the built-in compiler. Source/copy hashes and compilation are recorded in
  `results/paper_ii_referee_response/evidence/R4a/2026-10-28/figure1_adoption.json`.
- Preserved the original four-panel source bundle and the pre-adoption paper
  figure/manuscript under the evidence directory. Scientific data and the
  disclosed TP display normalization are unchanged; R4a remains in progress.

## 2026-10-28 — Figure 1 on page 3 and author caption edits

- Preserved the author's shortened Figure 1 caption verbatim. Removed the
  main-text sentence about 19 matched pairs and the TP display normalization,
  matching the author's removal of these details from the caption.
- Moved the double-column Figure 1 declaration to the beginning of the
  Observation section so that it floats to the top of page 3.
- The built-in compiler reports success. Rebuilt the same manuscript PDF
  in place and verified Figure 1 and its caption on PDF page 3 by text
  extraction and visual inspection. The figure data and artwork are unchanged.

## 2026-10-28 — remove redundant UV-coverage Appendix A

- Deleted the former UV-coverage Appendix A and its standalone figure from
  the working manuscript source. The other appendices renumber automatically;
  free-free emission is now Appendix A. Original figure assets are retained.
- Redirected all three former `fig:uv` references to Figure 1b. Reworded the
  TP reference to describe the actual inset beam-response guides, replacing
  the obsolete shaded-overlap description.
- Updated the response's Figure A.1 reference to Figure 1b and removed the
  fixed Appendix A assignment for the still-planned dirty-PSF appendix.
- Both sources compile with the built-in compiler. The manuscript PDF was
  rebuilt in place: Figure 1 remains on page 3, and the UV appendix is absent.

## 2026-10-29 — mosaic on-source and effective exposure times

- Measured unique observation/scan/field/time integrations from unflagged
  cross-correlations in the C5, C2 and ACA continuum measurement sets. Mean
  times per pointing are 36.288, 12.096 and 146.929 seconds, respectively.
- Evaluated Gaussian primary-beam-squared weights at the pointing centres
  (26 and 43 arcsec FWHM at 225 GHz). Median sensitivity-equivalent exposure
  is 85.34 s for combined 12 m and 286.36 s for ACA; approximation ignores
  antenna-count, baseline-flag and imaging-weight differences. Inputs and
  reproducible calculation are in `results/combination_diagnostic/reports/`
  as `mosaic_time_2026-10-29.{py,json}`.
- Added the on-source times and qualified overlap estimates to Section 2.1
  under `\revise`; the built-in compiler reports success.

## 2026-10-29 — R4a response

- Reorganized the response around Figure 1a--c and the measured 8.8--31.5
  kλ common radial range, noting the model dependence of the ACA PB translation.
- Marked R4a solved at the author's request; residual sky/PB model dependence
  remains a limitation on interpreting the apparent amplitude agreement.

## 2026-10-29 — R4b CASA/MIRIAD regional image comparison

- Added maintained `joint-images` stage and `joint_image_comparison` task
  settings. PB-corrected CASA MJy/sr is convolved to the MIRIAD 0.71 x 0.61
  arcsec restoring beam, regridded, and restricted to a shared valid PB>=0.2
  footprint. The 100-arcsec central/south cutouts use identical per-row asinh
  scales, with CASA left and MIRIAD right. No flux rescaling or background
  subtraction. Input hashes, matched FITS and display settings are recorded.
- Added the comparison to R4b in the response only, with a concise description
  of the stronger mottled positive/negative structures in the CASA reduction.
  Qualifies this as a comparison of existing reductions with differing settings.
- Evidence: `results/combination_diagnostic/joint_image_comparison/` including
  contract, exports, measurements, QA and response-copy provenance. Alignment
  PASS; minimum glyph 9 pt; collision zero FAIL/two reviewed label-edge WARNs.
- Response compiles with the built-in compiler; final comparison page inspected.
  The paper was not edited. R4b remains open for the separate dirty-PSF work
  being performed by the author and Claude; no PSF files/code were modified.

## 2026-10-29 — revise regional comparison palette and pointing

- Used the same truncated mako palette as the manuscript. Recentered the
  Central cutout lower-left to show both clusters, and South rightward to
  include more elongated filaments. Recorded centres in the task configuration.
- Regenerated the matched cutouts and PDF/SVG/PNG with the existing common
  resolution, valid-footprint rule and paired brightness scales. Repositioned
  beam markers/scale bars to avoid the blank masked edges.
- Inspected all panels; alignment PASS, glyph minimum 9 pt, collision zero
  FAIL/two reviewed readable label-edge WARNs. Updated the verified response
  figure copy and caption; no manuscript figure was changed.

## 2026-10-29 — match comparison colour scale to manuscript Figure 2

- Replaced the asinh/percentile scaling with Figure 2's full mako palette,
  square-root stretch and -20 to 500 MJy/sr limits, identical in all four
  panels. The renderer reads the original image_comparison_plot configuration.
- Regenerated and inspected PDF/SVG/PNG; alignment PASS, collision zero FAIL
  and two reviewed label-edge WARNs. Updated the response figure and caption.

## 2026-10-29 — evidence-based R4b image interpretation

- Inspected the actual delivered CASA/MIRIAD FITS histories; archived them
  with an imaging-settings comparison and assessment. CASA uses mosaic
  MT-MFS, auto-masking, threshold=300 microJy, niter=10000, nmajor=20.
  MIRIAD uses ACA-model synthetic data with adjusted weights and positive
  mossdi, cutoff=20 microJy, niters=10^7, plus a 0.1-arcsec taper.
- Expanded the response around the observed patchy positive/negative
  structures and these concrete processing differences. States incomplete
  deconvolution as a plausible contributor, not an established cause; notes
  that ACA model regularization can also smooth the MIRIAD background.
- CASA supports heterogeneous mosaics; the response no longer implies PB
  differences alone prevent joint imaging. Nominal cutoff/iteration settings
  do not establish achieved convergence; a matched restoring beam does not
  equalize dirty PSFs, masks, model inputs or weights. Image remains in the
  response only. Built-in compilation succeeds. PSF work remains separate.
- Evidence: `results/combination_diagnostic/joint_image_comparison/`
  `imaging_settings_assessment.md`, `casa_fits_history.txt`, and
  `miriad_fits_history.txt`. Main-text generalizations flagged for later
  synchronization in the retained author note.

## 2026-10-29 — focus R4b on the staged primary-beam treatment

- Reframed the image-comparison response around direct heterogeneous-PB joint
  reconstruction versus separately deconvolving ACA, applying the 12 m PB at
  each pointing, and generating synthetic short-spacing visibilities for the
  final MIRIAD combination. Clarifies that synthetic 12 m means a common PB
  convention/pointing geometry, not newly measured long baselines.
- Explains why broad ACA response to extended emission motivates this staging.
  Acknowledges CASA mosaic PB support and that the comparison cannot isolate
  PB treatment from model, weighting or cleaning differences. Removed the
  threshold/iteration details from the public response's main explanation.
- Preserved the image, caption and independently updated PSF paragraphs.
  The manuscript source and PSF products were not edited in this work unit.

## 2026-10-29 — concise Appendix A and direct PSF response

- Shortened Appendix A to the final mosaic PSF, its narrow core/broad shoulder,
  sidelobe level and adopted restoring-beam fit. Preserved Fig. A.1 and labels;
  shortened its caption while retaining the contour and annotation definitions.
- Removed the general CLEAN tutorial, alternate fitting-box experiments and
  detailed beam-volume discussion from the manuscript. Full diagnostics remain
  in the owning analysis products. Retained the distinction between the Gaussian
  restoring beam and the dirty-beam response of residual emission.
- Replaced the long PSF portion of R4b with a direct confirmation that Appendix A
  and Fig. A.1 provide the requested final combined PSF. Both working sources
  compile successfully with the built-in compiler.

## 2026-10-29 — support the Appendix A wing comparison

- Retained the 0.75-arcsec reference radius and explained its significance:
  the azimuthally averaged restoring Gaussian is below 5% (0.0335 in the
  radial-profile CSV), whereas the combined PSF is 0.1734.
- Replaced the specific C5/C2 causal attribution with the directly supported
  comparison: the peak-normalized 12 m-only PSF has mean 0.1403 in the same
  0.7--0.8 arcsec annulus, before synthetic 7 m is added. The stored combined
  and twelve_m_only robust0.5 PSF FITS share robust=0.5 and fwhm=0.1,0.1.
- Polished the two short appendix paragraphs; retained the figure and caption.

## Entry format for subsequent work

Use one dated entry per completed work unit:

```text
Subtask:
Question / scope:
Mode: diagnostic | maintained_update | manuscript
Inputs and configuration:
Command / verification:
Evidence and outcome:
Owning task/script changes:
Affected downstream products:
Clean-manuscript location / red revision:
Response location:
Remaining work / status:
```

## 2026-10-29 — R4b MIRIAD dirty PSF (Appendix A)

```text
Subtask: R4b dirty PSF of the 12m + 7m MIRIAD combination.
Question / scope: Show the combined PSF the referee asked about (narrow core on a broad component).
Mode: diagnostic + manuscript
Inputs and configuration: Shuting Lin's workflow from data/reduction/alma/7MTM2TM1/cont/script.zip
  (uvrandom synthetic 7m visibilities, Tsys 14500 K), rerun through the final invert in
  miriad_r4b/work_shuting/ with the miriad80 Docker image (MAXANT=80) and full-channel 12m UVFITS.
  The ms/combine_cont.sh copy (ACA field-50 tracks, Tsys 1800 K) and one-channel UVFITS do not
  reproduce the delivered beam (1.10" x 1.04").
Command / verification: combination_diagnostic --stage dirty-beam; restor beam at robust 0.5 is
  0.707" x 0.600" PA 47.2 vs 0.712" x 0.606" PA 47.4 in the paper image. Robust 1 and 2 are 4-5%
  narrower than the delivered versions (unexplained). tests/test_dirty_beam.py passes.
Evidence and outcome: results/combination_diagnostic/dirty_beam/ (figure, radial_profile.csv,
  summary.json). PSF shoulder: azimuthal mean 0.17 at 0.75", first null 2.0"; mostly present in
  the 12m-only PSF (C-5 + C-2); synthetic 7m adds ~0.03. Dirty/clean volume 2.1 within 2".
Owning task/script changes: combination_diagnostic dirty-beam stage, src/lancet/dirty_beam.py
  (commit e96356d on combination-diagnostic-dirty-beam).
Affected downstream products: none; diagnostic only.
Clean-manuscript location / red revision: new Appendix A and Fig. A.1 (placed first, restoring
  the v1 letters B-F for the other appendices); Sect. 3.1 sentence on an
  "approximately Gaussian dirty beam" replaced.
Response location: R4b response box, PSF paragraphs after the CASA-MIRIAD comparison.
Remaining work / status: check the residual level of the adopted image before keeping the
  residual-flux statement; author review of wording.
```

## 2026-10-29 — R4b residual check on the MIRIAD rerun (inconclusive)

```text
Subtask: R4b residual level behind the Appendix A residual-flux statement.
Mode: diagnostic
Inputs and configuration: final mossdi/restor of the work_shuting rerun (run_finaldeconv.sh;
  niters 1e7, cutoff 20 uJy, gain 0.05, positive). Stopped on niters (10,113,307), as did the
  delivered run (10,034,751); residual peak 0.41 mJy stays above the cutoff.
Evidence and outcome: work_shuting/residual_check.json and residual_check_sources.csv.
  Rerun residual: MAD rms 76 uJy, median +236 uJy, 31% of the footprint flux. Per source
  (231, r ~ 0.93"): residual fraction median 13%, 95th pct 40%; implied JvM bias median 3%,
  95th pct 15%. BUT the rerun does not reproduce the delivered image beyond the beam: its
  dirty map has a +360 uJy median pedestal and its clean image a +506 uJy median (delivered
  +129 uJy), with positive CLEAN speckle across the footprint; per-source aperture fluxes are
  0.85-1.56x delivered (16-84th pct). The short-spacing (synthetic 7m) flux differs, likely
  through the ACA input/model. The PSF is unaffected (it does not depend on model flux).
Remaining work / status: residual numbers describe the rerun only. Need the residual (and
  model) FITS of the adopted robust-0.5 image from Shuting Lin to quantify the adopted image.
Author decision (2026-10-29): stop here and adopt Shuting Lin's delivered image; R4b needs
  only the PSF. All rerun intermediates were deleted. Kept in miriad_r4b/work_shuting/: the
  beam cube, mospsf PSF, imfit log, 12m-only PSF, run scripts/log and Shuting's script copy.
```

## 2026-10-29 — single-column Appendix A figure typography

- Sized the PSF figure at 88 mm width for the existing one-column inclusion.
  Used 8 pt axis/panel labels and 7 pt legend/colourbar text at final size.
  Retained both panels, moved the colourbar above (a), and shortened labels.
- Regenerated PDF/SVG/PNG, inspected both panels, and verified alignment PASS,
  minimum glyph 7 pt, collision PASS (zero FAIL/WARN). Updated the verified
  paper figure copy and compiled the working manuscript successfully.
- No PSF analysis inputs or numerical methods changed. Evidence and copy
  provenance: dirty_beam/single_column_qa.md and single_column_adoption.json.

## 2026-10-29 — fix figure placement inside R4b response

- The author added a numbered `figure` float inside the breakable response
  box, which produced `Not in outer par mode`. Retained the author's wording,
  image, caption and label; used a non-floating minipage and `captionof` so
  the image and numbered caption remain together inside the box.
- Set the caption to small text to remain compatible with the local partial
  font cache. Built-in compilation and the existing local cached build succeed.
  Checked the resulting document: Figure 1 and its caption are together on
  page 5, and the R4b cross-reference resolves correctly.

## 2026-10-29 — zoom Appendix A PSF map

- Changed panel (a) from +/-3 to +/-1.5 arcsec, retaining the 0--3 arcsec
  radial profile in panel (b) and 88 mm single-column sizing. Removed the
  negative-contour caption because those levels lie outside the cropped map.
- Regenerated and visually inspected the figure; alignment/collision PASS,
  minimum glyph 7 pt. Radial profile and PSF summary hashes are unchanged.
  Updated the hash-verified manuscript PDF copy and its marked caption.

## 2026-10-29 — remove panel-label backgrounds from Fig. A.1

- Removed white backgrounds behind (a) and (b), then regenerated and visually
  inspected the figure. Alignment and collision audits PASS; the updated
  figure copy and its hash are recorded in single_column_adoption.json.

## 2026-10-29 — trim Appendix A figure margins

- Removed excess top and bottom whitespace from the single-column PSF figure,
  retaining its 88 mm width and the physical panel/text sizes. Inspected the
  output; alignment and collision PASS and minimum glyph remains 7 pt.
- Updated the manuscript figure copy and its recorded hash.

## 2026-10-29 — mark R4b solved

- Marked referee comment 4b solved at the author’s request and synchronized
  the task status. Existing response text and author notes were retained.

## 2026-10-29 — begin R4c source and combination audit

- Measured the existing paired catalogue: 169 sources pass both peak-S/N>=3
  selections; median change +1.49%, untrimmed SD 17.59%, clipped |change|<30%
  SD 8.17% (154 sources); 13 above +30% and 2 below -30%. The prior 10%
  statement describes rounded clipped scatter.
- Verified the nominal feather high-resolution input matches the delivered
  interferometer image pixel-for-pixel. Existing feather sdfactor=1 differs
  markedly from delivered J-comb on the common gain>0.2 footprint. At a
  common 18.1-arcsec beam the preliminary median J-comb / supplied TP model
  is 2.83, while nominal feather / TP is 1.28 (correlated sample pixels;
  conservative spatial masks). Exact J-comb inputs and scale remain unknown.
- Recorded full methods, selection, diagnostic limits and reproducible
  read-only script/JSON under diagnostics/R4c/2026-10-29. No physical scaling
  or manuscript/response change was made. The 2.7 multiplier in Fig. 1 is
  display-only and its current caption omits it; it must not be interpreted
  as absolute flux agreement. R4c remains in progress.

## 2026-10-29 — quantify test-field catalogue entries for R4c

- Inspected the retained 50-arcsec test-field catalogues. They contain 19,
  26, and 28 entries for 12m, 12m7m, and 12m7mTP. Each contains the same
  source near the field boundary. Requiring peak positions to be at least
  approximately 0.5 arcsec from the boundary gives 18, 25, and 27 entries.
- Replaced the qualitative count statement in the working manuscript with
  these explicitly selected counts and updated the corresponding response.
  The historical reason for the original 18 and 25 counts is not known; the
  stated boundary cut is an explicit, reproducible selection. See the R4c
  diagnostic report for inputs and details. R4c remains in progress.

## 2026-10-30 — clarify the R4c flux comparison

- At the author's request, rephrased and marked the existing qualitative
  comparison of integrated fluxes and peak intensities in Appendix D, and
  quoted the revised sentence in the referee response. Retained the existing
  test-field counts and the qualification on the clipped 10% flux scatter.
- Compiled both LaTeX sources and verified that the response quote matches
  the revised sentence. The comment's status remains unchanged.

## 2026-10-30 — clarify J-comb and feathering

- Added one marked sentence in Sect. 3.2 describing J-comb's weighting and
  tapering scheme and its intended nearly Gaussian beam response, following
  the existing citation to Jiao et al. (2022). The manuscript compiled and
  the manuscript/response consistency check passed.

## 2026-10-30 — mark R4c solved

- Marked R4c solved at the author's request and updated its task status.
  The diagnostic record retains the limitation that the original J-comb
  inputs and settings are unavailable.

## 2026-10-31 — SIMBA bandpass and bright-peak comparison

- Removed the obsolete bandpass reference and caption text, corrected the
  retained comparison caption's Jy/beam units, and explained in the response
  that the original bandpass data and response units are unavailable.
- The author confirmed that 0.20 is the observational SIMBA calibration
  factor. Added a reproducible 250 GHz, 24-arcsec bright-peak and local-
  background comparison under the existing continuum synthesis task, and
  added the comparison figure and bounded interpretation to the appendix.
- Verified the figure, matching manuscript-copy hash, and LaTeX compilation.
  Methods, selections, annulus sensitivity, limitations and provenance are
  recorded in `diagnostics/R10/2026-10-31/report.md`. R10 remains in progress.

## 2026-10-31 — SIMBA spatial-frequency comparison

- Replaced the peak figure in Appendix C with a reproducible Fourier-amplitude
  comparison on a common apodized footprint, and updated the response.
- Recorded the scale-dependent discrepancy and beam, noise and window limits
  in `diagnostics/R10/2026-10-31/fourier-report.md`. R10 remains in progress.

## 2026-10-31 — close SIMBA response

- Replaced the manuscript Fourier figure with the pixel comparison in
  Figure C1 and marked comment 10 solved at the author's request.

## 2026-10-31 — R6 footprint test

- Applied the adopted 12m+7m PB >= 0.2 centroid selection to Paper I
  fragments and updated the corresponding overlays and regional mass plot.
- Updated the manuscript and response with the measured selection and its
  limits. R6 remains in progress: a centroid cut does not measure completeness
  or clip boundary-crossing masses; cloud masses remain full-region totals.
- Inputs, selection, mass sums, figure provenance and baseline products are
  recorded under `diagnostics/R6/2026-10-31/`.

## 2026-10-31 — SIMBA comparison figure

- Added the SIMBA-selected pixel comparison to Figure C1 and removed the
  separate Fourier figure and its manuscript discussion at the author's request.
  The Fourier analysis remains available as a diagnostic. R10 remains in progress.

## 2026-10-31 — close R6

- Marked comment 6 solved at the author's request following the full-footprint
  selection, revised figures and mass fractions. Next discussion item: R7.

## 2026-10-31 — R7 smooth CORNISH diagnostic

- Integrated the optically thin 223.8 GHz extrapolation of the native
  2.5-arcsec CORNISH map in saved getsf 12m7m footprint ellipses.
- Recorded the 94 covered southern sources, signed ratios and the distinction
  between aperture template emission and background-subtracted source flux.
  See `diagnostics/R7/2026-10-31/report.md`. R7 remains in progress.

## 2026-10-31 — R7 radio noise selection

- Added an artifact-sensitive local radio noise map and native positive
  3-sigma pixel mask before core-aperture assignment. Updated the paired
  diagnostic and saved window-size sensitivity checks. See
  `diagnostics/R7/2026-10-31/noise-selection-report.md`.

## 2026-10-31 — R7 fixed threshold and original-map fluxes

- Applied native CORNISH >5 mJy/beam selection and original 12m7m aperture
  fluxes without background subtraction; updated both diagnostic panels.
- The adopted profile and results are in
  `diagnostics/R7/2026-10-31/fixed-threshold-original-map-report.md`.

## 2026-10-31 — R7 nondetection limits

- Adopted the author's uniform CORNISH rms of 1.5 mJy/beam, updated the
  contour/detection threshold and added aperture-integrated 3-sigma limits
  for the 67 covered radio nondetections. See
  `diagnostics/R7/2026-10-31/uniform-rms-upper-limits-report.md`.

- R7 follow-up: changed nondetection markers to nominal 1-sigma limits at the
  author's request; detections and radio contours retain the 3-sigma cutoff.

## 2026-10-31 — R7 mass-fraction sensitivity and response

- Extended the radio test to retained ACA leaf masks with beam matching and
  matching PB convention. Updated marked manuscript passages, Appendix B,
  Caveats and the referee response with both-reservoir sensitivity results.
- See `diagnostics/R7/2026-10-31/mass-efficiency-report.md`. R7 remains in
  progress pending author review; adopted mass catalogues are unchanged.

## 2026-10-31 — adopt CORNISH rms 1.0 mJy/beam

- Promoted the 1.0 mJy/beam profile to the default, regenerated core/fragment
  diagnostics, mass sensitivity and response figure, and synchronized the
  manuscript and response. Both documents compile.
- Evolved-region model mass changes are 9.9% for fragments and 3.3% for
  cores; fragment/cloud becomes 14.9%, core/cloud 3.00%, core/fragment 20.1%.

## 2026-10-31 — close R7

- Marked comment 7 solved at the author's request after the adopted
  1.0 mJy/beam CORNISH test, mass-fraction sensitivity and concise response.
- Next discussion item is R8, beginning with sizes and unresolved sources.

## 2026-11-01 — R8a deconvolution audit

- Verified that the maintained catalogue uses circularized quadratic beam
  subtraction with the 1.25 beam-area threshold and a half-beam upper size
  for unresolved sources. All stored size values and flags are reproduced.
- Added the explicit prescription to the marked manuscript and response,
  and clarified the interpretation of unresolved-source densities.
- See `diagnostics/R8a/2026-11-01/report.md`; statistical treatment of size
  limits in later fits remains to review.

## 2026-11-01 — R8d temperature stress test

- Tested all 231 cores from their adopted Tenv to author-selected 40/50/60 K
  endpoints for Young/Intermediate/Evolved, using full Planck functions and
  recomputed optical depths while retaining the environment at Tenv.
- Core sums fall by 74.4%, 68.4%, 64.6%. Both regional ratio trends persist
  along the coherent heating sequence, but arbitrary independent regional
  temperatures can change their ordering. Absolute masses are sensitive.
- Saved the per-core and regional scenarios, affected correlations and a
  verified two-panel diagnostic in `diagnostics/R8d/2026-11-01/`.
- The ALMA-IMF studies motivate heating scenarios but do not establish the
  chosen endpoints as measured or universal upper temperature limits.
  Adopted catalogues, manuscript, response and publication figures unchanged;
  R8 remains in progress for author review.

## 2026-11-01 — R8d temperature sampling

- At the author's request, extended the endpoint test to 100000 realizations
  of uniform core temperatures and an approved brightness-dependent beta
  heating prior, with the same 40/50/60 K regional bounds.
- Kept the environment at Tenv and recalculated core optical depths. Added
  a backward-compatible environment-temperature argument to the maintained
  opacity function; six physical-property tests pass.
- Brightness-model independent draws give core sums of 17.37, 119.63, 118.18
  solar masses (medians); both median regional ratio trends persist. Shared
  regional heating gives wider ranges and weaker ordering retention, which
  is recorded rather than treated as independent per-core uncertainty.
- Replotted the Figure-9-style hierarchy under both priors and both
  dependence choices. Reports, draws, intervals and verified figures are in
  `diagnostics/R8d/2026-11-01/temperature_sampling/`. Adopted catalogues,
  publication figures and TeX sources unchanged; R8 remains in progress.

## 2026-11-01 — Figure 9 temperature integration trials

- Created a compact overlay of the adopted hierarchy and brightness-prior
  median/16–84% range, and a two-panel alternative adding an explicit
  comparison of core/fragment fractions under both temperature priors.
- Both trials use existing verified sampling products, retain the original
  reservoir selection and keep counting and temperature errors distinct.
  Source/preflight and rendered checks are recorded with the trial captions.
- Products: `diagnostics/R8d/2026-11-01/temperature_sampling/figure9_trials/`.
  Manuscript Figure 9, TeX and adopted mass catalogues remain unchanged.

## 2026-11-01 — revise temperature endpoints to 30/40/50 K

- Replaced the diagnostic upper temperatures with 30 K Young, 40 K
  Intermediate, 50 K Evolved, at the author's request. Both the endpoint
  and sampling runners now take these defaults from the same configuration.
- Archived all preceding 40/50/60 K products and their task configuration in
  `diagnostics/R8d/2026-11-01_caps_40_50_60/` before regeneration.
- Recomputed deterministic endpoints/correlations, four 100000-realization
  prior/dependence cases, all five figure layouts, numerical reports and
  trial captions. Brightness-prior independent medians are 22.00,142.96,137.80
  solar masses; core/cloud fractions are 0.200%,0.841%,1.722%, and
  core/fragment fractions 28.03%,21.47%,10.42%.
- Both median trends persist. The core/fragment order occurs in 99.550% of
  independent uniform draws and all sampled independent brightness draws;
  the shared-regional alternatives give 71.153% and 83.498%, respectively.
  These fractions are conditional on the chosen priors, not confidence
  levels or unconditional guarantees.
- Six physics tests pass. All five regenerated PDFs pass text and rendered
  collision checks and all comparable panels pass alignment. Visually
  inspected each full layout and panel. Scientific input checksums remain
  unchanged. Main manuscript, response, publication Figure 9 and adopted
  physical-property catalogue remain untouched pending the layout decision.

## 2026-11-01 — incorporate temperature systematics in manuscript and reply

- Added the 30/40/50 K uniform/brightness-prior calculation and quantitative
  mass/fraction results to Sect. 5, including the wider shared-heating ranges.
  Updated the measurement-error description and compressed the temperature
  caveat. All new or replaced manuscript prose uses `revise`.
- Adopted the compact Figure 9 overlay, with source/copy provenance recorded.
  The maintained core-formation efficiency stage and trial plotter now share
  the renderer, with mass/denominator consistency checks. Catalogue masses
  remain the Tenv baseline.
- Response 8d now describes the completed analysis and verified locations,
  using the author's anti-defensive-writing preference. R8 stays in progress
  for its other subitems.
- Both working TeX sources compile with the built-in compiler; package
  consistency has no findings, final Figure 9 text/collision checks pass,
  and all 12 physical-property/core-formation tests pass. Evidence is in
  `diagnostics/R8d/2026-11-01/manuscript_update/`.

## 2026-11-01 — exclude unresolved size and density limits

- Confirmed that the 22 N-flagged sizes had been carried numerically into
  density summaries and correlations without a downstream resolution flag.
  Propagated the flag through physical/environment catalogues and restricted
  size/core-density summaries and affected pairs to the 209 resolved cores.
- Updated Sect. 4.3, bound markers in Tables 1–2, Figure 10, the correlation
  discussion/conclusion and response 8a. Corrected stale physical-property
  summary values using the current catalogue. All new manuscript text is red.
- Unchanged integrated masses retain all 231 cores; the main intercore fit
  retains its original 228 positions. Temperature density diagnostics and
  the radius-based BE task now use resolved cores; mass-function exports
  retain their mass sample and carry the flags.
- Verified 38 passing relevant tests and saved before/after products and
  invariants in `diagnostics/R8a/2026-11-01/unresolved_exclusion/`. The figure's
  three automatic mathtext bounding-box warnings match pre-update labels
  and were checked visually. R8 remains in progress for other mass-model
  issues; this exclusion policy is implemented.

## 2026-11-01 — correct the surface-density selection

- The author correctly noted that surface density is derived from peak
  intensity rather than deconvolved FWHM. Restored all 231 beam-averaged
  surface-density/optical-depth measurements and their applicable pairs.
  The earlier intrinsic-surface interpretation was inappropriate for this
  reported quantity; the prior audit is marked superseded on this point.
- Size and volume-density statistics and correlations retain the 209-source
  resolved selection. Mass/surface uses 231, surface/local-intercore 228,
  and pairs involving volume density use up to 209.
- Removed the spurious surface/opacity lower-limit markers and corrected
  Sect. 4.3, Figure 10, the discussion/conclusion and response 8a. Current
  evidence is in `unresolved_exclusion/surface_density_correction/`.

## 2026-11-01 — audit the apparent number-density median change

- Confirmed the pre-exclusion 231-core catalogue median was 6.4847e6 cm^-3,
  and the 209 resolved-core median is 5.2427e6, a 19.15% selection effect.
  The original prose value of 8e7 was already inconsistent with the retained
  full catalogue and its recovery checkpoint.
- Checked units, rho-to-n conversion, nominal mass/radius reconstruction,
  the initial committed radius prescription, and the retained manuscript ZIP
  archives. The FWHM-as-effective-radius convention was unchanged; a
  hypothetical half-FWHM radius gives an eightfold shift but still does not
  reproduce 8e7. No new source or method change was made during this audit.
- Evidence: `diagnostics/R8a/2026-11-01/number_density_median_audit/`.

## 2026-11-01 — close comment 8

- Revised the median statement to explicitly identify 5.2e6 cm^-3 for the
  209 resolved cores. The response records the 6.5e6 preselection catalogue
  median and the correction of the old 8e7 prose value.
- Rechecked aligned flags, units, rho-to-n conversion, per-quantity samples,
  optical-depth wording and manuscript/reply consistency. All 38 relevant
  tests pass and both working documents compile.
- Marked comment 8 Solved and R8 complete at the author's request. Evidence
  and final source snapshots are in `diagnostics/R8a/2026-11-01/closure/`.

## 2026-11-01 — start comment 9

- Next item is the environmental-estimator distinction (9a), followed by
  multiple velocity components in the CO moments (9b).
- Code inspection confirms that the local intercore background is sampled
  at the getsf peak pixel, while the environmental estimator is the median
  of valid positive map values within a 15-arcsec-radius aperture. The
  compact core's own surface density is a separate peak-opacity quantity.
  No R9 manuscript or response change has been made yet.

## 2026-11-01 — clarify the environment measurements (R9)

- Clarified in the working manuscript and the captions of Figs. 10 and 11
  that the local intercore quantity is derived at each core peak, whereas
  the environmental quantity is the median of positive finite map pixels
  in a 30-arcsec-diameter aperture. The scaled MAD describes variation within
  the aperture. Distinguished the peak-sampled temperature used in the core
  calculations from the environmental aperture-median temperature.
- Added the author's inspection finding of a coherent C18O component along
  the filament spine and described moment 2 as the observed line-of-sight
  velocity dispersion. This statement does not assert that the moment mask
  separates all velocity components.
- Aligned the abstract's slope with the retained fit and Fig. 11 (1.22 +/−
  0.07), and corrected the interpretation of its intercept. Replaced the
  speculative response with a concise account of the completed clarification;
  the author note remains. The response marker remains To be solved pending
  author review.
- Both documents compile with the desktop editor, and the response/manuscript
  package consistency check reports no findings. Source snapshots from before
  these changes are in `diagnostics/R9/2026-11-01/definition_clarification/`.

## 2026-11-01 — close R9 and begin R1

- Marked comment 9 Solved and R9 complete at the author's request; the
  response compiles with the desktop editor.
- Next is comment 1 on the physical scales and interpretation of mass
  ratios. The revised Introduction and Fig. 9 caption already give the
  approximately 4, 0.15 and 0.01 pc scales, while the abstract, Section 5
  and conclusions still use efficiency/transfer language to varying degrees.
- Checked the current core-formation summary: the cloud-to-fragment,
  cloud-to-core and fragment-to-core ratios reproduce the currently quoted
  regional values. Fragment selections use whole dendrogram footprints within
  the PB boundary; cloud masses refer to full subregions. No R1 manuscript or
  response prose has yet been changed in this step.

## 2026-11-01 — revise R1 terminology and scales

- Following the author's preference, the working manuscript now uses mass
  transfer efficiency for the observed regional mass ratios, states the
  characteristic 4, 0.15 and 0.01 pc scales, and defines the cloud-to-core,
  cloud-to-fragment and fragment-to-core ratios in Section 5. Cited the
  existing Motte2026 bibliography entry. Updated the abstract, Introduction,
  Figure 9 caption, conclusions and concise response; retained the author
  note and the To be solved marker for review.
- Rechecked the three ratios and their product using the unrounded current
  core-formation summary. The intermediate cloud-to-core ratio is 1.5051%,
  rounded to 1.51% in the manuscript (previously 1.50%). Figure 9 already
  displays mass/cloud mass and the three physical scales; no redraw required.
- Both TeX documents compile with the desktop editor and the package
  consistency check reports no findings. Before-edit snapshots are in
  `diagnostics/R1/2026-11-01/terminology_revision/`.

## 2026-11-01 — close R1 and address R2a

- Marked comment 1 Solved and R1 complete at the author's request.
- Added the observational cost of large-area, high-resolution mosaicking to
  the existing Introduction paragraph. The paragraph already acknowledges
  single-dish/interferometer combination literature, while the Introduction
  describes the common filament and its three environments.
- Replaced the R2a future-tense draft with a concise account of these changes.
  R2a remains To be solved pending author review. Both TeX documents compile
  and the manuscript/response consistency check reports no findings.
- Pre-edit sources are under `diagnostics/R2a/2026-11-01/coverage_motivation/`.

## 2026-11-01 — add large-mosaic literature to the Introduction

- Added the author's newly saved GM2013, Battersby2020, Ginsburg2026 and
  Beuther2022 references to the existing data-combination paragraph, alongside
  Diaz2023 and EMERGE-I. Continuum and molecular-line precedents are cited
  separately. All six bibliography keys are present exactly once.
- The working manuscript compiles; the manuscript/response consistency check
  reports no findings. Comment 2a remains To be solved for author review.

## 2026-11-01 — close R2a and begin R2b

- Marked comment 2a Solved and R2a complete at the author's request; the
  response compiles with the desktop editor. R2b is in progress.
- Rechecked the retained R4c evidence. In the 50-arcsec test field the counts
  are 18, 25 and 27 with the documented boundary cut. For 169 paired sources
  the median TP-minus-interferometer integrated-flux difference is +1.49%;
  the rounded 10% scatter describes the clipped central sample.
- Traced the reported physical properties: Ibkg2, the TP-containing
  source-removed background, enters the core optical-depth calculation, and
  the derived optical depth enters the mass correction. Consequently the
  flux comparison alone does not quantify the effect on the final masses.
  The core-formation task also uses Paper-I cloud mass denominators, and the
  intercore map comes from the 12m7m source-removed image. No R2b manuscript
  or response prose has been changed during this initial audit.

## 2026-11-01 — explain the scientific role of TP (R2b)

- Replaced the R2b response with the author's two main points: TP supplies
  large-scale ambient emission and the background used in the core
  radiative-transfer calculation, while the compact-source counts/fluxes
  change modestly. Preserved the author note and To be solved marker.
- Calculated nominal optical depths for all 231 catalogue sources at their
  adopted temperatures and peak intensities, comparing Ibkg2 with zero
  background. The optical depth increases by median 1.26%, maximum 3.15%.
  This is not an independent full reduction without TP. Distinguished this
  change from the catalogue's maximum 1.084 mass opacity-correction factor.
- At 20 K and 3000 MJy/sr background, the existing Appendix E model gives
  14.8–19.8% increases in core optical depth over the zero-background case
  for the 1–100 mJy/beam examples. Clarified this stronger-background example
  and added the catalogue comparison in Section 4.3.
- Reproducible calculation, table, summary and source snapshots are in
  `diagnostics/R2b/2026-11-01/background_effect/`. Both TeX files compile and
  the manuscript/response consistency check returns no findings.

## 2026-11-01 — emphasize TP's scientific role across the paper

- Traced the all-scale maps, source-removed backgrounds and environmental
  inputs. Strengthened the abstract, Section 3.2, source decomposition and
  conclusions to explain the recovered large-scale context and the background
  used for the core radiative-transfer calculation. The environmental
  aperture analysis uses the 12m7m intercore map and Paper-I SED surface
  density, so it is not represented as an all-scale TP measurement.
- Replaced a passage claiming that <10% of compact-core flux is recovered by
  interferometers with the cited paper's supported qualitative result:
  compact structure is preferentially recovered while extended emission is
  filtered. See Miao & Li, arXiv:2607.12396.
- Preserved the author's current referee response and its unresolved marker.
  The main manuscript compiles, and the package consistency check returns no
  findings. A pre-edit snapshot is in
  `diagnostics/R2b/2026-11-01/tp_narrative/`.

## 2026-11-01 — close R2b

- Marked comment 2b Solved and R2b complete at the author's request. The
  referee response compiles with the desktop editor.

## 2026-11-01 — prepare manuscript versions

- Preserved the marked manuscript as `lancet-II-v2-mark.tex` and removed
  revision markup from the working `lancet-II-v2.tex`. An exact snapshot of
  the original marked source is in `diagnostics/finalization/2026-11-01/`;
  the marked and clean files use TeX punctuation for compiler compatibility.
- Generated `lancet-II-v2-latexdiff.tex` directly from the submitted v1
  source and the clean v2. The deleted v1 figure references were rendered as
  their original figure number. All three compile with the built-in editor
  and the VS Code workspace's LaTeX Workshop recipe.
- The response calls the C18O line (1–0) in comment 9b, while the manuscript
  identifies the SEDIGISM data as (2–1); this requires author review.

## 2026-11-01 — update clean and selective latexdiff versions

- Used the author's updated `lancet-II-v2-mark.tex` as the source. The clean
  manuscript retains its content with `\revise` wrappers removed; the marked
  source was not modified.
- Generated a selective latexdiff from the submitted v1 and the updated v2.
  Highlighting is confined to the author's marked spans, checked against the
  referee topics; unmarked changes retain the revised text without markup.
  In particular, author-unmarked environmental estimator text remains
  unhighlighted. The selective diff is not a complete record of all v1–v2
  changes.
- The latexdiff `DIFNOMARKUP` implementation retained old deleted tokens in
  suppressed regions; patched the local generator to discard those tokens.
  A normalized source projection matches the clean manuscript exactly, and
  a markup scan reports no highlight commands outside selected spans.
  Scripts, inputs, the generator patch, and previous output snapshots are in
  `diagnostics/finalization/2026-11-01/selective_latexdiff/`.
- Both requested files compile with the built-in compiler and the configured
  VS Code LaTeX Workshop recipe. The package consistency check returns no
  findings.

## 2026-11-04 — adopt the multiwavelength strip as Figure 1

- Copied the author-approved figure from
  `results/multiwavelength_strip/lancet_multiwavelength_strip.pdf` to
  `papers/II/submission_a&a_v2/figures/Multiwavelength_strip.pdf` without
  changing its contents. Maintained sources are
  `tasks/multiwavelength_strip/task.toml`, `scripts/multiwavelength_strip.py`,
  and the VVV configuration in `tasks/vvv_mosaic/task.toml`. Adopted PDF
  SHA-256: `6a7c202535b4e5af4c4a50a5f19110aaf28a44f2803c8d7747d55384834bd66d`.
- Added Figure 1 to the clean manuscript, spanning both columns on a
  dedicated page 3. The image preserves its aspect ratio and occupies
  0.87 text height, with a caption identifying the five bands/composites,
  the evolutionary sequence, orientation, scale, and TP continuation
  outside the interferometric coverage. The Introduction first cites it
  while introducing the three regions.
- The combination diagnostic is now Figure 2. Updated the author-response
  references to the shifted manuscript figure numbers, preserving the
  referee's original comments. The marked and selective-diff variants
  were not regenerated; future regeneration must retain this addition
  from the working clean source.
- Both edited TeX documents compile with the built-in compiler. The
  manuscript also compiles with the configured LaTeX Workshop recipe;
  pages 2–4 were visually inspected, with no oversized-float or undefined-
  reference warnings. The package consistency check returns no findings.
  Pre-edit sources are in
  `diagnostics/figure1_multiwavelength/2026-11-04/`.

## 2026-11-05 — Henrik coauthor comments

Compared the annotated sent PDF with the current clean manuscript and applied
remaining textual clarifications: evolutionary directions, mosaic layouts,
unresolved density limits, explicit intercore measurement scales and their
interpretation, and selected wording changes. Missing citations were already
resolved. Per the author's instruction, the execution-block inset remains
unnumbered; its caption now explains the categorical rows. No figure changed.
The built-in editor compilation succeeded. Comment audit and pre-edit source:
`results/paper_ii_referee_response/diagnostics/henrik/2026-11-05/`.
Marked/latexdiff versions were not regenerated.

## 2026-11-05 — approximation notation

Standardized 24 separated approximation-number pairs in the working clean
manuscript, including captions and appendices: `$\sim$9` becomes `$\sim9$`.
Kept angular symbols and percentages within the same math expression and
used nonbreaking spaces before physical units. Range dashes remain in text
mode. Existing statistical distribution notation is unchanged. Built-in
compilation succeeded; backup and diff are in
`results/paper_ii_referee_response/diagnostics/latex_approximation/2026-11-05/`.

## 2026-11-05 — synchronize marked and selective latexdiff versions

On the author's instruction, synchronized the marked and selective latexdiff
versions with the current clean source, including the multiwavelength Figure 1,
Henrik edits and approximation notation. Mapped the 90 existing referee-related
revision spans to the current text; the author-added multiwavelength figure
remains unmarked. Used submitted v1 as the latexdiff baseline with markup
suppressed outside these selected regions. The accepted latexdiff projection
matches all 23,308 current body tokens; stripping the marked wrappers gives
an exact copy of the clean source. No highlight commands occur outside selected
regions. Both variants compile successfully with the built-in editor compiler.
Backups, selection mapping, reproducible sync script and verification:
`results/paper_ii_referee_response/diagnostics/finalization/2026-11-05/sync/`.
The local selective latexdiff generator and its suppression patch from the
1 November workflow were reused. Updated the Paper II workflow instructions
to preserve synchronization of these two variants in future manuscript edits.

## 2026-11-05 — remove the selective-diff blank page

The visible additions/deletions enlarged the A&A full-width title block beyond
textheight, producing an empty page before the Introduction. Reduced only the
latexdiff header-to-title skip from 49 pt to 24 pt and added a preamble hook
to start the Introduction on page 2. Body content and all diff markup are
byte-for-byte unchanged. Both compilers succeed; the too-tall warning is gone
and rendered page 2 contains the Introduction. Preserve these layout settings
when regenerating the selective diff. Backup, reusable preamble and checks:
`results/paper_ii_referee_response/diagnostics/finalization/2026-11-05/blank_page/`.

## 2026-11-05 — move interferometric comparison to Appendix A

Moved the former main-text Figure 3 (individual-array/combined continuum
comparison) and its detailed morphological description to Appendix A in clean,
marked and selective latexdiff sources. Appendix A is now titled
“12 m + 7 m interferometric combination”. Kept beam/noise specifications in
the main text with a concise appendix/figure cross-reference. Existing labels
provide automatic renumbering: comparison A.1, PSF A.2, all-scale image main
Figure 3. This author-directed relocation adds no new difference highlighting.
The accepted latexdiff text matches the clean source; all three built-in
compilations succeed. Pre-edit versions are saved in
`results/paper_ii_referee_response/diagnostics/finalization/2026-11-05/figure3_appendix/`.

Figure A.1 placement follow-up: use the float package and `[H]` for this
figure in all three versions, ensuring it follows the Appendix A heading.
All three compile successfully; clean Appendix A page visually inspected
and confirmed to show the heading immediately above the figure.

## 2026-11-05 — synchronize subsequent author edits

Updated marked and selective latexdiff sources from the latest clean text,
preserving the 90 existing referee-related selections, the title-page fix and
Appendix A figure placement. Both content projections match the clean text
and both variants compile successfully. Evidence and backups: `results/paper_ii_referee_response/diagnostics/finalization/2026-11-05/author_sync_195513`.
