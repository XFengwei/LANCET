# Referee-comment dependencies and proposed revision order

Assessment: 25 October 2026. This is a dependency review, not a record of
completed scientific validation. Comment IDs follow `task.toml`.

**Status update, 27 October 2026:** R5 and R3 are complete. The author
accepted the revised Introduction and the detailed Paper-I characterization
as the response to R3. The earlier references below to a later R3 boundary
perturbation describe proposed author-side work, not an outstanding referee
commitment or a prerequisite for moving on. R4a is the next discussion item,
followed by R4b and R4c; R10 can proceed alongside the imaging checks.

Follow-up evidence is recorded in
`results/paper_ii_referee_response/diagnostics/R5_R3/2026-10-25/report.md`
and `progress.md`. That audit establishes extraction-file relationships and
centroid footprint differences, but does not close R5/R3 or establish a
fully corrected mosaic response. The recommended sequence below remains
applicable; its original “no discrepancy established” statement describes
the starting assessment, before this follow-up audit.

## Recommended starting point

Start with **R5: primary-beam conventions and the common selection footprint**,
alongside the **definition/provenance part of R3: regional boundaries**.
These determine which sky, sources and background measurements the later
comparisons actually describe. R5 is an audit first: a discrepancy has not
yet been established, and a rerun is needed only if the adopted processing
requires correction.

The first deliverable should be a product/selection table and an overlay:

- Identify the 7m, 12m, 12m7m and 12m7mTP products, their units, beams,
  sky-brightness versus PB-attenuated status, effective gains and noise maps.
- Trace the combination, detection and photometry order; establish that PB
  correction is applied once and that the TP and interferometer inputs have
  compatible sky responses.
- Overlay the individual gains, adopted PB=0.2 selection, regional masks,
  cores and fragments. Record which objects and denominators are inside the
  common footprint. Keep regional boundaries distinct from sensitivity masks.
- Define whether the core/fragment ratio is a regional mass ratio or a
  parent-associated ratio. Do not silently exchange the two definitions.

This work supports R4, R6 and R2b directly, and identifies whether downstream
catalogues and maps need rebuilding. R3's boundary-perturbation calculation
comes later, once the selected masses can be recomputed consistently.

## Dependency table

“Needs before completion” describes dependencies of a final scientific
answer, not a prohibition on starting a diagnostic or writing definitions.
Conditional dependencies apply if the upstream work changes an adopted input.

| Item | Work that can begin now | Needs before completion | What it informs |
|---|---|---|---|
| R5: PB/selection | Trace products, gain, extraction and photometry; inspect noise and masks | Actual processing provenance; resolve any PB inconsistency | R4 comparisons, R6 recovery/footprint, R2b provenance; affected downstream maps/catalogues |
| R3: regions | Document coordinates, signposts, apertures and mask versions | For boundary sensitivity: R5 footprint, selected core/fragment/cloud masses and R7/R8 mass conventions | R6 selection; regional R1 and R2b results |
| R4a: overlap amplitudes | Inspect calibrated inputs and construct a common-response test | R5 response definition; frequency/unit/beam conventions | Confidence in adopted 12m+7m product; changes propagate to R4b/c, R6 and affected results |
| R4b: conventional joint imaging | Inspect existing CASA products, settings and dirty PSF availability | R5 common footprint and matched photometry; consider any R4a scale discrepancy | Validation or revision of the adopted interferometric image |
| R4c: J-comb/feather/TP effect | Inspect inputs and existing matched-source statistics | R5 conventions and a documented adopted interferometer input; repeat if R4a/b changes that input | R2b TP narrative; affected compact/background measurements in R1/R9 |
| R10: SIMBA benchmark | Locate bandpass data; check units; design matched comparison | Common frequency/beam/mask and documented dust+free–free reconstruction | Independent observational consistency check for TP-dependent claims in R2b/R4c |
| R6: completeness | Design and pilot flux/size/background/position recovery and parent matching | R5 selection, R3 boundaries, adopted extraction products; R8 sizes and R7/R8 assumptions for mass-space interpretation | Selected samples, mass-fraction sensitivity in R1; sample changes in R9 |
| R7: compact free–free | Check source-level radio/spectral constraints and meaningful limits | Adopted source positions/fluxes; propagate corrections or bounded scenarios consistently | Core masses, R1 regional sums and affected core-property relations in R9 |
| R8a: sizes/limits | Audit deconvolution, resolved flags and treatment of limits | Adopted catalogue/beam; R6 sample for final tables | Size/density distributions and correlations; size-dependent recovery interpretation |
| R8b–d: mass model | Audit opacity/gas-to-dust definitions; calculate brightness and temperature scenarios | Adopted fluxes/selection and R7 contamination handling for joint scenarios | R6 mass completeness, R1 mass fractions and affected R9 quantities |
| R9a: environmental estimators | Explain local versus aperture estimates and inspect implementation | Adopted maps/sample; propagate relevant R5–R8 changes before final fits | Final environmental interpretation and R2b assessment |
| R9b: velocity components | Inspect spectra/velocity interval and test component sensitivity | Its own component-aware or flagged moment product | Kinematic environmental relations; does not block continuum-only mass fractions |
| R1: scales/terminology/ratios | Agree scales and instantaneous-mass-fraction language | For final numbers: R3/R5/R6 selection plus R7/R8 masses; any adopted imaging correction | Abstract, mass-fraction figure, discussion, conclusions |
| R2a: motivation | Revise coverage/dynamic-range motivation | Verified observational description; no new analysis dependency | Introduction and abstract |
| R2b: role of TP | Trace which results use which image/background/denominator | R4c statistics and R10 benchmark; affected R1/R9 comparisons | Evidence-based claims about TP's added value |
| E1: final response/variants | Keep reply synchronized with completed units | All scientific answers and user-finalized clean manuscript | Submission deliverables; bold/diff remain frozen until then |

## Suggested sequence

1. **Establish definitions and inputs:** R5 audit plus R3 masks/provenance.
   Agree R1 terminology and R8 opacity conventions at this stage, without
   freezing numerical conclusions. R2a is an independent, low-risk text edit.
2. **Validate the images:** R4a and R4b share the response/footprint defined
   above. R4c uses the adopted interferometric input. Run R10 alongside these
   as a check of the reconstructed low-resolution input, rather than treating
   agreement between J-comb and feather as validation of that input.
   Existing-product inspection should precede any costly new CASA run.
3. **Assess selection and physical systematics:** R6 recovery/associations,
   R7 compact free–free and R8 size/mass/temperature checks. These are partly
   independent: flux-space recovery can run before temperatures are settled;
   mass-space recovery and summed-mass conclusions use the adopted scenarios.
   Begin these diagnostics earlier where their inputs are already stable.
4. **Recompute the conclusions:** finish R3 boundary sensitivity and R1
   unrounded mass accounting; update affected R9 measurements/fits. R9a
   definitions and R9b spectra can be examined independently much earlier.
   Finalize R2b using the measured compact-source and background differences.
5. **Synchronize manuscript and response:** verify all numbers and locations,
   then follow the existing finalization checklist. This review does not
   authorize generation of bold/diff versions.

The longest scientific chain is likely **PB/response → adopted image and
extraction → selection/recovery and physical assumptions → regional sums and
interpretation**. This is a dependency-based judgment, not an execution-time
estimate. We should not answer the report strictly in numerical order.

## Details verified in the current repository

- `source_catalogue/task.toml` reads the three `*_shrink.reg` masks;
  `column_density/task.toml` references the unshrunk `trichotomy_*.reg`
  masks. Both versions are present. Their different purposes and relationship
  to the comparison footprint need checking; their existence alone does not
  establish an error.
- `scripts/source_catalogue.py` selects and labels sources using those region
  masks before writing the adopted table. It reads an existing getsf
  extraction. Rebuilding this table is not itself a new extraction or a
  completeness test.
- `core_formation/task.toml` supplies fixed cloud and Paper-I fragment mass
  sums. `scripts/core_formation.py` summarizes the physical catalogue with
  these settings. Moving core boundaries alone will therefore not recompute
  all three mass reservoirs. Boundary/footprint sensitivity needs matching
  fragment selection and an auditable cloud-mass integration or a clearly
  bounded alternative.
- The same configuration already converts fragment masses from dust opacity
  0.9 to 0.5 cm²/g. The v1-based reply's numerical concerns must be reconciled
  with v2; this conversion must not be applied twice. The map configuration's
  opacity 0.005 and the core configuration's dust opacity 0.5 with gas-to-dust
  ratio 100 are compatible conventions if implemented as intended.
- The maintained chain is `source_catalogue → physical_properties →
  core_formation`, with `physical_properties + column_density + moment_maps
  → environment_measurements → correlation` as a second branch. Changes
  should be propagated along the affected branch, rather than rerunning all
  tasks indiscriminately.
- The intercore–intercore fit uses the source-removed 12m7m component; a change
  to TP alone does not automatically change it. Core-mass temperature/free–free
  corrections likewise need not change that specific fit, although they
  affect core-property relations. Check actual columns before choosing reruns.
- CASA joint-image and PB-corrected image paths declared in
  `prepare_continuum_images/task.toml` exist, as do the configured combined
  12m and ACA measurement-set paths. Their suitability, settings and PSFs have
  not yet been audited. A new joint-imaging run is not automatically necessary.
- `combination_diagnostic.py` currently exposes inspection, uv coverage,
  image comparison and beam-match stages. Existing coverage plots do not
  constitute the requested quantitative overlap or method-validation evidence.

Sources inspected: the full comment boxes in the response draft, the current
clean v2, and the relevant task configurations and runners. No scientific
pipeline was run and no manuscript/reply text was changed for this review.
Statistical covariance/spatial-dependence checks in R9a are author-side
robustness work supporting the environmental claim; the referee explicitly
asks for estimator clarity and velocity-component assessment.
