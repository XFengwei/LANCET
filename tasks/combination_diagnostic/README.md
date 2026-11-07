# Continuum combination diagnostics

This task owns continuum uv coverage, visibility amplitude/weight diagnostics,
and image-combination comparisons. It extends the former
`compare_interferometer_arrays` task; its canonical name, runner and result root
are now `combination_diagnostic`. Historical results were moved, not discarded.
An old-name result symlink preserves existing historical links.

## Read-only inputs

The four self-calibrated measurement sets, their original combination script,
and the TP model FITS image are declared in `task.toml`. The runner never
executes the historical shell script or modifies an MS. Its CASA tables are
opened with `nomodify=True`; automatic CASA data updates and file logs are
disabled. CASA can emit a sandbox gRPC-listener warning even when table reads
succeed; this is not evidence of a failed data read.

The configured prepared MJy/sr maps remain the inputs for the existing image
comparison and beam-matching stages. Input paths and output policies are checked
by the runner. All new products belong under `results/combination_diagnostic/`.

## Run

Use the project's `py312` environment, which includes `casatools`:

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage inspect
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage uv --refresh-uv --overwrite
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage visibilities --overwrite
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage total-power --overwrite
```

- `inspect`: report declared inputs and prepared-image metadata without writing.
- `uv`: reread the four MSs on request; retain cross-correlations not row-flagged
  and with at least one unflagged correlation/channel. Deduplicate positions to
  0.1 m, then plot every retained position and its conjugate. This is sampling
  geometry, not a density or imaging-weight plot. Full coverage and short
  baselines use the same array colours and equal square plot areas. The small
  square in the full view marks the zoom. No inferred TP disk or purported
  consistency region is drawn. Later runs can reuse the flag-aware cache;
  use `--refresh-uv` after input/flag changes.
- `visibilities`: read every DDID in bounded chunks. Exclude autocorrelations,
  `FLAG_ROW`, per-sample flags, nonfinite samples, and nonpositive/nonfinite
  weights. Use parallel-hand correlations only. Prefer `WEIGHT_SPECTRUM` when
  defined; otherwise broadcast `WEIGHT` over channels. These self-calibrated
  inputs contain `DATA` and no `CORRECTED_DATA`; `DATA` is therefore explicit.
  Vector-average channels and parallel hands within each row. Do not average
  across fields, times, baselines or DDIDs. Convert projected baseline to kλ
  using the row's weight-averaged frequency. Plot the weight-averaged scalar
  amplitude, mean row weight, summed weight and row count in 80 common
  logarithmic bins. Export full-range and short-baseline figures, global,
  per-field and per-DDID CSVs, additive NPZ statistics and full MS metadata.
  Empty-bin means are NaN; all usable rows must fall within the configured
  bin range. Saved amplitude standard deviations describe within-bin spread,
  not a calibration uncertainty. Weights are reported in native MS units,
  before Briggs/taper weighting or the synthetic-ACA reweighting in MIRIAD.
- `total-power`: convert the beam-convolved Jy/beam model to Jy/pixel and
  Fourier-transform it after a frame Hann window and a cosine-squared taper
  of the finite-data mask over two beam FWHM. Show complete Fourier annuli
  and the Gaussian beam transfer separately. Preserve NaN counts, taper,
  zero-spacing value and beam metadata in JSON. This is a finite/windowed
  model-image transform, not measured TP visibilities. The FITS file has no
  visibility-weight distribution. The header has no frequency axis; 225 GHz
  is explicitly the script-based reference frequency, not a FITS measurement.
- `compare-images`: retained vertical ACA/12 m/12 m+7 m comparison.
- `beam-match`: retained convolution/reprojection to the ACA beam, FITS
  residuals, summary statistics and diagnostic PNG.
- `all`: run every stage; use only when all image products are wanted too.

The combined 12 m MS is a bookkeeping control. It is not added to C2 and C5
as an independent observation. The control JSON compares its total weight,
row count and global scalar amplitude with the separate C2+C5 sets after
averaging; agreement checks the averaging implementation, not flux calibration.

## Interpretation and next R4a step

The current visibility profiles are initial diagnostics. ACA has 118 pointings,
12 m has 150; their primary beams, channel bandwidths and noise differ. Scalar
amplitudes have a positive noise bias and retain the individual mosaic PBs.
No spectral-index scaling or common-sky-response correction is applied. A raw
ACA/12 m amplitude ratio is therefore not an inter-array calibration factor.
The next scientific step is a common-response comparison over the actual
shared spatial frequencies, with a defensible uncertainty.

The supplied script records ACA-model-to-12 m synthetic visibility generation
and artificial reweighting. Its `robust=1.0` differs from the `robust0.5` names
of the adopted images, and its uv selection is commented. Treat it as supplied
provenance, not proof of the exact settings of every adopted image. No imaging
rerun, synthetic weight reconstruction or script correction is implied here.

## Controlled R4a overlap tests

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage overlap-tests --overwrite
```

This explicit stage writes `results/combination_diagnostic/overlap_tests/`.
It is separate from `all`, because it is a focused scientific test rather
than a routine figure refresh. It caches channel-averaged XX and YY rows
separately for the C5, C2 and ACA sets out to 40 kλ. The combined 12 m set is
not counted twice. Remove the three `*_rows.npz` and `*_metadata.json` caches
if the underlying data or flags change. No MS, weight or flag is altered.

The three tests are:

1. **Noise.** Isolate broad windows with total bandwidth at least 1 GHz;
   ACA then retains DDIDs 0/7. Estimate effective noise from non-overlapping
   adjacent integrations of the same observation, scan, field, DDID and
   baseline, with a maximum 30-second gap. Compare the same rows before/after
   coherent two-integration averaging and bound fringe smearing. Estimate
   real/imaginary variance via robust MAD of normalized time differences;
   compare XX-YY and XX/YY cross-power as independent-hand controls. Keep
   signed noise-subtracted power before taking a square root; do not clip
   individual negative powers. Signal RMS and mean scalar amplitude are
   different statistics and are labelled separately.
2. **Sky and sampling.** Nearest-first one-to-one field matching within 5″,
   with a 2″ sensitivity subset. Match broad sideband, radial bin and one of
   eight orientation sectors modulo 180 degrees. Retain cells with at least
   four rows from each array; use natural weights within cells and equal
   weights between common cells. This changes the estimand deliberately:
   both arrays now represent the same distribution of common sky/uv cells.
   The compared 12 m sample pools C2/C5. Shading is the 16th–84th percentile
   of 500 resamples of matched pointing pairs; it is a diagnostic scatter
   interval, not an independent-field calibration confidence interval.
3. **PB response.** Forward-transform two existing restored sky models
   (all-scale and 12 m only) through Gaussian PBs from the supplied script,
   26″/43″ at 225 GHz scaled inversely with frequency. Use each array's actual
   pointing and uv coordinates. Repeat with a common 26″ beam centred on
   each matched pair midpoint. Compare native/common model RMS ratios on
   exactly the same cells as the data. The conditional PB adjustment is the
   observed noise-corrected RMS ratio times model-common/model-native.
   It is an adjustment of an ensemble statistic, not division of individual
   extended-source visibilities by a scalar PB value. Models are binned by
   pixel area, sampled on 0.5″ tangent grids over 256″, FFT padded by two,
   and interpolated cubically. The w term is bounded at PB=0.01. A 0.25″
   sampling check uses every row of the nearest matched pointing pair.

The 12 m-only sky model provides a check independent of ACA imaging, but
both models are data-derived restored images. They are not an external
flux standard or an exact joint-deconvolution model. Gaussian PB shape,
residual sky-model errors, calibration and sample dependence remain relevant.
No scalar flux rescaling is adopted from these tests. Results and limitations
are recorded in the owning result report; R4a remains in progress.

## Simple ACA before/after PB curve

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage pb-curve --overwrite
```

This lightweight stage reuses the completed overlap-test row caches and saved
model predictions; it neither rereads the MSs nor reruns the sky transforms.
It checks that settings and selected cell keys match the completed run.
The full `overlap-tests` stage also generates this curve on subsequent runs.

`overlap_tests/aca_pb_amplitude_curve.{pdf,svg,png,csv}` shows native-PB ACA,
model-translated common-PB ACA, and common-PB 12 m as a reference. These are
**noise-corrected RMS amplitudes of the same matched sky/uv cells**, not the
original all-field mean scalar amplitudes. For each array separately:

    A_common = A_measured × A_model,common / A_model,native

Both arrays are expressed through the 26″ PB centred on the matched-pair
midpoint. The small 12 m adjustment includes this centring change. The
all-scale sky model supplies the main curves; the shaded ACA band spans
adjustments from the all-scale and 12 m-only sky models. It is a model range,
not a statistical confidence interval. No normalization is fitted to make
the two measured curves agree, and no observed visibilities are changed.

## Comprehensive Figure 1 candidate

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage full-sampling --overwrite
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage figure1 --overwrite
```

`full-sampling` reads the full 12 m baseline range for the same 19 matched
pointing pairs and reuses ACA's full-range row cache. It no longer restricts
both arrays to the common 8–32 kλ interval. Its source CSV, cell statistics
and metadata are in `results/combination_diagnostic/full_sampling/`. Clear
its row caches after input/flag changes. Rerun this stage after scientific
settings change; a figure-only redraw reuses its verified products.

Signal RMS uses independent-hand cross-power throughout, avoiding an
extrapolation of short-baseline noise subtraction to very faint long modes.
The auto-power-minus-noise estimate remains in the CSV as a check. Powers
are signed until averaged; non-detections use conservative limits based on
two propagated within-cell standard errors. These limits do not include
full calibration covariance or sky-model uncertainty. Cells receive natural
weights internally and equal weights between occupied cells. Within overlap,
both arrays use the same occupied cells; outside overlap each retains its
own cells, preserving complementary sampling.

ACA PB translation is model-based at every observed ACA radius; no factor
is extrapolated from the overlap. The target PB is now centred on the actual
matched 12 m pointing, so measured 12 m data provide the reference directly.
The RMS ratios in overlap remain near unity after PB matching, with residual
model/selection dependence retained. There is no fitted flux renormalization.

`figure1` assembles `figures/Figure1_combination_diagnostics.{pdf,svg,png}`:
(a) pointings without matched-pair boxes, (b) full amplitude ranges and observation
sampling, (c) full uv coverage, (d) short baselines plus TP zero spacing.
Only the panel letters and legends are internal. Panel b uses a logarithmic
uv-distance scale from 1 to 1000 kλ. It shows the ACA curve translated to the
12 m PB, separate C2/C5 references and the TP model, without the native-ACA curve
or shaded intervals. Below a divider, eleven observation tracks show sampling
in 45 logarithmic uv bins; dot area is proportional to the unweighted valid-row
count on a common scale. Eight tracks are ACA, one C2 and two C5. These use
all broad-window rows in the same 19 paired pointings before common-cell
selection; track height is an observation index, not amplitude. Counts are
exported in `figures/Figure1_observation_sampling.csv`, with all rows accounted
for and no thinning. Raw scalar amplitudes are not overplotted on the corrected
RMS axis. The strip heading is "Visibility counts"; its area/count key is
explained in the caption. `figure1` derives `full_sampling/separate_array_amplitudes.csv`
from the existing row caches. In each radial bin, it uses common cells across
all arrays present where possible, and each array's own cells otherwise,
requiring at least four rows per array/cell. Original pooled profiles remain
unchanged. C5/C2 has a median RMS ratio of 1.008 across 14 matched bins.

The PDF is exactly 184 mm wide, matching the local A&A class text width.
Text uses a 10 pt base and 11 pt panel labels; ordinary log-axis exponents
are 7 pt. Export preserves the canvas instead of applying a tight crop,
so fitting to text width preserves these sizes.
The source remains in this maintained task; the caption, contract and
manifest are saved beside the figure. No manuscript is renumbered by this
stage and no paper figure is automatically replaced.

The TP reference is evaluated after **multiplication by the same actual
12 m pointing PB**, using the pointing/sideband/orientation population of
ACA's shortest occupied radial bin. The 18.1″ TP restoring beam is retained,
and plotting stops at its 10% response below native pixel Nyquist. Equal
PB treatment establishes a common sky aperture; it does not remove native
TP smoothing or justify forcing the TP curve onto an interferometer curve.
TP source data are `Figure1_tp_reference.csv`. Panel d shows zero spacing
and 50%/10% beam-response radii (6.70/12.21 m at 225 GHz), not TP uv samples.

## TP restoring-beam response test

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage tp-beam-test --overwrite
```

This explicitly requested diagnostic uses the ACA row cache and both existing
sky models. It convolves each sky model to the 18.1″ TP beam using a normalized
Gaussian kernel after subtracting the native restoring-beam covariance. It
then applies the actual matched 12 m PB before Fourier transformation, rather
than multiplying an already PB-windowed visibility by a beam transfer.
Model target/native power ratios translate the measured ACA RMS conditionally.
TP is sampled at the exact same ACA uv coordinates, weights and occupied cells
in every bin. All inputs remain read-only and no relative flux scale is fitted.

Outputs and limitations are in `results/combination_diagnostic/tp_beam_test/`.
The all-scale and 12 m-only models give substantially different corrections:
beam-matched ACA/TP spans 2.11–2.58 and 0.63–0.81, respectively. Thus this test
does not establish a unique smooth join or justify a common vertical shift.
The all-scale model contains TP information; neither model is an independent
flux standard. C2/C5 have no cells with four rows below the TP 10% beam-response
limit, so a similarly controlled direct TP/12 m comparison is unsupported.
Figure 1 retains its physical amplitude scale pending scientific interpretation.

## Per-pointing amplitude profiles

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage pointing-profiles --overwrite
```

This separate stage produces a 19-panel overview and a four-page, 184 mm-wide
atlas with 10 pt text under `results/combination_diagnostic/pointing_profiles/`.
Each matched pointing pair has its own ACA/C2/C5/TP curves, with consistent
colours, axes and native field IDs. It reuses the row caches and existing
cross-power/cell estimator, retaining full baseline ranges and common cells
where possible. The original radial grid is split at 9.165 kλ for the TP
comparison. The main Figure 1 is unchanged.

For every pointing, TP is evaluated at exactly the measured ACA uv samples
with identical natural row weights and selected cells, up to the TP beam's
10% response limit. Its 18.1″ restoring beam remains; ACA is translated to
the 12 m PB using the all-scale sky model. No restoring-beam matching or
flux rescaling is applied here. Source CSVs, pointings/coordinates, NPZ cell
statistics, input hashes, selection counts and QA accompany the figures.
All 19 pointings and 76 TP/ACA comparison bins are verified. See the result
README for exclusions, interpretation limits and reproduction details.

The `pointing-range` stage adds a single-axis view with all 19 curves per
array overlaid. It reuses the per-pointing source CSV and retains upper
limits. Continuous TP curves extend to the 1 kλ left edge, evaluated directly
from the TP model through each pointing's actual matched 12 m PB with a fixed
first-bin sideband/orientation reference population. Grey diamonds retain
the exact-ACA-sampled comparisons. Zero-spacing values are in the source CSV;
zero cannot appear on a log axis. Export is 184 × 135 mm with 10 pt text.
The plot shows the actual pointing-to-pointing range without a fitted flux
factor or confidence shading. Existing atlases, Figure 1 and source profiles
are unchanged. Details and QA are in the pointing-profile result README.

## Adopted first manuscript figure

On 28 October 2026 the author clarified that the **four-panel** figure is the
manuscript's first figure: (a) pointings, (b) amplitudes and visibility counts,
(c) full uv coverage and (d) short baselines plus TP. The single-axis 19-pointing
overlay remains a separate diagnostic. The maintained
`pointing_range.tp_display_scale = 2.7` applies a common display normalization
to TP curves/diamonds only, disclosed in the manuscript caption. It is rounded
from the median of 74 detected matched ACA/TP comparisons and is not adopted
as a physical calibration correction.

`figure1.tp_display_scale = 2.7` applies this normalization in panel (b) of
`Figure1_combination_diagnostics.*`; its source CSV stays unscaled. This
184 mm-wide four-panel figure has tighter gutters and margins and the short
axis label "Amplitude (mJy)", with 10 pt text. Its original paper copy was
`papers/II/submission_a&a_v2/figures/Continuum_combination_diagnostics.pdf`.
`pointing-range` separately writes `pointing_profiles_overlay_display_scaled.*`,
with a 184 × 110 mm canvas and 10 pt text.
The original unscaled overlay and physical source tables are retained. The
RMS estimator remains defined in the caption. Verified PDF copy provenance
and successful manuscript/response compilation are recorded in the R4a
evidence directory. This figure placement does not close the unresolved
R4a scientific assessment or update any observed flux scale.

## Figure checks and referee links

The `figure1-three-panel` stage produces the adopted horizontal Figure 1,
`figures/Figure1_combination_diagnostics_three_panel.{pdf,svg,png}`. It keeps
(a) pointings, (b) full uv coverage and (c) amplitudes/counts in one row. The
former short-baseline panel becomes an unlabelled upper-right inset in (b),
covering −25 to +25 m in both axes and highlighting ACA/TP. Main uv axes
use km; inset coordinates use m. The combined legend and caption describe
both views. All eleven count tracks and the disclosed TP display factor
are retained. The figure is 184 × 60.454 mm with 10 pt text.
The inset has a white background, no axis ticks or labels, and two black
zoom connectors.

The inset is deliberately excluded from the equal-main-panel alignment
group, with a recorded reason. The combined caption, contract, manifest, source exports and QA accompany
the figure. The verified three-panel PDF replaces the manuscript Figure 1
copy on explicit user adoption (28 October 2026); provenance and compiler
success are recorded in the R4a evidence directory. The original four-panel
source bundle remains unchanged. Rendering this stage does not itself copy
outputs into the paper folder.

For render-time panel QA, put the Nature Figure skill's `scripts` directory on
`PYTHONPATH` alongside `src`, and set `LANCET_FIGURE_QA=1`. This activates
`require_matplotlib_panel_alignment()` before export. Run its source validator,
PDF-text and rendered-collision audits after the final render. Diagnostics are
exported as PDF/SVG/PNG. The run report and QA notes belong to this result root.

`tasks/paper_ii_referee_response/combination_diagnostic` links to this task.
`results/paper_ii_referee_response/combination_diagnostic` links to its outputs.
R4a remains in progress until the common-response scientific assessment is
complete. Figure 1 uses the explicitly adopted three-panel layout; no manuscript
figure is replaced merely by generating these diagnostics.

## R4b: matched CASA/MIRIAD image cutouts

Run `--stage joint-images --overwrite` with the existing Python figure-QA
environment. The `joint_image_comparison` configuration selects two 100-arcsec
central/south cutouts. CASA PB-corrected MJy/sr is convolved to the MIRIAD beam
and reprojected; both columns share the valid PB>=0.2 footprint and each row
shares an asinh scale. Products and input hashes are under
`results/combination_diagnostic/joint_image_comparison/`. The figure is
included only in the referee response. Dirty-PSF work remains separate.

## MIRIAD joint-imaging dirty beam (R4b)

```bash
PYTHONPATH=src conda run -n py312 python scripts/combination_diagnostic.py --stage dirty-beam --overwrite
```

The MIRIAD products come from a rerun of Shuting Lin's continuum workflow
(`data/reduction/alma/7MTM2TM1/cont/script.zip`) in
`data/reduction/alma/7MTM2TM1/cont/miriad_r4b/work_shuting/`, stopped after the
final `invert`. It runs in the `miriad80:latest` Docker image (MAXANT=80), built from
`~/Software/miriad-alma-maxant80-u2204.tar.gz`. The combined TM1TM2 MS has 71
antenna rows, more than the 64 the older CARMA build accepts. The 12m UVFITS
keep their native 30 channels per spw (`export_uvfits_fullchan.py`). Averaging
them to one channel removes the per-cell density that robust weighting acts on,
and broadens the beam. `mospsf` gives the effective mosaic PSF at the beam-cube
reference position.

At robust 0.5 the rerun's `restor` beam is 0.707″ × 0.600″ at PA 47.2°.
The paper image has 0.712″ × 0.606″ at PA 47.4°. Robust 1 and 2 come out about
4–5% narrower than Shuting's delivered versions. The figure compares the PSF
with the paper restoring beam and shows the cumulative dirty/clean beam
volume. Outputs are in `results/combination_diagnostic/dirty_beam/`.
