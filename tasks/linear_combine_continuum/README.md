# linear_combine_continuum

Linear uv-plane combination of the 12m+7m continuum mosaic with the SED-based total-power (TP)
models, by J-comb and by feather. It belongs to Paper II referee item R4c and is meant to replace
the unshared `allscale_combined_v3_*` images once a version is promoted.

**Status (2026-11-17).** The TP is combined as delivered (TP factor 1.0, no change of flux
scale). Every combined image is checked for conservation of the zero-spacing information: it is
smoothed to the TP beam and compared with the TP. No absolute calibration against Planck is
applied.

## Run

```bash
PYTHONPATH=src conda run -n py312 python scripts/linear_combine_continuum.py --stage all --overwrite
```

Stages, in order: `prepare`, `tp_factor`, `jcomb`, `feather`, `box_spectrum`, `planck_calibration`; `all` runs the six. The whole task
takes about 5 minutes, 7 GB of memory and 1.4 GB of disk. Existing outputs are replaced only with
`--overwrite`. All settings are in `task.toml`.

## Inputs

| Name | File | Notes |
|---|---|---|
| `interferometer` | `data/delivered/alma/7MTM2TM1/cont/miriad/IRDC316.8.7MTM2TM1.robust0.5.clean.fits` | MIRIAD 12m+7m joint mosaic, PB-corrected, beam 0.71″ × 0.61″, 0.1″ pixels |
| `interferometer_attenuated` | `…clean.pbant.fits` | the same image times the mosaic gain; gives the gain |
| `tp_dust` | `results/synthesize_continuum/model/dust_224GHz.fits` | Paper I SED dust model at 223.8 GHz (= `data/delivered/sedfitter/Dustmodel_223.8GHz.fits`), 18.1″ beam |
| `tp_dust_ff` | `results/synthesize_continuum/model/continuum_224GHz.fits` | dust plus free-free at 223.8 GHz |
| `core_catalogue` | `results/source_catalogue/tables/source_catalogue.ecsv` | core positions for the core-peak measurement |

The TP is an SED model, not a single-dish observation.

## Outputs (`results/linear_combine_continuum/`)

```
prepared/     inputs on the 12m+7m grid, one copy each
  interferometer_mjysr.fits      12m+7m in MJy/sr, zero where the mosaic gain is undefined (< 0.2)
  tp_dust_mjysr.fits, tp_dust_ff_mjysr.fits        TP models in MJy/sr (18.1″ beam)
  tp_dust_footprint.fits, tp_dust_ff_footprint.fits  1 where the TP model exists
  gain.fits                      mosaic gain; also the coverage map of the products
  input_check.json               beams, TP margin around the mosaic, regridding check
jcomb/        J-comb products
  jcomb_<tp>.fits                combined image, Jy per 12m+7m beam, float32
  jcomb_<tp>_at_tp_beam.fits     the combined image smoothed to the TP beam, on the TP's 6″ grid (MJy/sr)
  zero_spacing_ratio_<tp>.fits   that image divided by the TP
  zero_spacing_<tp>.{pdf,png}    TP, smoothed combined image, their ratio, pixel-by-pixel comparison
  comparison_<tp>.{pdf,png}      before (12m+7m, TP) and after (combined)
  power_spectrum_<tp>.{pdf,png}  power spectra of inputs and product; product amplitude relative to each
  measurements.json
feather/      the same for feather (feather_<tp>.fits, …)
tp_factor/    how a factor of 1.4 arises from matching power spectra (documented, not applied)
box_spectrum/ power spectra of the TP and the 12m+7m image inside boxes
planck_calibration/   TP model against Planck (absolute flux scale; nothing is applied)
manifest.json one record per stage run: settings, input checksums, outputs
```

`<tp>` is `dust` or `dust_ff`. Product headers carry the settings as `COMB_*` keywords
(method, thresholds, TP factor, padding, `COMB_ABSOLUTE_SCALE = pending`).

## Method

Both images are brought to the 12m+7m grid in surface brightness (MJy/sr) and Fourier
transformed. With T the TP transform, I the 12m+7m transform, G_low and G_high the two beams
and f the TP factor:

- **J-comb:** C = W · f · (G_high / G_low) · T + (1 − W) · I. W is 1 below k1 = 5 kλ, 0 above
  k2 = 9 kλ and a Butterworth taper (order 5) between. Dividing by G_low removes the TP beam, so
  the TP is used only where G_low ≥ 0.1. The result has the 12m+7m beam at every scale.
- **Feather:** C = f · T + (1 − G_low) · I. No thresholds; the result is not exactly at the
  12m+7m beam.

The sum is transformed back and converted to Jy per 12m+7m beam. Details:

- The TP is padded by 120″ with its mirror image, rolled off with a cosine; the 12m+7m image is
  padded with zeros.
- The 12m+7m image is not used where the mosaic gain is undefined (< 0.2).
- Products are kept wherever the TP model exists (`fill_outside_mosaic = true`). Outside the
  mosaic they are TP-only: f × TP for feather, f × TP sharpened to k2 for J-comb.
- The Butterworth taper is the well-behaved form of the J-comb paper's Eq. 19 (its printed form
  diverges at the mid-point), rescaled to run from 1 at k1 to 0 at k2.

## TP factor

`[tp_factor] value = 1.0` multiplies the TP in both methods: the TP is used as delivered. One
factor serves both TP models and is determined from the complete model, dust + free-free
(`determined_from`); the dust-only model has the free-free component removed, serves only for
dust column density, and takes the same factor.

A factor of 1.4 (`psd_matched`) had been proposed from matching power spectra. It is documented
by the `tp_factor` stage and not applied. The power-spectrum match works as follows: the 12m+7m image is smoothed to the TP beam, both spectra are taken over a common
region, and the TP is raised until its lowest point within 20–50″ touches the 12m+7m spectrum.
The `tp_factor` stage repeats this in two ways (`tp_factor.json`, dust + free-free TP):

| Variant | Region | Factor |
|---|---|---|
| `as_determined` | all non-zero 12m+7m pixels, including the rim where the gain is undefined; blank pixels count as zero in the smoothing; TP regridded bilinearly | **1.39** |
| `prepared` | the inputs actually combined: gain defined, edge-normalised smoothing, spline-regridded TP | **1.07** |

The first reproduces the procedure that gave 1.4. The difference between the two is the noisy
mosaic rim (gain < 0.2), which adds 12m+7m power (about × 1.19), and the bilinear regridding,
which lowers the TP power at these scales (about × 1.07).

A mock test with a known factor of 1.2 on the real mosaic geometry
(`results/paper_ii_referee_response/diagnostics/R4c/2026-11-17/estimator_mock/`) shows that this
power-spectrum estimate at the TP resolution is itself biased low (it returns 1.00 where the gain
is defined and 0.87 for gain ≥ 0.5), while the cross-spectrum at the 12m+7m resolution returns
1.20. On the real data the cross-spectrum gives 1.18. If the 12m+7m image is incomplete across
5–9 kλ, every estimator returns too low a value. The estimates of the factor therefore range from
about 1.0 to 1.4 with the method and the region; Planck in the brightest beam gives 0.97–1.02
(see "Absolute flux scale"). The TP is kept at 1.0.

### Boxes (`box_spectrum/`)

The same comparison inside boxes listed in `[box_spectrum] boxes` (`box_<name>.json`, `.png`).
The window is the part of the box where the 12m+7m gain is defined. TP factor that aligns the
dust + free-free TP with the 12m+7m image over 20–50″, and in brackets what a mock of true
factor 1.2 returns (mock 12m+7m without scales below 3 kλ):

| Estimate | middle, 120″ × 120″ at 14:44:57.89 −59:48:11.69 | south, 180″ × 90″ at 14:45:23.04 −59:49:22.49 |
|---|---|---|
| Power spectra at the TP resolution, lowest TP point | 0.99 (mock 1.03) | 1.20 (mock 1.09) |
| Amplitudes at the 12m+7m resolution, median | 1.08 (mock 1.24) | 1.21 (mock 1.03) |
| Cross-spectrum at the 12m+7m resolution | 0.91 (mock 1.16) | 1.08 (mock 1.01) |
| Corrected by the mock (× 1.2 / mock) | 0.94–1.15 | 1.29–1.41 |

The middle box aligns with a factor of about 1.0; the south box asks for about 1.1–1.4. The uv
bins are 1.7 kλ (middle) and 2.3 kλ (south) wide, so the matching range holds three bins, and the
mock correction depends on the assumed completeness of the 12m+7m image. The south box is 38″
from the edge of the image grid.

## Zero-spacing check

For every combined image: convert to surface brightness, smooth from the 12m+7m beam to the
18.1″ TP beam, sample on the TP's own 6″ grid, and compare with the TP it was combined with.
Only positions with a full TP beam inside the product are used. `measurements.json` gives, for
the mosaic (gain ≥ 0.5), the TP-only area and everything together: the ratio of summed brightness
(`flux_ratio`, the total-flux test), percentiles of the pixel ratio, and the rms difference.

Results (2026-11-17), flux ratio and pixel ratio median (16th–84th percentile):

| Image | Mosaic: flux ratio | Mosaic: pixel ratio | TP-only area: flux ratio | TP-only area: pixel ratio |
|---|---|---|---|---|
| J-comb dust+ff | 1.0001 | 1.000 (0.973–1.029) | 0.9997 | 1.000 (0.980–1.020) |
| J-comb dust | 1.0000 | 0.999 (0.967–1.034) | 0.9996 | 1.000 (0.979–1.020) |
| Feather dust+ff | 0.987 | 0.990 (0.953–1.038) | 1.006 | 1.009 (0.973–1.046) |
| Feather dust | 0.990 | 0.989 (0.945–1.042) | 1.005 | 1.007 (0.973–1.044) |

- J-comb keeps the TP's flux to 0.01% inside the mosaic and 0.04% outside it. The pixel scatter
  of 2–3% comes from the 12m+7m content between 5 and 9 kλ, which the 18.1″ beam still passes, and
  from pixel-scale structure in the TP model.
- Feather is 1.0–1.3% low inside the mosaic: it down-weights the TP by its beam at scales where
  the 12m+7m image does not have the full signal. Outside the mosaic the feather product is the TP
  itself, already at 18.1″; smoothing it again blurs it, which gives the 4% pixel scatter there
  while the flux stays within 0.6%. For feather the check is meaningful inside the mosaic.

## Other measurements (`measurements.json`)

All inside the mosaic where the gain is at least 0.5.

| Measurement | Meaning | J-comb dust+ff; dust | Feather dust+ff; dust |
|---|---|---|---|
| `high_k_join` | product / 12m+7m amplitude from 9 to 20 kλ | 1.001; 1.001 | 0.973; 0.973 (flagged) |
| `core_contrast` | background-subtracted core peaks, product / 12m+7m (211 cores) | 0.999; 0.999 | 0.999; 0.999 |
| `negative_fraction` | pixels below −3σ (12m+7m alone: 1.3%) | 0%; 0.56% | 0%; 0.21% |
| `mock_sky_test` (J-comb) | J-comb on a mock sky with a known answer | passes (flux 0.03%, rms 0.2% of peak) | — |
| `beam_volume_outside_high_beam` (feather) | part of feather's effective beam outside the 12m+7m beam | — | 0.05% |

`checks` flags a ratio that misses 1 by more than its tolerance in `[measure]`. Feather's
high-k join is 2.7% low because it still down-weights the 12m+7m image by 1 − G_low just above
9 kλ. The dust-only images keep negative patches in the south, around the H II region: the
12m+7m image contains free-free emission that the dust-only TP lacks.

## Absolute flux scale (pending)

No absolute calibration is applied; `[absolute_scale] status = "pending"`. The products are on
the scale of the 12m+7m calibration at small scales and of f × TP at large scales. The
`planck_calibration` stage provides the comparison the decision rests on. It follows the check of
2026-11-15 (`results/paper_ii_referee_response/diagnostics/R4c/2026-11-15/report.md`), reproduces
its numbers, and adds a measurement in the brightest Planck beam alone.

Method: Planck 217 GHz (K_CMB) minus CO(2-1) from the Type-1 CO map, converted to brightness at
224 GHz through the band transmission for a modified black body (β = 1.8, T = 20 K;
479.4 MJy/sr per K_CMB). The dust + free-free TP model is regridded flux-conservingly to 30″
pixels and smoothed to the 5.0′ Planck beam; the fraction of each beam inside the TP model
("coverage") is carried along. Planck's beam sees only the TP part of a combined image, so this
constrains the TP scale, not the 12m+7m scale.

Linear fit TP = gain × Planck + offset (`planck_calibration.json`):

| Selection | gain | offset (MJy/sr) | independent beams |
|---|---|---|---|
| coverage ≥ 0.9 | 1.15 | −4.7 | 1.6 |
| coverage ≥ 0.8 | 1.14 | −4.8 | 2.8 |
| coverage ≥ 0.7 | 1.18 | −5.7 | 3.8 |
| outside the model filled with scaled 870 µm | 1.07 | −2.8 | 5.6 |
| no CO removal, coverage ≥ 0.9 | 1.29 | −10.5 | 1.6 |

Brightest beam: a disc of 2.5′ radius around the peak of the smoothed TP (RA 221.308°,
Dec −59.815°, the southern and middle part of the filament; coverage 0.95 at the peak, at least
0.90 in the disc). TP / Planck there, and the TP factor that would match Planck:

| Treatment of the background | TP / Planck (disc mean; peak pixel) | Implied TP factor |
|---|---|---|
| none (direct ratio) | 0.97; 1.02 | 1.03; 0.98 |
| fitted uniform background (4.1 MJy/sr) removed from Planck | 1.17; 1.22 | 0.85; 0.82 |
| differential against the faintest 20% of the well-covered area | 1.29; 1.36 | 0.77; 0.74 |
| none, and no CO removal | 0.89; 0.93 | 1.12; 1.07 |

In the brightest beam the TP as delivered equals Planck within 3% when no background is removed.
The answer depends on the background: the more of Planck's signal is assigned to a background the
TP model lacks, the brighter the TP is relative to Planck. The brightest beam is one measurement
(the disc holds under one independent beam), and up to 5% of its beam weight falls outside the TP
model, which makes the TP read slightly low.

## Limitations

- The TP is used as delivered. Overlap estimates of a TP factor range from 1.0 (middle of the
  filament) to 1.1–1.4 (south); records of the earlier tests are in
  `results/paper_ii_referee_response/diagnostics/R4c/2026-11-17/`.
- The dust-only images are not free of free-free emission at small scales (see above).
- The 12m+7m image is assumed complete across 5–9 kλ; its shortest spacing is 6.6 kλ.
- Beyond the edge of the 12m+7m grid the TP is mirrored, not read from the wider TP map. The
  mosaic is at least 35″ from that edge.
- CASA 6.7.2 `feather` was tested and not used: on these images it returns the sum of the two
  inputs. The feather here is the task's own implementation.

## Promotion

Nothing is promoted. After the TP factor and the absolute scale are decided: copy the chosen
J-comb images to `data/delivered/`, repoint the downstream tasks, and state the Planck check and
the thresholds in Paper II.
