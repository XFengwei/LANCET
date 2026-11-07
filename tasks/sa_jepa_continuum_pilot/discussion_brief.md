# SA-JEPA × LANCET: discussion brief

Prepared for discussion with Guang-Xing Li. This note separates what the
method establishes, what the current LANCET pilot actually shows, and what
still needs controlled testing.

## The core idea

ScaleAware-JEPA is a label-free representation learner for continuous scalar
fields. It does not reconstruct missing pixels. Instead, an online encoder
receives a masked multiscale field and predicts the latent representation made
by an unmasked exponential-moving-average target encoder. The loss is evaluated
only at hidden 3-by-3 target patches. A weak channel-spread hinge discourages a
collapsed, nearly constant latent representation.

The distinctive step is that both the encoder input and the masking task use
the same physical scale hierarchy:

1. Constrained Diffusion Decomposition (CDD) separates the input into
   pixel-registered scale components. Unlike a wavelet decomposition, the
   nonlinear diffusion constraint suppresses ringing around sharp structures
   and gives non-negative component maps for a non-negative input.
2. Each CDD component is paired with a mask indicator and a scale code, passed
   through a shared per-scale adapter, and fused coarse-to-fine before a dense
   ConvNeXt encoder. The output is a 32-channel latent vector at every pixel.
3. The context footprint at CDD scale `sigma` is approximately
   `sigma * mask_multiplier + offset`, odd-rounded and hard-capped. Thus the
   model is asked to predict fine, intermediate, and coarse structure using
   context matched to each scale, rather than one arbitrary image patch.
4. At inference, the frozen EMA target encoder produces a dense latent atlas.
   PCA or UMAP is only a visualization of that atlas; scientific meaning comes
   from mapping selected latent neighborhoods back to the sky.

The 2026 paper demonstrates back-mappable organization in MHD density,
nighttime lights, and PHANGS-ALMA NGC 3627 CO emission. It explicitly treats
the latent coordinates as an inspection space, not supervised classes, a
unique phase diagram, or an exhaustive segmentation.

For the NGC 3627 example, the selected configuration used CDD scales
`[2, 4, 8, 16, 32]`, mask multiplier 1.6, a 35-pixel mask cap, ConvNeXt
dilations `[1, 1, 2, 4]`, batch size 1, gradient accumulation 2, ten epochs,
and targets restricted to emission above 3.5 times the RMS.

The 2022 CDD paper is also independently relevant to LANCET. Its component
maps sum back to the input and support a "scale spectrum": the amount of
emission or mass assigned to structures at each scale in a chosen region.
That is a quantitative, interpretable baseline against which the learned
SA-JEPA atlas should be compared.

## Post-discussion UMAP and S/N-mask test

Following Guang-Xing's advice, the current reference field has been replaced
by an exact 1500 by 1500 native-pixel cutout centred at ICRS
`14:44:58.79 -59:48:01.94`. The native sampling is 0.10 arcsec, giving a
150-arcsec (about 2.04-pc) field. The exact cutout is preserved as a new FITS
derivative; delivered data remain unchanged.

For a tractable CPU test, the cutout was resampled to 256 by 256 pixels
(0.586 arcsec or 0.00795 pc per pixel). The CDD coordinates therefore
correspond to approximately 0.016, 0.032, 0.064, 0.127, and 0.255 pc.

The full uncut dust image was used as context. An exploratory S/N target and
embedding mask adopts the project's provisional 3.5 MJy/sr continuum RMS,
estimates a 45.98 MJy/sr background from the crop's 10th percentile, and
requires 3.5-sigma excess. The resulting threshold is 58.23 MJy/sr and retains
81.1% of prepared pixels. This background convention is explicit and must be
replaced by a spatial noise map if one becomes available.

The masked model was trained for five epochs with 16 stochastic masks per
epoch. Mean loss fell from 24.29 to 6.53; final-epoch prediction MSE was 0.0316.
Mean final-epoch prediction/target similarity was 0.933.
The scale-drop fractions were 0.084, 0.114, 0.196, 0.053, and 0.553, making the
0.255-pc component dominant and the 0.064-pc component second. The dense
masked atlas contains 37,588 latent vectors and has effective rank 4.32.

PCA and TorchDR UMAP were fitted to exactly the same S/N-selected vectors.
The first three PCs explain 60.3%, 16.1%, and 8.6% of the variance. UMAP
produces a clearer nonlinear topology, including a distinct low-intensity
hook and a broad main manifold, and its spatial RGB back-map more sharply
outlines the ridge network. It does not remove the brightness sequence:
correlations with log intensity are `[-0.799, 0.062, 0.090]` for PCA and
`[0.887, -0.319, 0.239]` for UMAP. UMAP-1 is therefore more, not less,
intensity ordered; secondary UMAP coordinates are the candidates for testing
morphology beyond brightness.

The spread hinge remains near 0.989 and final context-channel standard
deviation is about 0.0039. The larger masked field provides roughly 19 target
patches per realization and improves sampling, but does not solve the spread
problem. The next controlled comparison should repeat this same trained atlas
with stronger analysis masks and then run a 10-epoch, multi-seed version only
after fixing the preferred mask/noise convention.

## Earlier pre-discussion run

- Input: delivered primary-beam-cut 12-m + 7-m + total-power dust continuum,
  `allscale_combined_v3_dust.MJysr.rp.pbcut.fits`.
- Inputs remained read-only. All derivatives and checkpoints are under
  `results/sa_jepa_continuum_pilot/`.
- The preferred experiment is a 112 by 112 pixel middle-field crop, selected
  from the 532 by 768 prepared overview as `y=98:210, x=330:442`. Its centre
  is approximately RA 14:44:54.99, Dec -59:48:01.9 (Galactic
  `l=316.76005, b=-0.01078`). It deliberately contains the dominant central
  condensation, secondary compact peaks, connecting ridges, and diffuse
  emission.
- At an assumed 2.8 kpc, the crop spans about 1.26 pc per side with a pixel
  scale of 0.01128 pc (0.831 arcsec). The five CDD coordinates therefore span
  approximately 0.023, 0.045, 0.090, 0.181, and 0.361 pc.
- The crop has 96.6% reliable resampling coverage. An 18-pixel erosion limits
  eligible training-target centres to a conservative interior while retaining
  all reliable continuum pixels as encoder context. Dense inference yields
  3,600 finite interior latent vectors after the network footprint is applied.
- Architecture and masking follow the paper's NGC-like setup. The crop run
  used 10 epochs and 32 stochastic masked samples per epoch, with gradient
  accumulation 2. It used upstream SA-JEPA commit
  `3ca0c252d28baee0e673b07154df571f9ae9dd78` on CPU.
- A lower-resolution, whole-cloud 10-epoch run with eight masks per epoch is
  retained as a feasibility comparison, not the preferred result.

## Earlier results and limits

The code path works end to end: CDD decomposition, scale-aware masking,
training, tiled dense inference, latent extraction, PCA back-mapping, and
scale-response diagnostics all completed. On the selected middle crop:

- Mean training loss fell from 21.16 to 4.33 over ten epochs. Final-epoch raw
  latent-prediction MSE is 0.0030 and mean prediction/target similarity is
  0.995. The corresponding overview values were 0.025 and 0.946, so the middle
  crop is a substantially cleaner prediction task.
- The interior 32-channel target atlas has effective rank 4.34. Its first
  three PCs explain 58.1%, 19.5%, and 7.3% of the variance. This is not a
  one-dimensional numerical collapse, but it remains strongly compressed.
- The scale-drop sensitivity fractions are 0.086, 0.111, 0.233, 0.076, and
  0.494. The broadest 32-pixel component (about 0.36 pc) dominates and the
  8-pixel component (about 0.09 pc) is second. This differs from the overview,
  where 8 pixels (about 0.18 pc at its coarser sampling) dominated.
- The focused PCA atlas contains coherent regions around the central
  condensation and surrounding ridge/diffuse transition. The six-group map is
  only a back-mappable inspection aid.
- PC1 correlates strongly with log intensity (`r = -0.779`); PC2 is nearly
  uncorrelated (`r = -0.007`). The dominant organization is therefore still
  partly an intensity sequence, while PC2 is a more interesting candidate for
  morphology tests.
- The spread hinge remains near 0.99 and contributes about 4.95 of the final
  active batches. Context channel standard deviation stays around 0.0022
  against a target of 1. This is the most important warning: increasing from
  80 overview masks to 320 crop masks greatly improves prediction but does not
  improve channel spread. Normalization, spread formulation, and seed/config
  controls therefore need attention.
- The conservative target mask removes primary-beam boundary contamination,
  but the reported latent atlas covers only the central interior. Additional
  overlapping, fully covered crops are preferable to relaxing this safeguard.

Bottom line: this is a successful feasibility run and a useful figure for the
discussion, but not evidence yet that SA-JEPA has discovered a new physical
coordinate in G316.

## Questions for Guang-Xing

1. What is the preferred treatment of an irregular mosaic footprint: zero
   fill plus a validity mask, smoothly tapered/inpainted context, overlapping
   rectangular crops, or training only on fully covered tiles?
2. For a total-power continuum map with a positive diffuse baseline, should
   targets be selected by local S/N, primary-beam response, intensity
   percentile, CDD energy, or a supplied science mask? The NGC 3.5-RMS rule is
   not directly transferable.
3. Since the spread hinge persists with 320 stochastic crop masks, does it
   point mainly to the single-field setup, input normalization/log floor, or a
   need to retune the spread target/weight for continuum data?
4. What sampling count did the reported NGC run use per epoch, beyond batch
   size 1 and accumulation 2? Is the repository default of 200 repeated
   samples the right baseline for a single map?
5. Should CDD scales be chosen in beam units or physical units? For LANCET,
   useful explicit targets are beam/core/clump/filament widths rather than
   powers of two inherited from a resampled overview.
6. Would he recommend training one shared model on matched dust-only,
   dust-plus-free-free, array-only, and all-scale images, or keeping them as
   separate controlled domain-shift tests?
7. Which ablation is essential for a convincing first result: CDD-only scale
   spectra, fixed-box JEPA, wavelet-JEPA, or the plain ConvNeXt control?
8. For compact-source science, should training use high-resolution crops and
   full-map tiled inference, or can a lower-resolution whole-field model be
   meaningfully combined with a second core-scale model?

## Recommended next experiment

1. Keep this middle crop as the reference case, then build several overlapping,
   fully covered, beam-aware crops from the dust-only continuum. Set scales
   from the measured 0.65-arcsec beam through chosen core/clump scales; do not
   let the primary-beam boundary enter a context box.
2. Run mask-multiplier controls at 0.8, 1.2, 1.6, and 2.0 plus one fixed-box or
   plain-ConvNeXt control, using at least three random seeds and enough repeated
   masks for the spread diagnostic to evolve.
3. Select on non-collapse and stability, not maximum effective rank. Compare
   seed-to-seed neighborhood consistency and whether similar morphology maps
   to similar latent neighborhoods across crops.
4. Only after freezing the representation, compare latent coordinates with
   source catalogue properties, column density, temperature, free-free
   fraction, C18O moments, and evolutionary regions. These are validation
   variables, not training labels.
5. Repeat the frozen-model inference on dust-plus-free-free and array-only
   continuum maps after matching grid, beam, units, and valid footprint. This
   directly tests whether latent changes trace astrophysics or observing/
   imaging transfer functions.
6. Move to a C18O moment-0 field next. Full spectral-cube training should wait
   because the repository describes 3D training as under active development.

## Sources

- [ScaleAware-JEPA paper](https://arxiv.org/abs/2606.29723)
- [SA-JEPA code](https://github.com/gxli/SA-JEPA)
- [CDD paper](https://doi.org/10.3847/1538-4365/ac4bc4)
