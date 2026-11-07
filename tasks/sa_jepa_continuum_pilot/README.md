# SA-JEPA continuum pilot

This is an exploratory, non-destructive application of Guang-Xing Li's
ScaleAware-JEPA to the delivered G316 all-scale dust-continuum image. It does
not alter any delivered FITS product. All generated inputs, sessions, logs,
and figures live under `results/sa_jepa_continuum_pilot/`.

The current post-discussion pilot is the exact 1500 by 1500 native-pixel field
centred at ICRS `14:44:58.79 -59:48:01.94`. It spans 150 arcsec, or about
2.04 pc per side at 2.8 kpc. A 256 by 256 training representation is used for
the CPU test. The pilot uses the NGC 3627-like configuration reported in Li
(2026): CDD
scales `[2, 4, 8, 16, 32]`, scale-aware mask multiplier 1.6, mask hard cap 35
pixels, ConvNeXt dilations `[1, 1, 2, 4]`, batch size 1, and gradient
accumulation 2. Low-S/N target and embedding pixels are removed with an
explicit mask, while the complete field remains available as encoder context.
The completed test uses five epochs and 16 stochastic masks per epoch; its
output must not be treated as a final scientific result.

Prepare the scalar field and explicit valid-target mask:

```bash
conda run -n py312 python scripts/sa_jepa_continuum_pilot.py prepare \
  --input data/delivered/alma/TP7MTM2TM1/cont/allscale_combined_v3_dust.MJysr.rp.fits \
  --center-ra 14:44:58.79 --center-dec=-59:48:01.94 \
  --native-crop-size 1500 --max-dimension 256 \
  --snr-min 3.5 --noise-rms 3.5 --background-percentile 10 \
  --output-root results/sa_jepa_continuum_pilot/center1500_snr_umap_v1
```

Train with an independently checked out copy of the upstream repository:

```bash
PYTHONPATH=/path/to/SA-JEPA \
  conda run -n py312 python /path/to/SA-JEPA/scripts/train.py \
  --config tasks/sa_jepa_continuum_pilot/center1500_snr_umap_5epochs.yaml \
  --name lancet_center1500_snr_5epoch \
  --sessions-dir results/sa_jepa_continuum_pilot/center1500_snr_umap_v1/sessions
```

Create matched PCA and TorchDR UMAP back-maps from the selected pixels:

```bash
conda run -n py312 python scripts/sa_jepa_continuum_pilot.py summarize \
  --output-root results/sa_jepa_continuum_pilot/center1500_snr_umap_v1 \
  --session-dir results/sa_jepa_continuum_pilot/center1500_snr_umap_v1/sessions/lancet_center1500_snr_5epoch \
  --summary-name center1500_snr_5epoch \
  --umap --umap-neighbors 50 --umap-min-dist 0.2 --umap-max-iter 500
```

For a same-day CPU smoke test, prepare the whole cloud at half the pilot
resolution and run `smoke.yaml` into a distinct session/output tree:

```bash
conda run -n py312 python scripts/sa_jepa_continuum_pilot.py prepare \
  --max-dimension 384 \
  --output-root results/sa_jepa_continuum_pilot/smoke
PYTHONPATH=/path/to/SA-JEPA \
  conda run -n py312 python /path/to/SA-JEPA/scripts/train.py \
  --config tasks/sa_jepa_continuum_pilot/smoke.yaml \
  --name lancet_continuum_cpu_smoke \
  --sessions-dir results/sa_jepa_continuum_pilot/smoke/sessions
conda run -n py312 python scripts/sa_jepa_continuum_pilot.py summarize \
  --output-root results/sa_jepa_continuum_pilot/smoke \
  --session-dir results/sa_jepa_continuum_pilot/smoke/sessions/lancet_continuum_cpu_smoke
```

Before scientific interpretation, repeat at 10 epochs and run controlled
comparisons against (1) a fixed-box or plain ConvNeXt control, (2) the
array-only continuum image, and (3) at least two random seeds. Compare latent
maps with the existing source catalogue, column density, temperature,
free-free model, and evolutionary-region masks. Unsupervised latent groups are
coordinates for inspection, not physical classes or a source catalogue.

The current result and historical crop/overview tests are compared in
[`discussion_brief.md`](discussion_brief.md). The exact native cutout is saved
as `data/native_crop.fits` inside the current result directory; delivered FITS
products remain unchanged.
